import os

import networkx as nx

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_core_gene_nodes(G, threshold, num_isolates):
    #Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        if float(G.nodes[node]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(alignments_dir, core_names, output_dir,
                                       output_prefix):
    #Open up each alignment that is assosciated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    #Read in all these alginemnts
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, 'fasta')
        gene_dict = {}
        for record in alignment:
            genome_id = record.id.split(";")[0]
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(
                        gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    #Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    #Write out the two output files
    SeqIO.write(isolate_aln, output_dir + output_prefix + '_alignment.aln', 'fasta')

    write_alignment_header(gene_alignments, output_dir, output_prefix)
    return core_filenames

def write_alignment_header(alignment_list, outdir, prefix):
    out_entries = []
    #Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        #Get length and name from one sequence in the alignment
        #Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        #Create the 3 line feature entry
        gene_entry1 = "FT   feature         " + str(gene_start) + ".." + str(
            gene_end) + '\n'
        gene_entry2 = "FT                   /label=" + gene_name + '\n'
        gene_entry3 = "FT                   /locus_tag=" + gene_name + '\n'
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        #Add it to the output list
        out_entries.append(gene_entry)
        #Alter the post-output variables
        gene_start += gene[2]
    #Create the header and footer
    header = ("ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key" +
              "             Location/Qualifiers\nFH\n")
    footer = ("XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" +
              " 0 other;\n//\n")
    #open file and output
    with open(outdir + prefix+"_alignment_header.embl", "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)
    return True



if __name__ == '__main__':
    import argparse
    description = 'Concatenate recombination-free core gene alignmets'
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory",
                         )
    parser.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)
    parser.add_argument("--filtered-core",
                        dest="prefiltered",
                        default=None,
                        help="""Provide a folder in the panaroo output 
                        containing a filtered set of core genes""")
    parser.add_argument("--prefix",
                        dest="prefix",
                        default="filtered_core_gene",
                        help="Provide outprefix")
    args = parser.parse_args()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    
    # Load isolate names
    seen = set()
    isolate_names = []
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            iso = line.split(",")[0]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)
    
    if type(args.prefiltered) == str:
        files = os.listdir(args.prefiltered)
        filt_names = [x.split('.')[0] for x in files]
        concatenate_core_genome_alignments(args.prefiltered, 
                                           filt_names, args.output_dir, args.prefix)
        
    else:
        #load graph
        G = nx.read_gml(args.output_dir + "final_graph.gml")
    
        #Identify core
        core_nodes = get_core_gene_nodes(G, args.core, len(isolate_names))
        core_names = [G.nodes[x]["name"] for x in core_nodes]
        concatenate_core_genome_alignments(args.output_dir + "recombination_free_aligned_genes/",
                                           core_names, args.output_dir,
                                           "recombination_free_core_gene")
