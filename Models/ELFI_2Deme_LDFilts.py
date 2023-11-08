import pandas as pd
import numpy as np
import allel
import msprime
import matplotlib.pyplot as plt
import elfi
import argparse
import sys
import pickle


# print('Parse Arguments')
parser=argparse.ArgumentParser(
description='''Running this code to test different GPSC and number of neutral genes impact on migration parameter estimates ''',
epilog="""Have fun!""")
parser.add_argument('--gpsc', type=int, default=8, help='Add a GPSC number from options 2,8,5,22,26,10 (default 8)')
parser.add_argument('--genes', type=int,default=81, help='Add number of genes either 81 or 355 (default 81)')
parser.add_argument('--true_data', type=str,default='simulated', help='Add either [simulated] or [true] to test with simulated data or get the truth (default simulated)')
parser.add_argument('--country1', type=str,default='sa', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default sa)')
parser.add_argument('--country2', type=str,default='mal', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default mal)')
parser.add_argument('--country3', type=str,default='ken', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (do not include in 2Deme)')
parser.add_argument('--country4', type=str,default='gam', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (do not include in 2Deme)')
parser.add_argument('--input_dir', type=str,default='/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs', help='Path just before GPSCX/ in which the core alignment file is (default local)')
parser.add_argument('--output_dir', type=str,default='/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs', help='Path just before 2Deme/outputs/ inclusive to write summary stats and plots (default local)')
parser.add_argument('--evidence', type=int, default=200, help='Number of evidence points for BOLFI (default 200)')
parser.add_argument('--sample', type=int, default=2000, help='Number of samples for BOLFI (default 2000)')
parser.add_argument('--bounds', type=int, default=1, help='Upper bound for parameter (default 1)')
parser.add_argument('--suffix', type=str, default="", help='Suffix for saved files (default "")')
parser.add_argument('--sampler', type=str, default='metropolis', help='Type of sampler. Either [metropolis] or [nuts] (default metropolis)')
parser.add_argument('--prior', type=str, default='exponential', help='Prior distribution. Either [exponential] or [uniform] (default exponential)')
parser.add_argument('--vcf_in', type=str, default='test.vcf', help='Full name of VCF file for input')


args = vars(parser.parse_args())

lins = args["gpsc"]
cnt = args["genes"]
data4yobs = args["true_data"]
input_dir1 = args["input_dir"]
output_dir1 = args["output_dir"]
country1_in = args["country1"]
country2_in = args["country2"]
country3_in = args["country3"]
country4_in = args["country4"]
evs = args["evidence"]
samps = args["sample"]
ub = args["bounds"]
suf= args["suffix"]
sampler = args["sampler"]
prior = args["prior"]
vcf_in = args["vcf_in"]
print('Prior Distribution:' +str(prior) )


# gpsc1=8
# genes1=81
# cat1='simulated'
# print('Define functions')
print('GPSC: '+str(lins))
print('Genes: '+str(cnt))
print('Data: '+str(data4yobs))
print('Input Directory:'+str(input_dir1))
print('Output Directory:'+str(output_dir1))
print(str(country1_in), str(country2_in), str(country3_in), str(country4_in))

# This was for when I was running the function within python rather than calling the python script
# def estimates_2demes(gpsc=gpsc1,genes=genes1,cat=cat1):
# lins=8 ## gpsc
# cnt=81 ### 355 or 81 genes
# data4yobs=true## 

## Define the initial population size depending on which countries being input
if country1_in=='sa':
    intpopA=6000
if country1_in=='mal':
    intpopA=2000
if country1_in=='ken':
    intpopA=5000
if country1_in=='gam':
    intpopA=1000

if country2_in=='sa':
    intpopB=6000
if country2_in=='mal':
    intpopB=2000
if country2_in=='ken':
    intpopB=5000
if country2_in=='gam':
    intpopB=1000

###Set up for true Fst
####Read in VCF and index files
# def get_callset_4demes(gpsc=lins,genes=cnt,c1='sa',c2='mal',c3='ken',c4='gam'):
def get_callset_4demes(gpsc=lins,genes=cnt,c1=country1_in,c2=country2_in,c3=country3_in,c4=country4_in):
    header=genes##gene number
    var=gpsc###GPSC
    # callset = allel.read_vcf('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/core_'+str(header)+'_gpsc'+str(var)+'_alignment.snp.aln.biallelic.vcf')
    callset = allel.read_vcf(str(input_dir1)+'/GPSC'+str(var)+'/' + str(vcf_in))
    print('Read in alignment')
#     # define empty lists
#     c1_list=[]
#     c2_list=[]
#     c3_list=[]
#     c4_list=[]
    # open file and read the content in a list
    ## south africa
    with open(str(input_dir1)+'/GPSC'+str(var)+'/index/'+str(c1)+'_'+str(header)+'_index.csv', 'r') as filehandle:
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/index/'+str(c1)+'_'+str(header)+'_index.csv', 'r') as filehandle:
        c1_list = [current_place.rstrip() for current_place in filehandle.readlines()]
    c1_list=[int(x)-1 for x in c1_list]#convert strings to integers for each element
    ##malawi
    with open(str(input_dir1)+'/GPSC'+str(var)+'/index/'+str(c2)+'_'+str(header)+'_index.csv', 'r') as filehandle:
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/index/'+str(c2)+'_'+str(header)+'_index.csv', 'r') as filehandle:
        c2_list = [current_place.rstrip() for current_place in filehandle.readlines()]
    c2_list=[int(x)-1 for x in c2_list]#convert strings to integers for each element
    ## kenya
    with open(str(input_dir1)+'/GPSC'+str(var)+'/index/'+str(c3)+'_'+str(header)+'_index.csv', 'r') as filehandle:
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/index/'+str(c3)+'_'+str(header)+'_index.csv', 'r') as filehandle:
        c3_list = [current_place.rstrip() for current_place in filehandle.readlines()]
    c3_list=[int(x)-1 for x in c3_list]#convert strings to integers for each element
    ## the gambia
    with open(str(input_dir1)+'/GPSC'+str(var)+'/index/'+str(c4)+'_'+str(header)+'_index.csv', 'r') as filehandle:
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/index/'+str(c4)+'_'+str(header)+'_index.csv', 'r') as filehandle:
        c4_list = [current_place.rstrip() for current_place in filehandle.readlines()]
    c4_list=[int(x)-1 for x in c4_list]#convert strings to integers for each element
#     if sample==True:
#         sa1_sp=list(np.random.choice(c1_list,50,replace=True))
#         mal1_sp=list(np.random.choice(c2_list,50,replace=True))
#         ken1_sp=list(np.random.choice(c3_list,50,replace=True))
#         gam1_sp=list(np.random.choice(c4_list,50,replace=True))
#     else:
#         sa1_sp=c1_list
#         mal1_sp=c2_list
#         ken1_sp=c3_list
#         gam1_sp=c4_list
    return callset,c1_list,c2_list,c3_list,c4_list

###set up function for True Fst
def true_fst(callset,country1,country2,country3=None,country4=None,sample=True): 
    def pairwise_fst(country1,country2,callset):
        country1_sp=list(np.random.choice(country1,50,replace=True))
        country2_sp=list(np.random.choice(country2,50,replace=True))
        gt = allel.GenotypeArray(callset['calldata/GT'])
        gt_sma = gt[:,country1_sp+country2_sp]
        subpops_sma=[list(range(0,len(country1_sp))),list(range(len(country1_sp),len(country1_sp)+len(country2_sp)))]  
        ac1 = gt_sma.count_alleles(subpop=subpops_sma[0])
        ac2 = gt_sma.count_alleles(subpop=subpops_sma[1])
        num, den = allel.hudson_fst(ac1, ac2)
        fst=num/den
        obsa=np.nanmean(fst)
        obsb=np.nanstd(fst)
        return obsa
    if sample==True:
        reps=500
        obs12_mean=np.zeros(reps)
        obs13_mean=np.zeros(reps)
        obs14_mean=np.zeros(reps)
        obs23_mean=np.zeros(reps)
        obs24_mean=np.zeros(reps)
        obs34_mean=np.zeros(reps)
        for i in np.arange(reps):
            obs12_mean[i] = pairwise_fst(country1,country2,callset)
            if country3 is not None:
                obs13_mean[i] = pairwise_fst(country1,country3,callset)
                obs23_mean[i] = pairwise_fst(country2,country3,callset)
            if country4 is not None:
                obs14_mean[i]= pairwise_fst(country1,country4,callset)
                obs24_mean[i]= pairwise_fst(country2,country4,callset)
                obs34_mean[i]= pairwise_fst(country3,country4,callset)
        # obs12fst=np.nanmean(obs12_mean)
        obs12fst=np.quantile(obs12_mean,np.arange(0.1,1,0.1),axis=None)
        obs12fst_std=np.nanstd(obs12_mean)
        obs13fst=np.nanmean(obs13_mean)
        obs14fst=np.nanmean(obs14_mean)
        obs23fst=np.nanmean(obs23_mean)
        obs24fst=np.nanmean(obs24_mean)
        obs34fst=np.nanmean(obs34_mean)

        if country3 is None:
            # return np.array([obs12fst,obs12fst_std]).reshape(1, -1)
            return np.array([obs12fst]).reshape(1, -1)
        if country3 is not None and country4 is None:
            return np.array([obs12fst,obs13fst,obs23fst]).reshape(1, -1)
        if country3 is not None and country4 is not None:
            return np.array([obs12fst,obs13fst,obs14fst,obs23fst,obs24fst,obs34fst]).reshape(1, -1)
        
    else:
        gt = allel.GenotypeArray(callset['calldata/GT'])
        gt_sma = gt[:,country1+country2]
        subpops_sma=[list(range(0,len(country1))),list(range(len(country1),len(country1)+len(country2)))]  
        ac1 = gt_sma.count_alleles(subpop=subpops_sma[0])
        ac2 = gt_sma.count_alleles(subpop=subpops_sma[1])
        num, den = allel.hudson_fst(ac1, ac2)
        fst=num/den
        obsa=np.nanmean(fst)
        obsb=np.nanstd(fst)
# ### Test it works
# callset,country1,country2,country3,country4 = get_callset_4demes()
# true_fst(callset,country1,country2,sample=True)

####Function for Simulated Migration Rate
def simulate_migration(mig_ab,mig_ba, batch_size=1, random_state=None):
    mig_ab = mig_ab / 5000
    mig_ba = mig_ba / 5000
    # mig_ab = mig_ab / 1000
    # mig_ba = mig_ba / 1000
    theta=2.5e-5
    num_replicates=500
    seq_length=500
    n_sampA=600
    n_sampB=600
    
    demography = msprime.Demography()
    demography.add_population(name="A",initial_size=intpopA,
                            description="South Africa")
    demography.add_population(name="B",initial_size=intpopB,
                            description="Malawi")

    demography.set_migration_rate(source="A", dest="B", rate=mig_ab)
    demography.set_migration_rate(source="B", dest="A", rate=mig_ba)

    def sim_replicates(sample_size, num_replicates, demography, theta=0):
        ancestry_reps = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=seq_length,
            ploidy = 1,
            num_replicates=num_replicates,
            record_migrations=True)
        
        for ts in ancestry_reps:
            mutated_ts = msprime.sim_mutations(ts, rate=theta, discrete_genome=False)
            yield mutated_ts
    
    def get_hudsons_fst(reps):
        data1 = np.zeros(num_replicates)
        data2 = np.zeros(num_replicates)
        for rep_index, ts in enumerate(reps):
            data1[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=0),ts.samples(population=1)], windows = "sites"))
            data2[rep_index]=np.nanstd(ts.Fst(sample_sets=[ts.samples(population=0),ts.samples(population=1)], windows = "sites"))
        return data1, data2

    truth1, truth_social = get_hudsons_fst(sim_replicates(sample_size={"A":n_sampA, "B":n_sampB}, 
                                        num_replicates=num_replicates, 
                                        demography=demography, 
                                        theta=theta))

    # truemean1=np.nanmean(truth1,axis=None)
    truemean1=np.quantile(truth1,np.arange(0.1,1,0.1),axis=None)

    # truemean2=np.nanmean(truth_social,axis=None)
    return np.array((truemean1)).reshape(1,-1)
    # return np.array((truemean1,truemean2)).reshape(1,-1)
print('Set Data: '+str(data4yobs))
if data4yobs=='simulated':
    ## Simulated
    mig_ab_true=0.1
    mig_ba_true=0.6
    yobs = simulate_migration(mig_ab=mig_ab_true,mig_ba=mig_ba_true)
if data4yobs=='true':            
    ## TRUE
    callset,country1,country2,country3,country4 = get_callset_4demes(lins,cnt)
    yobs=true_fst(callset,country1,country2,sample=True)

print(yobs)
model = elfi.new_model()
if prior=='exponential':
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_ba', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_ab', model=model)
if prior=='uniform':
    elfi.Prior('uniform',0, ub, name='mig_ba', model=model)
    elfi.Prior('uniform',0, ub, name='mig_ab', model=model)
# elfi.Prior('exponential',0.000001, 1.0000, name='mig_ba', model=model)
# elfi.Prior('exponential',0.000001, 1.0000, name='mig_ab', model=model)
# elfi.Prior('exponential',1,10, name='mig_ba', model=model)
# elfi.Prior('exponential',1,10, name='mig_ab', model=model)
# print('Set priors')

elfi.Simulator(simulate_migration, model['mig_ab'],model['mig_ba'], name='SFS', observed=yobs)

elfi.Distance('euclidean', model['SFS'], name='d')
bolfi = elfi.BOLFI(model['d'], initial_evidence=50, bounds={'mig_ab': (0.0, ub),'mig_ba': (0.0, ub)}, acq_noise_var=0.1)

print('Fitting BOLFI with '+str(evs)+' evidence')
bolfi.fit(n_evidence=evs)

if data4yobs=='simulated':
    print('Sampling '+str(samps))
    sys.stdout = open(str(output_dir1)+'/2Deme/sims/samples_stdout'+str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.txt', "w")
    # sys.stdout = open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/sims/samples_stdout.txt', "w")
    if (sampler=='metropolis'):
        samples = bolfi.sample(samps,algorithm='metropolis',sigma_proposals={'mig_ab': 0.05,'mig_ba': 0.05})
    if (sampler=='nuts'):
        samples = bolfi.sample(samps)
    

if data4yobs=='true':
    # print('Sampling')
    # sys.stdout = open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/truth/samples_stdout.txt', "w")
    sys.stdout = open(str(output_dir1)+'/2Deme/GPSC'+str(lins)+'/truth/samples_stdout'+str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.txt', "w")
    # samples = bolfi.sample(samps,algorithm='metropolis',sigma_proposals={'mig_ab': 0.05,'mig_ba': 0.05})
    if (sampler=='metropolis'):
        samples = bolfi.sample(samps,algorithm='metropolis',sigma_proposals={'mig_ab': 0.05,'mig_ba': 0.05})
    if (sampler=='nuts'):
        samples = bolfi.sample(samps)

if data4yobs=='simulated':
    bolfi.plot_gp(true_params={'mig_ab': mig_ab_true, 'mig_ba': mig_ba_true});
    plt.savefig(str(output_dir1)+'/2Deme/sims/contour2DemeGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.pdf' )
    # plt.savefig('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/sims/contour2DemeGPSC'+str(lins)+'_'+ str(cnt)+'.pdf' )
    # ###asymmetric rate
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    plt.hist2d(samples.samples['mig_ab'],samples.samples['mig_ba'], 20)
    fig.suptitle("GPSC"+str(lins)+" with "+str(cnt)+" genes - "+str(data4yobs))
    plt.savefig(str(output_dir1)+'/2Deme/sims/hist2DemeGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.pdf' )
    # plt.savefig('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/sims/hist2DemeGPSC'+str(lins)+'_'+ str(cnt)+'.pdf' )

if data4yobs=='true':
    # ###asymmetric rate
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    plt.hist2d(samples.samples['mig_ab'],samples.samples['mig_ba'], 20)
    fig.suptitle("GPSC"+str(lins)+" with "+str(cnt)+" genes - "+str(data4yobs))
    # plt.savefig('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/truth/hist2DemeGPSC'+str(lins)+'_'+ str(cnt)+'.pdf' )
    plt.savefig(str(output_dir1)+'/2Deme/GPSC'+str(lins)+'/truth/hist2DemeGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.pdf' )


if data4yobs=='simulated':
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/sims/samples_means_summary.txt', 'w') as f:
    with open(str(output_dir1)+'/2Deme/sims/samples_means_summary'+str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.txt', 'w') as f:
        print(samples.sample_means_summary)
        f.write(str(samples.sample_means_summary))
if data4yobs=='true':
    # with open('/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs/2Deme/GPSC'+str(lins)+'/truth/samples_means_summary.txt', 'w') as f:
    with open(str(output_dir1)+'/2Deme/GPSC'+str(lins)+'/truth/samples_means_summary'+str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.txt', 'w') as f:
        print(samples.sample_means_summary)
        f.write(str(samples.sample_means_summary))

if data4yobs=='true':
    with open(str(output_dir1)+'/2Deme/GPSC'+str(lins)+'/truth/samples_posteriors'+str(cnt)+'_'+str(country1_in)+'_'+str(country2_in)+str(ub)+str(suf)+'.pickle', 'wb') as samples_file:
        pickle.dump(samples, samples_file)
if data4yobs=='simulated':
    with open(str(output_dir1)+'/2Deme/sims/samples_posteriors'+str(cnt)+'_'+str(ub)+str(suf)+'.pickle', 'wb') as samples_file:
        pickle.dump(samples, samples_file)   


print('Complete')


# estimates_2demes(gpsc=8,genes=81,cat='simulated')






# %%
