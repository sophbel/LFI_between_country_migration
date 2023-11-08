###This is the same ELFI 4 Deme script but with a second stat (SS) of the Ystatistic and symmetric migration between 4 demes

import pandas as pd
import numpy as np
import allel
import msprime
import matplotlib.pyplot as plt
import elfi
import argparse
import sys
import seaborn as sns
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
parser.add_argument('--country3', type=str,default='ken', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default ken)')
parser.add_argument('--country4', type=str,default='gam', help='Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default gam)')
parser.add_argument('--input_dir', type=str,default='/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs', help='Path just before GPSCX/ in which the core alignment file is (default local)')
parser.add_argument('--output_dir', type=str,default='/Users/sb62/Documents/Migration/DemographyModel/RunModel/outputs', help='Path just before 2Deme/outputs/ inclusive to write summary stats and plots (default local)')
parser.add_argument('--evidence', type=int, default=4000, help='Number of evidence points for BOLFI (default 4000)')
parser.add_argument('--sample', type=int, default=10000, help='Number of samples for BOLFI (default 10000)')
parser.add_argument('--bounds', type=int, default=3, help='Upper bound for parameter (default 3)')
parser.add_argument('--suffix', type=str, default="", help='Suffix string (default "")')
parser.add_argument('--sampler', type=str, default='metropolis', help='Type of sampler. Either [metropolis] or [nuts] (default metropolis)')
parser.add_argument('--prior', type=str, default='exponential', help='Prior distribution. Either [exponential] or [uniform] (default exponential)')
parser.add_argument('--epistasis', type=bool, default=False, help='Indicate whether epistatic sites within 1kb with >0.5r2 have been removed')

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
suf = args["suffix"]
sampler = args["sampler"]
prior = args["prior"]
epistasis = args["epistasis"]



print('GPSC: '+str(lins))
print('Genes: '+str(cnt))
print('Data: '+str(data4yobs))
print('Input Directory:'+str(input_dir1))
print('Output Directory:'+str(output_dir1))
print(str(country1_in), str(country2_in), str(country3_in), str(country4_in))
print('Prior Distribution:' +str(prior) )


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

if country3_in=='sa':
    intpopC=6000
if country3_in=='mal':
    intpopC=2000
if country3_in=='ken':
    intpopC=5000
if country3_in=='gam':
    intpopC=1000

if country4_in=='sa':
    intpopD=6000
if country4_in=='mal':
    intpopD=2000
if country4_in=='ken':
    intpopD=5000
if country4_in=='gam':
    intpopD=1000

# This was for when I was running the function within python rather than calling the python script
# def estimates_2demes(gpsc=gpsc1,genes=genes1,cat=cat1):
#   lins=gpsc1 ## gpsc
#   cnt=genes1 ### 355 or 81 genes
#   data4yobs=cat1## 

###Set up for true Fst
####Read in VCF and index files
# def get_callset_4demes(gpsc=lins,genes=cnt,c1='sa',c2='mal',c3='ken',c4='gam'):
def get_callset_4demes(gpsc=lins,genes=cnt,c1=country1_in,c2=country2_in,c3=country3_in,c4=country4_in):
    header=genes##gene number
    var=gpsc###GPSC
    # callset = allel.read_vcf('/Users/sb62/Documents/Migration/DemographyModel/RunModel/inputs/GPSC'+str(var)+'/core_'+str(header)+'_gpsc'+str(var)+'_alignment.snp.aln.biallelic.vcf')
    if epistasis==False:
        callset = allel.read_vcf(str(input_dir1)+'/GPSC'+str(var)+'/core_'+str(header)+'_gpsc'+str(var)+'_alignment.snp.aln.biallelic.vcf')
    if epistasis==True:
        callset = allel.read_vcf(str(input_dir1)+'/GPSC'+str(var)+'/test.filter_r2_5_1kb.vcf')

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

    return callset,c1_list,c2_list,c3_list,c4_list

###set up function for True Fst
def true_fst(callset,country1,country2,country3=None,country4=None,sample=True):
    def pairwise_fst(country1,country2,callset):
        country1_sp=list(np.random.choice(country1,30,replace=True))
        country2_sp=list(np.random.choice(country2,30,replace=True))
        gt = allel.GenotypeArray(callset['calldata/GT'])
        gt_sma = gt[:,country1_sp+country2_sp]
        subpops_sma=[list(range(0,len(country1_sp))),list(range(len(country1_sp),len(country1_sp)+len(country2_sp)))]  
        ac1 = gt_sma.count_alleles(subpop=subpops_sma[0])
        ac2 = gt_sma.count_alleles(subpop=subpops_sma[1])
        num, den = allel.hudson_fst(ac1, ac2)
        fst=num/den
        obsa=np.nanmean(fst)
        return obsa 
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

    obs12fst=np.quantile(obs12_mean,np.arange(0.1,1,0.1),axis=None)
    obs13fst=np.quantile(obs13_mean,np.arange(0.1,1,0.1),axis=None)
    obs14fst=np.quantile(obs14_mean,np.arange(0.1,1,0.1),axis=None)
    obs23fst=np.quantile(obs23_mean,np.arange(0.1,1,0.1),axis=None)
    obs24fst=np.quantile(obs24_mean,np.arange(0.1,1,0.1),axis=None)
    obs34fst=np.quantile(obs34_mean,np.arange(0.1,1,0.1),axis=None)

    # obs12fst=np.nanmean(obs12_mean)
    # obs13fst=np.nanmean(obs13_mean)
    # obs14fst=np.nanmean(obs14_mean)
    # obs23fst=np.nanmean(obs23_mean)
    # obs24fst=np.nanmean(obs24_mean)
    # obs34fst=np.nanmean(obs34_mean)

    if country3 is None:
        return np.array(obs12fst).reshape(1, -1)
    if country3 is not None and country4 is None:
        return np.array([obs12fst,obs13fst,obs23fst]).reshape(1, -1)
    if country3 is not None and country4 is not None:
        return np.array([obs12fst,obs13fst,obs14fst,obs23fst,obs24fst,obs34fst]).reshape(1, -1)
        

### Test it works
# callset,country1,country2,country3,country4 = get_callset_4demes()
# true_fst(callset,country1,country2,sample=True)

####Function for Simulated Migration Rate
def simulate_migration(mig_ab, 
                       mig_ac,
                       mig_ad,
                       mig_bc,
                       mig_bd,
                       mig_dc,
                       batch_size=1, random_state=None):
    mig_ab = mig_ab / 5000
    mig_ac = mig_ac / 5000
    mig_ad = mig_ad / 5000
    mig_bc = mig_bc / 5000
    mig_bd = mig_bd / 5000
    mig_dc = mig_dc / 5000
    
    theta=2.5e-5
    num_replicates=800
    seq_length=500
    n_sampA=600
    n_sampB=600
    n_sampC=600
    n_sampD=600

    
    demography = msprime.Demography()
    demography.add_population(name="A",initial_size=intpopA,
                             description="South Africa")
    demography.add_population(name="B",initial_size=intpopB,
                             description="Malawi")
    demography.add_population(name="C",initial_size=intpopC,
                             description="Kenya")
    demography.add_population(name="D",initial_size=intpopD,
                             description="The Gambia")

    demography.set_symmetric_migration_rate(["A", "B"], rate=mig_ab)

    demography.set_symmetric_migration_rate(["A", "C"], rate=mig_ac)
    
    demography.set_symmetric_migration_rate(["A", "D"], rate=mig_ad)
    
    demography.set_symmetric_migration_rate(["B", "C"], rate=mig_bc)
        
    demography.set_symmetric_migration_rate(["B", "D"], rate=mig_bd)
        
    demography.set_symmetric_migration_rate(["D", "C"], rate=mig_dc)


    def sim_replicates(sample_size, num_replicates, demography, theta):
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
        data3 = np.zeros(num_replicates)
        data4 = np.zeros(num_replicates)
        data5 = np.zeros(num_replicates)
        data6 = np.zeros(num_replicates)

        for rep_index, ts in enumerate(reps):
            data1[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=0),ts.samples(population=1)]))
            data2[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=0),ts.samples(population=2)]))
            data3[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=0),ts.samples(population=3)]))
            data4[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=1),ts.samples(population=2)]))
            data5[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=1),ts.samples(population=3)]))
            data6[rep_index]=np.nanmean(ts.Fst(sample_sets=[ts.samples(population=2),ts.samples(population=3)]))

        return data1, data2, data3, data4, data5, data6

    truth1, truth2, truth3, truth4, truth5, truth6 = get_hudsons_fst(sim_replicates(sample_size={"A":n_sampA, "B":n_sampB,"C":n_sampC, "D":n_sampD}, 
                                           num_replicates=num_replicates, 
                                           demography=demography, 
                                           theta=theta))
    truemean1=np.quantile(truth1,np.arange(0.1,1,0.1),axis=None)
    truemean2=np.quantile(truth2,np.arange(0.1,1,0.1),axis=None)
    truemean3=np.quantile(truth3,np.arange(0.1,1,0.1),axis=None)
    truemean4=np.quantile(truth4,np.arange(0.1,1,0.1),axis=None)
    truemean5=np.quantile(truth5,np.arange(0.1,1,0.1),axis=None)
    truemean6=np.quantile(truth6,np.arange(0.1,1,0.1),axis=None)





    return np.array((truemean1, truemean2, truemean3, truemean4,truemean5,truemean6)).reshape(1,-1)


print('Set Data: '+str(data4yobs))
if data4yobs=='simulated':
    
    ## Simulated - testing low, med, high
    mig_ab_true = 2.5

    mig_ac_true = 2.5

    mig_ad_true = 2.5

    mig_bc_true = 1.5

    mig_bd_true = 1.5

    mig_cd_true = 0.5

 
    yobs = simulate_migration(mig_ab_true,
                        mig_ac_true,
                        mig_ad_true,
                        mig_bc_true,
                        mig_bd_true,
                        mig_cd_true)

if data4yobs=='true':            
    ## TRUE
    callset,country1,country2,country3,country4 = get_callset_4demes(lins,cnt)
    yobs=true_fst(callset,country1,country2,country3,country4,sample=True)

print(yobs)
model = elfi.new_model()

if prior=='exponential':
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_ab', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_ac', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_ad', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_bc', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_bd', model=model)
    elfi.Prior('exponential',0.000001, 1.0000, name='mig_cd', model=model)

if prior=='uniform':
    elfi.Prior('uniform',0, ub, name='mig_ab', model=model)
    elfi.Prior('uniform',0, ub, name='mig_bc', model=model)
    elfi.Prior('uniform',0, ub, name='mig_ac', model=model)

    elfi.Prior('uniform',0, ub, name='mig_bd', model=model)
    elfi.Prior('uniform',0, ub, name='mig_cd', model=model)
    elfi.Prior('uniform',0, ub, name='mig_ad', model=model)

elfi.Simulator(simulate_migration, model['mig_ab'],
               model['mig_bc'],
               model['mig_ac'],
               model['mig_ad'],
               model['mig_bd'],
               model['mig_cd'], name='SFS', observed=yobs)

elfi.Distance('euclidean', model['SFS'], name='d')

bolfi = elfi.BOLFI(model['d'], initial_evidence=200, bounds={'mig_ab': (0.0, ub),
                                                           'mig_bc': (0.0, ub),
                                                           'mig_ac': (0.0, ub),
                                                           'mig_ad': (0.0, ub),
                                                           'mig_bd': (0.0, ub),
                                                           'mig_cd': (0.0, ub)}, acq_noise_var=0.1)


print('Fitting BOLFI with '+str(evs)+' evidence')
bolfi.fit(n_evidence=evs)

if data4yobs=='simulated':
    print('Sampling '+str(samps))
    sys.stdout = open(str(output_dir1)+'/4Demes/sims/samples_stdout'+str(cnt)+'_'+str(ub)+str(suf)+'.txt', "w")
    if (sampler=='metropolis'):
        samples = bolfi.sample(samps,algorithm='metropolis',sigma_proposals={'mig_ab': 0.05,
                                                                'mig_bc': 0.05,
                                                                'mig_ac': 0.05,
                                                                'mig_ad': 0.05,
                                                                'mig_bd': 0.05,
                                                                'mig_cd': 0.05}, n_chains=10)
    if (sampler=='nuts'):
        samples = bolfi.sample(samps)

if data4yobs=='true':
    # print('Sampling')
    sys.stdout = open(str(output_dir1)+'/4Demes/GPSC'+str(lins)+'/truth/samples_stdout'+str(cnt)+'_'+str(ub)+str(suf)+'.txt', "w")
    if (sampler=='metropolis'):
        samples = bolfi.sample(samps,algorithm='metropolis',sigma_proposals={'mig_ab': 0.05,
                                                                'mig_bc': 0.05,
                                                                'mig_ac': 0.05,
                                                                'mig_ad': 0.05,
                                                                'mig_bd': 0.05,
                                                                'mig_cd': 0.05}, n_chains=10)
    if (sampler=='nuts'):
        samples = bolfi.sample(samps)

if data4yobs=='simulated':
    bolfi.plot_gp(true_params={'mig_ab': mig_ab_true,
                          'mig_bc': mig_bc_true,
                          'mig_ac': mig_ac_true,
                          'mig_ad': mig_ad_true,
                          'mig_bd': mig_bd_true,
                          'mig_cd': mig_cd_true});
    plt.savefig(str(output_dir1)+'/4Demes/sims/contour4DemesGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(ub)+str(suf)+'.pdf' )
    # ###asymmetric rate
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    # sns.kdeplot(samples.samples['mig_ba'])
    sns.kdeplot(samples.samples['mig_ab'],alpha=0.5)
    sns.kdeplot(samples.samples['mig_ac'])
    sns.kdeplot(samples.samples['mig_ad'])
    sns.kdeplot(samples.samples['mig_bc'])
    sns.kdeplot(samples.samples['mig_bd'])
    sns.kdeplot(samples.samples['mig_cd'])
    plt.legend(labels=[ 'mig_ab','mig_ac','mig_ad','mig_bc','mig_bd','mig_cd'])    
    fig.suptitle("GPSC"+str(lins)+" with "+str(cnt)+" genes - "+str(data4yobs))
    plt.savefig(str(output_dir1)+'/4Demes/sims/hist4DemesGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(ub)+str(suf)+'.pdf' )

if data4yobs=='true':
    # ###asymmetric rate
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    sns.kdeplot(samples.samples['mig_ab'],alpha=0.5)
    sns.kdeplot(samples.samples['mig_ac'])
    sns.kdeplot(samples.samples['mig_ad'])
    sns.kdeplot(samples.samples['mig_bc'])
    sns.kdeplot(samples.samples['mig_bd'])
    sns.kdeplot(samples.samples['mig_cd'])
    plt.legend(labels=[ 'mig_ab','mig_ac','mig_ad','mig_bc','mig_bd','mig_cd'])    
    fig.suptitle("GPSC"+str(lins)+" with "+str(cnt)+" genes - "+str(data4yobs))
    plt.savefig(str(output_dir1)+'/4Demes/GPSC'+str(lins)+'/truth/hist4DemesGPSC'+str(lins)+'_'+ str(cnt)+'_'+str(ub)+str(suf)+'.pdf' )


if data4yobs=='simulated':
    with open(str(output_dir1)+'/4Demes/sims/samples_means_summary'+str(cnt)+'_'+str(ub)+str(suf)+'.txt', 'w') as f:
        print(samples.sample_means_summary)
        f.write(str(samples.sample_means_summary))
    print('mig_ab: '+str(mig_ab_true))
    # print('mig_ba: '+str(mig_ba_true))
    print('mig_ac: '+str(mig_ac_true))
    # print('mig_ca: '+str(mig_ca_true))
    print('mig_ad: '+str(mig_ad_true))
    # print('mig_da: '+str(mig_da_true))
    print('mig_bc: '+str(mig_bc_true))
    # print('mig_cb: '+str(mig_cb_true))
    print('mig_bd: '+str(mig_bd_true))
    # print('mig_db: '+str(mig_db_true))
    print('mig_cd: '+str(mig_cd_true))
    # print('mig_dc: '+str(mig_dc_true))

if data4yobs=='true':
    with open(str(output_dir1)+'/4Demes/GPSC'+str(lins)+'/truth/samples_means_summary'+str(cnt)+'_'+str(ub)+str(suf)+'.txt', 'w') as f:
        print(samples.sample_means_summary)
        f.write(str(samples.sample_means_summary))
if data4yobs=='true':
    with open(str(output_dir1)+'/4Demes/GPSC'+str(lins)+'/truth/samples_posteriors'+str(cnt)+'_'+str(ub)+str(suf)+'.pickle', 'wb') as samples_file:
        pickle.dump(samples, samples_file)
if data4yobs=='simulated':
    with open(str(output_dir1)+'/4Demes/sims/samples_posteriors'+str(cnt)+'_'+str(ub)+str(suf)+'.pickle', 'wb') as samples_file:
        pickle.dump(samples, samples_file)       
if data4yobs=='true':
    with open(str(output_dir1)+'/4Demes/GPSC'+str(lins)+'/truth/bolfi_posteriors'+str(cnt)+'_'+str(ub)+str(suf)+'.pickle', 'wb') as bolfi_file:
        pickle.dump(bolfi, bolfi_file)
if data4yobs=='simulated':
    with open(str(output_dir1)+'/4Demes/sims/bolfi_posteriors'+str(cnt)+'_'+str(ub)+str(suf)+'.pickle', 'wb') as bolfi_file:
        pickle.dump(bolfi, bolfi_file)


print('Complete')








# %%
