import pandas as pd
import numpy as np
import pickle
import elfi
import matplotlib.pyplot as plt
gpscs = [10,2,22,26,5,8]
# gpscs=[8]
deme_vec=['sa_mal','sa_ken','sa_gam','mal_ken','mal_gam','gam_ken']
gpsc_meanlist=[]
for var in enumerate(gpscs):
    print(var[1])
    deme_list = []
    sym_list = []
    for d in enumerate(deme_vec):
        print('1) Concatenating posteriors to estimate mean migration by GPSC for GPSC'+str(var[1])+str(d[1]))
        ## chain1
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(var[1])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(var[1])+"_chain1.pickle",'rb') ## uniform
        object_file_chain1 = pickle.load(file)
        file.close()
        ## chain2
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(var[1])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(var[1])+"_chain2.pickle",'rb') ## uniform
        object_file_chain2 = pickle.load(file)
        file.close()
        ## chain3
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(var[1])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(var[1])+"_chain3.pickle",'rb') ## uniform
        object_file_chain3 = pickle.load(file)
        file.close()
        ### concatenate chains
        out_ab=np.concatenate((object_file_chain1.outputs['mig_ab'],object_file_chain2.outputs['mig_ab'],object_file_chain3.outputs['mig_ab']),axis=0)
        out_ba=np.concatenate((object_file_chain1.outputs['mig_ab'],object_file_chain2.outputs['mig_ba'],object_file_chain3.outputs['mig_ba']),axis=0)
        # out_ab=np.concatenate((object_file_chain2.outputs['mig_ab'],object_file_chain3.outputs['mig_ab']),axis=0)
        # out_ba=np.concatenate((object_file_chain2.outputs['mig_ba'],object_file_chain3.outputs['mig_ba']),axis=0)
        ## probability that mig_ab is > mig_ba
        # symvec=np.mean(object_file.outputs['mig_ab']>object_file.outputs['mig_ba'])
        symvec=np.mean(out_ab>out_ba)
        sym_list.append(symvec)
        ## concatenate all posteriors to calculate mean for each GPSC
        vec_outarr=np.concatenate((out_ab,out_ba),axis=0)
        gpsc_meanlist.append(np.mean(vec_outarr))
        np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/PostProcessing/outputs/posteriors/out_ab_"+str(var[1])+d[1]+".txt",out_ab)
        np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/PostProcessing/outputs/posteriors/out_ba_"+str(var[1])+d[1]+".txt",out_ba)
    ### mean migrationfor each GPSC
    mean_mig_byGPSC=np.array(gpsc_meanlist)
    ### directional probability for each GPSc
    gpsc_symmat=np.array(sym_list)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/PostProcessing/outputs/directional_migs/prob_ab_"+str(var[1])+".txt",gpsc_symmat)

    ##### Find mean migration
    gpscs = [10,2,22,26,5,8]
    # gpscs=[8]
for var in range(6):
    print('2) Calculating relative migration to calculate relative migration rates for GPSC'+str(gpscs[var]))
    # print(var)
    deme_list = []
    quants_mig_list_ba= []
    quants_mig_list_ab= []
    for d in enumerate(deme_vec):
        ## chain1
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(gpscs[var])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(gpscs[var])+"_chain1.pickle",'rb') ## uniform
        object_file_chain1 = pickle.load(file)
        file.close()
        ## chain2
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(gpscs[var])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(gpscs[var])+"_chain2.pickle",'rb') ## uniform
        object_file_chain2 = pickle.load(file)
        file.close()
        ## chain3
        file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/2Deme/GPSC"+str(gpscs[var])+"/samples_posteriors81_"+str(d[1])+"3LD_0.5_1000_"+str(gpscs[var])+"_chain3.pickle",'rb') ## uniform
        object_file_chain3 = pickle.load(file)
        file.close()
        ### concatenate chains
        out_ab=np.concatenate((object_file_chain1.outputs['mig_ab'],object_file_chain2.outputs['mig_ab'],object_file_chain3.outputs['mig_ab']),axis=0)
        out_ba=np.concatenate((object_file_chain1.outputs['mig_ab'],object_file_chain2.outputs['mig_ba'],object_file_chain3.outputs['mig_ba']),axis=0)
        # out_ab=np.concatenate((object_file_chain2.outputs['mig_ab'],object_file_chain3.outputs['mig_ab']),axis=0)
        # out_ba=np.concatenate((object_file_chain2.outputs['mig_ba'],object_file_chain3.outputs['mig_ba']),axis=0)
        #### 
        rel_posts_ab=out_ab/mean_mig_byGPSC[var]
        rel_posts_ba=out_ba/mean_mig_byGPSC[var]
        # print(i[1])
        ### mean relative different
        quants_mig_ab=np.quantile(rel_posts_ab,[0.025,0.5,0.975],axis=None)
        quants_mig_ba=np.quantile(rel_posts_ba,[0.025,0.5,0.975],axis=None)
        quants_mig_list_ab.append(quants_mig_ab)
        quants_mig_list_ba.append(quants_mig_ba)
    quants_mig_arr_ab=np.array(quants_mig_list_ab)
    quants_mig_arr_ba=np.array(quants_mig_list_ba)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/PostProcessing/outputs/relative_migs/relative_migquants_0.5_1000_"+str(gpscs[var])+"_ab.txt",quants_mig_arr_ab)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/PostProcessing/outputs/relative_migs/relative_migquants_0.5_1000_"+str(gpscs[var])+"_ba.txt",quants_mig_arr_ba)