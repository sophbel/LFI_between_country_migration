import pandas as pd
import numpy as np
import pickle
import elfi
import matplotlib.pyplot as plt
gpscs = [10,2,22,26,5,8]
gpsc_meanlist=[]
for var in enumerate(gpscs):
    print(var[1])
    deme_list = []
    sym_list = []
    ## chain
    file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/samples_posteriors81_5_symmetric_9stat_GPSC"+str(var[1])+"_LD.pickle",'rb') ## uniform
    object_file = pickle.load(file)
    file.close()
    ###extract chains
    mig_ab=object_file.outputs['mig_ab']
    mig_ac=object_file.outputs['mig_ac']
    mig_ad=object_file.outputs['mig_ad']
    mig_bc=object_file.outputs['mig_bc']
    mig_bd=object_file.outputs['mig_bd']
    mig_cd=object_file.outputs['mig_cd']
    ## concatenate all posteriors to calculate mean for each GPSC
    vec_outarr=np.concatenate((mig_ab,mig_ac,mig_ad,mig_bc,mig_bd,mig_cd),axis=0)
    gpsc_meanlist.append(np.mean(vec_outarr))
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_ab_"+str(var[1])+".txt",mig_ab)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_ac_"+str(var[1])+".txt",mig_ac)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_ad_"+str(var[1])+".txt",mig_ad)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_bc_"+str(var[1])+".txt",mig_bc)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_bd_"+str(var[1])+".txt",mig_bd)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(var[1])+"/mig_cd_"+str(var[1])+".txt",mig_cd)
    ### mean migrationfor each GPSC
    mean_mig_byGPSC=np.array(gpsc_meanlist)


    gpscs = [10,2,22,26,5,8]
for var in range(6):
    print(var)
    deme_list = []
    quants_mig_list_ab= []
    quants_mig_list_ac= []
    quants_mig_list_ad= []
    quants_mig_list_bc= []
    quants_mig_list_bd= []
    quants_mig_list_cd= []

     ## chain
    file = open("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/samples_posteriors81_5_symmetric_9stat_GPSC"+str(gpscs[var])+"_LD.pickle",'rb') ## uniform
    object_file = pickle.load(file)
    file.close()
    ###extract chains
    mig_ab=object_file.outputs['mig_ab']
    mig_ac=object_file.outputs['mig_ac']
    mig_ad=object_file.outputs['mig_ad']
    mig_bc=object_file.outputs['mig_bc']
    mig_bd=object_file.outputs['mig_bd']
    mig_cd=object_file.outputs['mig_cd']
    #### calculate relative migraiton
    rel_posts_ab=mig_ab/mean_mig_byGPSC[var]
    rel_posts_ac=mig_ac/mean_mig_byGPSC[var]
    rel_posts_ad=mig_ad/mean_mig_byGPSC[var]
    rel_posts_bc=mig_bc/mean_mig_byGPSC[var]
    rel_posts_bd=mig_bd/mean_mig_byGPSC[var]
    rel_posts_cd=mig_cd/mean_mig_byGPSC[var]
    # print(i[1])
    ### mean relative different
    quants_mig_ab=np.quantile(rel_posts_ab,[0.025,0.5,0.975],axis=None)
    quants_mig_ac=np.quantile(rel_posts_ac,[0.025,0.5,0.975],axis=None)
    quants_mig_ad=np.quantile(rel_posts_ad,[0.025,0.5,0.975],axis=None)
    quants_mig_bc=np.quantile(rel_posts_bc,[0.025,0.5,0.975],axis=None)
    quants_mig_bd=np.quantile(rel_posts_bd,[0.025,0.5,0.975],axis=None)
    quants_mig_cd=np.quantile(rel_posts_cd,[0.025,0.5,0.975],axis=None)

    quants_mig_list_ab.append(quants_mig_ab)
    quants_mig_list_ac.append(quants_mig_ac)
    quants_mig_list_ad.append(quants_mig_ad)
    quants_mig_list_bc.append(quants_mig_bc)
    quants_mig_list_bd.append(quants_mig_bd)
    quants_mig_list_cd.append(quants_mig_cd)

    quants_mig_arr_ab=np.array(quants_mig_list_ab)
    quants_mig_arr_ac=np.array(quants_mig_list_ac)
    quants_mig_arr_ad=np.array(quants_mig_list_ad)
    quants_mig_arr_bc=np.array(quants_mig_list_bc)
    quants_mig_arr_bd=np.array(quants_mig_list_bd)
    quants_mig_arr_cd=np.array(quants_mig_list_cd)

    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_ab.txt",quants_mig_arr_ab)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_ac.txt",quants_mig_arr_ac)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_ad.txt",quants_mig_arr_ad)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_bc.txt",quants_mig_arr_bc)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_bd.txt",quants_mig_arr_bd)
    np.savetxt("/Users/sb62/Documents/Migration/DemographyModel/DemographyAnalysis/RunModel/outputs/4Deme/GPSC"+str(gpscs[var])+"/relative_migquants_0.5_1000"+str(gpscs[var])+"_cd.txt",quants_mig_arr_cd)