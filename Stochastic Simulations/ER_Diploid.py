#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 22:56:03 2025

@author: puneeth
"""

import numpy as np
import multiprocessing as mp 
import pandas as pd 
import sys 
def prechange_diploid_det(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,timeseries=False, NoG=10**5):
    wAA,waa,wBB,wbb,wAa,wAB,wAb,waB,wab,wBb = 1-np.array(hpre_list)*spre 
    
    pA = 1
    pB = 0 
    pa = 0
    pb = 0 
    
    gen = 0
    
    if timeseries: 
        pAlist = []
        pBlist = []
        palist = []
        pblist = [] 
        
    
    while True: 
        #print(pA,pa,pB,pb)
        pAs = (wAA*pA+wAa*pa+wAB*pB+wAb*pb)*pA
        pas = (wAa*pA+waa*pa+waB*pB+wab*pb)*pa
        pBs = (wAB*pA+waB*pa+wBB*pB+wBb*pb)*pB
        pbs = (wAb*pA+wab*pa+wBb*pB+wbb*pb)*pb
        
        pAmut = (1-mApre)*(1-mu)*pAs + tApre*(1-nu)*pBs + mu*pas 
        pamut = (1-mapre)*(1-mu)*pas + tapre*(1-nu)*pbs + mu*pAs 
        pBmut = (1-tApre)*(1-nu)*pBs + mApre*(1-mu)*pAs + nu*pbs
        pbmut = (1-tapre)*(1-nu)*pbs + mapre*(1-mu)*pas + nu*pBs
        
        wavg = pAmut + pamut + pBmut + pbmut    
        
        pAnext = pAmut/wavg
        pBnext = pBmut/wavg
        panext = pamut/wavg
        pbnext = 1-(pAnext+pBnext+panext)
        
        if pAnext + pBnext + panext + pbnext > 1:
            print("Negatives imminent", pAnext + pBnext + panext, pbnext)     
        
        gen += 1 
        
        if gen > NoG: 
            break
        pA = pAnext 
        pB = pBnext 
        pa = panext 
        pb = pbnext 
        
        if pA < 0 or pB < 0 or pa <0 or pb <0 : 
            raise TypeError("Negativesss", [pAlist[-1], pAs, pAmut, pA], [pBlist[-1], pBs, pBmut, pB], [palist[-1], pas, pamut, pa], [pblist[-1], pbs, pbmut, pb])
        
        if timeseries: 
            pAlist += [pA]
            pBlist += [pB]
            palist += [pa]
            pblist += [pb]
            
    if timeseries: 
        return [pAlist, palist, pBlist, pblist, pA, pa, pB, pb]
    
    return [pAnext,panext,pBnext,pbnext]

   
def postchange_diploid_stoch(pA,pB,pa,pb,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,N_start,mu,nu,timeseries=False) :
    """
    
    Parameters
    ----------
    spre : float (positive)
    spost : float (positive)
    rpost : float (positive)
    hpre_list : list of length 
    hpost_list : list of length
    mapre : float (between 0 and 1)
    mApre : float (between 0 and 1)
    tapre : float (between 0 and 1)
    tApre : float (between 0 and 1)
    mapost : float (between 0 and 1)
    mApost : float (between 0 and 1)
    tapost : float (between 0 and 1)
    tApost : float (between 0 and 1)
        
    Returns
    -------
    Ext : int (0 or 1)
        Is 0 is the population survives and 1 if the population went extinct
    Xseries : list of lists
        The timeseries data of Xm and Xw. 
    T : int (postive)
        The time to first mutation in generations. 
    """
    
    np.random.seed()
    
    lamAA,lamaa,lamBB,lambb,lamAa,lamAB,lamAb,lamaB,lamab,lamBb = 1 - rpost + np.array(hpost_list)*spost
    
    
    Ext = 0 
    
    NA = pA*N_start
    NB = pB*N_start
    Na = pa*N_start
    Nb = pb*N_start
    
    NAlist = [NA]
    NBlist = [NB]
    Nalist = [Na]
    Nblist = [Nb]
    NumSGV = Na+Nb
    NumMutants = 0
    while True: 
        N_tot = NA + NB + Na + Nb     
        
        if N_tot == 0:
            Ext = 1
            break
        elif N_tot > 2*max(N_start, 10**4): 
            break
        
        pA = NA/float(N_tot)
        pB = NB/float(N_tot)
        pa = Na/float(N_tot)
        pb = Nb/float(N_tot)
        
        NAs = (lamAA*pA+lamAa*pa+lamAB*pB+lamAb*pb)*NA
        Nas = (lamAa*pA+lamaa*pa+lamaB*pB+lamab*pb)*Na
        NBs = (lamAB*pA+lamaB*pa+lamBB*pB+lamBb*pb)*NB
        Nbs = (lamAb*pA+lamab*pa+lamBb*pB+lambb*pb)*Nb
        
        #print(NAs,NBs,Nas,Nbs)
        
        NAmut = (1-mApost)*(1-mu)*NAs + tApost*(1-nu)*NBs + mu*Nas 
        Namut = (1-mapost)*(1-mu)*Nas + tapost*(1-nu)*Nbs + mu*NAs 
        NBmut = (1-tApost)*(1-nu)*NBs + mApost*(1-mu)*NAs + nu*Nbs
        Nbmut = (1-tapost)*(1-nu)*Nbs + mapost*(1-mu)*Nas + nu*NBs
        #print((1-tapost)*(1-nu)*Nbs,  mapost*(1-mu)*Nas, nu*NBs, Nbmut)
        
        NumMutants += mu*NAs + nu*NBs
        
        NA = np.random.poisson(NAmut)
        NB = np.random.poisson(NBmut)
        Na = np.random.poisson(Namut)
        Nb = np.random.poisson(Nbmut)
        
        if timeseries: 
            NAlist += [NA]
            NBlist += [NB]
            Nalist += [Na]
            Nblist += [Nb]
        
    if timeseries: 
        return [NAlist, Nalist, NBlist, Nblist, Ext]
    
    return Ext, NumMutants, NumSGV



#wAA,waa,wBB,wbb,wAa,wAB,wAb,waB,wab,wBb = 1-np.array(hpre_list)*spre 
#wAA = 0 wAB = 0.1 wBB = 0.2 wAa = 0.5 wAb = 0.5 waB = 0.5 wBb = 0.5 wbb = 0.7 wab = 0.8 waa = 1

#hpre_list = [ 0,1,0.2,0.7,0.5,0.1,0.5,0.5,0.8,0.5 ] 

hpre_list = [0,1,0.1,0.3,0.7,0.3,0.3,0.3,0.2,0.6]
hpost_list = hpre_list

spre = 0.01 
spost = 0.04
rpost = 0.01

N_start = 10**4
NoR = 10**5
mu = 10**(-6)
nu = 10**(-6)

y_range = np.arange(0,1.01,0.05) 
#y_range = [1]
a = mp.cpu_count() - 2


ExtProb = [] 


x = sys.argv[1]
x = round(float(x), 5)
mApre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
tApre = (1-x)*mApre/x if x != 0 else np.random.random()
mapre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
tapre = (1-x)*mapre/x if x != 0 else np.random.random()

pA, pB, pa, pb = prechange_diploid_det(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,NoG = 10**5, timeseries = False)


Ext_Row = []


output_writer = pd.ExcelWriter('Diploid_ER_'+str(x)+'.xlsx')

for y in y_range: 
    y = round(y,4)
    mApost = np.random.random()*min(y/(1-y),1) if y != 1 else np.random.random()
    tApost = (1-y)*mApost/y if y != 0 else np.random.random()
    mapost = np.random.random()*min(y/(1-y),1) if y != 1 else np.random.random()
    tapost = (1-y)*mapost/y if y != 0 else np.random.random()
    pool = mp.Pool(a)
    results = pool.starmap(postchange_diploid_stoch, [(pA,pB,pa,pb,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,N_start,mu,nu) for rep in range(NoR)])
    pool.close()
    pool.join()
    Ext = np.mean(results,axis=0)
    Ext_Row += [Ext]
    
    print(Ext,flush=True)

ExtProb = np.transpose(Ext_Row)

ExtDF = pd.DataFrame(ExtProb, columns = y_range, index = [x,'NumMutant','NumSGV'])
ExtDF.to_excel(output_writer,sheet_name = str(x))
output_writer.save()

    
        
        