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

def Diploid_PreChange(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,NoG = 10**5, timeseries = False):
    wAA,waa,wBB,wAa,wAB,waB = 1-np.array(hpre_list)*spre
    
    pA = 0
    pa = 1
    pB = 0 
    
    gen = 0
    
    if timeseries: 
        pAlist = []
        palist = []
        pBlist = [] 
        
    
    while True: 
        #print(pA,pa,pB,pb)
        pas = (waa*pa+wAa*pA+waB*pB)*pa
        pAs = (wAa*pa+wAA*pA+wAB*pB)*pA
        pBs = (waB*pa+wAB*pA+wBB*pB)*pB
        
        pamut = (1-mu)*pas + mu*pAs + nu*pBs 
        pAmut = (1-mu)*pAs - mApre*(1-mu)*pA*(wAa*pa+wAB*pB) + tApre*(1-nu)*pBs + mu*pBs 
        pBmut = (1-tApre)*(1-nu)*pBs + mApre*(1-mu)*pA*(wAa*pa+wAB*pB)
        
        wavg = pamut + pAmut + pBmut    
        
        pAnext = pAmut/wavg
        panext = pamut/wavg
        pBnext = 1-(pAnext+panext)
        
        if pAnext + panext + pBnext > 1:
            print("Negatives imminent")     
        
        gen += 1 
        
        if gen > NoG: 
            break
        pA = pAnext 
        pa = panext 
        pB = pBnext 
        
        if pA < 0 or pa <0 or pB <0 : 
            raise TypeError("Negativesss", [pAlist[-1], pAs, pAmut, pA], [palist[-1], pas, pamut, pa], [pBlist[-1], pBs, pBmut, pB])
        
        if timeseries: 
            pAlist += [pA*N_start]
            palist += [pa*N_start]
            pBlist += [pB*N_start]
            
    return pA,pa,pB

def Diploid(pA,pa,pB,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,mu,nu,N_start,timeseries = False) :
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
    
    
    
    
    
    lamAA,lamaa,lamBB,lamAa,lamAB,lamaB = 1 - rpost + np.array(hpost_list)*spost
    
    Ext = 0 
    
    NA = round(pA*N_start)
    Na = round(pa*N_start)
    NB = round(pB*N_start)
    
    NAlist = [NA]
    Nalist = [Na]
    NBlist = [NB]
    
    while True: 
        N_tot = NA + Na + NB
        
        if N_tot == 0:
            Ext = 1
            break
        elif N_tot > 2*max(N_start, 10**4): 
            break
        
        pA = NA/float(N_tot)
        pa = Na/float(N_tot)
        pB = NB/float(N_tot)
        
        Nas = (lamaa*pa+lamAa*pA+lamaB*pB)*Na
        NAs = (lamAa*pa+lamAA*pA+lamAB*pB)*NA
        NBs = (lamaB*pa+lamAB*pA+lamBB*pB)*NB
        
        Namut = (1-mu)*Nas + mu*NAs + nu*NBs 
        NAmut = (1-mu)*NAs - mApost*(1-mu)*NA*(lamAa*pa+lamAB*pB) + tApost*(1-nu)*NBs + mu*NBs 
        NBmut = (1-tApost)*(1-nu)*NBs + mApost*(1-mu)*NA*(lamAa*pa+lamAB*pB)
        
        
        NA = np.random.poisson(NAmut)
        Na = np.random.poisson(Namut)
        NB = np.random.poisson(NBmut)
        
        if timeseries: 
            NAlist += [NA]
            Nalist += [Na]
            NBlist += [NB]
        
    if timeseries: 
        return [NAlist, Nalist, NBlist, Ext]
    
    return Ext



#wAA,waa,wbb,wAa,wAb,wab = 1-np.array(hpre_list)*spre 

h = 0.5
hpre_list = [1,0,1,h,1,h] 
hpost_list = hpre_list

spre = 0.01 
spost = 0.04
rpost = 0.01

N_start = 10**4
NoR = 10**5
mu = 10**(-6)
nu = 10**(-6)

x_range = np.arange(0,1,0.05) 
y_range = np.arange(0,1,0.05) 

a = mp.cpu_count() - 2

ExtProb = [] 

x = sys.argv[1]
x = round(float(x),4)
mApre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
tApre = (1-x)*mApre/x if x != 0 else np.random.random()
mapre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
tapre = (1-x)*mapre/x if x != 0 else np.random.random()

pA,pa,pB=Diploid_PreChange(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,NoG = 10**5, timeseries = False)

output_writer = pd.ExcelWriter('ParamutableWildtype_ER_'+str(x)+'.xlsx')

Ext_Row = []
for y in y_range: 
    y = round(y,4)
    mApost = np.random.random()*min(y/(1-y),1) if y != 1 else np.random.random()
    tApost = (1-y)*mApost/y if y != 0 else np.random.random()
    mapost = np.random.random()*min(y/(1-y),1) if y != 1 else np.random.random()
    tapost = (1-y)*mapost/x if x != 0 else np.random.random()
    
    pool = mp.Pool(a)
    results = pool.starmap(Diploid, [(pA,pa,pB,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,mu,nu,N_start) for rep in range(NoR)])
    pool.close()
    pool.join()
    Ext = np.mean(results)
    Ext_Row += [Ext]
    
    print(x,y,Ext,flush=True)
ExtProb += [Ext_Row]

ExtDF = pd.DataFrame(ExtProb, columns = y_range, index = [x])
ExtDF.to_excel(output_writer,sheet_name = 'ExtProb')
output_writer.save()

    
        
        