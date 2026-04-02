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

def prechange_paramutableMt_det(spre,hpre_list,mapre,tapre,mu,nu,timeseries=False, NoG=10**5):
    wAA,waa,wbb,wAa,wAb,wab = 1-np.array(hpre_list)*spre 
    #print(wAA,waa,wbb,wAa,wAb,wab)
    pA = 1
    pa = 0
    pb = 0 
    
    gen = 0
    
    if timeseries: 
        pAlist = []
        pBlist = []
        palist = []
        pblist = [] 
        
    
    while True: 
        
        pAs = (wAA*pA+wAa*pa+wAb*pb)*pA
        pas = (wAa*pA+waa*pa+wab*pb)*pa
        pbs = (wAb*pA+wab*pa+wbb*pb)*pb
        
        pAmut = (1-mu)*pAs + mu*pas + nu*pbs 
        pamut = (1-mu)*pas - mapre*(1-mu)*pa*(wAa*pA+wab*pb) + tapre*(1-nu)*pbs + mu*pAs 
        pbmut = (1-tapre)*(1-nu)*pbs + mapre*(1-mu)*pa*(wAa*pA+wab*pb)
        
        wavg = pAmut + pamut + pbmut    
        
        pAnext = pAmut/wavg
        panext = pamut/wavg
        pbnext = 1-round(pAnext+panext,5)
        
        if round(pAnext + panext + pbnext,5) > 1:
            print("Negatives imminent", pAnext + panext+ pbnext)     
        gen += 1 
        
        if gen > NoG: 
            break
        pA = pAnext 
        pa = panext 
        pb = pbnext 
        
        if pA < 0 or pa <0 or pb <0 : 
            raise TypeError("Negativesss", [pAs, pAmut, pA], [pas, pamut, pa], [pbs, pbmut, pb])
        
        if timeseries: 
            pAlist += [pA]
            palist += [pa]
            pblist += [pb]
            
    if timeseries: 
        return [pAlist, palist, pblist, pA, pa, pb]
    
    return [pAnext,panext,pbnext]

    
def postchange_paramutableMt_stoch(pA,pa,pb,spost,rpost,hpost_list,mapost,tapost,mu,nu,N_start, timeseries = False) :
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
      
    lamAA,lamaa,lambb,lamAa,lamAb,lamab = 1 - rpost + np.array(hpost_list)*spost
    Ext = 0 
    
    NA = pA*N_start
    Na = pa*N_start
    Nb = pb*N_start
    
    Nsgv = Na + Nb
    
    NAlist = [NA]
    Nalist = [Na]
    Nblist = [Nb]
    
    Nmut = 0
    while True: 
        N_tot = NA + Na + Nb     
        
        if N_tot == 0:
            Ext = 1
            break
        elif N_tot > 2*max(N_start, 10**4): 
            break
        
        pA = NA/float(N_tot)
        pa = Na/float(N_tot)
        pb = Nb/float(N_tot)
        
        NAs = (lamAA*pA+lamAa*pa+lamAb*pb)*NA
        Nas = (lamAa*pA+lamaa*pa+lamab*pb)*Na
        Nbs = (lamAb*pA+lamab*pa+lambb*pb)*Nb
        
        NAmut = (1-mu)*NAs + mu*Nas + nu*Nbs
        Namut = (1-mu)*Nas - mapost*(1-mu)*Na*(lamAa*pA+lamab*pb) + tapost*(1-nu)*Nbs + mu*NAs 
        Nbmut = (1-tapost)*(1-nu)*Nbs + mapost*(1-mu)*Na*(lamAa*pA+lamab*pb)
        
        Nmut += mu*NAs
        
        NA = np.random.poisson(NAmut)
        Na = np.random.poisson(Namut)
        Nb = np.random.poisson(Nbmut)
        
        if timeseries: 
            NAlist += [NA]
            Nalist += [Na]
            Nblist += [Nb]
        
    if timeseries: 
        return NAlist, Nalist, Nblist, Ext
    
    return Ext, Nsgv, Nmut


def postchange_WF_stoch(pA,pa,pb,spost,rpost,hpost_list,mapost,tapost,mu,nu,N_start, timeseries = False) :
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
      
    lamAA,lamaa,lambb,lamAa,lamAb,lamab = 1 - rpost + np.array(hpost_list)*spost
    Ext = 0 
    
    NA = pA*N_start
    Na = pa*N_start
    
    
    Nsgv = Na
    
    NAlist = [NA]
    Nalist = [Na]
    
    Nmut = 0
    while True: 
        N_tot = NA + Na
        
        if N_tot == 0:
            Ext = 1
            break
        elif N_tot > 2*max(N_start, 10**4): 
            break
        
        pA = NA/float(N_tot)
        pa = Na/float(N_tot)
        
        NAs = (lamAA*pA+lamAa*pa)*NA
        Nas = (lamAa*pA+lamaa*pa)*Na
        
        NAmut = (1-mu)*NAs + mu*Nas 
        Namut = (1-mu)*Nas + mu*NAs 
        
        Nmut += mu*NAs
        
        NA = np.random.poisson(NAmut)
        Na = np.random.poisson(Namut)
        
        if timeseries: 
            NAlist += [NA]
            Nalist += [Na]
            
    if timeseries: 
        return NAlist, Nalist, Ext
    
    return Ext, Nsgv, Nmut
#wAA,waa,wbb,wAa,wAb,wab = 1-np.array(hpre_list)*spre 
spre = 0.01 
spost = 1.08 #0.04
rpost = 1 #0.01
mu = 10**(-4)
nu = 10**(-4)
N_start = 10**4


hpre_list = [0, 1, 1, 0.8, 0.7, 1] #[0, 1, 0.48, 0.47, 0.4, 0.27]
hpost_list = hpre_list

NoR = 10**5


x_range = np.arange(0,1.01,0.2) 
y_range = np.arange(0,1.01,0.1) 

a = mp.cpu_count() - 2


ExtProb = [] 
Nsgv = []
Nmut = []

for x in x_range: 
    x = round(x,4)
    mapre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
    tapre = (1-x)*mapre/x if x != 0 else np.random.random()
    pA, pa,pb = prechange_paramutableMt_det(spre,hpre_list,mapre,tapre,mu,nu,NoG = 10**5,timeseries=False)
    output_writer = pd.ExcelWriter('ParamutableMutant_ER_1.xlsx')
    
    Ext_Row = []
    Nsgv_Row = []
    Nmut_Row = []
    for y in y_range: 
        y = round(y,4)
        mapost = np.random.random()*min(y/(1-y),1) if y != 1 else np.random.random()
        tapost = (1-y)*mapost/y if y != 0 else np.random.random()
        
        pool = mp.Pool(a)
        results = pool.starmap(postchange_WF_stoch, [(pA,pa,pb,spost,rpost,hpost_list,mapost,tapost,mu,nu,N_start) for rep in range(NoR)])
        pool.close()
        pool.join()
        Ext = np.mean([replicate[0] for replicate in results])
        Nsgvavg = np.mean([replicate[1] for replicate in results])
        Nmutavg = np.mean([replicate[2] for replicate in results])
        
        Ext_Row += [Ext]
        Nsgv_Row += [Nsgvavg]
        Nmut_Row += [Nmutavg]
        
        print(x,y,Ext,flush=True)
    ExtProb += [Ext_Row]
    Nsgv += [Nsgv_Row]
    Nmut += [Nmut_Row]

ExtDF = pd.DataFrame(ExtProb, columns = y_range, index = x_range)
ExtDF.to_excel(output_writer,sheet_name = 'RescueProb')

NsgvDF = pd.DataFrame(Nsgv, columns = y_range, index = x_range)
NsgvDF.to_excel(output_writer,sheet_name = 'SGV')

NmutDF = pd.DataFrame(Nmut, columns = y_range, index = x_range)
NmutDF.to_excel(output_writer,sheet_name = 'DNM')


output_writer.save()

    
        
        