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

def Haploid_PreChange(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,NoG = 10**5, timeseries = False): 
    
    wA,wa,wB,wb = 1-np.array(hpre_list)*spre 
    
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
        pAs = wA*pA
        pas = wa*pa
        pBs = wB*pB
        pbs = wb*pb
        
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
    
    return pA, pB, pa, pb
    
def Haploid(pA,pB,pa,pb,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,N_start,mu,nu,timeseries=False) :
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
    
    lamA,lama,lamB,lamb = 1 - rpost + np.array(hpost_list)*spost
    #print(lamA, lama, lamB, lamb)
    
    
    Ext = 0 
    
    NA = round(pA*N_start)
    NB = round(pB*N_start)
    Na = round(pa*N_start)
    Nb = round(pb*N_start)
    
    NAlist = [NA]  
    NBlist = [NB]
    Nalist = [Na]
    Nblist = [Nb]
    
    
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
        
        NAs = lamA*NA
        Nas = lama*Na
        NBs = lamB*NB
        Nbs = lamb*Nb
        
        #print(NAs,NBs,Nas,Nbs)
        
        NAmut = (1-mApost)*(1-mu)*NAs + tApost*(1-nu)*NBs + mu*Nas 
        Namut = (1-mapost)*(1-mu)*Nas + tapost*(1-nu)*Nbs + mu*NAs 
        NBmut = (1-tApost)*(1-nu)*NBs + mApost*(1-mu)*NAs + nu*Nbs
        Nbmut = (1-tapost)*(1-nu)*Nbs + mapost*(1-mu)*Nas + nu*NBs
        #print((1-tapost)*(1-nu)*Nbs,  mapost*(1-mu)*Nas, nu*NBs, Nbmut)
        
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
    
    return Ext


ind = sys.argv[1]
ind = round(float(ind))
xpre = sys.argv[2]
xtot = round(float(xpre),5)
#wAA,waa,wBB,wbb,wAa,wAB,wAb,waB,wab,wBb = 1-np.array(hpre_list)*spre 
#wAA = 0 wAB = 0.1 wBB = 0.2 wAa = 0.3 wAb = 0.4 waB = 0.5 wBb = 0.6 wbb = 0.7 wab = 0.8 waa = 1
rpost=0.01
spre=0.01
spost=0.02
hbpre=0.8
hBpre=0.1
#hbpost=0.625
#hBpost=0.3
mu=0
nu=0
N_start=10**4
NoR = 10**5


#xpre_range = np.arange(0,1.1,0.2)
xpost_range = np.arange(0,1.05,0.05)

hbpost, hBpost = 0,1
hpost_list = [0,1,hBpost,1-hbpost]

a = mp.cpu_count() - 2

ExtProb = [] 

output_writer = pd.ExcelWriter('Haploid_ER_PE'+str(ind)+'_'+str(xtot)+'.xlsx')

Ext_row = []

mApre = 0.8
tApre = 0.8
mapre = 0
tapre = 1

hpre_list = [0,1,hBpre,1-hbpre]
pA,pB,pa,pb = Haploid_PreChange(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,NoG = 10**5, timeseries = False)


for xpost in xpost_range: 
    y = round(xpost,5)
    
    mApost = xtot*xpost
    tApost = (1-xtot)*xpost
    mapost = 0 # np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
    tapost = 1 #(1-x)*mapost/x if x != 0 else np.random.random()
    
    #Presval = Pres(r,spre,spost,hbpre,hBpre,hbpost,hBpost,mapre,mApre,mapost,mApost,mu,nu,N)
    pool = mp.Pool(a)
    results = pool.starmap(Haploid, [(pA,pB,pa,pb,spost,rpost,hpost_list,mapost,mApost,tapost,tApost,N_start,mu,nu,) for rep in range(NoR)])
    pool.close()
    pool.join()
    Ext = np.mean(results)
    print(Ext)
    Ext_row += [Ext]    
ExtProb += [Ext_row]

ExtDF = pd.DataFrame(ExtProb, columns = xpost_range, index = [xtot])
ExtDF.to_excel(output_writer,sheet_name = str(xtot))
output_writer.close()

    
        
        
