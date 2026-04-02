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
def prechange_diploid_det(spre,hpre_list,mu,timeseries=False, NoG=10**5):
    wAA,waa,wAa = 1-np.array(hpre_list)*spre 
    
    pA = 1
    pa = 0
    
    gen = 0
    
    if timeseries: 
        pAlist = []
        palist = []
        
    
    while True: 
        pAs = (wAA*pA+wAa*pa)*pA
        pas = (wAa*pA+waa*pa)*pa
        
        pAmut = (1-mu)*pAs + mu*pas 
        pamut = (1-mu)*pas + mu*pAs 
        
        wavg = pAmut + pamut
        
        pAnext = pAmut/wavg
        panext = pamut/wavg
        
        if pAnext + panext > 1:
            print("Negatives imminent", pAnext, panext)     
        
        gen += 1 
        
        if gen > NoG: 
            break
        pA = pAnext 
        pa = panext 
        
        if pA < 0 or pa <0 : 
            raise TypeError("Negativesss", [pAlist[-1], pAs, pAmut, pA], [palist[-1], pas, pamut, pa])
        
        if timeseries: 
            pAlist += [pA]
            palist += [pa]
            
    if timeseries: 
        return [pAlist, palist, pA, pa]
    
    return [pAnext,panext]

   
def postchange_diploid_stoch(pA,pa,spost,rpost,hpost_list,N_start,mu,timeseries=False) :
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
    
    lamAA,lamaa,lamAa = 1 - rpost + np.array(hpost_list)*spost
    
    
    Ext = 0 
    
    NA = pA*N_start
    Na = pa*N_start
    
    NAlist = [NA]
    Nalist = [Na]
    
    NumSGV = Na
    NumMutants = 0
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
        
        #print(NAs,NBs,Nas,Nbs)
        
        NAmut = (1-mu)*NAs + mu*Nas 
        Namut = (1-mu)*Nas + mu*NAs 
        #print((1-tapost)*(1-nu)*Nbs,  mapost*(1-mu)*Nas, nu*NBs, Nbmut)
        
        NumMutants += mu*NAs
        
        NA = np.random.poisson(NAmut)
        Na = np.random.poisson(Namut)
        
        if timeseries: 
            NAlist += [NA]
            Nalist += [Na]
        
    if timeseries: 
        return [NAlist, Nalist, Ext]
    
    return Ext, NumMutants, NumSGV

def Pres_theoretical(spre,spost,rpost,N_start,hpre_list,hpost_list,mu):
    wAA,waa,wAa = 1- np.array(hpre_list)*spre
    lamAA,lamaa,lamAa = 1-rpost + np.array(hpost_list)*spost
    saaeff = lamaa - 1
    sAaeff = lamAa - 1
    pasgv = mu/(wAa-waa) 
    pAasgv = 2*pasgv
    
    return 1 - np.exp(-2*saaeff*( 0.5*pAasgv*N_start*(0.25+2*mu)/(0.25+2*mu-sAaeff) + pasgv**2 ))

hpre_list = [0,1,0.8]
hpost_list = hpre_list

spre = 0.01 
spost = 1.02
rpost = 1

N_start = 10**4
NoR = 10**5
mu = 10**(-6)
nu = 10**(-6)

y_range = np.arange(0,1.01,0.05) 
a = mp.cpu_count() - 2


pA, pa = prechange_diploid_det(spre,hpre_list,mu,timeseries=False, NoG=10**5)

pool = mp.Pool(a)
results = pool.starmap(postchange_diploid_stoch, [(pA,pa, spost,rpost,hpost_list,N_start,mu) for rep in range(NoR)])
pool.close()
pool.join()
Ext = np.mean(results,axis=0)

print(Ext,flush=True)

    
        
        