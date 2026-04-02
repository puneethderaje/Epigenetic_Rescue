#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 16:29:29 2025

@author: puneeth
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 22:56:03 2025

@author: puneeth
"""

import numpy as np
import multiprocessing as mp 
import pandas as pd 

import time 

def Diploid(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,N_start,NoG = 10**5, timeseries = False) :
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
        
        pAexpected = pAmut/wavg
        pBexpected = pBmut/wavg
        paexpected = pamut/wavg
        pbexpected = 1-(pAexpected+pBexpected+paexpected)
        
        pAnext, pBnext, panext, pbnext = np.around( np.random.multinomial(N_start,[pAexpected,pBexpected,paexpected,pbexpected])/float(N_start),10)
        
        if np.around(pAnext + pBnext + panext + pbnext,10) > 1:
            print("Negatives imminent", pAnext, pBnext, panext, pbnext, pAnext + pBnext + panext + pbnext)     
        
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
            pAlist += [pA*N_start]
            pBlist += [pB*N_start]
            palist += [pa*N_start]
            pblist += [pb*N_start]
    
    return [pAnext, pBnext, panext, pbnext]



#wAA,waa,wBB,wbb,wAa,wAB,wAb,waB,wab,wBb = 1-np.array(hpre_list)*spre 

hpre_list = [ 0,1,0.2,0.7,0.3,0.1,0.5,0.4,0.8,0.6 ] 

spre = 0.01 

N_start = 10**4
NoR = 10**4
mu = 10**(-6)
nu = 10**(-6)

x_range = np.arange(0,1,0.05) 

a = mp.cpu_count() - 2

output_writer = pd.ExcelWriter('Prechange_Stochastic_Diploid.xlsx')

for x in x_range: 
    x = round(x,4)
    print(x)
    mApre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
    tApre = (1-x)*mApre/x if x != 0 else np.random.random()
    mapre = np.random.random()*min(x/(1-x),1) if x != 1 else np.random.random()
    tapre = (1-x)*mapre/x if x != 0 else np.random.random()
    
    Ext_Row = []
    st_time = time.time()
    pool = mp.Pool(a)
    results = pool.starmap(Diploid, [(spre,hpre_list,mapre,mApre,tapre,tApre,mu,nu,N_start) for rep in range(NoR)])
    pool.close()
    pool.join()
    end_time = time.time()
    print(end_time-st_time)

    
    EquiDF = pd.DataFrame(results, columns = ['A','B','a','b'], index = range(NoR))
    EquiDF.to_excel(output_writer,sheet_name = str(x))
output_writer.save()
    
        
        