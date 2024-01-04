#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 1 17:55:46 2021

@author: Ginny Wei
"""
from hierarchy_package import *
from hierarchy_plot_package import *
from scipy.stats import entropy
import numpy as np
import matplotlib.pyplot as plt
import statistics
import pickle
'''
Global constants
'''
n        = 5    #group size
b        = 1     #contributing a benefit of magnitude b to a common pool
alpha    = 0.5   #(g/h) ratio
ini_pC   = 0.9   #initial proportion of cooperators/population
maxstep  = 500   #times of iteration
test_num = 200   #repeating times
data     = {}    #for storing indivisuals' levels and status(C,D)
'''
Variables
'''
#c_set    = [0.31, 0.3162500000, 0.31625000001, 0.32]
#G_set    = [0.75]

##c_set    = [0.47, 0.4741499999, 0.474150, 0.48] ; G_set = [0.95]
##c_set    = [0.4741499999, 0.474150] ; G_set = [0.95]

##c_set    = [0.31, 0.3162500000, 0.31625000001, 0.32] ; G_set = [0.75]
##c_set    = [0.3162500000, 0.31625000001] ; G_set = [0.75]

c_set    = list(np.linspace(0,1,201))   #keeping a benefit of magnitude c for defector
G_set    = list(np.linspace(0,1,201))   #value of Gini Coefficient (0 ~ 0.6?)

#c_set    = list(np.linspace(0,1,401))   #keeping a benefit of magnitude c for defector
#G_set    = list(np.linspace(0,1,401))   #value of Gini Coefficient (0 ~ 0.6?)

def main(c, G):
    name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(c/b,11))+',G='+str(G)
    data[name]={}
    print(name)

    data[name]['final_L_i'] = []
    data[name]['CD_ratio'] = []
    data[name]['num_levels'] = []
    data[name]['times_overMaxG'] = []
    data[name]['steps'] = []
    data[name]['entropy'] = []
    
    data[name]['C_final_L_i'] = []
    data[name]['C_num_levels'] = []
    data[name]['C_times_overMaxG'] = []
    data[name]['C_steps'] = []
    data[name]['C_entropy'] = []
    data[name]['C_end'] = []
    
    data[name]['D_final_L_i'] = []
    data[name]['D_num_levels'] = []
    data[name]['D_times_overMaxG'] = []
    data[name]['D_steps'] = []
    data[name]['D_entropy'] = []
    data[name]['D_end'] = []
    
    for i in range(test_num):
        step = 0 ; overMaxG = 0 ; C_end = 0 ; D_end = 0
        pC = ini_pC
        S_i = states(ini_pC, n)  
        L_i = ini_levels(n)      
        distr = L_distribution(L_i, S_i, n) #distribution of levels, format: ['level','#total','#C','#D']
        H = hierarchicalness(distr,n)   
        #levelstates_hist(distr, step)
        L_i = level_next(L_i, S_i, distr, n, G, H, alpha)
        distr = L_distribution(L_i, S_i, n)
        H = hierarchicalness(distr,n)  
        income = income_CD(distr, G, n, b, c, H)
        income_i = income[0]
        mark = income[1]
        if mark == 1:
            overMaxG += 1
        W_C = income[2] ; W_D = income[3]
    
        while len(set(S_i)) > 1 and step <= maxstep-1:
            step += 1
            pC = pC_Next(pC, W_C, W_D)
            #print('pC',pC,W_C,W_D)
            L_i = level_next(L_i, S_i, distr, n, G, H, alpha)
            S_i = states(pC, n)
            distr = L_distribution(L_i, S_i, n) 
            H = hierarchicalness(distr,n)  
            #levelstates_hist(distr, step)
            income = income_CD(distr, G, n, b, c, H)
            income_i = income[0]
            mark = income[1]
            if mark == 1:
                overMaxG += 1
            W_C = income[2] ; W_D = income[3]
    
        data[name]['final_L_i'].append(L_i)
        data[name]['CD_ratio'].append(round(S_i.count('C')/len(S_i), 4))
        data[name]['num_levels'].append(len(list(set(L_i))))
        data[name]['steps'].append(step)
        data[name]['times_overMaxG'].append(overMaxG/step)
        data[name]['entropy'].append(entropy(L_i, ini_levels(n)))
        
        if S_i.count('D') == 0:
            C_end += 1
            data[name]['C_final_L_i'].append(L_i)
            data[name]['C_num_levels'].append(len(list(set(L_i))))
            data[name]['C_times_overMaxG'].append(overMaxG/step)
            data[name]['C_steps'].append(step)
            data[name]['C_entropy'].append(entropy(L_i, ini_levels(n)))
        
        if S_i.count('C') == 0:
            D_end += 1
            data[name]['D_final_L_i'].append(L_i)
            data[name]['D_num_levels'].append(len(list(set(L_i))))
            data[name]['D_times_overMaxG'].append(overMaxG/step)
            data[name]['D_steps'].append(step)
            data[name]['D_entropy'].append(entropy(L_i, ini_levels(n)))
    
    data[name]['ave_overMaxG'] = np.mean(data[name]['times_overMaxG'])
    data[name]['ave_CD_ratio'] = round(np.mean(data[name]['CD_ratio']), 3)
    data[name]['ave_num_levels'] = np.mean(data[name]['num_levels'])
    data[name]['ave_step'] = np.mean(data[name]['steps'])
    data[name]['ave_entropy'] = np.mean(data[name]['entropy'])
    
    data[name]['C_ave_overMaxG'] = np.mean(data[name]['C_times_overMaxG'])
    data[name]['C_ave_num_levels'] = np.mean(data[name]['C_num_levels'])
    data[name]['C_ave_step'] = np.mean(data[name]['C_steps'])
    data[name]['C_ave_entropy'] = np.mean(data[name]['C_entropy'])
    data[name]['C_end'].append(C_end)
    
    data[name]['D_ave_overMaxG'] = np.mean(data[name]['D_times_overMaxG'])
    data[name]['D_ave_num_levels'] = np.mean(data[name]['D_num_levels'])
    data[name]['D_ave_step'] = np.mean(data[name]['D_steps'])
    data[name]['D_ave_entropy'] = np.mean(data[name]['D_entropy'])
    data[name]['D_end'].append(D_end)

    print('ave_step:', data[name]['ave_step'], 'ave_CD:', data[name]['ave_CD_ratio'], 'ave_num_levels:', data[name]['ave_num_levels'], 'OverMaxG:', data[name]['ave_overMaxG'])
    print(' ')
    
if __name__ == '__main__':
    data['title'] = [n, alpha, ini_pC, b, c_set, G_set]
    for c in c_set:
        for G in G_set:
            main(c, G)
    with open('alpha_'+str(alpha)+'inipC_'+str(ini_pC)+'_201x201_CD'+'.pickle', 'wb') as file:
        pickle.dump(data, file)