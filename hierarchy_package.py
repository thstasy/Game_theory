#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Ginny Wei
"""
import random
import statistics
import numpy as np
from hierarchy_plot_package import *
from scipy.stats import entropy

'''
S_i : states (C,D) of indivisuals
'''
def states(pC, n):
    states = []
    num_C = round(n*pC)     #how many C are in the population?
    for i in range(1,n+1):  #make C or D states according to pC
        if i <= num_C:
            states.append('C')
        else:
            states.append('D')
    
    random.shuffle(states)
    return states

'''
L_i(i) - Initial levels (all = 1) of indivisuals
'''
def ini_levels(n):
    Level = []
    for i in range(n):
        Level.append(1)
    return Level

'''
Distrubution of levels  (level,number)
'''
def L_distribution(L_i, S_i, n):
    Level_num = []         #distribution of levels, form: ['level','#total','#C','#D']
    
    sim_L = sorted(list(set(L_i)))
    for k in range(len(sim_L)):
        Level_num.append([sim_L[k], 0, 0, 0])  
        
    sort_L = sorted(L_i) ; k=0
    for i in range(len(L_i)):
        if sim_L[k] == sort_L[i]:
            Level_num[k][1] += 1
        else: 
            k += 1
            Level_num[k][1] += 1
            
        if S_i[i] == 'C':
            Level_num[k][2] += 1
        elif S_i[i] == 'D':
            Level_num[k][3] += 1
        else:
            raise ValueError('State must be "C" or "D"')
        
    return Level_num

'''
H : hierarchicalness, also written as GCR
'''
def hierarchicalness(distr,n):
    x = distr[-1][1]   #x: how many nodes are on the highest level
    C_R_max = (n-x)/(n-1)
    
    howmanylower = 0  #for a given node, how many nodes are at the lower level
    GCR_terms = []; k=0
    num = distr[k][1]
    for i in range(n):
        if i+1 > num:
            howmanylower += distr[k][1]
            k += 1
            num += distr[k][1]
        
        C_R_i = (howmanylower/(n-1))   
        GCR_terms.append(C_R_max - C_R_i)

    H = round(sum(GCR_terms)/(n-1),3)
    
    return H

'''
Max possible G (Gini Coefficient) with given level distribution and NO negative incomes
'''
def nn_maxG(distr, n):
    if len(distr) == 1:
        maxG = 0
    else:
        f_last = distr[-1][1]
        maxG = round(1 - f_last/n, 3)
        
    return maxG

'''
Income pool (each 'C' has a probability H to contribute b into the pool)
'''
def income_pool(distr, b, H):
    pool = 0
    for i in range(len(distr)):
        #pool += distr[i][2]
        for j in range(distr[i][2]):
            r = random.random()
            if r <= H:
                pool += b
    return pool

'''
Divide the pool into TWO levels (NO negative incomes) *Readme +Picture
'''
def twolevel_div(distr, G, n, pool):
    f1 = distr[0][1]
    f2 = distr[1][1]
    h1 = 1-G -f2/n
    
    div_i = []
    for i in range(n):
        if i+1 <= f1:
            div_i.append(round(pool*h1/f1, 4))
        else:
            div_i.append(round(pool*(1-h1)/f2, 4))
    return div_i

'''
Lorenz Curve - area check, No negative incomes (area of 'B')
'''
def G_LCarea(X,Y,pool): #(number of levels >= 3)
    subarea=[]
    for i in range(len(Y)-1):
        subarea.append((Y[i]/pool)*(X[i+1]-X[i]))
        subarea.append((Y[i]/pool)*(X[i+2]-X[i+1]))
    subarea.append((Y[-1]/pool)*(X[-1]-X[-2]))
    return 1-sum(subarea)

'''
Numerical Div into 3 (parabola fit)
'''
def find_para(upper_a, lower_a, xdata, distr, G, pool, n):
    if pool == 0:
        guess_a = -1; guess_G = -1
        div_income = []
        for i in range(n):
            div_income.append(0) 
    else:
        err = 1
        while abs(err) > 0.0001:
            guess_a = (upper_a+lower_a)/2
            div_income = []
            cumu_income = []
            for i in range(len(distr)):
                cumu_income.append(round((xdata[i+1]**guess_a)*pool,5))
                for j in range(distr[i][1]):
                    div_income.append(round((xdata[i+1]**guess_a-xdata[i]**guess_a)*pool/distr[i][1], 5))
            guess_G = G_LCarea(xdata, cumu_income, pool) 
            #print('div_income',div_income)
            #print('cumu_income',cumu_income)
            #print('guess_a:', guess_a,', guess_G:', guess_G, ', G:', G)      
    
            err = G - guess_G
            if err > 0:
                lower_a = guess_a
            elif err < 0:
                upper_a = guess_a

    return guess_a, guess_G, div_income

'''
Divide the pool into 3 or more levels by the approach of parabola
'''
def parabola_div(distr, G, n, pool):  # y = x^a
    #print('pool:',pool)
    xdata = [0]
    for i in range(len(distr)):
        xdata.append(round((distr[i][1]+xdata[i]*10)/n,3))
    #print('xdata', xdata)
    
    upper_a = 100 ; lower_a = 1
    findpara = find_para(upper_a, lower_a, xdata, distr, G, pool, n)
    sol_a = findpara[0]
    sol_G = findpara[1]
    div_income = findpara[2]
    
    return div_income

'''
Given Gini Coefficient -> individual income
'''
def div_pool(distr, G, n, b, c, H):
    pool = income_pool(distr, b, H)
    MaxG = nn_maxG(distr,n)
    #print('* Max G without negative incomes:', MaxG)
    
    divpool_i = []
    if len(distr) == 1: #Only '1' level exists -> evenly divide the pool
        for i in range(n):
            divpool_i.append(round(pool/n, 4))
            mark = 2 # 2: G --> 1only 1 level exists
    
    else:
        if G > MaxG:
            #G_input = MaxG
            G_input = MaxG*1
            mark = 1 # 1: G --> MaxG
        else:
            G_input = G
            mark = 0 # 0: keep the orinigal G
        
        if len(distr) == 2:
            #print('* Only 2 levels')
            divpool_i = twolevel_div(distr, G_input, n, pool)
        else:
            divpool_i = parabola_div(distr, G_input, n, pool)

    return divpool_i, mark

'''
Income after b,"c"(D keep for theirselves) effects
'''
def income_CD(distr, G, n, b, c, H):
    div = div_pool(distr, G, n, b, c, H)
    divpool_i = div[0]
    mark = div[1]
    #print('mark:', mark)
    #print('div_i',divpool_i, ', sum:', round(sum(divpool_i),5))
    #Lorenz_curve(sorted(divpool_i),G)
    
    income_i = []; num = 0
    sumC = 0 ; sumD = 0 ; nC = 0 ; nD = 0
    for i in range(len(distr)):
        for j in range(distr[i][2]): #number of C
            income_i.append(divpool_i[num])
            sumC = sumC + divpool_i[num]
            nC +=1 ; num +=1
        for k in range(distr[i][3]): #number of D
            income_i.append(divpool_i[num]+c)
            sumD = sumD + (divpool_i[num]+c)
            nD +=1 ; num +=1
    #print('num:', num, ' income_i ', income_i)
    
    if nC > 0:
        W_C = round(sumC/nC,4)
    else:
        W_C = 0
    if nD > 0:
        W_D = round(sumD/nD,4)
    else:
        W_D = 0
    
    return income_i, mark, W_C, W_D

'''
Iteration of pC : proportion of C next step
'''
def pC_Next(pC, W_C, W_D): # W(C) / (W(C)+W(D)
    pD = 1-pC
    if W_C == 0 and W_D == 0:
        return 0
    else:
        return round( (pC * W_C)/((pC * W_C)+(pD * W_D)), 4)

'''
h : probability of promotion related to hierarchicalness
g : probability of promotion related to Gini coefficient
'''
def gh(n, G, H, L_i):
    h_i = []
    g_i = []
    pPR = []
    ### kind of "PR" value ###
    sigma = statistics.stdev(L_i)
    ave_L = sum(L_i)/len(L_i)
    range_PR = 3
    
    for i in range(n):
        if sigma == 0:
            psuedo_PR = 0
        else:
            psuedo_PR = (L_i[i]-ave_L)/sigma

        if psuedo_PR >= range_PR: # set the upper bound of psuedo_PR = 3
            psuedo_PR = range_PR
        elif psuedo_PR <= -range_PR: # set the lower bound of psuedo_PR = -3
            psuedo_PR = -range_PR
            
        pPR.append(psuedo_PR)
            
        ### h ; Ax+By+C = 0 ###
        h_x0 = 0.5 - (psuedo_PR/range_PR)*0.5
        #h_value = h_x0 * (1-H)
        h_value = h_x0 * (1-H*0.9) # Give some chances when H=1
        h_i.append(h_value)
        
        ### g ; Ax+By+C = 0 ###
        g_x1 = 0.5 - (psuedo_PR/range_PR)*0.5
        g_value = g_x1 * G
        g_i.append(g_value)
        
    return h_i, g_i   #, pPR

'''
Iteration of L_i (levels)
'''
def level_next(L_i, S_i, distr, n, G, H, alpha):
    L_next = []
    #print('len(distr) :', len(distr))
    
    g_and_h = gh(n, G, H, L_i) # calculate g & h
    h = g_and_h[0] # probability of promotion related to hierarchicalness
    g = g_and_h[1] # probability of promotion related to Gini coefficient
    #print('h:', h); print('g:', g)
    
    for i in range(n):
    
        ### generate the probability of promotion ###
        if len(distr) == 1: # All nodes are in the sane level
            pp = 1/n # probability of promotion
        else:
            pp = alpha*g[i] + (1-alpha)*h[i] # combinition of g and h (0 <= pp <= 1)
    
        ### roll the die ###
        if S_i[i] == 'C':  # C gets promoted
            r = random.random()
            #print('r,pp:',r,pp)
            if r < pp:
                L_next.append(L_i[i]+1) # get promoted (level+1)
            else:
                L_next.append(L_i[i]) # doesn't get promoted
         
        else:  # D doesn't get promoted
            #print('r,pp:','--',0)
            L_next.append(L_i[i]) # doesn't get promoted

    return L_next