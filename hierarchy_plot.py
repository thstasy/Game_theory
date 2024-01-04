#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 17:32:09 2021

@author: weiyujie
"""
from hierarchy_plot_package import *
from hierarchy_package import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
import statistics


with open('alpha_0.5inipC_0.9_201x201_CD.pickle', 'rb') as handle:
    data = pickle.load(handle)

n = data['title'][0]
alpha = data['title'][1]
ini_pC = data['title'][2]
b = data['title'][3]
c_set = data['title'][4]
G_set = data['title'][5]

#print(c_set)
#print(G_set)
item = 'D_final_L_i'

#name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(0.2/b,4))+',G='+str(0.7000000000000001)
#v1 = data[name][item]
#name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(0.3/b,4))+',G='+str(0.6)
#v2 = data[name][item]
#name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(0.4/b,4))+',G='+str(0.5)
#v3 = data[name][item]
name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(0.5/b,4))+',G='+str(0.4)
#v4 = data[name][item]

v1 = data[name][item]

#print(v1, v2, v3, v4)
#print(v1)
p1 = []
for i in v1:
    p1.append(sorted(i))
sp1 = sorted(p1)

uni_L = []
for j in sp1:
    if j not in uni_L:
        uni_L.append(j)
print(uni_L)

finals = []; sumn = 0
for k in uni_L:
    num = sp1.count(k)
    finals.append([num,k])
    sumn += num
sfinals = sorted(finals)
print(sfinals)
print(sumn)






# Plot_aveCD(n, alpha, ini_pC, b, c_set, G_set, data)
# Plot_aveStep(n, alpha, ini_pC, b, c_set, G_set, data)
# Plot_aveLevels(n, alpha, ini_pC, b, c_set, G_set, data)
# Plot_aveOverMaxG(n, alpha, ini_pC, b, c_set, G_set, data) 
# Plot_aveEntropy(n, alpha, ini_pC, b, c_set, G_set, data)



Z = []
for G in G_set:
    z=[]
    for c in c_set:
        name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(c/b,4))+',G='+str(G)
        z.append(data[name]['C_ave_step'])
    Z.append(np.log(z))
#print(Z)


#Plot_C_aveStep(n, alpha, ini_pC, b, c_set, G_set, data)
#Plot_C_aveLevels(n, alpha, ini_pC, b, c_set, G_set, data)
#Plot_C_aveOverMaxG(n, alpha, ini_pC, b, c_set, G_set, data) 
#Plot_C_aveEntropy(n, alpha, ini_pC, b, c_set, G_set, data)

# Plot_D_aveStep(n, alpha, ini_pC, b, c_set, G_set, data)
# Plot_D_aveLevels(n, alpha, ini_pC, b, c_set, G_set, data)
# Plot_D_aveOverMaxG(n, alpha, ini_pC, b, c_set, G_set, data) 
# Plot_D_aveEntropy(n, alpha, ini_pC, b, c_set, G_set, data)


'''
time_shapes = [[1, (17, 21, 22, 22, 23, 23, 23, 24, 24, 25)], [1, (12, 14, 14, 14, 15, 15, 15, 15, 15, 16)], [1, (3, 4, 5, 5, 5, 5, 5, 5, 6, 7)], [1, (11, 11, 11, 11, 11, 12, 13, 14, 14, 14)]]
print('1',time_shapes)
time_shapes.sort()
print('2',time_shapes)
'''
'''
data2 = {}

Z = [] ; S_i = ['C','C','C','C','C','C','C','C','C','C'] ; TimeShapes = [] ; Majority = []
for G in G_set:
    z=[] ; most_shape = [] ; Time_shapes = [] ; majority = []
    for c in c_set:
        name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(c/b,4))+',G='+str(G)
        z.append(data[name]['C_ave_step'])
        
        if data[name]['C_final_L_i'] == []:
            most_shape.append([])
            majority.append([])
        else:
            sortL_set = []
            #print(data[name]['C_final_L_i'])
            for i in range(len(data[name]['C_final_L_i'])):
                #print(data[name]['C_final_L_i'][i])
                a = np.sort(data[name]['C_final_L_i'][i])
                #print('a',a)
                sortL_set.append(tuple(a.tolist()))
            
            time_shapes = []
            #print('sortL_set',sortL_set)
            shapes = set(sortL_set)
            #print('shapes',shapes)
            for j in shapes:
                #print('j',j)
                keer = sortL_set.count(j)
                time_shapes.append([keer, j])
            
            #print('!!!!time_shapes',time_shapes)
        
            time_shapes.sort()
            #print('time_shapes',time_shapes)
            majority.append(time_shapes[-1])
            
            Time_shapes.append(time_shapes)
        
    Z.append(z)
    TimeShapes.append(Time_shapes)
    Majority.append(majority)
print(len(Z))
print('TimeShapes',TimeShapes)
print('Majority',Majority)

data2['C_TimeShapes'] = TimeShapes
data2['C_Majority'] = Majority

###--------------------------------------------------#

Z = [] ; S_i = ['D','D','D','D','D','D','D','D','D','D'] ; TimeShapes = [] ; Majority = []
for G in G_set:
    z=[] ; most_shape = [] ; Time_shapes = [] ; majority = []
    for c in c_set:
        name = 'n='+str(n)+',alpha='+str(alpha)+',ini_pC='+str(ini_pC)+',c/b='+str(round(c/b,4))+',G='+str(G)
        z.append(data[name]['D_ave_step'])
        
        if data[name]['D_final_L_i'] == []:
            most_shape.append([])
            majority.append([])
        else:
            sortL_set = []
            #print(data[name]['C_final_L_i'])
            for i in range(len(data[name]['D_final_L_i'])):
                #print(data[name]['C_final_L_i'][i])
                a = np.sort(data[name]['D_final_L_i'][i])
                #print('a',a)
                sortL_set.append(tuple(a.tolist()))
            
            time_shapes = []
            #print('sortL_set',sortL_set)
            shapes = set(sortL_set)
            #print('shapes',shapes)
            for j in shapes:
                #print('j',j)
                keer = sortL_set.count(j)
                time_shapes.append([keer, j])
            
            #print('!!!!time_shapes',time_shapes)
        
            time_shapes.sort()
            #print('time_shapes',time_shapes)
            majority.append(time_shapes[-1])
            
            Time_shapes.append(time_shapes)
        
    Z.append(z)
    TimeShapes.append(Time_shapes)
    Majority.append(majority)
print(len(Z))
print('TimeShapes',TimeShapes)
print('Majority',Majority)

data2['D_TimeShapes'] = TimeShapes
data2['D_Majority'] = Majority

with open('CD_shapes.pickle', 'wb') as file:
        pickle.dump(data2, file)
'''