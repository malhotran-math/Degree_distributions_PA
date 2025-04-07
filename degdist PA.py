#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:41:49 2022

@author: malhotran
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time 


#parameters for linear g(x) = ax+b and fitness f(x)=0
n = 100
itern = 1

#function to simulate Preferential attachment tree and return degrees 
def degreedist1(n): 
    nodes = [0]  #initial vertex set
    indegree = [0] #initial indegree
    P=[1] 
    p=[1] #probability of attachment 
    A = np.array(np.zeros((n,n)))   #adjacency
    counter = 1
    for k in range(1,n):
        counter += 1
        C = np.random.choice(nodes,1,p=p) #choosing a vertex to attach to with law p
        c = int(C)
        A[k,c]+=1 #attach new vertex to choses old vertex
        nodes.append(k) #add new vertex to tree
        indegree.append(0)
        P = np.zeros(counter)
        p = P 
        for i in range(counter): 
            indegree[i]=A[:,i].sum() 
        for j in range(counter): 
            P[j] = indegree[j] + j**2 #law of attachment
        SUM = (sum(P[i] for i in range(counter))) 
        for l in range(counter): 
            p[l] = P[l]/SUM     #scaling to get probability 
    IN = indegree           
    deglist = np.arange(n)
    count = []
    for i in deglist:
        c = IN.count(i) 
        count.append(c)
    return deglist, count


def analyticalsup():
    AN = []
    for k in range(n):
        den = 2**(k+1)
        frac = 1/den 
        AN.append(frac) 
    return AN

deglist, count = degreedist1(n)
avgcount = np.array(count) 
for i in range(itern-1):
    d, c = degreedist1(n)
    C = np.array(c)
    avgcount = avgcount + c
avgcount = np.multiply(avgcount,1/(itern*n)) 
#averaging the frequency of degrees and scaling it by number of vertices to obtain probabilities

AVC1 = [] 
#obtaining CCDF 
for K in range(len(avgcount)):
    temp = float(sum(avgcount[i] for i in range(K,len(avgcount))))
    AVC1.append(temp)
    


    

    
AN1 = analyticalsup() #analytical expression
AVC2 = [] 
for K in range(len(avgcount)):
    temp = float(sum(AN1[i] for i in range(K,len(avgcount))))
    AVC2.append(temp)
 
plt.scatter(np.log(deglist + 1),np.log(AVC1),color='green',label='Simulated for f(x)=x^3')  
plt.plot(np.log(deglist[:50] + 1),np.log(AVC2[:50]),label='Analytical') 

X= np.log(AVC1)


x = np.log(deglist + 1)
y = np.log(AVC1)
z = np.polyfit(x,y,3) 
p = np.poly1d(z) 


plt.xlabel('log(k+1)')
plt.ylabel('log(F(k)')
plt.legend()
plt.show()