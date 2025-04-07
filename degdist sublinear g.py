#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:22:42 2022

@author: malhotran
"""


#SUBLINEAR CASE
#PART 1: Execute the program and save the data to a text file

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse



#parameters for linear g(x) = ax+b and fitness f(x)=0
a=2
b=3
itern = 1 #number of iterations 
n =3000 #number of vertices
lam = 5

#function to simulate Preferential attachment tree and return degrees 
#the same code can simulate a sublinear g, by putting an np.log() for the indegree function
def degreedist1(n): 
    nodes = [0] 
    indegree = [0]
    P=[1]
    p=[1]
    A = np.array(np.zeros((n,n)))  
    counter = 1
    for k in range(1,n):
        counter += 1
        C = np.random.choice(nodes,1,p=p)
        c = int(C)
        A[k,c]+=1
        nodes.append(k)
        indegree.append(0)
        P = np.zeros(counter)
        p = P 
        for i in range(counter): 
            indegree[i]=A[:,i].sum() 
        for j in range(counter): 
            P[j] = 2*indegree[j]+3 # + np.random.poisson(10) (for linear fitness)
        SUM = (sum(P[i] for i in range(counter))) 
        for l in range(counter): 
            p[l] = P[l]/SUM     
    IN = indegree          
    deglist = np.arange(n)
    count = []
    for i in deglist:
        c = IN.count(i) 
        count.append(c)
    return deglist, count
    
'''def degreedist2(n): 
    nodes = [0] 
    indegree = [0]
    P=[1]
    p=[1]
    A = np.array(np.zeros((n,n)))  
    counter = 1
    for k in range(1,n):
        counter += 1
        C = np.random.choice(nodes,1,p=p)
        c = int(C)
        A[k,c]+=1
        nodes.append(k)
        indegree.append(0)
        P = np.zeros(counter)
        p = P 
        for i in range(counter): 
            indegree[i]=A[:,i].sum() 
        for j in range(counter): 
            P[j] = a*indegree[j] + b + j**11
        SUM = (sum(P[i] for i in range(counter))) 
        for l in range(counter): 
            p[l] = P[l]/SUM     
    IN = indegree          
    deglist = np.arange(n)
    count = []
    for i in deglist:
        c = IN.count(i) 
        count.append(c)
    return deglist, count
    
'''
#analytical expression for "a+b" model (sublinear case) where f(u) = 0 and g(k) = log(ak+b)
def analyticalsub():
    AN = []
    for k in range(n): 
        prod = 1
        for i in range(k):
            frac = (np.log(2*i + 3)+5)/(10+np.log(2*i + 3)) 
            prod = prod*frac
        prod = prod*5/(10+np.log(2*k + 3))  
        AN.append(prod)  
    return AN

#analytical function for g(x) = ax+b and linear/sublinear fitness. Set CL =0 for sublinear fitness
def analyticalab():
    AN = []
    CL = 0
    for k in range(n): 
        prod = 1
        for i in range(k):
            frac = (a*i + b + CL)/(a*(i+1)+2*b+2*CL)  
            prod = prod*frac
        prod = prod*(a+b+CL)/(a*(k+1)+2*b+2*CL)
        AN.append(prod)  
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
    




'''deglist2, count2 = degreedist2(n)
avgcount2 = np.array(count2) 
for i in range(itern-1):
    d, c = degreedist2(n)
    C = np.array(c)
    avgcount2 = avgcount2 + c
avgcount2 = np.multiply(avgcount2,1/(itern*n)) 
#averaging the frequency of degrees and scaling it by number of vertices to obtain probabilities
AVC11 = [] 
#obtaining CCDF 
for K in range(len(avgcount2)):
    temp = float(sum(avgcount2[i] for i in range(K,len(avgcount2))))
    AVC11.append(temp)
    
plt.scatter(np.log(deglist2 + 1),np.log(AVC11),color='red',label='Simulated for g(x)=x^11') 

    '''
    

    
AN1 = analyticalab() #analytical expression
AVC2 = [] 
for K in range(len(avgcount)):
    temp = float(sum(AN1[i] for i in range(K,len(avgcount))))
    AVC2.append(temp)
plt.scatter(np.log(deglist[:50] + 1),np.log(AVC2[:50]),label='analytical')     
plt.scatter(np.log(deglist + 1),np.log(AVC1),color='yellowgreen',label='simulated')  


plt.xlabel('log(k+1)')
plt.ylabel('log(F(k)')
plt.title('')
plt.legend()
plt.savefig("Linear g with Poi(10) fitness.jpg") 
  
