#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 14:12:13 2021

@author: malhotran
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 15:09:39 2021

@author: Nandan
"""

#SUBLINEAR CASE
#PART 1: Execute the program and save the data to a text file

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time 

start = time.time()


#parameters for linear g(x) = ax+b and fitness f(x)=0
a=2
b=3
itern = 1 #number of iterations 
n =50 #number of vertices
lam = 5

#function to simulate Preferential attachment tree and return degrees 
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
            P[j] = a*indegree[j] + b + j**3
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
    
def degreedist2(n): 
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
    

#analytical expression for "a+b" model (sublinear case) where f(u) = 0 and g(k) = ak+b
def analyticalab():
    AN = []
    mu =a+b
    for k in range(n): 
        prod = 1
        for i in range(k):
            frac = (a*i+b)/(mu+a*i+b) 
            prod = prod*frac
        prod = prod*mu/(mu+a*k+b) 
        AN.append(prod)  
    return AN

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
    

plt.scatter(np.log(deglist + 1),np.log(AVC1),color='green',label='Simulated for g(x)=x^3') 

def analyticalsup():
    AN = []
    for k in range(n):
        den = 2**(k+1)
        frac = 1/den 
        AN.append(frac) 
    return AN

deglist2, count2 = degreedist2(n)
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

    
    

    
AN1 = analyticalsup() #analytical expression
AVC2 = [] 
for K in range(len(avgcount)):
    temp = float(sum(AN1[i] for i in range(K,len(avgcount))))
    AVC2.append(temp)
 
plt.scatter(np.log(deglist + 1),np.log(AVC1),color='green',label='Simulated for f(x)=x^3') 
plt.scatter(np.log(deglist2 + 1),np.log(AVC11),color='red',label='Simulated for f(x)=x^11') 
plt.plot(np.log(deglist[:50] + 1),np.log(AVC2[:50]),color='blue',label='Analytical') 

plt.xlabel('log(k+1)')
plt.ylabel('log(F(k)')
plt.legend()
plt.show()
  

'''  
textfile = open("CCDF_Sublinear_simulated.txt", "w")
for element in AVC1:
    textfile.write(str(element) + "\n")
textfile.close()
   
AN1 = analyticalab() #analytical expression
AVC2 = [] 
for K in range(len(avgcount)):
    temp = float(sum(AN1[i] for i in range(K,len(avgcount))))
    AVC2.append(temp)

textfile = open("CCDF_Sublinear_analytical.txt", "w")
for element in AVC2:
    textfile.write(str(element) + "\n")
textfile.close()

textfile = open("X_Axis_Sublinear.txt", "w")
for element in deglist:
    textfile.write(str(element) + "\n")
textfile.close()

#PART 2: Plotting the data 

simulated = []
analytical =[]
deg =[]
plt.figure()
with open("CCDF_Sublinear_simulated.txt") as f:
    for y in f.readlines():
        simulated.append(float(y))
with open("CCDF_Sublinear_analytical.txt") as f:
    for y in f.readlines():
        analytical.append(float(y))
with open("X_Axis_sublinear.txt") as f:
    for y in f.readlines():
        deg.append(float(y)+1)

plt.plot(np.log(deg),np.log(analytical),label='Analytical CCDF')
plt.scatter(np.log(deg),np.log(simulated),color='red',label='Simulated CCDF')
plt.xlabel('log(k+1)')
plt.ylabel('F(k)')
plt.title('Sublinear Regime')
plt.legend() 
plt.show() '''
print(time.time()-start) 
