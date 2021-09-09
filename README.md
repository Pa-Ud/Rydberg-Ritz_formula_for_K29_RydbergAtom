# Rydberg-Ritz_formula_for_K29_RydbergAtom
# This is python code
# Authour Mr Parinya Udommai, Date Sept 2021
#####################################################################

#Use Rydberg-Ritz formula to calculate 
#https://iopscience.iop.org/article/10.1088/0031-8949/27/4/012/pdf
import math
import numpy as np
import plotly as py
import plotly.graph_objs as go
from plotly import tools
from plotly.subplots import make_subplots
py.offline.init_notebook_mode(connected=True)
from array import array
import csv
from scipy import optimize
from scipy.optimize import curve_fit
import os
import plotly.express as px
import time
from scipy import stats
import time
#start = time.time()
#end = time.time()

def readParameters(filename):
    clight = 299792458
    f = open(filename, "r")
    # first line is a header so ignore it
    f.readline()
    L,j,a,b,c,d,e,Ei = [],[],[],[],[],[],[],[]
    for line in f:
        L = np.append(L,line.split(',')[0])
        j = np.append(j,str(int(2*float(line.split(',')[1])))+'/2')
        a = np.append(a,float(line.split(',')[2]))
        b = np.append(b,float(line.split(',')[3]))
        c = np.append(c,float(line.split(',')[4]))
        d = np.append(d,float(line.split(',')[5]))
        e = np.append(e,float(line.split(',')[6]))
        #Ei= np.append(Ei,29979.2458*float(line.split(',')[7])) #MHz ionization
        Ei= np.append(Ei,1e-6*clight*(1e2*float(line.split(',')[7])))
    return L,j,a,b,c,d,e,Ei

def RRF(n,L,j,a,b,c,d,e,Ei):
    R39K = 3289820654.552509
    #Find the data row to use
    if L == 'S':
        r = [0] #r stands for row.
    elif L == 'D':
        r = [1,2]
    elif L == 'P':
        r = [3,4]
    EnL = [] #Energies at n and L
    #d stands for devider
    d0 = n
    for i in r:
        d1 = -a[i]
        dd = d0 + d1
        d2 = -b[i]/dd**2
        d3 = -c[i]/dd**4
        d4 = -d[i]/dd**6
        d5 = -e[i]/dd**8
        divider = (d0+d1+d2+d3+d4+d5)**2
        Ei0 = 1049568101.92 #1048371288.1797894+288.6 
        EnLr = 0*Ei0 + 1*(Ei[i]+288.6) - R39K/divider
        EnL.append(EnLr) 
    return EnL
######################################################################################################
filename = 'Rydberg-RitzFormulaParameters.csv'
clight = 299792458
L,j,a,b,c,d,e,Ei = readParameters(filename)
a4PhalfF2 = 389286058.716+288.6+20.8 #MHz
Lstate = ['S','D','P']
jstate = [[0.5], [1.5,2.5], [0.5,1.5]]
Rydstate, Rydfreq, Rydwave = [],[],[]
nStates = np.linspace(55,75,21)
nStates = np.append(nStates, np.linspace(115,130,16) )
print(nStates)
for v in range(0,len(Lstate)):
    lstate = Lstate[v]
    for q in nStates:#range(50,130):
        q = int(q)
        Ans = RRF(q,lstate,j,a,b,c,d,e,Ei)
        for u in range(0,len(Ans)):
            freq = (Ans[u] - a4PhalfF2)*1e-6   
            #stringOUT = str(q)+lstate+'(j='+str(int(2*jstate[v][u]))+'/2) '+'%10.7f' %energy 
            #print(stringOUT) #Use Rydberg-Ritz formula to calculate 
            Rydstate.append(str(q)+lstate+' (j='+str(int(2*jstate[v][u]))+'/2)')
            bluefreq = '%10.7f' %freq
            Rydfreq.append(bluefreq)
            bluewave = '%10.7f' %(1e-3*clight/freq)
            Rydwave.append(bluewave)
            
sortA, Rydstate = zip(*sorted(zip(Rydfreq, Rydstate)))
Rydfreq, Rydwave = zip(*sorted(zip(Rydfreq, Rydwave)))
#for w in range(0,len(Rydfreq)):
    #print(Rydstate[w],Rydwave[w],Rydfreq[w])
with open('RydbergRitzCalculation.csv','w') as f1:
    writer=csv.writer(f1, delimiter=',',lineterminator='\n',)
    for w in range(0,len(Rydfreq)):
        row = Rydstate[w],Rydwave[w],Rydfreq[w]
        writer.writerow(row)    
 ##################################################################################################################
