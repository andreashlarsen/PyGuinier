"""
PyGuninier: make Guinier analysis for SAS Data
"""

## import python packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## import data (has to be adjusted to your data format)
q,I,dI = np.genfromtxt('dat1.dat',skip_header=3,usecols=[0,1,2],unpack=True)

## change 0 to 1 to plot data (to check import)
if 0:
    plt.errorbar(q,I,yerr=dI,linestyle='none',marker='.',color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('q')
    plt.ylabel('Intensity')
    plt.show()


## define linear function
def linfunc(x,a,b):
    y = a*x+b
    return y

## make scan over various qmax values for the Guinier fit
Rg,qRg = [],[]

# adjust range for your data
for idx_qmax in range(30,50,1):
    
    # define relevant quantities: q-squared and ln(I)
    idx = range(idx_qmax)
    qmax = q[idx_qmax]
    q2 = q[idx]**2
    lnI = np.log(I[idx])
    dlnI = dI[idx]/I[idx] # propagated uncertainty

    # linear fit (Guinier fit)
    popt,pcov = curve_fit(linfunc,q2,lnI,sigma=dlnI)
    a = np.sqrt(-3*popt[0])
    Rg.append(a)
    qRg.append(qmax*a)
    fit = linfunc(q2,*popt)
    
    # change 0 to 1 to plot the Guinier fit for this range
    if 0:
        plt.errorbar(q2,lnI,yerr=dlnI,linestyle='none',marker='.',color='red',zorder=0)
        plt.plot(q2,fit,color='black',zorder=1)
        plt.xlabel('q^2')
        plt.ylabel('ln(I)')
        plt.show()

# change 0 to 1 to plot the scan
if 0:
    plt.plot(qRg,Rg,linestyle='none',marker='.',color='black')
    plt.xlabel('qmax*Rg')
    plt.ylabel('Rg')

# choose range, assessed from the scan
# this is one way of estimating Rg and uncertainty on Rg:
# by finding a range of values and taking mean and standard deviation

if 0:
    # change these nubmers for your data and Guinier fit
    sel_min,sel_max = 2,9
    Rg_sel = Rg[sel_min:sel_max] 
    qRg_sel = qRg[sel_min:sel_max]

    plt.plot(qRg_sel,Rg_sel,linestyle='none',marker='o',fillstyle='none',color='green')
    plt.show()

    #estimate Rg error from selected range
    mRg = np.mean(Rg_sel)
    dRg = np.std(Rg_sel)/np.sqrt(len(Rg_sel))

    print('Rg from range:')
    print('Rg=%f +/- %f' % (mRg,dRg))

# choose single value (can be selected using the scan) 
# this is another way (more standard/classic) way of estimating Rg and uncertainties:
# by selecting one range, and propagate uncertainties from fitting uncertainties

# range need to be adjusted for your data! 
idx_sel = range(37)
qmax = q[37]

# define relevant quantities
q2 = q[idx_sel]**2
lnI = np.log(I[idx_sel])
dlnI = dI[idx_sel]/I[idx_sel]

# linear fitting (Guinier fit)
popt,pcov = curve_fit(linfunc,q2,lnI,sigma=dlnI)
Rg_sel = np.sqrt(-3*popt[0])
qRg_sel = qmax*Rg_sel
fit = linfunc(q2,*popt)
I0_sel = np.exp(popt[1])

# assess fit by chi2r
R = (lnI-fit)/dlnI
chi2_0 = np.sum(R**2) # chi-square
df = (len(lnI)-len(popt)) # degrees of freedom
chi2r = chi2_0/df # reduced chi-square (often just called chi-square)

# get chi2 p-value
from scipy.stats import chi2
x=np.linspace(0,500,1000)
dist=chi2.pdf(x,df)
idx_chi2 = np.where(x>chi2_0)
p_chi2r = np.sum(dist[idx_chi2])/np.sum(dist)

print('chi2r = %f, p-value = %f' % (chi2r,p_chi2r))

# change 0 to 1 to illustrate the chi2 p-value
if 0:
    plt.plot(x,dist)
    plt.fill_between(x[idx_chi2],dist[idx_chi2],y2=0)
    plt.show()

# plot the Guinier fit
if 1:
    plt.errorbar(q2,lnI,yerr=dlnI,linestyle='none',marker='.',color='red',zorder=0)
    plt.plot(q2,fit,color='black',zorder=1)
    plt.xlabel('q^2')
    plt.ylabel('ln(I)')
    plt.show()

# propagate fit uncertainties and print results
var_a=pcov[0,0] # fit uncertainty, from cov matrix
da = np.sqrt(var_a)
dRg_sel = abs(-3/(2*Rg_sel)*da) # error propagation
print('Rg = %f +/- %f' % (Rg_sel,dRg_sel))

qRgmin = q[0]*Rg_sel
print('qRg min = %f, qRg max = %f' % (qRgmin,qRg_sel))

var_b=pcov[1,1] # fit uncertainty, from cov matrix
db = np.sqrt(var_b)
dI0_sel = I0_sel*db # error propagation
print('I(0) = %f +/- %f' % (I0_sel,dI0_sel))

