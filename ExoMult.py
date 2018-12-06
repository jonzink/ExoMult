#!/usr/bin/python

### E x o M u l t -  Created by Jon Zink ###
### Devolped on R 3.4.3 ###
### Devolped on python 3.7.1 ###

### If you make use of this code, please cite: ###
### J. K. Zink, J. L. Christiansen, and  B. M. S. Hansen 2018, MNRAS




import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
import numpy as np

pandas2ri.activate()
robjects.r('''
       source('ExoMult.R')
''')

def ExoMult(rMin=0.5,rMax=16,alpha_1=-1.65,rad_break=2.66,alpha_2=-4.35,pMin=.5,pMax=500,beta_1=0.76,per_break=7.09,beta_2=-0.64,
 		mut_Ray=1,ecc_alpha=0,ecc_beta=1,frac_m1=.72,frac_m2=.68,frac_m3=.66,frac_m4=.63,frac_m5=.60,frac_m6=.55,frac_m7=.40,export_csv=True):
		
		ExoMultR=robjects.globalenv['ExoMult']
		
		df=ExoMultR(rMin,rMax,alpha_1,rad_break,alpha_2,pMin,pMax,beta_1,per_break,beta_2,
 		mut_Ray,ecc_alpha,ecc_beta,frac_m1,frac_m2,frac_m3,frac_m4,frac_m5,frac_m6,frac_m7,export_csv)
		
		return_df=pandas2ri.ri2py(df)
		return_df.index += 1
		
		return(return_df)
		
def ExoMult_Prob(radius,period,eccentricity,mut_Ray=0):
		
		ExoMultR=robjects.globalenv['ExoMult.Prob']
		
		df=ExoMultR(radius,period,eccentricity,mut_Ray)
		
		return_df=pandas2ri.ri2py(df)
		return_df.index += 1 
		
		return(return_df)		


