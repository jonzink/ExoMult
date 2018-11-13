# ExoMult

This program will simulate the transiting exoplanet population found Kepler from an underlaying power-law distribution. Multiplicity and its effects on detection efficiency are also considered here.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run the software begin by install R. The latest version of R can be downloaded from:
```
https://www.r-project.org
```

## Running ExoMult

Before starting, make sure that "kepler_solar-like_stellar_parameters.csv" and "ExoMult.R" are located in the same directory. Start by opeaning R in this directory.

```
R
```
Now load ExoMult into R
```
source("ExoMult")
```
Now you can execute the ExoMult function
```
ExoMult()
```

ExoMult will return a data frame of the detected planets and export the results to a csv file "simulated_detect_pop.csv"

### The ExoMult Function

The above opperation will run ExoMult with the default setting.
```
ExoMult(rMin=0.5,rMax=16,alpha_1=-1.76,rad_break=2.66,alpha_2=-4.39,pMin=.5,pMax=500,beta_1=0.79,per_break=7.025,beta_2=-0.61,frac_m1=.74,frac_m2=.71,frac_m3=.68,frac_m4=.66,frac_m5=.64,frac_m6=.60,frac_m7=.46,export_csv=TRUE)
```
Arguments

rMin, rMax = The minimum/maximum values of radius considered in this population simulation. The units are provided in Earth radii.
pMin, pMax = The minimum/maximum values of period considered in this population simulation. The units are provided in days.

alpha_1, alpha_2 = The population parameter for the planet radius distribution. The values correspond to exponent of the power-law.
beta_1, beta_2 = The population parameter for the planet period distribution. The values correspond to exponent of the power-law.

rad_break = The radius value where the power-law distribution changes exponents.  
per_break = The radius value where the power-law distribution changes exponents.  

frac_m1,frac_m2,frac_m3,frac_m4,frac_m5,frac_m6,frac_m7 = the fraction of the stellar population with at least m planets.

export_csv = Tell ExoMult whether it should or should not print the results to a csv file.

## Recommendations

To utilize this code for an alternative stellar population, be sure that your stellar sample file has the same format as "kepler_solar-like_stellar_parameters.csv". 

## Please cite as
J. K. Zink, J. L. Christiansen, and  B. M. S. Hansen 2018, MNRAS 


