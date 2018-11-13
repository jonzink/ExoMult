#!/usr/bin/Rscript

### E x o M u l t -  Created by Jon Zink ###
### Devolped on R 3.4.3 ###

### If you make use of this code, please cite: ###
### J. K. Zink, J. L. Christiansen, and  B. M. S. Hansen 2018, MNRAS



ExoMult <- function(rMin=0.5,rMax=16,alpha_1=-1.76,rad_break=2.66,alpha_2=-4.39,pMin=.5,pMax=500,beta_1=0.79,per_break=7.025,beta_2=-0.61,frac_m1=.74,frac_m2=.71,frac_m3=.68,frac_m4=.66,frac_m5=.64,frac_m6=.60,frac_m7=.46,export_csv=TRUE){

	###Input Errors
	if(frac_m1<frac_m2 | frac_m2<frac_m3 | frac_m3<frac_m4 | frac_m4<frac_m5 | frac_m5<frac_m6 | frac_m6<frac_m7){
		print("frac values must be in descending order")
		return()
		}
	
	if(rMin>rMax | pMin>pMax){
		print("Population limits are defined incorrectly")
		return()
		}	

	if(rad_break<rMin | rad_break>pMax){
		print("rad_break must be within the population limits")
		return()
		}	

	if(per_break<rMin | per_break>pMax){
		print("per_break must be within the population limits")
		return()
		}			




	print("Loading required libraries")
	tryCatch(library(VGAM),
		error=function(e) {
			install.packages("VGAM", dependencies=TRUE)
			library(VGAM)})

	print("Importing stellar population parameters")
	###Import the stellar population
	tryCatch({star=read.table("kepler_solar-like_stellar_parameters.csv", header=FALSE, sep=",",skip = 2, 
		col.names=c("KIC","T_eff","Teff_errp","Teff_errm","logg","logg_errp","logg_errm","Rad","Rad_errp","Rad_errm","Mass","Mass_errp","Mass_errm",
		"duty","d_span","CDPP1.5","CDPP2","CDPP2.5","CDPP3","CDPP3.5","CDPP4.5","CDPP5","CDPP6","CDPP7.5","CDPP9","CDPP10.5","CDPP12","CDPP12.5","CDPP15"))},
		error=function(e) {
			print("Error loading stellar data")
			return()})

	number_stars=nrow(star)
	gravConstant=6.674e-11*1.98855e30*(86400)^2/(6.957e8)^3


	####Read CDPP Data into a matrix
	CDPPx=c(1.5,2,2.5,3,3.5,4.5,5,6,7.5,9,10.5,12,12.5,15)
	CDPPy=matrix(data=0, nrow=14,ncol=number_stars)

	CDPPy[1,]=star$CDPP1.5
	CDPPy[2,]=star$CDPP2
	CDPPy[3,]=star$CDPP2.5
	CDPPy[4,]=star$CDPP3
	CDPPy[5,]=star$CDPP3.5
	CDPPy[6,]=star$CDPP4.5
	CDPPy[7,]=star$CDPP5
	CDPPy[8,]=star$CDPP6
	CDPPy[9,]=star$CDPP7.5
	CDPPy[10,]=star$CDPP9
	CDPPy[11,]=star$CDPP10.5
	CDPPy[12,]=star$CDPP12
	CDPPy[13,]=star$CDPP12.5
	CDPPy[14,]=star$CDPP15


	######Determine Stellar Limb Darkening parameters
	limb_u1=-.000193*star$T_eff+1.5169
	limb_u2=.000125*star$T_eff-0.4601



	####Broken Powerlaw Distribution Draw####
	brokenPow <- function(n,aCon,bCon,xmix,xmin,xmax){
	 	constant1=((xmix^(aCon+1)-xmin^(aCon+1))/(aCon+1)+(xmix^aCon/xmix^bCon)*(xmax^(bCon+1)-xmix^(bCon+1))/(bCon+1))^-1
	 	constant2=constant1*(xmix^aCon/xmix^bCon)
		turn=constant1*(xmix^(aCon+1)-xmin^(aCon+1))/(aCon+1)
		rValue=runif(n)
		results=ifelse(rValue<=turn,(rValue*(aCon+1)/constant1+xmin^(aCon+1))^(1/(aCon+1)),((rValue-turn)*(bCon+1)/constant2+xmix^(bCon+1))^(1/(bCon+1)))
	  	return(results)
	  	}



	######Function for detection probability#################
	probability_detection<-function(star, period, radius,b,cdpp,u1,u2,plNum){
	
		semiMajor=(period^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
		transit_chord=(1-b^2)^.5
		transit_duration=star$Rad*period/(semiMajor*pi)*transit_chord*24

		number_transits=star$d_span/period
	
		######Window Function#########
		probabilty_window=1-(1-star$duty)^number_transits-number_transits*star$duty*(1-star$duty)^(number_transits-1)-number_transits*(number_transits-1)/2*star$duty^2*(1-star$duty)^(number_transits-2)
		probabilty_window=ifelse(number_transits<1,0,probabilty_window)
		probabilty_window=ifelse(probabilty_window<0,0,probabilty_window)
		probabilty_window=ifelse(probabilty_window>1,1,probabilty_window)


		######Calculate Transit Depth####
		krp=radius/star$Rad*0.00916399
		c0  = 1.0 - (u1+u2)
		omega = c0/4+0/5+(u1+2*u2)/6+0/7-u2/8
		ksq = 1.0 - krp^2
		tmp0 = c0/4 * ksq
		tmp2 = (u1+2*u2)/6* ksq^(3.0/2.0)
		tmp4 = -u2/8 * ksq^2
		depth=(1.0 - (tmp0 + tmp2 + tmp4)/omega) * 1.0e6

		######Calculate MES####
		MES=depth/(cdpp)*1.003*number_transits^.5

		###Calulate the fraction detected at this MES
		if(plNum==1){
			fdet=ifelse(period>200,pgamma(MES-1.0984,shape=18.4119, scale=0.3959)*0.9051,pgamma(MES-0.0102,shape=29.3363, scale=0.2856)*0.9825)
		}	else{	
			fdet=ifelse(period>200,pgamma(MES-2.9774,shape=5.5213, scale=1.2307)*0.7456,pgamma(MES-0.0093,shape=21.3265, scale=0.4203)*0.9276)
		}
		return(fdet*probabilty_window)
		}
	
	
	
	impactDraw <- function(n,star,inclination,average_mutInc,per,int_w,plNum){
		
		if(plNum==1){
			mutual_inclination=0
		}	else{
			mutual_inclination=rrayleigh(n,average_mutInc)*(2*pi)/360
		}
	
		###Draw argument of periapsis
		w=asin(runif(n,0,1))
	
		####Calculate Planet Inclination
		if(plNum==1){
			inclination_m=acos(cos(inclination)*cos(mutual_inclination)+sin(inclination)*sin(mutual_inclination)*cos(int_w))
		}	else{
			inclination_m=acos(cos(inclination)*cos(mutual_inclination)+sin(inclination)*sin(mutual_inclination)*cos(w))
		}
	
		###Planet Semi-major axis
		semiMajor=(per^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	
		return(abs(semiMajor/star$Rad*(sin(inclination_m)*sin(int_w)-cos(inclination_m)*cos(int_w))))
		}



	########Create Blank Arrays

	cdpp_m1=numeric(number_stars)
	cdpp_m2=numeric(number_stars)
	cdpp_m3=numeric(number_stars)
	cdpp_m4=numeric(number_stars)
	cdpp_m5=numeric(number_stars)
	cdpp_m6=numeric(number_stars)
	cdpp_m7=numeric(number_stars)

	rad_sorted_m1=numeric(number_stars)
	rad_sorted_m2=numeric(number_stars)
	rad_sorted_m3=numeric(number_stars)
	rad_sorted_m4=numeric(number_stars)
	rad_sorted_m5=numeric(number_stars)
	rad_sorted_m6=numeric(number_stars)
	rad_sorted_m7=numeric(number_stars)

	per_sorted_m1=numeric(number_stars)
	per_sorted_m2=numeric(number_stars)
	per_sorted_m3=numeric(number_stars)
	per_sorted_m4=numeric(number_stars)
	per_sorted_m5=numeric(number_stars)
	per_sorted_m6=numeric(number_stars)
	per_sorted_m7=numeric(number_stars)

	impact_sorted_m1=numeric(number_stars)
	impact_sorted_m2=numeric(number_stars)
	impact_sorted_m3=numeric(number_stars)
	impact_sorted_m4=numeric(number_stars)
	impact_sorted_m5=numeric(number_stars)
	impact_sorted_m6=numeric(number_stars)
	impact_sorted_m7=numeric(number_stars)

	print("Drawing exoplanet population from broken power-law distribution")
	#######Draw Planet Population
	rad_m1=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m2=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m3=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m4=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m5=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m6=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)
	rad_m7=brokenPow(number_stars,alpha_1,alpha_2,rad_break,rMin,rMax)

	per_m1=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m2=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m3=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m4=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m5=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m6=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)
	per_m7=brokenPow(number_stars,beta_1,beta_2,per_break,pMin,pMax)

	#Remove planets to correct for population fractions
	##In oder to maintain arrays of similar size we set these discared planets to parameter at undetectable limits
	truePl=runif(number_stars,0,1)

	per_m1=ifelse(truePl<frac_m1,per_m1,Inf)
	per_m2=ifelse(truePl<frac_m2,per_m2,Inf)
	per_m3=ifelse(truePl<frac_m3,per_m3,Inf)
	per_m4=ifelse(truePl<frac_m4,per_m4,Inf)
	per_m5=ifelse(truePl<frac_m5,per_m5,Inf)
	per_m6=ifelse(truePl<frac_m6,per_m6,Inf)
	per_m7=ifelse(truePl<frac_m7,per_m7,Inf)

	rad_m1=ifelse(truePl<frac_m1,rad_m1,0)
	rad_m2=ifelse(truePl<frac_m2,rad_m2,0)
	rad_m3=ifelse(truePl<frac_m3,rad_m3,0)
	rad_m4=ifelse(truePl<frac_m4,rad_m4,0)
	rad_m5=ifelse(truePl<frac_m5,rad_m5,0)
	rad_m6=ifelse(truePl<frac_m6,rad_m6,0)
	rad_m7=ifelse(truePl<frac_m7,rad_m7,0)


	###Draw Inclination Parameters
	system_inclination=asin(runif(number_stars,0,1))
	system_average_mutInc=rrayleigh(number_stars,1)
	intial_w=asin(runif(number_stars,0,1))

	###Calculate Impact parameters
	impact_m1=impactDraw(number_stars,star,system_inclination,0,per_m1,intial_w,1)
	impact_m2=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m2,intial_w,2)
	impact_m3=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m3,intial_w,3)
	impact_m4=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m4,intial_w,4)
	impact_m5=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m5,intial_w,5)
	impact_m6=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m6,intial_w,6)
	impact_m7=impactDraw(number_stars,star,system_inclination,system_average_mutInc,per_m7,intial_w,7)


	####Discard planets with impact parameters grreater than 1
	rad_m1=ifelse(impact_m1>1,0,rad_m1)
	rad_m2=ifelse(impact_m2>1,0,rad_m2)
	rad_m3=ifelse(impact_m3>1,0,rad_m3)
	rad_m4=ifelse(impact_m4>1,0,rad_m4)
	rad_m5=ifelse(impact_m5>1,0,rad_m5)
	rad_m6=ifelse(impact_m6>1,0,rad_m6)
	rad_m7=ifelse(impact_m7>1,0,rad_m7)

	### Set an upper limit of 1 for impact parameters to avoid imaginary transit chord values
	impact_m1=ifelse(impact_m1<=1,impact_m1,1)
	impact_m2=ifelse(impact_m2<=1,impact_m2,1)
	impact_m3=ifelse(impact_m3<=1,impact_m3,1)
	impact_m4=ifelse(impact_m4<=1,impact_m4,1)
	impact_m5=ifelse(impact_m5<=1,impact_m5,1)
	impact_m6=ifelse(impact_m6<=1,impact_m6,1)
	impact_m7=ifelse(impact_m7<=1,impact_m7,1)

	####Calculate signal scale
	MES_m1=(rad_m1)^2*per_m1^(-1/3)*(1-impact_m1^2)^(1/4)
	MES_m2=(rad_m2)^2*per_m2^(-1/3)*(1-impact_m2^2)^(1/4)
	MES_m3=(rad_m3)^2*per_m3^(-1/3)*(1-impact_m3^2)^(1/4)
	MES_m4=(rad_m4)^2*per_m4^(-1/3)*(1-impact_m4^2)^(1/4)
	MES_m5=(rad_m5)^2*per_m5^(-1/3)*(1-impact_m5^2)^(1/4)
	MES_m6=(rad_m6)^2*per_m6^(-1/3)*(1-impact_m6^2)^(1/4)
	MES_m7=(rad_m7)^2*per_m7^(-1/3)*(1-impact_m7^2)^(1/4)


	print("Sorting planets by MES values")
	#####Sort MES Values
	for (i in 1:number_stars){
		#####Sorting MES Values
		sortOrder=sort(c(MES_m1[i],MES_m2[i],MES_m3[i],MES_m4[i],MES_m5[i],MES_m6[i],MES_m7[i]))
	
		####Sorting Period Values
		per_sorted_m1[i]=ifelse(MES_m1[i]==sortOrder[7],per_m1[i],ifelse(MES_m2[i]==sortOrder[7],per_m2[i],ifelse(MES_m3[i]==sortOrder[7],per_m3[i],ifelse(MES_m4[i]==sortOrder[7],per_m4[i],ifelse(MES_m5[i]==sortOrder[7],per_m5[i],ifelse(MES_m6[i]==sortOrder[7],per_m6[i],per_m7[i]))))))
		per_sorted_m2[i]=ifelse(MES_m1[i]==sortOrder[6],per_m1[i],ifelse(MES_m2[i]==sortOrder[6],per_m2[i],ifelse(MES_m3[i]==sortOrder[6],per_m3[i],ifelse(MES_m4[i]==sortOrder[6],per_m4[i],ifelse(MES_m5[i]==sortOrder[6],per_m5[i],ifelse(MES_m6[i]==sortOrder[6],per_m6[i],per_m7[i]))))))
		per_sorted_m3[i]=ifelse(MES_m1[i]==sortOrder[5],per_m1[i],ifelse(MES_m2[i]==sortOrder[5],per_m2[i],ifelse(MES_m3[i]==sortOrder[5],per_m3[i],ifelse(MES_m4[i]==sortOrder[5],per_m4[i],ifelse(MES_m5[i]==sortOrder[5],per_m5[i],ifelse(MES_m6[i]==sortOrder[5],per_m6[i],per_m7[i]))))))
		per_sorted_m4[i]=ifelse(MES_m1[i]==sortOrder[4],per_m1[i],ifelse(MES_m2[i]==sortOrder[4],per_m2[i],ifelse(MES_m3[i]==sortOrder[4],per_m3[i],ifelse(MES_m4[i]==sortOrder[4],per_m4[i],ifelse(MES_m5[i]==sortOrder[4],per_m5[i],ifelse(MES_m6[i]==sortOrder[4],per_m6[i],per_m7[i]))))))
		per_sorted_m5[i]=ifelse(MES_m1[i]==sortOrder[3],per_m1[i],ifelse(MES_m2[i]==sortOrder[3],per_m2[i],ifelse(MES_m3[i]==sortOrder[3],per_m3[i],ifelse(MES_m4[i]==sortOrder[3],per_m4[i],ifelse(MES_m5[i]==sortOrder[3],per_m5[i],ifelse(MES_m6[i]==sortOrder[3],per_m6[i],per_m7[i]))))))
		per_sorted_m6[i]=ifelse(MES_m1[i]==sortOrder[2],per_m1[i],ifelse(MES_m2[i]==sortOrder[2],per_m2[i],ifelse(MES_m3[i]==sortOrder[2],per_m3[i],ifelse(MES_m4[i]==sortOrder[2],per_m4[i],ifelse(MES_m5[i]==sortOrder[2],per_m5[i],ifelse(MES_m6[i]==sortOrder[2],per_m6[i],per_m7[i]))))))
		per_sorted_m7[i]=ifelse(MES_m1[i]==sortOrder[1],per_m1[i],ifelse(MES_m2[i]==sortOrder[1],per_m2[i],ifelse(MES_m3[i]==sortOrder[1],per_m3[i],ifelse(MES_m4[i]==sortOrder[1],per_m4[i],ifelse(MES_m5[i]==sortOrder[1],per_m5[i],ifelse(MES_m6[i]==sortOrder[1],per_m6[i],per_m7[i]))))))
	
		####Sorting Radius Values
		rad_sorted_m1[i]=ifelse(MES_m1[i]==sortOrder[7],rad_m1[i],ifelse(MES_m2[i]==sortOrder[7],rad_m2[i],ifelse(MES_m3[i]==sortOrder[7],rad_m3[i],ifelse(MES_m4[i]==sortOrder[7],rad_m4[i],ifelse(MES_m5[i]==sortOrder[7],rad_m5[i],ifelse(MES_m6[i]==sortOrder[7],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m2[i]=ifelse(MES_m1[i]==sortOrder[6],rad_m1[i],ifelse(MES_m2[i]==sortOrder[6],rad_m2[i],ifelse(MES_m3[i]==sortOrder[6],rad_m3[i],ifelse(MES_m4[i]==sortOrder[6],rad_m4[i],ifelse(MES_m5[i]==sortOrder[6],rad_m5[i],ifelse(MES_m6[i]==sortOrder[6],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m3[i]=ifelse(MES_m1[i]==sortOrder[5],rad_m1[i],ifelse(MES_m2[i]==sortOrder[5],rad_m2[i],ifelse(MES_m3[i]==sortOrder[5],rad_m3[i],ifelse(MES_m4[i]==sortOrder[5],rad_m4[i],ifelse(MES_m5[i]==sortOrder[5],rad_m5[i],ifelse(MES_m6[i]==sortOrder[5],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m4[i]=ifelse(MES_m1[i]==sortOrder[4],rad_m1[i],ifelse(MES_m2[i]==sortOrder[4],rad_m2[i],ifelse(MES_m3[i]==sortOrder[4],rad_m3[i],ifelse(MES_m4[i]==sortOrder[4],rad_m4[i],ifelse(MES_m5[i]==sortOrder[4],rad_m5[i],ifelse(MES_m6[i]==sortOrder[4],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m5[i]=ifelse(MES_m1[i]==sortOrder[3],rad_m1[i],ifelse(MES_m2[i]==sortOrder[3],rad_m2[i],ifelse(MES_m3[i]==sortOrder[3],rad_m3[i],ifelse(MES_m4[i]==sortOrder[3],rad_m4[i],ifelse(MES_m5[i]==sortOrder[3],rad_m5[i],ifelse(MES_m6[i]==sortOrder[3],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m6[i]=ifelse(MES_m1[i]==sortOrder[2],rad_m1[i],ifelse(MES_m2[i]==sortOrder[2],rad_m2[i],ifelse(MES_m3[i]==sortOrder[2],rad_m3[i],ifelse(MES_m4[i]==sortOrder[2],rad_m4[i],ifelse(MES_m5[i]==sortOrder[2],rad_m5[i],ifelse(MES_m6[i]==sortOrder[2],rad_m6[i],rad_m7[i]))))))
		rad_sorted_m7[i]=ifelse(MES_m1[i]==sortOrder[1],rad_m1[i],ifelse(MES_m2[i]==sortOrder[1],rad_m2[i],ifelse(MES_m3[i]==sortOrder[1],rad_m3[i],ifelse(MES_m4[i]==sortOrder[1],rad_m4[i],ifelse(MES_m5[i]==sortOrder[1],rad_m5[i],ifelse(MES_m6[i]==sortOrder[1],rad_m6[i],rad_m7[i]))))))
	
		####Sorting Impact Values
		impact_sorted_m1[i]=ifelse(MES_m1[i]==sortOrder[7],impact_m1[i],ifelse(MES_m2[i]==sortOrder[7],impact_m2[i],ifelse(MES_m3[i]==sortOrder[7],impact_m3[i],ifelse(MES_m4[i]==sortOrder[7],impact_m4[i],ifelse(MES_m5[i]==sortOrder[7],impact_m5[i],ifelse(MES_m6[i]==sortOrder[7],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m2[i]=ifelse(MES_m1[i]==sortOrder[6],impact_m1[i],ifelse(MES_m2[i]==sortOrder[6],impact_m2[i],ifelse(MES_m3[i]==sortOrder[6],impact_m3[i],ifelse(MES_m4[i]==sortOrder[6],impact_m4[i],ifelse(MES_m5[i]==sortOrder[6],impact_m5[i],ifelse(MES_m6[i]==sortOrder[6],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m3[i]=ifelse(MES_m1[i]==sortOrder[5],impact_m1[i],ifelse(MES_m2[i]==sortOrder[5],impact_m2[i],ifelse(MES_m3[i]==sortOrder[5],impact_m3[i],ifelse(MES_m4[i]==sortOrder[5],impact_m4[i],ifelse(MES_m5[i]==sortOrder[5],impact_m5[i],ifelse(MES_m6[i]==sortOrder[5],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m4[i]=ifelse(MES_m1[i]==sortOrder[4],impact_m1[i],ifelse(MES_m2[i]==sortOrder[4],impact_m2[i],ifelse(MES_m3[i]==sortOrder[4],impact_m3[i],ifelse(MES_m4[i]==sortOrder[4],impact_m4[i],ifelse(MES_m5[i]==sortOrder[4],impact_m5[i],ifelse(MES_m6[i]==sortOrder[4],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m5[i]=ifelse(MES_m1[i]==sortOrder[3],impact_m1[i],ifelse(MES_m2[i]==sortOrder[3],impact_m2[i],ifelse(MES_m3[i]==sortOrder[3],impact_m3[i],ifelse(MES_m4[i]==sortOrder[3],impact_m4[i],ifelse(MES_m5[i]==sortOrder[3],impact_m5[i],ifelse(MES_m6[i]==sortOrder[3],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m6[i]=ifelse(MES_m1[i]==sortOrder[2],impact_m1[i],ifelse(MES_m2[i]==sortOrder[2],impact_m2[i],ifelse(MES_m3[i]==sortOrder[2],impact_m3[i],ifelse(MES_m4[i]==sortOrder[2],impact_m4[i],ifelse(MES_m5[i]==sortOrder[2],impact_m5[i],ifelse(MES_m6[i]==sortOrder[2],impact_m6[i],impact_m7[i]))))))
		impact_sorted_m7[i]=ifelse(MES_m1[i]==sortOrder[1],impact_m1[i],ifelse(MES_m2[i]==sortOrder[1],impact_m2[i],ifelse(MES_m3[i]==sortOrder[1],impact_m3[i],ifelse(MES_m4[i]==sortOrder[1],impact_m4[i],ifelse(MES_m5[i]==sortOrder[1],impact_m5[i],ifelse(MES_m6[i]==sortOrder[1],impact_m6[i],impact_m7[i]))))))
		}

	####Calculate Semi-Major axis
	semiMajor_sorted_m1=(per_sorted_m1^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m2=(per_sorted_m2^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m3=(per_sorted_m3^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m4=(per_sorted_m4^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m5=(per_sorted_m5^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m6=(per_sorted_m6^2*star$Mass*gravConstant/(4*pi^2))^(1/3)
	semiMajor_sorted_m7=(per_sorted_m7^2*star$Mass*gravConstant/(4*pi^2))^(1/3)


	###Calculate Transit Chords
	tran_chord_m1=(1-impact_sorted_m1^2)^.5
	tran_chord_m2=(1-impact_sorted_m2^2)^.5
	tran_chord_m3=(1-impact_sorted_m3^2)^.5
	tran_chord_m4=(1-impact_sorted_m4^2)^.5
	tran_chord_m5=(1-impact_sorted_m5^2)^.5
	tran_chord_m6=(1-impact_sorted_m6^2)^.5
	tran_chord_m7=(1-impact_sorted_m7^2)^.5


	#Calculate Transit Duration
	tran_duration_m1=star$Rad*per_sorted_m1/(semiMajor_sorted_m1*pi)*tran_chord_m1*24
	tran_duration_m2=star$Rad*per_sorted_m2/(semiMajor_sorted_m2*pi)*tran_chord_m2*24
	tran_duration_m3=star$Rad*per_sorted_m3/(semiMajor_sorted_m3*pi)*tran_chord_m3*24
	tran_duration_m4=star$Rad*per_sorted_m4/(semiMajor_sorted_m4*pi)*tran_chord_m4*24
	tran_duration_m5=star$Rad*per_sorted_m5/(semiMajor_sorted_m5*pi)*tran_chord_m5*24
	tran_duration_m6=star$Rad*per_sorted_m6/(semiMajor_sorted_m6*pi)*tran_chord_m6*24
	tran_duration_m7=star$Rad*per_sorted_m7/(semiMajor_sorted_m7*pi)*tran_chord_m7*24

	#Remove NaN from transit durations
	tran_duration_m1=ifelse(is.nan(tran_duration_m1),Inf,tran_duration_m1)
	tran_duration_m2=ifelse(is.nan(tran_duration_m2),Inf,tran_duration_m2)
	tran_duration_m3=ifelse(is.nan(tran_duration_m3),Inf,tran_duration_m3)
	tran_duration_m4=ifelse(is.nan(tran_duration_m4),Inf,tran_duration_m4)
	tran_duration_m5=ifelse(is.nan(tran_duration_m5),Inf,tran_duration_m5)
	tran_duration_m6=ifelse(is.nan(tran_duration_m6),Inf,tran_duration_m6)
	tran_duration_m7=ifelse(is.nan(tran_duration_m7),Inf,tran_duration_m7)

	print("Interpolating stellar fluctuation values for each transit duration")
	###Interpolate the CDPP values for each star
	for (w in  1:number_stars) {
		cdpp_m1[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m1[w], rule=2)$y)
		cdpp_m2[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m2[w], rule=2)$y)
		cdpp_m3[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m3[w], rule=2)$y)
		cdpp_m4[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m4[w], rule=2)$y)
		cdpp_m5[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m5[w], rule=2)$y)
		cdpp_m6[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m6[w], rule=2)$y)
		cdpp_m7[w]=as.vector(approx(x=CDPPx,y=CDPPy[,w], xout=tran_duration_m7[w], rule=2)$y)
		}


	#Calculate the probabilty of detection for each planet
	prob_detection_m1=probability_detection(star, per_sorted_m1, rad_sorted_m1, impact_sorted_m1,cdpp_m1,limb_u1,limb_u2,1)
	prob_detection_m2=probability_detection(star, per_sorted_m2, rad_sorted_m2, impact_sorted_m2,cdpp_m2,limb_u1,limb_u2,2)
	prob_detection_m3=probability_detection(star, per_sorted_m3, rad_sorted_m3, impact_sorted_m3,cdpp_m3,limb_u1,limb_u2,3)
	prob_detection_m4=probability_detection(star, per_sorted_m4, rad_sorted_m4, impact_sorted_m4,cdpp_m4,limb_u1,limb_u2,4)
	prob_detection_m5=probability_detection(star, per_sorted_m5, rad_sorted_m5, impact_sorted_m5,cdpp_m5,limb_u1,limb_u2,5)
	prob_detection_m6=probability_detection(star, per_sorted_m6, rad_sorted_m6, impact_sorted_m6,cdpp_m6,limb_u1,limb_u2,6)
	prob_detection_m7=probability_detection(star, per_sorted_m7, rad_sorted_m7, impact_sorted_m7,cdpp_m7,limb_u1,limb_u2,7)



	#Remove planets based on their detection probabilty
	random_detection_threshold=runif(number_stars,0,1)
	rad_sorted_m1=ifelse(prob_detection_m1>random_detection_threshold,rad_sorted_m1,NA)
	rad_sorted_m2=ifelse(prob_detection_m2>random_detection_threshold,rad_sorted_m2,NA)
	rad_sorted_m3=ifelse(prob_detection_m3>random_detection_threshold,rad_sorted_m3,NA)
	rad_sorted_m4=ifelse(prob_detection_m4>random_detection_threshold,rad_sorted_m4,NA)
	rad_sorted_m5=ifelse(prob_detection_m5>random_detection_threshold,rad_sorted_m5,NA)
	rad_sorted_m6=ifelse(prob_detection_m6>random_detection_threshold,rad_sorted_m6,NA)
	rad_sorted_m7=ifelse(prob_detection_m7>random_detection_threshold,rad_sorted_m7,NA)

	per_sorted_m1=ifelse(prob_detection_m1>random_detection_threshold,per_sorted_m1,NA)
	per_sorted_m2=ifelse(prob_detection_m2>random_detection_threshold,per_sorted_m2,NA)
	per_sorted_m3=ifelse(prob_detection_m3>random_detection_threshold,per_sorted_m3,NA)
	per_sorted_m4=ifelse(prob_detection_m4>random_detection_threshold,per_sorted_m4,NA)
	per_sorted_m5=ifelse(prob_detection_m5>random_detection_threshold,per_sorted_m5,NA)
	per_sorted_m6=ifelse(prob_detection_m6>random_detection_threshold,per_sorted_m6,NA)
	per_sorted_m7=ifelse(prob_detection_m7>random_detection_threshold,per_sorted_m7,NA)
	
	print("Removing undetected systems")
	###Move data to date frame and remove non-detected systems
	return_df=data.frame(star$KIC,rad_sorted_m1,per_sorted_m1,rad_sorted_m2,per_sorted_m2,rad_sorted_m3,per_sorted_m3,rad_sorted_m4,per_sorted_m4,rad_sorted_m5,per_sorted_m5,rad_sorted_m6,per_sorted_m6,rad_sorted_m7,per_sorted_m7)
	return_df=return_df[!is.na(return_df["rad_sorted_m1"]),]

	if(export_csv){
		#Print to csv
		print("Printing detected population to CSV file...simulated_detect_pop.csv")
		write.csv(return_df, file = "Simulated_detected_population.csv",row.names=FALSE, na="",quote = FALSE)
		}

	return(return_df)
	}