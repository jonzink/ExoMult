#!/usr/bin/python3

### E x o M u l t -  Created by Jon Zink ###


### If you make use of this code, please cite: ###
### J. K. Zink, J. L. Christiansen, and  B. M. S. Hansen 2018, MNRAS



import numpy as np
import pandas as pd
#import math
from scipy.stats import gamma
from scipy.special import gammaincc
from scipy.special import gamma as gammaFunc


# d_span=80 ##days###
gravConstant=6.674e-11*1.98855e30*(86400)**2/(6.957e8)**3


class PopParams:
    """Initialize Population Model.
    """
    
    def __init__(self, stellar=None, prettyPrint=False, mission="Kepler", indep=False, export_csv=True, campaignCom=False, fast=False):
        """Help on ExoMult:
    
            Parameters
            ----------
            stellar : panda dataframe
                stellar sample which must include  Rad, Mass, limb1, limb2, d_span, dutycycle,CDPP values,CDPP slope, and CDPP int.
                Alongside the stars, the planets in your simulation populations must also be appended to the stellar dataframe.
            prettyPrint : bool, optional
                prints out summary stats.
            mission : string, optional
                Set mission completeness, "Kepler" or "K2"
            indep : bool, optional
                Ignore correlations assoicated with multiplicity. 
            export_csv : bool, optional
                Prints the output planet sample to a CVS
            campaignCom : bool, optional
                Use completeness on a campain by campain bases for K2. This should only be used the the stellar
                sample is heavily baises towards one campaign, as the completeness uncertainty it much higher here.
            fast : bool, optional
                Runs the simulations faster but loses some of the precision.
        
        Notes
        
        running:
         
        params=ExoMult.PopParms(stellar)
        outputPopulation=SimulateDetection(params)
        
        will return a observationally baises sample of planets and their host stars.
        """

        self.mission=mission
        self.indep=indep
        self.prettyPrint=prettyPrint
        self.campaignCom=campaignCom
        self.fast=fast
        self.export_csv=export_csv

        ###Input Errors
            
        if((mission!="Kepler") & (mission!="K2")):
            print("ERROR: Mission name must either be Kepler or K2")
            
        
        ###Import the stellar population
        if isinstance(stellar, pd.DataFrame):
            try:
                if mission=="Kepler":
                    cdx=np.array([1.5,2,2.5,3,3.5,4.5,5,6,7.5,9,10.5,12,12.5,15])
                    # stellar["slope1"]=(stellar.CDPP2-np.array(stellar.CDPP1_5,dtype=float))/(cdx[1]-cdx[0])
                    # stellar["int1"]=stellar.CDPP2-stellar.slope1*cdx[1]
                    # stellar["slope2"]=(stellar.CDPP2_5-stellar.CDPP2)/(cdx[2]-cdx[1])
                    # stellar["int2"]=stellar.CDPP2_5-stellar.slope2*cdx[2]
                    # stellar["slope3"]=(stellar.CDPP3-stellar.CDPP2_5)/(cdx[3]-cdx[2])
                    # stellar["int3"]=stellar.CDPP3-stellar.slope3*cdx[3]
                    # stellar["slope4"]=(stellar.CDPP3_5-stellar.CDPP3)/(cdx[4]-cdx[3])
                    # stellar["int4"]=stellar.CDPP3_5-stellar.slope4*cdx[4]
                    # stellar["slope5"]=(stellar.CDPP4_5-stellar.CDPP3_5)/(cdx[5]-cdx[4])
                    # stellar["int5"]=stellar.CDPP4_5-stellar.slope5*cdx[5]
                    # stellar["slope6"]=(stellar.CDPP5-stellar.CDPP4_5)/(cdx[6]-cdx[5])
                    # stellar["int6"]=stellar.CDPP5-stellar.slope6*cdx[6]
                    # stellar["slope7"]=(stellar.CDPP6-stellar.CDPP5)/(cdx[7]-cdx[6])
                    # stellar["int7"]=stellar.CDPP6-stellar.slope7*cdx[7]
                    # stellar["slope8"]=(stellar.CDPP7_5-stellar.CDPP6)/(cdx[8]-cdx[7])
                    # stellar["int8"]=stellar.CDPP7_5-stellar.slope8*cdx[8]
                    # stellar["slope9"]=(stellar.CDPP9-stellar.CDPP7_5)/(cdx[9]-cdx[8])
                    # stellar["int9"]=stellar.CDPP9-stellar.slope9*cdx[9]
                    # stellar["slope10"]=(stellar.CDPP10_5-stellar.CDPP9)/(cdx[10]-cdx[9])
                    # stellar["int10"]=stellar.CDPP10_5-stellar.slope10*cdx[10]
                    # stellar["slope11"]=(stellar.CDPP12-stellar.CDPP10_5)/(cdx[11]-cdx[10])
                    # stellar["int11"]=stellar.CDPP12-stellar.slope11*cdx[11]
                    # stellar["slope12"]=(stellar.CDPP12_5-stellar.CDPP12)/(cdx[12]-cdx[11])
                    # stellar["int12"]=stellar.CDPP12_5-stellar.slope12*cdx[12]
                    # stellar["slope13"]=(stellar.CDPP15-stellar.CDPP12_5)/(cdx[13]-cdx[12])
                    # stellar["int13"]=stellar.CDPP15-stellar.slope13*cdx[13]
                    self.stellar=stellar
                if mission=="K2":
                     cdx=np.array([1,1.5,2,2.5,3,4,5,6,7,8,9,10])
                     # stellar["slope1"]=(stellar.CDPP1_5-stellar.CDPP1)/(cdx[1]-cdx[0])
                     # stellar["int1"]=stellar.CDPP1_5-stellar.slope1*cdx[1]
                     # stellar["slope2"]=(stellar.CDPP2-stellar.CDPP1_5)/(cdx[2]-cdx[1])
                     # stellar["int2"]=stellar.CDPP2-stellar.slope2*cdx[2]
                     # stellar["slope3"]=(stellar.CDPP2_5-stellar.CDPP2)/(cdx[3]-cdx[2])
                     # stellar["int3"]=stellar.CDPP2_5-stellar.slope3*cdx[3]
                     # stellar["slope4"]=(stellar.CDPP3-stellar.CDPP2_5)/(cdx[4]-cdx[3])
                     # stellar["int4"]=stellar.CDPP3-stellar.slope4*cdx[4]
                     # stellar["slope5"]=(stellar.CDPP4-stellar.CDPP3)/(cdx[5]-cdx[4])
                     # stellar["int5"]=stellar.CDPP4-stellar.slope5*cdx[5]
                     # stellar["slope6"]=(stellar.CDPP5-stellar.CDPP4)/(cdx[6]-cdx[5])
                     # stellar["int6"]=stellar.CDPP5-stellar.slope6*cdx[6]
                     # stellar["slope7"]=(stellar.CDPP6-stellar.CDPP5)/(cdx[7]-cdx[6])
                     # stellar["int7"]=stellar.CDPP6-stellar.slope7*cdx[7]
                     # stellar["slope8"]=(stellar.CDPP7-stellar.CDPP6)/(cdx[8]-cdx[7])
                     # stellar["int8"]=stellar.CDPP7-stellar.slope8*cdx[8]
                     # stellar["slope9"]=(stellar.CDPP8-stellar.CDPP7)/(cdx[9]-cdx[8])
                     # stellar["int9"]=stellar.CDPP8-stellar.slope9*cdx[9]
                     # stellar["slope10"]=(stellar.CDPP9-stellar.CDPP8)/(cdx[10]-cdx[9])
                     # stellar["int10"]=stellar.CDPP9-stellar.slope10*cdx[10]
                     # stellar["slope11"]=(stellar.CDPP10-stellar.CDPP9)/(cdx[11]-cdx[10])
                     # stellar["int11"]=stellar.CDPP10-stellar.slope11*cdx[11]
                     self.stellar=stellar


            except:
                print("Error loading stellar data")

                
def mass(r):
    mass=np.where(r<1.5,((2.43+3.39*r)/5.51)*r**3,np.where(r<4,2.69*r**.93,np.where(r<=9,.86*r**1.89,318)))
    return mass 
    
def metalFunction(x,beta,mmax,mmin):
    const=np.max([10**(beta*mmax),10**(beta*mmin)])
    return(10**(beta*x)/const)    
    
def TeffFunction(x,alpha,beta,mmax,mmin,mid=None):
    if mid!=None:
        const1=10**(alpha*mid)
        const2=10**(beta*mid)
        stat=np.where(x<mid,10**(alpha*x),const1/const2*10**(beta*x))
        const=np.max([const1/const2*10**(beta*mmax),10**(alpha*mmin),10**(alpha*mid)])
        return(stat/const)
    else:
        const=np.max([10**(alpha*mmax),10**(alpha*mmin)])
        return(10**(alpha*x)/const)
        

def convert2flux(mass,teff,radius,period):
    #mass = [Solar mass]
    #teff = [k]
    #radius = [Solar Radius]
    #period= [days]
    #return flux in earth insolation units
    
    # Define physical constants
    sigma = 5.67e-8 # Stefan-Boltzmann constant in W m^-2 K^-4
    au = 1.496e11 # 1 astronomical unit in meters
    
    # Define input parameters
    Mstar = mass # stellar mass in solar masses
    Teff = teff # effective temperature in Kelvin
    Rstar = radius # stellar radius in solar radii
    P = period # planet period in days
    a = ((6.6743e-11*Mstar*1.989e30)*(P*24*3600)**2/(4*math.pi**2))**(1/3) # planet semi-major axis in meters
    #print(a/1e11)
    
    # Calculate insolation flux
    F = (sigma / ((a**2))) * (Teff)**4 * (Rstar*695510000)**2
    
    # Convert to Earth insolation flux
    F_Earth = F / 1361 # Earth's solar constant in W/m^2
    return(F_Earth)




def triPower(n,xmin,xmax,a,b,c,p1,p2):
    c1=(((p1**(a+1)-(xmin**(a+1)))/(a+1))
        +(p1**a/p1**b)*((p2**(b+1)-(p1**(b+1)))/(b+1))
        +(p1**a/p1**b)*(p2**b/p2**c)*(xmax**(c+1)-(p2**(c+1)))/(c+1))**(-1)
    c2=c1*p1**a/p1**b
    c3=c2*p2**b/p2**c
    
    turn1=c1*(p1**(a+1)-xmin**(a+1))/(a+1)
    turn2=c2*(p2**(b+1)-p1**(b+1))/(b+1)
    turn3=c3*(xmax**(c+1)-p2**(c+1))/(c+1)
    print(turn1)
    print(turn2)
    print(turn3)
    
    cdf=np.random.random_sample(n)
    
    m=np.where(cdf>(turn2+turn1),((cdf-turn2-turn1)*(c+1)/c3+p2**(c+1))**(1/(c+1)),
               np.where(cdf>turn1,((cdf-turn1)*(b+1)/c2+p1**(b+1))**(1/(b+1)),
                        (cdf*(a+1)/c1+xmin**(a+1))**(1/(a+1))))
    
    return(m)

def hillRad(starMass,a,a2,mass1,mass2):
	return np.where(mass1==0,200,np.where(mass2==0,200,abs(a-a2)/(((mass1+mass2)*0.000003003/(3*starMass))**(1/3)*(a+a2)/2)))

def brokenPow(n,xmin,xmax,aCon,bCon=None,xmix=None):
    foo=False
    if xmix==None:
        bCon=aCon
        xmix=xmax
        foo=True
    constant1=((xmix**(aCon+1)-xmin**(aCon+1))/(aCon+1)+(xmix**aCon/xmix**bCon)*(xmax**(bCon+1)-xmix**(bCon+1))/(bCon+1))**-1
    if foo:
        rValue=np.random.uniform(size=n)
        results=(rValue*(aCon+1)/constant1+xmin**(aCon+1))**(1/(aCon+1))
        return(results)
    else:    
        constant2=constant1*(xmix**aCon/xmix**bCon)
        turn=constant1*(xmix**(aCon+1)-xmin**(aCon+1))/(aCon+1)
        rValue=np.random.uniform(size=n)
        results=np.where(rValue<=turn,(rValue*(aCon+1)/constant1+xmin**(aCon+1))**(1/(aCon+1)),((rValue-turn)*(bCon+1)/constant2+xmix**(bCon+1))**(1/(bCon+1)))
        return(results)

def powerlawKnee(n,aCon,bCon,xmix,xmin,xmax):
    def gammaHere(x,y):
        return(gammaFunc(x)*gammaincc(x,y))

    def integral(x,aCon,bCon,xmix):
        g=aCon-bCon
        guy=(x**(1 + bCon)*(g + ((1 + bCon)*gammaHere((1 + bCon)/g, (x/xmix)**g))/((x/xmix)**g)**((1 + bCon)/g)))/((1 + bCon)*g)
        return(guy)
    
    constant=integral(xmax,aCon,bCon,xmix)-integral(xmin,aCon,bCon,xmix)
    rand=np.linspace(xmin,xmax,10000)    
    cdf=1/constant*(integral(rand,aCon,bCon,xmix)-integral(xmin,aCon,bCon,xmix))
    f = interp.interp1d(cdf, rand)
    rValue=np.random.uniform(size=n)
    outPut=f(rValue)
    return(outPut)

def rrayleigh(n,sigma):
    yValue=np.random.uniform(size=n)
    invCDF=np.sqrt(-2*sigma**2*np.log(1-yValue))
    return(invCDF)
    
def splitNorm(u,sigma1,sigma2,n):
    yValue=np.random.uniform(size=n)
    dist1=abs(np.random.normal(0,sigma1,n))+u
    dist2=-abs(np.random.normal(0,sigma2,n))+u
    return(np.where(yValue>.5,dist1,dist2))    
      		
def MES_calc(star, period, radius,cdpp, b, mission):
	######Calculate Transit Depth####
    MES=np.zeros(len(star))
    limb_u1=star.limb1
    limb_u2=star.limb2
    krp=radius/star.Rad*0.00916399
    c0  = 1.0 - (limb_u1+limb_u2)
    omega = c0/4+0/5+(limb_u1+2*limb_u2)/6+0/7-limb_u2/8
    ksq = 1.0 - krp**2
    tmp0 = c0/4 * ksq
    tmp2 = (limb_u1+2*limb_u2)/6* ksq**(3.0/2.0)
    tmp4 = -limb_u2/8 * ksq**2
    depth=(1.0 - (tmp0 + tmp2 + tmp4)/omega) * 1.0e6
        
######Calculate MES####
    number_transits=star.d_span/period
    if mission=="Kepler":
        MES=depth/(cdpp)*1.003*(number_transits*star.dutycycle)**.5
    if mission=="K2":
        MES=depth/(cdpp)*0.9488*np.floor(number_transits)**.5   
        #0.9488     
        #*.8735
    return(MES)
    
    
    
######Function for detection probability#################
def probability_detection(star, period, radius,cdpp,plNum, b, mission,campCon):
    
    ######Calculate MES####
    MES=MES_calc(star, period, radius,cdpp, b,mission)
    
    ######Window Function#########
    number_transits=star.d_span/period
    Ntr=number_transits
    
    ###Calulate the fraction detected at this MES
    if mission=="Kepler":
        
        #fdet=gamma.cdf(MES,29.3363,loc=0.0102,scale=0.2856)*0.9825
        Ntr=np.floor(Ntr*star.dutycycle)
        number_transits=number_transits*star.dutycycle
        probabilty_window=np.ones(len(star))
        probabilty_window=np.where(number_transits<=3,number_transits-2,probabilty_window)
        probabilty_window=np.where(number_transits>3,1,probabilty_window)
        probabilty_window=np.where(number_transits<2,0,probabilty_window)
        
        fdet=np.where(Ntr==3,gamma.cdf(MES,a=33.3884,loc=0,scale=0.264472)*0.699093,
                                      np.where(Ntr==4,gamma.cdf(MES,a=32.8860,loc=0,scale=0.269577)*0.768366,
                                      np.where(Ntr==5,gamma.cdf(MES,a=31.5196,loc=-0,scale=0.282741)*0.833673,
                                      np.where(Ntr==6,gamma.cdf(MES,a=30.9919,loc=-0,scale=0.286979)*0.859865,
                                      np.where((Ntr>=7) & (Ntr<=9) ,gamma.cdf(MES,a=30.1906,loc=-0,scale=0.294688)*0.875042,
                                      np.where((Ntr>=10)& (Ntr<=18) ,gamma.cdf(MES,a=31.6342,loc=0,scale=0.279425)*0.886144,
                                      np.where((Ntr>=19)& (Ntr<=36) ,gamma.cdf(MES,a=32.6448,loc=-0,scale=0.268898)*0.889724,
                                      np.where(Ntr>=37,gamma.cdf(MES,a=27.8185,loc= 0,scale=0.32432)*0.945075,0))))))))

    if mission=="K2":
        
        probabilty_window=np.ones(len(star))
        probabilty_window=np.where(number_transits<=3,number_transits-2,probabilty_window)
        probabilty_window=np.where(number_transits>3,1,probabilty_window)
        probabilty_window=np.where(number_transits<2,0,probabilty_window)
        
        if campCon==False:
            fdet=np.where(period>26,0.4635/(1+np.exp(-.6607*(MES-10.5441))),0.6619/(1+np.exp(-.6231*(MES-10.9072))))
            fdet=np.where(MES==0,0,fdet)
        else:
            fdet=np.where(star.Campaign==1,0.3923/(1+np.exp(-.7654*(MES-11.3914))),
            np.where(star.Campaign==2,0.6430/(1+np.exp(-.7173*(MES-10.8544))),
            np.where(star.Campaign==3,0.7462/(1+np.exp(-.6689*(MES-10.5701))),
            np.where(star.Campaign==4,0.6734/(1+np.exp(-.6344*(MES-11.1443))),
            np.where(star.Campaign==5,0.4425/(1+np.exp(-.5923*(MES-11.3923))),
            np.where(star.Campaign==6,0.7654/(1+np.exp(-.5759*(MES-10.8772))),
            np.where(star.Campaign==7,0.3941/(1+np.exp(-.6052*(MES-11.7002))),
            np.where(star.Campaign==8,0.6669/(1+np.exp(-.5726*(MES-10.0560))),
            np.where(star.Campaign==10,0.5572/(1+np.exp(-.6469*(MES-10.0056))),
            np.where(star.Campaign==11,0.2171/(1+np.exp(-.4759*(MES-12.3882))),
            np.where(star.Campaign==12,0.6192/(1+np.exp(-.7341*(MES-10.6272))),
            np.where(star.Campaign==13,0.6853/(1+np.exp(-.5698*(MES-11.3878))),
            np.where(star.Campaign==14,0.7505/(1+np.exp(-.6596*(MES-10.9776))),
            np.where(star.Campaign==15,0.6067/(1+np.exp(-.6480*(MES-10.4673))),
            np.where(star.Campaign==16,0.6809/(1+np.exp(-.7256*(MES-10.5453))),
            np.where(star.Campaign==17,0.5848/(1+np.exp(-.6633*(MES-10.3635))),
            np.where(star.Campaign==18,0.6116/(1+np.exp(-.4676*(MES-11.5783))),
            0)))))))))))))))))
            
            fdet=np.where(period>26,fdet*0.4635/0.6619,fdet)
        
    return(fdet*probabilty_window)
	
def impactDraw(n,star,inclination,average_mutInc,per,node,ecc,arg_peri,plNum):
    # if(plNum==1):
    #     mutual_inclination=0
    # else:
	mutual_inclination=rrayleigh(n,average_mutInc)*(2*np.pi)/360

	###Draw argument of periapsis
	w=np.arcsin(np.random.uniform(size=n))

	####Calculate Planet Inclination
	inclination_m=np.cos(inclination)*np.cos(mutual_inclination)+np.sin(inclination)*np.sin(mutual_inclination)*np.cos(w)
    
	###Planet Semi-major axis
	semiMajor=(per**2*star.Mass*gravConstant/(4*np.pi**2))**(1/3)

	return(abs(semiMajor/star.Rad*(1-ecc**2)/(1+ecc*np.sin(arg_peri))*inclination_m))


def SimulateDetection(params):
    star=params.stellar
    multLen=np.max(star.m)+1

    star=star.reset_index(drop=True)
    number_stars=len(star)

    gravConstant=6.674e-11*1.98855e30*(86400)**2/(6.957e8)**3
    
    if (params.prettyPrint) & (params.mission=="K2"):
        print("Applying population completeness corrections for K2")
    if (params.prettyPrint) & (params.mission=="Kepler"):
        print("Applying population completeness corrections for Kepler")    
    
    ###Draw Inclination Parameters
    if params.indep==False:
        system_inclination=np.arccos(np.random.uniform(size=number_stars))
    
    ### mutual incliation parameter
    mutIncParam=0    
    system_average_mutInc=rrayleigh(number_stars,mutIncParam)
    intial_w=1#np.arcsin(np.random.uniform(size=number_stars))
    
    star['dummyVariable']=0
    #print(len(star))
    for m in range(1,multLen):
        if params.indep==True:
            system_inclination=np.arccos(np.random.uniform(size=number_stars))
            #system_average_mutInc=180
            
        exec(
"semiMajor_m{}=(star.per_m{}**2*star.Mass*gravConstant/(4*np.pi**2))**(1/3);"
"eccentricity_m{}=0;"
"argPericenter_m{}=0;"
"star['impact_m{}']=impactDraw(number_stars,star,system_inclination,system_average_mutInc,star.per_m{},intial_w,eccentricity_m{},argPericenter_m{},{}); "
"star.impact_m{}=np.where(star.impact_m{}>.9,2,star.impact_m{});"
"tran_chord_m{}=(1-star.impact_m{}**2)**.5;"
"star['tran_duration_m{}']=star.Rad*star.per_m{}/(semiMajor_m{}*np.pi)*tran_chord_m{}*(1-eccentricity_m{}**2)**.5/(1+eccentricity_m{}*np.sin(argPericenter_m{}))*24;"
"star.tran_duration_m{}=np.where(np.isnan(star.tran_duration_m{}),0,star.tran_duration_m{});"
"star.tran_duration_m{}=np.where(star.rad_m{}==0,0,star.tran_duration_m{});"
"star['dummyVariable']=star.dummyVariable+star.tran_duration_m{}".format(m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m)
)

    star=star[star.dummyVariable>0]
#    print(len(star))
    
    if params.prettyPrint:
        number_stars2=len(star)
        print("- "+str(number_stars-number_stars2) + " non-transiting exoplanets removed.") 
        
    number_stars=len(star)
    star=star.reset_index(drop=True)
    star.dummyVariable=0
    
    random_detection_threshold=np.random.uniform(size=number_stars)
    if params.mission=="Kepler":
        cdx=np.array([1.5,2,2.5,3,3.5,4.5,5,6,7.5,9,10.5,12,12.5,15])
    if params.mission=="K2":
        cdx=np.array([1,1.5,2,2.5,3,4,5,6,7,8,9,10])   
    
    for m in range(1,multLen):
    ###Interpolate the CDPP values for each star
        if params.mission=="Kepler":
            if params.fast==False:
                comCode='''
cdppCat=pd.cut(star.tran_duration_m{},[0, 1.5,2,2.5,3,3.5,4.5,5,6,7.5,9,10.5,12,12.5,15,np.inf],labels=["low",'med1','med2','med3','med4','med5','med6','med7','med8','med9','med10','med11','med12','med13','high'])
cdpp_m{}=np.zeros(number_stars);
cdpp_m{}=np.where(cdppCat=="low",star.CDPP1_5,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med1",star.tran_duration_m{}*star.slope1+star.int1,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med2",star.tran_duration_m{}*star.slope2+star.int2,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med3",star.tran_duration_m{}*star.slope3+star.int3,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med4",star.tran_duration_m{}*star.slope4+star.int4,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med5",star.tran_duration_m{}*star.slope5+star.int5,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med6",star.tran_duration_m{}*star.slope6+star.int6,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med7",star.tran_duration_m{}*star.slope7+star.int7,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med8",star.tran_duration_m{}*star.slope8+star.int8,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med9",star.tran_duration_m{}*star.slope9+star.int9,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med10",star.tran_duration_m{}*star.slope10+star.int10,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med11",star.tran_duration_m{}*star.slope11+star.int11,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med12",star.tran_duration_m{}*star.slope12+star.int12,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med13",star.tran_duration_m{}*star.slope13+star.int13,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="high",star.CDPP15,cdpp_m{})'''.format(m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m)
            
            else:
                comCode='''
cdpp_m{}=np.zeros(number_stars);
cdpp_m{}=star.CDPP7_5*(star.tran_duration_m{}/7.5)**(-.5)
'''.format(m,m,m)
                


        if params.mission=="K2":
            if params.fast==False:
                comCode='''
cdppCat=pd.cut(star.tran_duration_m{},[0, 1,1.5,2,2.5,3,4,5,6,7,8,9,10,np.inf],labels=["low",'med1','med2','med3','med4','med5','med6','med7','med8','med9','med10','med11','high'])
cdpp_m{}=np.zeros(number_stars);
cdpp_m{}=np.where(cdppCat=="low",star.CDPP1,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med1",star.tran_duration_m{}*star.slope1+star.int1,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med2",star.tran_duration_m{}*star.slope2+star.int2,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med3",star.tran_duration_m{}*star.slope3+star.int3,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med4",star.tran_duration_m{}*star.slope4+star.int4,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med5",star.tran_duration_m{}*star.slope5+star.int5,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med6",star.tran_duration_m{}*star.slope6+star.int6,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med7",star.tran_duration_m{}*star.slope7+star.int7,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med8",star.tran_duration_m{}*star.slope8+star.int8,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med9",star.tran_duration_m{}*star.slope9+star.int9,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med10",star.tran_duration_m{}*star.slope10+star.int10,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="med11",star.tran_duration_m{}*star.slope11+star.int11,cdpp_m{});
cdpp_m{}=np.where(cdppCat=="high",star.CDPP10,cdpp_m{})'''.format(m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m,m)


            else:
                comCode='''
cdpp_m{}=np.zeros(number_stars);
cdpp_m{}=star.CDPP8*(star.tran_duration_m{}/8)**(-.5)
'''.format(m,m,m)

            
        code=compile(comCode, '<string>', 'exec')
        exec(code)
	#Calculate the probabilty of detection for each planet
        exec(
"prob_detection_m{}=probability_detection(star, star.per_m{}, star.rad_m{}, cdpp_m{},{}, star.impact_m{},params.mission,params.campaignCom);"
"star.rad_m{}=np.where(np.where(np.isnan(prob_detection_m{}),0,prob_detection_m{})>random_detection_threshold,star.rad_m{},0);"
"star.dummyVariable=star.dummyVariable+star.rad_m{}".format(m,m,m,m,m,m,m,m,m,m,m)
)
        
    star=star[star.dummyVariable>0]
    
    if params.prettyPrint:
        number_stars2=len(star)
        print("- "+str(number_stars-number_stars2) + " undetected exoplanets removed.")
    
    number_stars=len(star)
    star=star.reset_index(drop=True)
    star.dummyVariable=0
    star["m"]=0
    
    
    if params.prettyPrint:
        print("Drawing planet radii from stellar uncertainties")
        
    #####Add Rad uncertainity
    star_error=star.Rad
    #abs(splitNorm(star.Rad,abs(star.U_Rad),abs(star.u_Rad),len(star)))
    #radTotal=np.zeros(number_stars)

  
    for m in range(1,multLen):
        exec(
"star['dummyVariable']=star['dummyVariable']+star.rad_m{};"
"star['m']=star['m']+np.where(star.rad_m{}>0,1,0)".format(m,m,m,m,m,m,m,m,m,m,m,m,m,m,m)
)


    star=star[star.dummyVariable>0]
    if params.prettyPrint:
        number_stars2=len(star)
        print("- "+str(number_stars-number_stars2) + " planets now exceed thresholds and have been removed.")

    star=star.reset_index(drop=True)
    number_stars=len(star)
    del star['dummyVariable']
    
    return(star)
