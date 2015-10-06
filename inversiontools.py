import os,sys
import numpy as np
import pandas as pan
from copy import deepcopy
import random
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
font = {'family' : 'serif',
        'weight' : 'medium',
        'size'   : 22}

plt.rc('font', **font)
plt.rcParams['lines.linewidth']=2.0
from itertools import groupby
import operator
import copy
import pickle
from itertools import cycle

def std_atm(z):
    ########################################################################
    # Program that defines temperature, pressure, and density profiles for
    # a US Standard Atmosphere (1976) given a height
    # in meters [AMSL] up to 47,000 m
    # Returns:
    # T = temperature [K]
    # P = pressure [Pa]
    # d = density [kg/m^3]
    # 
    ########################################################################

    #constants

    g = 9.80665  #gravitational constant [m/s^2]
    M = 0.028964 #molar density of air [kg/mol]
    R = 8.31432  #universal gas constant [N*m/mol/K]

    #Define breakpoints for US Standard atmosphere, with associated altitude,
    #pressure, density, temperature, and dry adiabatic lapse rates

    alt = [0, 11.000, 20.000, 32.000, 47.000]
    press = [101325,22632.1, 5474.89,868.019,110.906]
    dense = [1.225, 0.36391, 0.08803, 0.01322, 0.00143]
    temp = [288.15, 216.65, 216.65, 228.65, 270.65]
    lapse = [-0.0065, 0, 0.001, 0.0028, 0]

    #fisrt test to make sure no altitudes exceed the maximum

    if z > max(alt):
        raise ValueError('Sorry, all altitudes must be below %d' %max(alt))

    #start by determining temprature through linear interpolation

    T = np.interp(z,alt,temp)

    #now determine Pressure and density using different functions based on
    #atmospheric region
    
    if alt[0] <= z <= alt[1]:
        P = press[0]*(temp[0]/T)**(g*M/(R*lapse[0]))
        d = dense[0]*(T/temp[0])**((-g*M/(R*lapse[0]))-1)
    elif alt[1] < z <= alt[2]:
        P = press[1]*np.exp(-g*M*(z-alt[1])/(R*temp[1]))
        d = dense[1]*np.exp(-g*M*(z-alt[1])/(R*temp[1]))
    elif alt[2] < z <= alt[3]:
        P = press[2]*(temp[2]/T)**(g*M/(R*lapse[2]))
        d = dense[2]*(T/temp[2])**((-g*M/(R*lapse[2]))-1)
    elif alt[3] < z <= alt[4]:
        P = press[3]*(temp[3]/T)**(g*M/(R*lapse[3]))
        d = dense[3]*(T/temp[3])**((-g*M(R*lapse[3]))-1)

    return T,P,d


def molecular(z,wave=532.0):
    """
    Function for generating molecular scattering and extinction coefficients based
    on an altitude and a laser wavelength.  Ozone absorption is ignored.

    Inputs:

    z = altitude [m]
    wave = wavelength [nm]

    Outputs:

    beta = backscatter coefficients [1/km*sr]
    sigma = extinction coefficients [1/km]
    """

    #calculate temperature and pressure profiles using US Standard atmosphere

    [T,P,d] = std_atm(z)

    #determine backscatter and extinction coefficients as per Rayleigh scattering
    #equations and constats from Kovalev pp 33-36

    T_s = 288.15  #[K] reference temperature
    P_s = 101325.0  #[Pa] reference pressure
    N_s = 2.547e25  #[1/km^3]  reference number concentration
    gamma = 0.0071 #[unitless] depolarization factor (for narrow band fillter, from Young[1980])

    #calculate reference index of refraction using a polynomial approximation

    nu = 1000/wave #frequency in 1/um
    m_s = 1 + 1e-8*(8342.13+(2406030/(130-nu**2))+(15997/(38.9-nu**2)))

    #now calculate index of refraction at altitude as function of temperature and pressure

    m = 1+(m_s-1)*((1+0.00367*T_s)/(1+0.00367*T))*(P/P_s)

    #convert air mass density to number density

    N_a = 6.02214e23 #Avogadro's number [#/mol]
    M_air = 0.028964 #molar density of air [kg/mol]

    N = N_a*d/M_air


    #without absorption, extinction is equal to total scattering

    sigma = (8*np.pi**3*(m**2-1)**2*N/(3*N_s**2*(wave*1e-9)**4))*((6+3*gamma)/(6-7*gamma))* \
            (P/P_s)*(T_s/T)

    sigma=sigma*1000.0 #convert to 1/km
    #For Rayleigh scattering the extinction to backscatter ratio is 8*pi/3

    beta = 3.0*sigma/(8.0*np.pi)

    return T,P,d,beta,sigma

def molprof(z,wave=532.0,T0=1.0,E0=1.0,C=1.0):
    """
    Function for generating a theoretical profile of normalized attenuated
    backscatter.  In other words, this provides
    an array that can be multiplied by lidar output power and system constant
    to provide a lidar response profile

    Inputs:

    z = an array of altitudes [m]
    wave = lidar wavelength [nm]
    T0 = transmissivity to the first z value [unitless], defaults to 1.0 (perfect transmition)
    E0 = output laser pulse energy in mJ, defaults to 1.0 (normalized placeholder value)
    C = calibrated conversion constant between energy at the apertrue and detected counts/us [counts/us/mJ]
    
    Outputs:
    P_out = a pandas dataframe with columns representing altitudes and the following index values
    vals = vector of normalized attenuated backscatter [unitless]
    beta_R = vector of molecular backscatter coefficients [1/m*sr]
    sigma_R = vector of molecular extinction coefficients [1/m]
    beta_p = vector of particulate backscatter coefficients (assumed to be zero here)
    sigma_p = vector of particulate extinction coefficients (assumed to be zero here)
    beta_t = vector of total backscatter coefficients [1/m*sr]
    sigma_t = vector of total extinction coefficients [1/m]
    Temp = vector of temperatures [K]
    Press = vector of pressures [Pa]
    Density = vector of atmospheric density [kg/m^3]
    """
    
    T = pan.Series(0,index=z,dtype=float)
    P = pan.Series(0,index=z,dtype=float)
    d = pan.Series(0,index=z,dtype=float)
    beta_R = pan.Series(0,index=z,dtype=float)
    sigma_R = pan.Series(0,index=z,dtype=float)
    P_z = pan.Series(0,index=z,dtype=float)
    beta_p = pan.Series(0,index=z,dtype=float)
    sigma_p = pan.Series(0,index=z,dtype=float)

    for alt in z:
        [T.loc[alt],P.loc[alt],d.loc[alt],beta_R.loc[alt],sigma_R.loc[alt]] = molecular(alt, wave)
    
    T_total = T0
    oldalt=z[0]
    if oldalt==0.0:
        P_z.loc[oldalt] = E0*C*beta_R.iloc[0]*T_total
    else:
       P_z.loc[oldalt] = E0*C*oldalt**-2*beta_R.iloc[0]*T_total
    
    for alt in z[1:]:
        T_step = np.exp(-2*sigma_R.loc[alt]*(alt-oldalt))
        T_total = T_total*T_step
        P_z.loc[alt] = (E0*C*beta_R.loc[alt]*T_total)*alt**-2.0
        oldalt=alt
        
    rsq = [val*alt**2.0 for val,alt in zip(P_z,z)]
    NRB = [val/E0 for val in rsq]
        
    
    beta_t = beta_R+beta_p
    sigma_t = sigma_R+sigma_p
    keys = ['vals','rsq','NRB','Temp','Press','Density','beta_R','sigma_R','beta_p','sigma_p','beta_t','sigma_t']
    vals = [P_z,rsq,NRB,T,P,d,beta_R,sigma_R,beta_p,sigma_p,beta_t,sigma_t]
    P_out = pan.DataFrame.from_dict(dict(zip(keys,vals)))
    
    return P_out
    

def calcprof(betaprof,sigmaprof,E0=1.0,C = 1.0,T_0=1.0):
    """
    Function for generating a profile of attenuated
    backscatter based on input backscatter and extinction coefficients.

    Inputs:

    betaprof [m-1*sr-1] = profile of backscatter coefficients
    sigmaprof [m-1] = profile of extinction coefficients
    E0 [mJ] = measured laser output, defaults to 1.0 for normalized profile
    C = conversion fro mJ to counts/us, defaults to 1.0 
    T_0 = transmissivity to the first z value [unitless], defaults to 1.0 (perfect transmition)

    Outputs:
    P_z = a profile of at-aperture signal values in mJ as a function of altitudw
    """
    
    alts=betaprof.index
    P_z=pan.Series(index=alts)    
    T_total = T_0
    oldz=alts[0]
    P_z.ix[oldz] = E0*C*oldz**-2*betaprof.ix[oldz]*T_total**2
    for z in alts[1:]:
        T_step = np.exp(-sigmaprof.ix[z]*(z-oldz))
        T_total = T_total*T_step
        P_z.ix[z] = E0*C*z**-2.0*betaprof.ix[z]*T_total**2
        oldz=z
    
    return P_z

def addlayer(P_in, beta, lrat, bottom=None, top=None, inplace=True):
    """
    Function that adds a layer of known backscatter coefficient
    and lidar ratio onto an existing lidar response profile

    Inputs:

    P_in = a pandas dataframe object containing the input profile with indexed values defined in mplprof

    beta = a pandas series object showing backscatter coefficients [1/m*sr] with altitude 
    lrat = lidar (extinction to backscater) ratio of the added layer [1/sr]

    Outputs:
    P_out = dataframe object like P_in with layer added

    """    
    if inplace:
        P_out=P_in
    else:
        P_out = deepcopy(P_in)
    if type(beta)!=float:
        z_in = np.array(beta.index.values, dtype='float64')
    else:
        if bottom is not None and top is not None:
            z_in = np.array(P_out.index[(P_out.index>=bottom)&(P_out.index<=top)],dtype='float64')
        else:
            print "You must define a layer top and bottom"
            return
            
    z_min = min(z_in)
    z_max = max(z_in)
    
    z_old = P_out.index[0]
    for z in P_out.index:
        if z < z_min:
            z_old=z
            oldval=P_out.loc[z,'vals']
            oldbeta=P_out.loc[z,'beta_t']
        elif z <= z_max:
            if type(beta)!=float:
                P_out.loc[z,'beta_p'] += np.interp(z,z_in,beta)
            else:
                P_out.loc[z,'beta_p'] += beta
            P_out.loc[z,'sigma_p'] = P_out.loc[z,'beta_p']*lrat
            P_out.loc[z,'beta_t'] = P_out.loc[z,'beta_p'] + P_out.loc[z,'beta_R']
            P_out.loc[z,'sigma_t'] = P_out.loc[z,'sigma_p'] + P_out.loc[z,'sigma_R']
            T_step = np.exp(-2.0*P_out.loc[z,'sigma_t']*(z-z_old))
            P_out.loc[z,'vals'] = oldval*(P_out.loc[z,'beta_t']/oldbeta)*(z_old**2/z**2)*T_step        
            z_old=z
            oldval=P_out.loc[z,'vals']
            oldbeta=P_out.loc[z,'beta_t']
        else:
            T_step = np.exp(-2.0*P_out.loc[z,'sigma_t']*(z-z_old))
            P_out.loc[z,'vals'] = oldval*(P_out.loc[z,'beta_t']/oldbeta)*(z_old**2/z**2)*T_step        
            z_old=z
            oldval=P_out.loc[z,'vals']
            oldbeta=P_out.loc[z,'beta_t']
    
    return P_out

def backandnoise(P_in,background = None,stdev = None,inplace=True):
    """
    Adds gaussian random noise and background signal to any profile
    
    Inputs:
    P_in = a pandas series with altitude index
    background = a float defining the backgorund signal to be applied(defaults to 0)
    stdev = the standard deviation of the gaussian noise component, can be a single
                value or a Pandas series
    
    if no standard deviation is defined, the noise added is standard
    shot noise - poisson distribution approximated by a gaussian with
    std = sqrt(signal)
    
    Outputs:
    P_out = a copy of P_in with background signal and noise applied
    """
    
    if inplace:
        P_out = P_in
    else:
        P_out=deepcopy(P_in)
    
    if background is not None:
        P_out['vals']=P_out['vals']+background
    
    if stdev is not None:
        if type(stdev)==pan.core.series.Series:
            P_out['vals'] = [v+random.gauss(0.0,s) for v,s in zip(P_out['vals'],stdev.values)]
        elif stdev:
            P_out['vals'] = [v+random.gauss(0.0,stdev) for v in P_out['vals']]
        elif stdev=='shot':
            P_out['vals'] = [v+random.gauss(background,s) for v,s in zip(P_out['vals'],np.sqrt(P_out['vals']))]
       
    return P_out

def background_subtract(P_in,back_avg=None,z_min=None,inplace=True):
    #subtracts background from signal, without background
    #takes advantage of inverse square law to calculate background signal for
    #a profile

    #start by selecting a region of sufficient altitude that the r squared law
    #will make the ratio of signal to background sufficiently low
    
    if inplace:
        P_out=P_in
    else:
        P_out = deepcopy(P_in)
    
    if back_avg is not None:
        P_out['vals']=P_out['vals']-back_avg
    else:    
        #select data from altitudes higher than z_min and muliply full signal by
        #range squared, if this gives less than 500 values, take uppermost 500
        if z_min is not None:
            z=P_out.index[P_out.index>=z_min]
        else:
            z=P_out.index[-100:]
    
        temprsq = P_out['vals'].loc[z]*z**2
    
        #since background is constant while signal is reduced as z^2, first
        #coefficient is equal to background
        
        coeffs = np.polyfit(z,temprsq,2,full=False)
        background = coeffs[0]    
        P_out['vals'] = P_out['vals']-background
        
    return P_out

def calcNRB(P_in,E0=1.0,background=None,inplace=True):
    P_out=background_subtract(P_in,back_avg=background,inplace=inplace)
    P_out['rsq'] = P_out['vals']*(P_in.index.values)**2.0
    P_out['NRB'] = P_out['rsq']/E0
    
    return P_out

def calc_slope(prof, winsize = 10):
    import pandas as pan
    import numpy as np
    """
    Calculates slope of data for a single profile using a smoothing window of
    predetermined size
    
    inputs:
    prof:  a pandas series where index is altitude
    n:  number of consecutive values to average
    
    output:
    slopeout: output series,same size as input,with profile slopes
    """
    data = prof.values
    altrange = np.asarray(prof.index.values,dtype='float')
    
    #Step 1: pad dataset to allow averaging
    
    leftpad = np.int(np.floor(winsize/2))
    rightpad = winsize-leftpad
      
    #Step 2: Calculate a linear fit to the data in the window
    
    slopes = np.empty(len(data)-winsize)
    for n in range(len(slopes)):       
        x = altrange[n:n+winsize]
        y = data[n:n+winsize]
        
        coeffs = np.polyfit(x,y,1,full=False)
        slopes[n] = coeffs[0]
        
    
    slopes = np.pad(slopes,(leftpad,rightpad),'edge')
    
    slope_out = pan.Series(slopes, index=altrange)
    
    
    return slope_out

def profgen(z,**kwargs):
    """
    generates a pandas dataframe containing simulated profiles of backscatter and extinction coefficients as well as
    relevant atmospheric parameters (see molprof)    
    
    Inputs:
    z - range of altitudes to calculate profiles for [m]
    
    Kwargs:
    wave - wavelength of scattered light [nm] Defaults to 532nm.
    E0 - laser pulse output energy in uJ.  Defaults to 1.0
    C - calibration coefficient for converting emitted energy to detected counts/us.  Defaults to 1.0
    layers - list of dict items containing layer properties.  Key-value pairs are:
            beta_p - particulate backscatter coefficients.  Can be either a 
                        single value or a Pandas series of values with altitude index
            lrat - layer lidar ratio [1/sr]
            bot - altitude of layer bottom.  Only used if beta_p is a single value
            top - altitude of layer top.  Only used if beta_p is a single value
    background - background signal level in counts/us.  Default is []
    noise - std of noise to add to profile.  Can be either a single value or a pandas
                series of values with altitude index.  Default is 0.0    
    donorm - boolean to determine whether values are normalized.  Default is False
    
    Outputs:
    P_out - Pandas dataframe containing the following profiles
        vals = vector of attenuated backscatter [counts/us]
        norm = vector of normalized attenuated backscatter [unitless] (not ot be confused with NRB)
        rsq = range squared corrected profile [counts*km^2/us]
        NRB = vactor of normalized relative backscatter [counts*km^2/us/uJ]
        beta_R = vector of molecular backscatter coefficients [1/m*sr]
        sigma_R = vector of molecular extinction coefficients [1/m]
        beta_p = vector of particulate backscatter coefficients [1/m*sr]
        sigma_p = vector of particulate extinction coefficients [1/m]
        beta_t = vector of total backscatter coefficients [1/m*sr]
        sigma_t = vector of total extinction coefficients [1/m]       
        Temp = vector of temperatures [K]
        Press = vector of pressures [Pa]
        Density = vector of atmospheric density [kg/m^3]
    """

    wave = kwargs.get('wave',532.0)  #nm
    E0 = kwargs.get('E0',1.0)
    C = kwargs.get('C',1.0)
    layers = kwargs.get('layers',None)
    background = kwargs.get('background',0.0) #background signal level
    noise = kwargs.get('noise',0.0) #use defined noise level, defaults to 0
    donorm = kwargs.get('donorm',False) #if true, profile is normalized by first value
    
    #start by creating a purely molecular profile
    P_out = molprof(z,wave,E0,C)
    
    #add particulate layers, if any
    
    if layers is not None:
        for n in range(len(layers)):
            templayer=layers[n]
            if type(templayer['beta_p'])==float:
                P_out = addlayer(P_out,beta=templayer['beta_p'],lrat=templayer['lrat'],
                                     bottom=templayer['bot'],top=templayer['top'],
                                        inplace=True)
            else:
                P_out = addlayer(P_out,beta=templayer['beta_p'],lrat=templayer['lrat'],
                                        inplace=True)
    
    P0 = P_out['vals'].iloc[0]
    P_out=backandnoise(P_out,background=background,stdev=noise,inplace=True)
    
    P_out=calcNRB(P_out,E0=E0,background=background)
    if donorm:
        P_out['vals']=P_out['vals']/P0
            
    return P_out
    

def fernald(P_in, lrat, wave = 532.0, E = 1.0, calrange = None):
    """
    Inputs:
        P_in: a pandas series depicting a 1-D profile of NRB from an MPL class object. 
        Indices are altitudes.  Description can be found in MPLtools.py
        lrat:  a pre-defined lidar ratio
        wave:  lidar wavelength for Rayleight scattering calculations
        E:  laser output energy in mJ (default to 1.0 for normalized calculations)
        calrange: a minimum and maximum altitude for the calibration range
                    (assumed to be clear air)
        
    
    Outputs:
        beta_out: a pandas series containing backscatter coefficients in units of [1/m*sr]
        sigma_out: a pandas series containing extinction coefficients in units of [1/m]
    
    In the langauge of Fernald's original paper, the lidar equation is reduced to the form:
    P(Z) = E*C*Z^-2*[beta_R(Z)+beta_P(Z)]*T_R(Z)^2*T_P(Z)^2
    Where:
    P(Z)=Power at the receiver
    E=laser output power
    C=system constant
    beta_R & beta_P=backscatter ceofficients for Rayleigh and particulate scattering,respectively
    T_R & T_P = integrated extinction (optical depth) from Rayleigh and particulate scattering
    and  
    T_R = exp[-int_0-Z(sigma_R(z)dz)]
    In these terms, the Normalized Realtive Backscatter of the MPL
    in units of [counts*km^2/us/mJ] is equivalent to 
    NRB = P(Z)*Z^2/E = C*[beta_R(Z)+beta_P(Z)]*T_R(Z)^2*T_P(Z)^2
    
    """
    
    #First calculate Pure Rayleigh scattering coefficients
    
    altitudes = P_in.index.values
    
    Rayleigh_coeffs = molprof(altitudes,wave)
    Rayleigh_lrat = 8.0*np.pi/3.0
    
    #Start with initial calculation of system constant C
    
#    tau_R = 0
#    T_Rstar = P_in[-1]/P_in[0]
    
    beta_total = pan.Series(index=altitudes)
    
        
    if calrange is None:
        for alt in reversed(altitudes):
            if alt == altitudes[-1]:
                #at calibration altitude assume no aerosols
                beta_total.loc[alt] = Rayleigh_coeffs['beta_R'].loc[alt]
                oldalt = alt
            
            else:
                X1 = P_in.loc[oldalt]*E
                X0 = P_in.loc[alt]*E
                beta_R1 = Rayleigh_coeffs['beta_R'].loc[oldalt]
                beta_R0 = Rayleigh_coeffs['beta_R'].loc[alt]
                delta_Z = oldalt-alt
                A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                
                beta_total.loc[alt] = X0*np.exp(A)/((X1/beta_total[oldalt])+lrat* \
                np.abs(X1+X0*np.exp(A))*delta_Z)
                
                oldalt = alt
    else:
        #first assign Rayleigh coefficients to calibration range
        
        minalt = calrange[0]
        maxalt = calrange[1]
        calalts = [z for z in altitudes if minalt <= z <= maxalt]
        
        for alt in calalts:
            beta_total.loc[alt] = Rayleigh_coeffs['beta_R'].loc[alt]
        
        #next calculate coefficients for range below 
        if len(calalts)==0 or len(altitudes)==0:
            print "trouble!"
        if calalts[0] > altitudes[0]:
            oldalt = calalts[0]
            below = [z for z in altitudes if z <= minalt]
            
            for alt in reversed(below):
                if np.isnan(P_in.loc[alt]):
                    beta_total.loc[:alt]=np.nan
                    break
                else:
                    X1 = P_in.loc[oldalt]*E
                    X0 = P_in.loc[alt]*E
                    beta_R1 = Rayleigh_coeffs['beta_R'].loc[oldalt]
                    beta_R0 = Rayleigh_coeffs['beta_R'].loc[alt]
                    delta_Z = oldalt-alt
                    A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                    
                    beta_total.loc[alt] = X0*np.exp(A)/((X1/beta_total[oldalt])+lrat* \
                    np.abs(X1+X0*np.exp(A))*delta_Z)
                    
                    oldalt = alt
        
        #then  calculate for altitudes above maxalt
        if calalts[-1] < altitudes[-1]:
            oldalt = calalts[-1]
            above = [z for z in altitudes if z >= maxalt]
            
            for alt in above:
                if np.isnan(P_in.loc[alt]):
                    beta_total.loc[alt:]=np.nan
                    break
                else:
                    X1 = P_in.loc[oldalt]*E
                    X0 = P_in.loc[alt]*E
                    beta_R1 = Rayleigh_coeffs['beta_R'].loc[oldalt]
                    beta_R0 = Rayleigh_coeffs['beta_R'].loc[alt]
                    delta_Z = oldalt-alt
                    A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                    
                    beta_total.loc[alt] = X0*np.exp(-A)/((X1/beta_total[oldalt])-lrat* \
                    np.abs(X1+X0*np.exp(-A))*delta_Z)
                    
                    oldalt = alt
        
    sigma_total=beta_total*lrat
    
    return beta_total,sigma_total


def klett(P_in,lrat_in,r_m=None,k=1,wave=532.0):
    """
    Function that calculates backscatter and extinction coefficients based on
    variable lidar ratios using the Klett algorithm 
    
    Inputs:
    P_in = a pandas series with values of signal strength and altitude index
    r_m = the reference altitude - maximum altitude for which calculations are done
            and the point at which the extinction coefficeint is assumed ot be known
    lrat_in = a Pandas series with values of gross lidar ratio and altitude index
    sigma_m = the extinction coefficient at altitude r_m
    k = the power coefficient in the power law relationship bwetween backscatter
        and extinction (defaults to 1)
        
    Outputs:
    beta = pandas series of backscatter coefficients
    sigma = pandas series of extinction coefficients
    
    """

    
    P_trans = P_in+0.001-np.min(P_in.values)
    S=np.log(P_trans)
    lrat = 1.0/lrat_in  #Klett definition of lidar ratio is backscatter/extintion not the other way round

    altitudes = S.index.values
    if r_m is None:
        for a in altitudes[::-1]:
            if not np.isnan(P_in.ix[a]):
                r_m=a
                break
    
    sigma = pan.Series(index=altitudes,dtype='float')
    beta = pan.Series(index=altitudes,dtype='float')
    
    for alt in reversed(altitudes):
       if alt > r_m:
           sigma.loc[alt] = np.nan
           beta.loc[alt] = np.nan      
       else:
           S_new = S.loc[:alt]
           break
    
    newalts = S_new.index.values
    sigma_m=molprof(z=newalts,wave=wave)['sigma_R'].loc[r_m]
    for alt in reversed(newalts):
        if alt == newalts[-1]:            
            sigma.loc[alt] = sigma_m
            beta.loc[alt] = lrat.loc[alt]*sigma_m**k
            oldalt = alt
            Xint=0.0
        
        else:
            X1 = (lrat.loc[r_m]/lrat.loc[alt])**(1/k)
            X2 = np.exp((S_new.loc[alt]-S_new.loc[r_m])/k)
            Xint += X1*X2*(oldalt-alt)
            sigma.loc[alt] = X1*X2/(sigma_m**-1+(2/k)*Xint)
            
            beta.loc[alt] = lrat.loc[alt]*sigma.loc[alt]**k
            oldalt = alt
        
    return beta, sigma  
    
def klett2(P_in,lrat_in,**kwargs):
    """
    Function that calculates backscatter and extinction coefficients based on
    variable lidar ratios using the Klett algorithm .  Note: requires lrat_in not
    contain zero values except for at r_m altitude
    
    Inputs:
    P_in = a pandas series with values proportional to range squared corrected signal strength and altitude index
    r_m = the reference altitude - maximum altitude for which calculations are done
            and the point at which the backscatter coefficeint is assumed ot be known
    lrat_in = a Pandas series with values of particulate only lidar ratio and altitude index
    beta_m = the backscatter coefficient at altitude r_m
        
    Outputs:
    beta = pandas series of backscatter coefficients
    sigma = pandas series of extinction coefficients
    
    """
    r_m=kwargs.get('r_m',None)
    wave=kwargs.get('wave',532.0)
    
#    lrat=pan.Series(data=0.0,index=lrat_in.index)
#    
#    grouped=lrat_in.groupby(lrat_in)
#    
#    #Klett definition of lidar ratio is backscatter/extintion not the other way round
#    for name,group in grouped:
#        if name==0.0:
#            continue
#        else:
#            lrat.ix[group.index]=1.0/group
            

    lrat_R=3.0/(8.0*np.pi)
    altitudes = P_in.index.values
    lrat_new=lrat_in.replace(0.0,(8.0*np.pi)/3.0,inplace=False)
    lrat=1.0/lrat_new
    if r_m is None:
        for a in altitudes[::-1]:
            if not np.isnan(P_in.ix[a]):
                r_m=a
                break
    beta = pan.Series(index=altitudes,dtype='float')
    sigma = pan.Series(index=altitudes,dtype='float')
    for alt in reversed(altitudes):
       if alt > r_m:
           beta.loc[alt] = np.nan      
       else:
           P_new = P_in.loc[:alt]
           break
    
    newalts = P_new.index.values
    P_mol=molprof(z=newalts,wave=wave)
    beta_R=P_mol['beta_R']
    beta_m=beta_R.loc[r_m]
    
    sigma_R=P_mol['sigma_R']
    sigma_m=sigma_R.loc[r_m]
    #add arbitrary translation before transformation to avoid negative numbers in log value
    #because inversion is based on the differential of dS/dr, this changes nothing in solution    
    if np.min(P_new.values)<=0:
        P_trans=P_new+(1e-2)-np.min(P_new.values)
    else:
        P_trans=P_new
        
    S=np.log(P_trans).fillna(method='pad')
    
    for alt in reversed(newalts):
        if alt == r_m: 
            beta.ix[alt]=beta_m
            sigma.ix[alt]=sigma_m
            S_m=S.ix[alt]
            beta_int1=0
            beta_int2=0
            S_int=0
            oldalt = alt      
        else:
            delta_r=oldalt-alt 
            if oldalt==r_m:
                #in this case, use simple reimann integration
                beta_int1 += (2.0/lrat_R)*(beta_R.ix[alt])*delta_r
                beta_int2 += 2.0*(beta_R.ix[alt]/lrat.ix[alt])*delta_r
                delta_Sprime = S.ix[alt]-S_m+beta_int1-beta_int2
                S_int+=2.0*(np.exp(delta_Sprime)/lrat.ix[alt])*delta_r
            else:
                #otherwise use trapezoidal integration
                beta_int1+=0.5*((2.0/lrat_R)*(beta_R.ix[alt]+beta_R.ix[oldalt]))*delta_r   
                beta_int2+=0.5*(2.0*((beta_R.ix[alt]/lrat.ix[alt])+(beta_R.ix[oldalt]/lrat.ix[oldalt])))*delta_r            
                delta_Sprime=S.ix[alt]-S_m+beta_int1-beta_int2                
                S_int+=0.5*(2.0*((np.exp(delta_Sprime)/lrat.ix[alt])+(np.exp(old_delta_Sprime)/lrat.ix[oldalt])))*delta_r
                
            beta.ix[alt]=np.exp(delta_Sprime)/((1.0/beta_m)+S_int)
                
#            if beta.ix[alt]<=beta_R.ix[alt]:
#                beta_p=0.0
#            else:
#                beta_p=beta.ix[alt]-beta_R.ix[alt]
#            sigma.ix[alt]=(beta_p/lrat.ix[alt])+beta_R.ix[alt]/lrat_R
            sigma.ix[alt]=beta.ix[alt]/lrat.ix[alt]
            oldalt = alt
            old_delta_Sprime=delta_Sprime
    return beta,sigma

#def iterative_klett(P_in,lrat_p,**kwargs):
#    wave=kwargs.get('wave',532.0)
#    k=kwargs.get('k',1.0)
#    r_m=kwargs.get('r_m',None)
#    maxiter=kwargs.get('maxiter',20)
#    deltathresh=kwargs.get('deltathresh',0.1)
#    lrat_min=kwargs.get('lrat_min',0.0)
#    lrat_max=kwargs.get('lrat_max',100.0)
#    verbose=kwargs.get('verbose',False)
#    
#    alts=P_in.index
#    Pmol=molprof(z=alts,wave=wave)
#    beta_m=Pmol.beta_R
#    lrat_m=8.0*np.pi/3.0
#    
#    lrat_old=lrat_p.replace(0.0,lrat_m,inplace=False)
#    back_old,ext_old=klett(P_in=P_in,lrat_in=lrat_old,r_m=r_m,k=k,wave=wave)
#    n=0
#    while True:     
#        n+=1
#        
#        if verbose:
#            print "Iteration #{0}".format(n)
#        lrat_new=(lrat_m*beta_m+lrat_p*(back_old-beta_m)).div(back_old)
#        for i in lrat_new.index:
#            if lrat_new.ix[i]<lrat_min:
#                lrat_new.ix[i]=lrat_old.ix[i]
#            elif lrat_new.ix[i]>lrat_max:
#                lrat_new.ix[i]=lrat_max
#            elif lrat_old.ix[i]==lrat_m:
#                lrat_new.ix[i]=lrat_m
#        
#        back_new,ext_new=klett(P_in=P_in,lrat_in=lrat_new,r_m=r_m,k=k,wave=wave)
#        
#        delta_ext=abs(100.0*(ext_new-ext_old).div(ext_old))
#        back_old=back_new
#        ext_old=ext_new
#        lrat_old=lrat_new
#        
#        if n>=maxiter:
#            if verbose:
#                print 'Failed to converge after {0} Iterations'.format(n)
#            back_out,ext_out=klett(P_in=P_in,lrat_in=lrat_old,r_m=r_m,k=k,wave=wave)
#            lrat_out=lrat_old
#            break
#        if delta_ext.max() <= deltathresh:
#            lrat_out=lrat_new
#            for i in lrat_new.index:
#                if lrat_new.ix[i]<lrat_min:
#                    lrat_out.ix[i]=lrat_old.ix[i]
#                elif lrat_new.ix[i]>lrat_max:
#                    lrat_out.ix[i]=lrat_max
#                elif lrat_old.ix[i]==lrat_m:
#                    lrat_out.ix[i]=lrat_m
#            back_out,ext_out=klett(P_in=P_in,lrat_in=lrat_out,r_m=r_m,k=k,wave=wave)
#                    
#            break
#    
#    return back_out,ext_out,lrat_out

#def invert_profile(profin,lratin,**kwargs):
#    method=kwargs.get('method','klett2')
#    refalt=kwargs.get('refalt',profin.index[-1])
#    backscatter=pan.Series(np.nan,index=profin.index)
#    extinction=pan.Series(np.nan,index=profin.index)
#    
#    mollayers=lratin[lratin==0.0]
#    moledges=[]
#    altstep=lratin.index[1]-lratin.index[0]
##    molcount=[int(round((x-mollayers.index[0])/altstep)) for x in mollayers.index]
#    
#    for key,alt in groupby(enumerate(mollayers.index),lambda (i,x):i-int(round((x-mollayers.index[0])/altstep))):
#        temprange=map(operator.itemgetter(1),alt)
##        tempalts=[x*altstep+mollayers.index[0] for x in temprange]
#        moledges.append((temprange[0],temprange[-1]))
#    
#    alt=profin.index[-1]
#    if len(moledges)==0:
#        if methos=='klett2':
#            backscatter,extinction=klett2(profin,lratin,r_m=refalt)
#        return backscatter,extinction
#    else:
#        tempedges=moledges.pop()        
#        while True:        
#            if tempedges[0]<alt<=tempedges[1]: 
#                layeralts=profin.ix[tempedges[0]:tempedges[1]].index.values
#                tempmol=molprof(z=layeralts)
#                backscatter.ix[layeralts]=tempmol['beta_R'].values
#                extinction.ix[layeralts]=tempmol['sigma_R'].values
#                alt=layeralts[0]
#            else:
#                try:
#                    tempedges=moledges.pop()
#                    layeralts=profin.ix[tempedges[1]:alt].index[1:]
#                    alt=tempedges[1]
#                except IndexError:
#                    layeralts=profin.ix[:alt].index
#                    alt=layeralts[0]
#                    
#                layerprof=profin.ix[layeralts]
#                layerlrat=lratin.ix[layeralts]
#    
#                if method=='klett2':
#                    tempback,tempext=klett2(layerprof,layerlrat,r_m=layeralts[-1])
#        
#                backscatter.ix[layeralts]=tempback.values
#                extinction.ix[layeralts]=tempext.values
#                
#            if alt-profin.index[0]<=altstep:
#                break
#        
#        return backscatter,extinction

#def kovalev(P_in,lrat_in,**kwargs):
#    wave=kwargs.get('wave',532.0)
#    threshval=kwargs.get('threshval',0.10)
#    iterthresh=kwargs.get('iterthresh',20)
#    divergethresh=kwargs.get('divergethresh',5)
#    verbose=kwargs.get('verbose',False)
#    
#    lrat_m=8.0*np.pi/3.0
#    alts=P_in.index
#    tempmol = molprof(z=alts,wave=wave)
#    sigma_m = tempmol['sigma_R']
#    
#    n=0
#    divergeflag=0
#    P_old=P_in
#    P_vals=pan.DataFrame(index=alts)
#    delta_vals=pan.DataFrame(index=alts)
#    sigma_vals=pan.DataFrame(index=alts)
#    while True:
#        #step 1: calculate I_r = integral of P_in from r_0 - r
#        oldalt=alts[0]
#        I_r = pan.Series(index=alts)
#        I_r.ix[oldalt]=0.0
#        for newalt in P_old.index[1:]:
#            I_r.ix[newalt]=I_r.ix[oldalt]+0.5*(P_old.ix[newalt]+P_old.ix[oldalt])*(newalt-oldalt)                        
#            oldalt=newalt
#        
#        I_max=I_r.iloc[-1]    
#        #step 2: assume R_b=0 and use eqn. 22 to find gamma_r    
#        gamma_r = 1.0-(2.0*I_max/(P_old.div(sigma_m)+2.0*I_r))
#    
#        #step 3: calculcate gamma_min from gamma_r and use eqn. 23 to calculate sigma_p
#        
#        gamma_min = gamma_r.min()
#        gamma_idxmin = gamma_r.idxmin() 
#        
#        denom=(I_max/(1-gamma_min))-I_r
#        sigma_p = 0.5*P_old.div(denom) - sigma_m
#        
#        #step 4: calculate Y_r using equation 14
#        
#        Y_r_num = sigma_m+sigma_p
#        Y_r_denom = sigma_m + lrat_in.div(lrat_m)*sigma_p
#        
#        Y_r = Y_r_num.div(Y_r_denom)
#        #step 5: nmormalize P_in by Y_r using eqn. 15
#        
#        P_new = P_in*Y_r
#        #step 6: recalculate fom step 1 until no difrferences between steps        
#        
#        delta_P = abs(100.0*(P_new-P_old).div(P_old))  #delta in percentage units
#        
#        
#        if delta_P.max() <= threshval:
#            extinction=sigma_p+sigma_m
#            backscatter=extinction.div(lrat_in)
#            flag=0
#            break
#        
#        if n>=1 and delta_P.mean()>= delta_P_old.mean():
#            divergeflag+=1
#        
#            if divergeflag>=divergethresh:
#                if verbose:
#                    print "Kovalev diverging error"
#                extinction=pan.Series(data=np.nan,index=P_in.index)
#                backscatter=pan.Series(data=np.nan,index=P_in.index)
#                flag=2
#                break
#            
#        P_old=P_new
#        delta_P_old=delta_P
#        n+=1  
#        
#        P_vals[n]=P_old
#        delta_vals[n]=delta_P_old
#        sigma_vals[n]=sigma_p
#        
#        if n>=iterthresh:
#            if verbose:
#                print 'Failed to converge after {0} iterations'.format(n)
#            extinction=sigma_p+sigma_m
#            backscatter=extinction.div(lrat_in)
#            flag=1
#            break
#
#    P_vals.plot()
#    delta_vals.plot()
#    sigma_vals.plot()    
#    return backscatter,extinction,flag

def lrat_tester_noise(z,**kwargs):
    background=kwargs.get('background',0.0)
    betalist=kwargs.get('betalist',[1e-6])
    altlist=kwargs.get('altlist',[10])
    lratlist=kwargs.get('lratlist',[30])
    layerwidth=kwargs.get('layerwidth',1.0)
    noise=kwargs.get('noise',[0.0])
    wave=kwargs.get('wave',532.0)
    E0=kwargs.get('E0',1.0)
    lrat_testrange=kwargs.get('lrat_testrange',np.arange(-.50,1.55,.05))
    r_m=kwargs.get('r_m',None)
    k=kwargs.get('k',1.0)
    verbose=kwargs.get('verbose',True)
    layeronly=kwargs.get('layeronly',True)
    plotall=kwargs.get('plotall',True)
    figlib=kwargs.get('figlib','./Figures')
    saveall=kwargs.get('saveall',False)
    proclib=kwargs.get('proclib','./Processed')
    
    lrat_m = 8.0*np.pi/3.0
                                                           
    repnum=1 
    indexlist=[]
    namelist=[]
    for vals,names in zip([noise,lrat_testrange],['Noise','Lidar Ratio Error']):
        if len(vals)>1:
            indexlist.append(vals)
            namelist.append(names)
    
    if not indexlist:
        lratcalc=pan.Series(index=z)
        betacalc=pan.Series(index=z)
        sigmacalc=pan.Series(index=z)
        deltalrat=pan.Series(index=z)
        deltabeta=pan.Series(index=z)
        deltasigma=pan.Series(index=z)
    else:
        index=pan.MultiIndex.from_product(indexlist,names=namelist)
        lratcalc=pan.DataFrame(index=z,columns=index)
        betacalc=pan.DataFrame(index=z,columns=index)
        sigmacalc=pan.DataFrame(index=z,columns=index)
        deltalrat=pan.DataFrame(index=z,columns=index)
        deltabeta=pan.DataFrame(index=z,columns=index)
        deltasigma=pan.DataFrame(index=z,columns=index)
    totreps=len(noise)*len(lrat_testrange)  
    repnum=1 
                                         
    for nval in noise:
        testlayer=[]
        for beta,lrat,alt in zip(betalist,lratlist,altlist):
            templayer=pan.Series(data=beta,index=[alt,alt+layerwidth])
            templayer={'beta_p':templayer,'lrat':lrat}
            testlayer.append(templayer)
        P_0=profgen(z,layers=testlayer,background=background,noise=nval,wave=wave,E0=E0)
        P_temp=P_0['NRB']
        lratprof=P_0['sigma_t'].div(P_0['beta_t'])
        lrat0=pan.Series(data=lratprof,index=z)
        beta0=pan.Series(data=P_0['beta_t'],index=z)
        sigma0=pan.Series(data=P_0['sigma_t'],index=z)
                    
        for lrat_multiplier in lrat_testrange: 
            if verbose:
                print 'Processing step {0} out of {1}'.format(repnum,totreps)
            repnum+=1
            if layeronly:
                lrattemp=pan.Series(data=lrat0.values,index=z)
                for alt in altlist:
                    lrattemp.ix[alt:alt+layerwidth]=lrattemp.ix[alt:alt+layerwidth]*(1.0+lrat_multiplier)
            else:
                lrattemp=lrat0*(1.0+lrat_multiplier)
            betatemp,sigmatemp=klett2(P_in=P_temp,lrat_in=lrattemp,r_m=r_m,k=k)
            deltabetatemp=100.0*(betatemp-P_0['beta_t'])/P_0['beta_t']
            deltasigmatemp=100.0*(sigmatemp-P_0['sigma_t'])/P_0['sigma_t']
            deltalrattemp=100.0*(lrattemp-lratprof)/lratprof
            
            indexer=tuple([ival for ival,iist in zip([nval,lrat_multiplier],[noise,lrat_testrange]) if len(indexlist)>1])
            
            if not indexer:
                lratcalc.loc[:]=lrattemp
                betacalc.loc[:]=betatemp
                sigmacalc.loc[:]=sigmatemp
                deltalrat.loc[:]=deltalrattemp
                deltabeta.loc[:]=deltabetatemp
                deltasigma.loc[:]=deltasigmatemp                
            else:
                lratcalc.loc[:,indexer]=lrattemp
                betacalc.loc[:,indexer]=betatemp
                sigmacalc.loc[:,indexer]=sigmatemp
                deltalrat.loc[:,indexer]=deltalrattemp
                deltabeta.loc[:,indexer]=deltabetatemp
                deltasigma.loc[:,indexer]=deltasigmatemp
                
            if plotall:
                fig=plt.figure()
                ax=fig.add_subplot(111)
                plt.xlabel('Altitude [km]')
                plt.ylabel('Percent Error')
                fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
                yticks = mtick.FormatStrFormatter(fmt)
                ax.yaxis.set_major_formatter(yticks)
                deltasigma.columns.name='Lidar Ratio Error'
                for col in deltasigma.columns:
                    ax.plot(deltasigma.index,deltasigma[col],label='{0}%'.format(col*100))
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Lidar Ratio Error')
                fig.savefig(os.path.join(figlib,'Deltasigma{0}{1}{2}.png'.format(np.log10(beta),alt,lrat)),
                            bbox_inches='tight')
                plt.close()
                
                fig=plt.figure()
                ax=fig.add_subplot(111)
                plt.xlabel('Altitude [km]')
                plt.ylabel('Percent Error')
                fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
                yticks = mtick.FormatStrFormatter(fmt)
                ax.yaxis.set_major_formatter(yticks)
                deltasigma.columns.name='Lidar Ratio Error'
                for col in deltabeta.columns:
                    ax.plot(deltabeta.index,deltabeta[col],label='{0}%'.format(col*100))
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Lidar Ratio Error')
                fig.savefig(os.path.join(figlib,'Deltabeta{0}{1}{2}.png'.format(np.log10(beta),alt,lrat)),
                            bbox_inches='tight')
                plt.close()
                                           
    beta_comp={'beta0':beta0,'betacalc':betacalc,'deltabeta':deltabeta}
    sigma_comp={'sigma0':sigma0,'sigmacalc':sigmacalc,'deltasigma':deltasigma}
    lrat_comp={'lrat0':lrat0,'lratcalc':lratcalc,'deltalrat':deltalrat}
    if saveall:
        pickle.dump(beta_comp,open(os.path.join(proclib,'testdat_beta.p','wb')))
        pickle.dump(sigma_comp,open(os.path.join(proclib,'testdat_sigma.p','wb')))
        pickle.dump(lrat_comp,open(os.path.join(proclib,'testdat_lrat.p','wb')))   
    
    return beta_comp,sigma_comp,lrat_comp
    
def lrat_tester_full(z,beta_list,alt_list,lrat_list,**kwargs):
    background=kwargs.get('background',0.0)
    layerwidth=kwargs.get('layerwidth',1.0)
    noise=kwargs.get('noise',0.0)
    wave=kwargs.get('wave',532.0)
    E0=kwargs.get('E0',1.0)
    method=kwargs.get('method','range')
    lrat_testrange=kwargs.get('lrat_testrange',np.arange(-.50,1.55,.05))
    numiter=kwargs.get('numiter',10)
    r_m=kwargs.get('r_m',None)
    k=kwargs.get('k',1.0)
    verbose=kwargs.get('verbose',True)
    layeronly=kwargs.get('layeronly',True)
    plotall=kwargs.get('plotall',True)
    figlib=kwargs.get('figlib','./Figures')
    saveall=kwargs.get('saveall',False)
    proclib=kwargs.get('proclib','./Processed')
    
    lrat_m = 8.0*np.pi/3.0
                                                           
    repnum=1 
    beta0_all=[]#np.empty([len(beta_list),len(alt_list),len(lrat_list)]) 
    sigma0_all=[]#np.empty_like(beta0_all)
    lrat0_all=[]#np.empty_like(beta0_all)
    beta_calc_all=[]#np.empty_like(beta0_all)
    sigma_calc_all=[]#np.empty_like(beta0_all)
    lrat_calc_all=[]#np.empty_like(beta0_all)
    delta_beta_all=[]#np.empty_like(beta0_all)
    delta_sigma_all=[]#np.empty_like(beta0_all)
    delta_lrat_all=[]#np.empty_like(beta0_all)  
    

    totreps=len(beta_list)*len(alt_list)*len(lrat_list)*len(lrat_testrange)  
    repnum=1                                          
    for beta in beta_list:
        outaltbeta0=[]
        outaltsigma0=[]
        outaltlrat0=[]
        outaltbetacalc=[]
        outaltsigmacalc=[]
        outaltlratcalc=[]
        outaltdeltabeta=[]
        outaltdeltasigma=[]
        outaltdeltalrat=[]
        for alt in alt_list:
            outlratbeta0=[]
            outlratsigma0=[]
            outlratlrat0=[]
            outlratbetacalc=[]
            outlratsigmacalc=[]
            outlratlratcalc=[]
            outlratdeltabeta=[]
            outlratdeltasigma=[]
            outlratdeltalrat=[]
            for lrat in lrat_list:
                layer=pan.Series(data=beta,index=[alt,alt+layerwidth])
                testlayer={'beta_p':layer,'lrat':lrat}
                P_0=profgen(z,layers=[testlayer],background=background,noise=noise,wave=wave,E0=E0)
                P_temp=P_0['NRB']
                lratprof=P_0['sigma_t'].div(P_0['beta_t'])
                lrat0=pan.Series(data=lratprof,index=z)
                beta0=pan.Series(data=P_0['beta_t'],index=z)
                sigma0=pan.Series(data=P_0['sigma_t'],index=z)
                if method=='range':
                    lratcalc=pan.DataFrame(index=z,columns=lrat_testrange)
                    betacalc=pan.DataFrame(index=z,columns=lrat_testrange)
                    sigmacalc=pan.DataFrame(index=z,columns=lrat_testrange)
                    deltalrat=pan.DataFrame(index=z,columns=lrat_testrange)
                    deltabeta=pan.DataFrame(index=z,columns=lrat_testrange)
                    deltasigma=pan.DataFrame(index=z,columns=lrat_testrange)
                    
                    for lrat_multiplier in lrat_testrange: 
                        if verbose:
                            print 'Processing step {0} out of {1}'.format(repnum,totreps)
                        repnum+=1
                        if layeronly:
                            lrattemp=pan.Series(data=lrat0.values,index=z)
                            lrattemp.ix[alt:alt+layerwidth]=lrattemp.ix[alt:alt+layerwidth]*(1.0+lrat_multiplier)
                        else:
                            lrattemp=lrat0*(1.0+lrat_multiplier)
                        betatemp,sigmatemp=klett2(P_in=P_temp,lrat_in=lrattemp,r_m=r_m,k=k)
                        deltabetatemp=100.0*(betatemp-P_0['beta_t'])/P_0['beta_t']
                        deltasigmatemp=100.0*(sigmatemp-P_0['sigma_t'])/P_0['sigma_t']
                        deltalrattemp=100.0*(lrattemp-lratprof)/lratprof
                        
                        lratcalc.loc[:,lrat_multiplier]=lrattemp
                        betacalc.loc[:,lrat_multiplier]=betatemp
                        sigmacalc.loc[:,lrat_multiplier]=sigmatemp
                        deltalrat.loc[:,lrat_multiplier]=deltalrattemp
                        deltabeta.loc[:,lrat_multiplier]=deltabetatemp
                        deltasigma.loc[:,lrat_multiplier]=deltasigmatemp
                        
                elif method=='assigned':
                    if verbose:
                        print 'Processing step {0} out of {1}'.format(repnum,totreps)
                    repnum+=1
                    lratcalc=pan.Series(data=lrat_m,index=z)
                    lratcalc.ix[alt:alt+layerwidth]=lrat
                    betacalc,sigmacalc=klett2(P_in=P_temp,lrat_in=lratcalc,r_m=r_m,k=k)
                    deltabeta=100.0*(betacalc-P_0['beta_t'])/P_0['beta_t']
                    deltasigma=100.0*(sigmacalc-P_0['sigma_t'])/P_0['sigma_t']
                    deltalrat=100.0*(lratcalc-lratprof)/lratprof
                
                if plotall:
                    fig=plt.figure()
                    ax=fig.add_subplot(111)
                    plt.xlabel('Altitude [km]')
                    plt.ylabel('Percent Error')
                    fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
                    yticks = mtick.FormatStrFormatter(fmt)
                    ax.yaxis.set_major_formatter(yticks)
                    deltasigma.columns.name='Lidar Ratio Error'
                    for col in deltasigma.columns:
                        ax.plot(deltasigma.index,deltasigma[col],label='{0}%'.format(col*100))
                    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Lidar Ratio Error')
                    fig.savefig(os.path.join(figlib,'Deltasigma{0}{1}{2}.png'.format(np.log10(beta),alt,lrat)),
                                bbox_inches='tight')
                    plt.close()
                    
                    fig=plt.figure()
                    ax=fig.add_subplot(111)
                    plt.xlabel('Altitude [km]')
                    plt.ylabel('Percent Error')
                    fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
                    yticks = mtick.FormatStrFormatter(fmt)
                    ax.yaxis.set_major_formatter(yticks)
                    deltasigma.columns.name='Lidar Ratio Error'
                    for col in deltabeta.columns:
                        ax.plot(deltabeta.index,deltabeta[col],label='{0}%'.format(col*100))
                    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Lidar Ratio Error')
                    fig.savefig(os.path.join(figlib,'Deltabeta{0}{1}{2}.png'.format(np.log10(beta),alt,lrat)),
                                bbox_inches='tight')
                    plt.close()
                                   
                if saveall:
                    outlratbeta0.append(beta0)
                    outlratsigma0.append(sigma0)
                    outlratlrat0.append(lrat0)
                    outlratbetacalc.append(betacalc)
                    outlratsigmacalc.append(sigmacalc)
                    outlratlratcalc.append(lratcalc)
                    outlratdeltabeta.append(deltabeta)
                    outlratdeltasigma.append(deltasigma)
                    outlratdeltalrat.append(deltalrat)
            if saveall:
                outaltbeta0.append(outlratbeta0)
                outaltsigma0.append(outlratsigma0)
                outaltlrat0.append(outlratlrat0)
                outaltbetacalc.append(outlratbetacalc)
                outaltsigmacalc.append(outlratsigmacalc)
                outaltlratcalc.append(outlratlratcalc)
                outaltdeltabeta.append(outlratdeltabeta)
                outaltdeltasigma.append(outlratdeltasigma)
                outaltdeltalrat.append(outlratdeltalrat)
        if saveall:
            beta0_all.append(outaltbeta0)
            sigma0_all.append(outaltsigma0)
            lrat0_all.append(outaltlrat0)
            beta_calc_all.append(outaltbetacalc)
            sigma_calc_all.append(outaltsigmacalc)
            lrat_calc_all.append(outaltlratcalc)
            delta_beta_all.append(outaltdeltabeta)
            delta_sigma_all.append(outaltdeltasigma)
            delta_lrat_all.append(outaltdeltalrat)            
    
    if saveall:
        beta_comp=[beta0_all,beta_calc_all,delta_beta_all]
        sigma_comp=[sigma0_all,sigma_calc_all,delta_sigma_all]
        lrat_comp=[lrat0_all,lrat_calc_all,delta_lrat_all]
        
        pickle.dump(beta_comp,open(os.path.join(proclib,'testdat_beta.p','wb')))
        pickle.dump(sigma_comp,open(os.path.join(proclib,'testdat_sigma.p','wb')))
        pickle.dump(lrat_comp,open(os.path.join(proclib,'testdat_lrat.p','wb')))
    #may expand to other methods at some future date
    #    if method=='fernald':
#        beta_out=pan.DataFrame(index=P_0.index,columns=lrat_fern)
#        sigma_out=pan.DataFrame(index=P_0.index,columns=lrat_fern)
#        deltabeta=pan.DataFrame(index=P_0.index,columns=lrat_fern)
#        deltasigma=pan.DataFrame(index=P_0.index,columns=lrat_fern)
#        for lrat in lrat_fern:
#            beta_out[lrat],sigma_out[lrat]=fernald(P_0['NRB'], lrat, wave=wave, E=E0, calrange=calrange_fern)
#            deltabeta[lrat]=100.0*(P_0['beta_t']-beta_out[lrat])/P_0['beta_t']
#            deltasigma[lrat]=100.0*(P_0['sigma_t']-sigma_out[lrat])/P_0['sigma_t']
#    elif method=='klett':
#        beta_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
#        sigma_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
#        deltabeta=pan.DataFrame(index=P_0.index,columns=lrat_klett)
#        deltasigma=pan.DataFrame(index=P_0.index,columns=lrat_klett)
#        lratprof=P_0['sigma_t']/P_0['beta_t'] 
#        lratprof.fillna(value=8.0*np.pi/3.0,inplace=True)
#        if not r_m:
#            r_m=P_0.index.values[-1]    
#        for tempdelt in lrat_klett:
#            lrattemp=lratprof*tempdelt
#            beta_out[tempdelt],sigma_out[tempdelt]=klett(P_in=P_0['NRB'],lrat_in=lrattemp,
#                                                            r_m=r_m,k=k,wave=wave)
#            deltabeta[tempdelt]=100.0*(P_0['beta_t']-beta_out[tempdelt])/P_0['beta_t']
#            deltasigma[tempdelt]=100.0*(P_0['sigma_t']-sigma_out[tempdelt])/P_0['sigma_t']
#    return data_out

def lrat_tester_profile(P_0,lrat_rangein,rangetype='assigned',**kwargs):
    wave=kwargs.get('wave',532.0)
    lrat0=kwargs.get('lrat0',None)
    E0=kwargs.get('E0',1.0)
    method=kwargs.get('method','klett2')
    calrange_fern=kwargs.get('calrange_fern',None)
    r_m=kwargs.get('r_m',None)
    ylim=kwargs.get('ylim',None)
    doplots=kwargs.get('doplots',True)
    saveplots=kwargs.get('saveplots',True)
    layerinfo=kwargs.get('layerinfo',None)

    
    NRB=P_0['NRB']
    beta0=P_0['beta_t']
    sigma0=P_0['sigma_t']
    
    z=NRB.index

    if method=='fernald':
        if rangetype=='assigned' or rangetype=='assigned_layer':
            lrat_range=lrat_rangein
            indexvals=["{:.1f}".format(v) for v in lrat_range]
            lratindex=pan.Index(indexvals,name='Lidar Ratio [sr]')
        elif rangetype=='percent' or rangetype=='percent_layer':
            lrat_range=[lrat0*(1.0+v) for v in lrat_rangein]
            indexvals=["{:+.0%}".format(v) for v in lrat_range]
            lratindex=pan.Index(indexvals,name='Lidar Ratio Error')
        
        df_sigmaout=pan.DataFrame(index=z,columns=lratindex)
        df_deltasigma=pan.DataFrame(index=z,columns=lratindex)
        for lrat,ival in zip(lrat_range,indexvals):
            beta_temp,sigma_temp=fernald(NRB, lrat, wave=wave, E=E0, calrange=calrange_fern)
            deltasigma_temp=100.0*(sigma0-sigma_temp)/sigma0
            
            df_sigmaout.loc[:,ival]=sigma_temp
            df_deltasigma.loc[:,ival]=deltasigma_temp
    elif method=='klett2':
        if rangetype=='assigned':
            indexvals=["{:.1f}".format(v) for v in lrat_rangein]
            lratindex=pan.Index(indexvals,name='Lidar Ratio [sr]')
            lrat_profs=[]
            for lrat in lrat_rangein:
                tempprof=pan.Series(data=lrat,index=z)
                lrat_profs.append(tempprof)
        elif rangetype=='assigned_layer':
            indexvals=["{:.1f}".format(v) for v in lrat_rangein]
            lratindex=pan.Index(indexvals,name='Layer Lidar Ratio [sr]')
            lrat_profs=[]
            for lrat in lrat_rangein:
                tempprof=lrat*(P_0['sigma_p']!=0)
                lrat_profs.append(tempprof)
        elif rangetype=='percent':
            lrat0=P_0['sigma_t'].div(P_0['beta_t']).fillna(0)
            indexvals=["{:+.0%}".format(v) for v in lrat_rangein]
            lratindex=pan.Index(indexvals,name='Lidar Ratio Error')
            lrat_profs=[]
            for lrat in lrat_rangein:
                tempprof=lrat0*(1.0+lrat)
                lrat_profs.append(tempprof)
        elif rangetype=='percent_layer':
            lrat0=P_0['sigma_p'].div(P_0['beta_p']).fillna(0)
            indexvals=["{:+.0%}".format(v) for v in lrat_rangein]
            lratindex=pan.Index(indexvals,name='Layer Lidar Ratio Error')
            lrat_profs=[]
            for lrat in lrat_rangein:
                tempprof=lrat0*(1.0+lrat)
                lrat_profs.append(tempprof)
             
        df_sigmaout=pan.DataFrame(index=z,columns=lratindex)
        df_deltasigma=pan.DataFrame(index=z,columns=lratindex)
        
        if not r_m:
            r_m=z.values[-1]
        
        beta_m=beta0.loc[r_m] 

        for lprof,ival in zip(lrat_profs,indexvals):            
            beta_temp,sigma_temp=klett2(P_in=NRB,lrat_in=lprof,r_m=r_m,beta_m=beta_m)
            deltasigma_temp=100.0*(sigma0-sigma_temp)/sigma0
            
            df_sigmaout.loc[:,ival]=sigma_temp
            df_deltasigma.loc[:,ival]=deltasigma_temp
    
    if doplots:
        fig1=plt.figure(figsize=(20, 5), dpi=80)
        ax1=fig1.add_subplot(111)
        sigma0.plot(ax=ax1,style='-k',linewidth=1.0,label='')
        lines = ["--","-.",":"]
        linecycler = cycle(lines)
        
        for col in df_sigmaout.columns:
            ax1.plot(df_sigmaout.index,df_sigmaout.loc[:,col],next(linecycler),linewidth=4.0,label=col)
#        df_sigmaout.plot(ax=ax1)
        plt.xlabel('Altitude [km]')
        plt.ylabel('Extionction [$km^{-1}$]')
        plt.legend(title=df_sigmaout.columns.name)
        if ylim is not None:
            plt.ylim(ylim)
        if saveplots:
            if layerinfo is not None:
                tempbeta=[]
                templrat=[]                
                for layer in layerinfo:
                    tempbeta.append('{:.0e}'.format(layer['beta_p'].max()))
                    templrat.append('{:.0f}'.format(layer['lrat']))
                
                if len(tempbeta)==1:
                    betaname=tempbeta
                    lratname=templrat
                else:
                    betaname='{0}-{1}'.format(tempbeta[0],tempbeta[-1])
                    lratname='{0}-{1}'.format(templrat[0],templrat[-1])
                    
            savename='Sigmaprofs_{}_{}-beta{}_lrat{}-range{:.0f}-{:.0f}.png'.format(method,rangetype,betaname,lratname,lrat_rangein[0],lrat_rangein[-1])
            plt.savefig(savename,bbox_inches='tight')
        
#        fig2=plt.figure()
#        ax2=fig2.add_subplot(111)
#        df_deltasigma.plot(ax=ax2)
#        if saveplots:
#            savename='Deltasigmaprofs_{0}-{1}.png'.format(lrat_range[0],lrat_range[-1])
#            plt.savefig(savename)
#        

    return df_sigmaout,df_deltasigma

def profile_input_tester(beta_list=[],alt_list=[],lrat_list=[]):
    
    betapanel=pan.Panel(items=beta_list,major_axis=alt_list,minor_axis=lrat_list)
    sigmapanel=pan.Panel(items=beta_list,major_axis=alt_list,minor_axis=lrat_list)
    total=len(beta_list)*len(alt_list)*len(lrat_list)
    n=0
    for beta in beta_list:
        for alt in alt_list:
            for lrat in lrat_list:
                templayer=[{'beta_p':beta,'bot':alt-500.0,'top':alt+500.0,'lrat':lrat}]
                
                P_1 = profgen(z,E0=E0,layers=templayer,donorm=True,background=background)
                betapoint_temp,sigmapoint_temp=lrat_tester_quick(P_1,method='klett2',lrat_klett=lrat_testrange)
                betapanel.loc[beta,alt,lrat]=betapoint_temp[0]
                sigmapanel.loc[beta,alt,lrat]=betapoint_temp[0]
                n+=1
                print "Completed run {0} out of {1}".format(n,total)
    
    return betapanel,sigmapanel

def MCtest(z,numruns=100,testmethod=lrat_tester_noise,**runkwargs):
    
    betacalcout=pan.DataFrame(columns=range(numruns))
    deltabetaout=pan.DataFrame(columns=range(numruns))
    sigmacalcout=pan.DataFrame(columns=range(numruns))
    deltasigmaout=pan.DataFrame(columns=range(numruns))
    
    for r in range(numruns):
        print "Calculating run {0} of {1}".format(r+1,numruns)
        beta_comp,sigma_comp,lrat_comp=testmethod(z,**runkwargs)
        betacalcout.loc[:,r]=beta_comp['betacalc']
        deltabetaout.loc[:,r]=beta_comp['deltabeta']
        
        sigmacalcout.loc[:,r]=sigma_comp['sigmacalc']
        deltasigmaout.loc[:,r]=sigma_comp['deltasigma']
        
    betameanout=betacalcout.mean(1)
    betastdout=betacalcout.std(1)
    deltabetameanout=deltabetaout.mean(1)
    deltabetastdout=deltabetaout.std(1)

    sigmameanout=sigmacalcout.mean(1)
    sigmastdout=sigmacalcout.std(1) 
    deltasigmameanout=deltasigmaout.mean(1)
    deltasigmastdout=deltasigmaout.std(1)
   
    betastats={'betamean':betameanout,'betastd':betastdout,
               'deltabetamean':deltabetameanout,'deltabetastd':deltabetastdout}     
    sigmastats={'sigmamean':sigmameanout,'sigmastd':sigmastdout,
                'deltasigmamean':deltasigmameanout,'deltasigmastd':deltasigmastdout} 
    
    return betastats,sigmastats
    
if __name__ == '__main__':
    os.chdir('C:\Users\dashamstyr\Dropbox\PhD-General\Thesis\Thesis Sandbox')
    z = np.arange(0.150,15.000,0.15,dtype='float')
    E0=1.0 
    background = 0
    numruns=500
#    SNRrange=[0.1,0.5,1.0,1.1,1.5,1.8,2.0,4.0,5.0,10.0,20.0]
    SNRrange=[1.0,2.0,5.0,20.0]
    layerwidth=3.0
    layeronly=True
    betalist=[1e-3]#[1.0*10**-exp for exp in np.arange(3,8,1)]
    altlist=[7.5]#np.arange(11.750,13.750,1.000)
    lratlist=[65]#np.arange(15.0,65.0,10.0)
    
#    lrat_testrange=[8.0*np.pi/3.0,35,65,85]
    lrat_testrange=[-0.5]
    interval=np.arange(1.5,15,1.5)
    
    beta_layer1=pan.Series(data=1e-4, index=[2.500,4.500])
    beta_layer2=pan.Series(data=1e-4, index=[7.500,9.500])
    layer1={'beta_p':beta_layer1,'lrat':65}
    layer2={'beta_p':beta_layer2,'lrat':65}
    layerinfo=[layer2]
#    P_0=profgen(z,layers=layerinfo,background=background)
#    sigma,deltasigma=lrat_tester_profile(P_0,lrat_testrange,method='fernald',rangetype='assigned',lrat0=35,
#                                         layerinfo=layerinfo,ylim=[0.0,0.12])
    
    
    
    
    SNRindex=pan.Index(SNRrange,name='Mean SNR for topmost 1 km')
    
    df_sigmamean=pan.DataFrame(index=z,columns=SNRindex)
    df_sigmastd=pan.DataFrame(index=z,columns=SNRindex)
    df_sigmasnr=pan.DataFrame(index=z,columns=SNRindex)
    snrvals=pan.Series(index=SNRindex)


    templayer=pan.Series(data=betalist[0],index=[altlist[0],altlist[0]+layerwidth])
    templayer={'beta_p':templayer,'lrat':lratlist[0]}

    P_0=profgen(z,layers=[templayer],background=background,noise=0.0,E0=E0)
    P_temp=P_0['NRB']
    lratprof=P_0['sigma_t'].div(P_0['beta_t'])
    lrat0=pan.Series(data=lratprof,index=z)
    beta0=pan.Series(data=P_0['beta_t'],index=z)
    sigma0=pan.Series(data=P_0['sigma_t'],index=z)
 
    
    for SNRval in SNRrange:
#        for betaval in betalist:
#            for altval in altlist:
#                for deltalrat in deltalratlist:
        
        tempmean=P_0['vals'].loc[(z[-1]-1):].mean()
        tempstd=tempmean/SNRval
        
        runkwargs={'background':background,
                   'layerwidth':layerwidth,
                   'betalist':betalist,
                   'altlist':altlist,
                   'lratlist':lratlist,
                   'noise':[tempstd],
                   'lrat_testrange':lrat_testrange,
                   'plotall':False,
                   'verbose':False}
                   
        betastats,sigmastats=MCtest(z=z,numruns=numruns,**runkwargs)
        tempsigmamean=sigmastats['sigmamean']
        tempsigmastd=sigmastats['sigmastd']        

        df_sigmamean.loc[:,SNRval]=tempsigmamean
        df_sigmastd.loc[:,SNRval]=tempsigmastd       
    
        fig1=plt.figure(figsize=(20, 5), dpi=80)
        ax1=fig1.add_subplot(111)
        ax1=sigma0.plot(ax=ax1,style='--k',linewidth=4.0,label='',logy=False)
        ax1.plot(tempsigmamean.index,tempsigmamean,'-',linewidth=2.0,color='b',label='{:.2f}'.format(SNRval))
        ax1.fill_between(tempsigmamean.index,(tempsigmamean+tempsigmastd),(tempsigmamean-tempsigmastd),color='b',alpha=0.2)
        ax1.fill_between(tempsigmamean.index,(tempsigmamean+2*tempsigmastd),(tempsigmamean-2*tempsigmastd),color='g',alpha=0.2)
        tempmean=tempsigmamean.loc[interval]        
        ax1.errorbar(interval,tempmean,yerr=tempsigmastd.loc[interval],fmt='ro')  
        plt.xlabel('Altitude [km]')
        plt.ylabel('Extinction [$km^{-1}$]')
        plt.legend(title='Mean SNR for topmost 1 km',loc=2)
        plt.savefig('SNR-{:.1f}_lrat-{:.0f}_deltalrat={:.1f}_layerbeta-{:.0e}_alt-{:.0f}-{:.0f}_numruns-{:.0f}.png'.format(SNRval,lratlist[0],lrat_testrange[0],betalist[0],altlist[0],(altlist[0]+layerwidth),numruns),
                    bbox_inches='tight')
    
    
#    df_betamean.plot(logy=False)    
#    beta0.plot(style='k')
#    for c in df_betamean.columns:
#        tempmean=df_betamean.loc[interval,c]        
#        plt.errorbar(interval,tempmean,yerr=df_betastd.loc[:,c],fmt='ro')
#    plt.savefig('Beta_{:.0e}.png'.format(noiserange[0]))
#    
#    df_deltabetamean.plot(logy=False)
#    beta0.plot(style='k')
#    for c in df_deltabetamean.columns:
#        tempmean=df_deltabetamean.loc[interval,c]        
#        plt.errorbar(interval,tempmean,yerr=df_deltabetastd.loc[:,c],fmt='ro')
#    plt.savefig('Delta-Beta_{:.0e}.png'.format(noiserange[0]))
#    
#    
#    df_deltasigmamean.plot(logy=False)
#    sigma0.plot(style='k')
#    for c in df_deltasigmamean.columns:
#        tempmean=df_deltasigmamean.loc[interval,c]        
#        plt.errorbar(interval,tempmean,yerr=df_deltasigmastd.loc[:,c],fmt='ro')    
#    plt.savefig('Delta-Sigma_{:.0e}.png'.format(noiserange[0]))
    
    
                                        
#    lrat_tester_full(z=z,beta_list=beta_list,alt_list=alt_list,lrat_list=lrat_list,
#                     lrat_testrange=lrat_testrange,noise=noise,background=background,
#                     layerwidth=layerwidth,method='range',layeronly=layeronly)
    
    
        
