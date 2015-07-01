import numpy as np
import pandas as pan
from copy import deepcopy
import random
import matplotlib.pyplot as plt
from itertools import groupby
import operator

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

    alt = [0, 11000, 20000, 32000, 47000]
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

    beta = backscatter coefficients [1/m*sr]
    sigma = extinction coefficients [1/m]
    """

    #calculate temperature and pressure profiles using US Standard atmosphere

    [T,P,d] = std_atm(z)

    #determine backscatter and extinction coefficients as per Rayleigh scattering
    #equations and constats from Kovalev pp 33-36

    T_s = 288.15  #[K] reference temperature
    P_s = 101325.0  #[Pa] reference pressure
    N_s = 2.547e25  #[1/m^3]  reference number concentration
    gamma = 0.0279 #[unitless] depolarization factor (from Kovalev pg.35)

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


    #For Rayleigh scattering the extinction to backscatter ratio is 8*pi/3

    beta = 3*sigma/(8*np.pi)

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
       P_z.loc[oldalt] = E0*C*(oldalt/1000.0)**-2*beta_R.iloc[0]*T_total
    
    for alt in z[1:]:
        T_step = np.exp(-2*sigma_R.loc[alt]*(alt-oldalt))
        T_total = T_total*T_step
        P_z.loc[alt] = (E0*C*beta_R.loc[alt]*T_total)*(alt/1000.0)**-2.0
        oldalt=alt
        
    rsq = [val*(alt/1000.0)**2.0 for val,alt in zip(P_z,z)]
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

def backandnoise(P_in,background = 0.0,stdev = 0.0,inplace=True):
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
    
    if type(stdev)==float:
        P_out['vals'] = [v+random.gauss(background,stdev) for v in P_out['vals']]
    elif stdev:
        P_out['vals'] = [v+random.gauss(background,s) for v,s in zip(P_out['vals'],stdev.values)]
    else:
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
    P_out['rsq'] = P_out['vals']*(P_in.index.values/1000.0)**2.0
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

    
    S=np.log(P_in).fillna(method='pad',inplace=True)
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
        
        else:
            X1 = (lrat.loc[oldalt]/lrat.loc[alt])**(1/k)   
            sigma.loc[alt] = X1*np.exp((S_new.loc[alt]-S_new.loc[oldalt])/k)/ \
            (sigma.loc[oldalt]**-1+(1.0/k)*(1+X1*np.exp((S_new.loc[alt]-S_new.loc[oldalt])/k))*(oldalt-alt))
            
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
            and the point at which the extinction coefficeint is assumed ot be known
    lrat_in = a Pandas series with values of particulate only lidar ratio and altitude index
    beta_m = the backscatter coefficient at altitude r_m
        
    Outputs:
    beta = pandas series of backscatter coefficients
    sigma = pandas series of extinction coefficients
    
    """
    r_m=kwargs.get('r_m',None)
    wave=kwargs.get('wave',532.0)
    
    lrat = (1.0/lrat_in).replace(np.inf,0.0) #Klett definition of lidar ratio is backscatter/extintion not the other way round
    lrat_R=3.0/(8.0*np.pi)
    altitudes = P_in.index.values
    
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
    P_trans = P_new+0.001-np.min(P_new.values)
    S=np.log(P_trans)
    
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
                
            if beta.ix[alt]<=beta_R.ix[alt]:
                beta_p=0.0
            else:
                beta_p=beta.ix[alt]-beta_R.ix[alt]
            sigma.ix[alt]=(beta_p/lrat.ix[alt])+beta_R.ix[alt]/lrat_R
            oldalt = alt
            old_delta_Sprime=delta_Sprime
    return beta,sigma

def invert_profile(profin,lratin,**kwargs):
    method=kwargs.get('method','klett2')
    refalt=kwargs.get('refalt',profin.index[-1])
    backscatter=pan.Series(np.nan,index=profin.index)
    extinction=pan.Series(np.nan,index=profin.index)
    
    mollayers=lratin[lratin==0.0]
    moledges=[]
    altstep=lratin.index[1]-lratin.index[0]
    
    for key,alt in groupby(enumerate(mollayers.index),lambda (i,x):i-(x-mollayers.index[0])/altstep):
        temprange=map(operator.itemgetter(1),alt)
        moledges.append((temprange[0],temprange[-1]))
    
    alt=profin.index[-1]
    tempedges=moledges.pop()
    while True:        
        if tempedges[0]<alt<=tempedges[1]: 
            layeralts=profin.ix[tempedges[0]:tempedges[1]].index.values
            tempmol=molprof(z=layeralts)
            backscatter.ix[layeralts]=tempmol['beta_R'].values
            extinction.ix[layeralts]=tempmol['sigma_R'].values
            alt=layeralts[0]
        else:
            try:
                tempedges=moledges.pop()
                layeralts=profin.ix[tempedges[1]:alt].index[1:]
                alt=tempedges[1]
            except IndexError:
                layeralts=profin.ix[:alt].index
                alt=layeralts[0]
                
            layerprof=profin.ix[layeralts]
            layerlrat=lratin.ix[layeralts]

            if method=='klett2':
                tempback,tempext=klett2(layerprof,layerlrat,r_m=layeralts[-1])
    
            backscatter.ix[layeralts]=tempback.values
            extinction.ix[layeralts]=tempext.values
            
        if alt==profin.index[0]:
            break
        
    return backscatter,extinction
            
    #Step 2: Use Fernald algorithm for molecular sections    
    #Step 3: Use klett2 for layers

def lrat_tester_full(P_0,**kwargs):
    wave=kwargs.get('wave',532.0)
    E0=kwargs.get('E0',1.0)
    method=kwargs.get('method','klett2')
    lrat_klett=kwargs.get('lrat_klett',np.arange(.50,1.55,.5))
    lrat_fern=kwargs.get('lrat_fern',np.arange(15,80))
    calrange_fern=kwargs.get('calrange_fern',None)
    r_m=kwargs.get('r_m',[])
    k=kwargs.get('k',1.0)
    
    if method=='fernald':
        beta_out=pan.DataFrame(index=P_0.index,columns=lrat_fern)
        sigma_out=pan.DataFrame(index=P_0.index,columns=lrat_fern)
        deltabeta=pan.DataFrame(index=P_0.index,columns=lrat_fern)
        deltasigma=pan.DataFrame(index=P_0.index,columns=lrat_fern)
        for lrat in lrat_fern:
            beta_out[lrat],sigma_out[lrat]=fernald(P_0['NRB'], lrat, wave=wave, E=E0, calrange=calrange_fern)
            deltabeta[lrat]=100.0*(P_0['beta_t']-beta_out[lrat])/P_0['beta_t']
            deltasigma[lrat]=100.0*(P_0['sigma_t']-sigma_out[lrat])/P_0['sigma_t']
    elif method=='klett':
        beta_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        sigma_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        deltabeta=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        deltasigma=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        lratprof=P_0['sigma_t']/P_0['beta_t'] 
        lratprof.fillna(value=8.0*np.pi/3.0,inplace=True)
        if not r_m:
            r_m=P_0.index.values[-1]    
        for tempdelt in lrat_klett:
            lrattemp=lratprof*tempdelt
            beta_out[tempdelt],sigma_out[tempdelt]=klett(P_in=P_0['NRB'],lrat_in=lrattemp,
                                                            r_m=r_m,k=k,wave=wave)
            deltabeta[tempdelt]=100.0*(P_0['beta_t']-beta_out[tempdelt])/P_0['beta_t']
            deltasigma[tempdelt]=100.0*(P_0['sigma_t']-sigma_out[tempdelt])/P_0['sigma_t']
    elif method=='klett2':
        beta_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        sigma_out=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        deltabeta=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        deltasigma=pan.DataFrame(index=P_0.index,columns=lrat_klett)
        lratprof=P_0['sigma_p']/P_0['beta_p']         
        lratprof.fillna(value=0.0,inplace=True)
        if not r_m:
            r_m=P_0.index.values[-1]

        for tempdelt in lrat_klett:
            lrattemp=lratprof*tempdelt
            beta_out[tempdelt],sigma_out[tempdelt]=klett2(P_in=P_0['NRB'],lrat_in=lrattemp,r_m=r_m)
            deltabeta[tempdelt]=100.0*(P_0['beta_t']-beta_out[tempdelt])/P_0['beta_t']
            deltasigma[tempdelt]=100.0*(P_0['sigma_t']-sigma_out[tempdelt])/P_0['sigma_t']
    
    panelout=pan.Panel({'beta_calc':beta_out,'sigma_calc':sigma_out,'delta_beta':deltabeta,
                        'delta_sigma':deltasigma})
    
    return panelout

def lrat_tester_quick(P_0,**kwargs):
    wave=kwargs.get('wave',532.0)
    E0=kwargs.get('E0',1.0)
    method=kwargs.get('method','klett2')
    lrat_klett=kwargs.get('lrat_klett',np.arange(-.50,.55,.5))
    lrat_fern=kwargs.get('lrat_fern',np.arange(15,80))
    calrange_fern=kwargs.get('calrange_fern',None)
    r_m=kwargs.get('r_m',None)
    k=kwargs.get('k',1.0)
    lrat_type=kwargs.get('lrat_type','part')
    
    old_deltabeta=0.0
    old_deltasigma=0.0
    
    if method=='fernald':

        for lrat in lrat_fern:
            beta_temp,sigma_temp=fernald(P_0['NRB'], lrat, wave=wave, E=E0, calrange=calrange_fern)
            deltabeta_temp=100.0*(P_0['beta_t']-beta_temp)/P_0['beta_t']
            deltasigma_temp=100.0*(P_0['sigma_t']-sigma_temp)/P_0['sigma_t']
            
            if max(abs(deltabeta_temp)) > old_deltabeta:
                loctemp=abs(deltabeta_temp).idxmax()                
                betapoint=(deltabeta_temp.loc[loctemp],loctemp,lrat)
            if max(abs(deltasigma_temp)) > old_deltasigma:
                loctemp=abs(deltasigma_temp).idxmax()                
                sigmapoint=(deltasigma_temp.loc[loctemp],loctemp,lrat)
                
    elif method=='klett':

        lratprof=P_0['sigma_t']/P_0['beta_t'] 
        lratprof.fillna(value=0,inplace=True)
        if not r_m:
            r_m=P_0.index.values[-1]
        
        sigma_m=P_0['sigma_t'].loc[r_m]        
        for tempdelt in lrat_klett:
            lrattemp=lratprof*tempdelt
            beta_temp,sigma_temp=klett(P_in=P_0['NRB'],lrat_in=lrattemp,
                                                            r_m=r_m,sigma_m=sigma_m,k=k)
            deltabeta_temp=100.0*(P_0['beta_t']-beta_temp)/P_0['beta_t']
            deltasigma_temp=100.0*(P_0['sigma_t']-sigma_temp)/P_0['sigma_t']
            
            if max(abs(deltabeta_temp)) > old_deltabeta:
                loctemp=abs(deltabeta_temp).idxmax()                
                betapoint=(deltabeta_temp.loc[loctemp],loctemp,tempdelt)
            if max(abs(deltasigma_temp)) > old_deltasigma:
                loctemp=abs(deltasigma_temp).idxmax()                
                sigmapoint=(deltasigma_temp.loc[loctemp],loctemp,tempdelt  )
    elif method=='klett2':

        lratprof=P_0['sigma_p']/P_0['beta_p']         
        lratprof.fillna(value=0,inplace=True)
        if not r_m:
            r_m=P_0.index.values[-1]
        
        beta_m=P_0['beta_t'].loc[r_m] 
        for tempdelt in lrat_klett:
            lrattemp=lratprof*tempdelt
            beta_temp,sigma_temp=klett2(P_in=P_0['NRB'],lrat_in=lrattemp,
                                                            r_m=r_m,beta_m=beta_m)
            deltabeta_temp=100.0*(P_0['beta_t']-beta_temp)/P_0['beta_t']
            deltasigma_temp=100.0*(P_0['sigma_t']-sigma_temp)/P_0['sigma_t']
            
            if max(abs(deltabeta_temp)) > old_deltabeta:
                loctemp=abs(deltabeta_temp).idxmax()                
                betapoint=(deltabeta_temp.loc[loctemp],loctemp,tempdelt)
            if max(abs(deltasigma_temp)) > old_deltasigma:
                loctemp=abs(deltasigma_temp).idxmax()                
                sigmapoint=(deltasigma_temp.loc[loctemp],loctemp,tempdelt)

    return betapoint,sigmapoint

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
    
if __name__ == '__main__':

    z = np.arange(150,15000,3,dtype='float')
    E0=1.0 
    background = 1e-6
    noise=0.0
    beta_list=[1.0*10**-exp for exp in range(0,8)]
    alt_list=np.arange(750.0,13750.0,500.0)
    lrat_list=np.arange(15.0,80.0,5.0)
    lrat_testrange=[0.5]
    
    beta_layer=pan.Series(data=[1e-5,1e-5], index=[12500,14500])
    layer1={'beta_p':beta_layer,'lrat':15}
    P_0=profgen(z,layers=[layer1],background=background,noise=noise)

    testpan=lrat_tester_full(P_0,lrat_klett=[1.0])
                
#    beta_fern = fernald(P_1,30,wave,1.0)
#    sigma_fern=beta_fern*lrat_p1
#    
#    #Klett's lrat is the inverse of Fernald's
#    
#    lrat = P_1.loc['sigma_t']/P_1.loc['beta_t']
#    
#    lrat_klett=P_1.loc['sigma_p']/P_1.loc['beta_p'] #klett2 takes only particulte lrat
#    lrat_klett.fillna(0,inplace=True)
#    
#    r_m = z[-2]
#    
#    sigma_m = P_mol.loc['sigma_R'].loc[r_m]
#    beta_klett, sigma_klett = klett(p_norm1,lrat,r_m,sigma_m)
#    beta_klett2 = klett2(p_norm1,lrat_klett,r_m=r_m)
#    sigma_klett2 = beta_klett2*lrat_klett    
#    
#    delta_fern=(beta_fern-P_1.loc['beta_t'])/P_1.loc['beta_t']
#    delta_klett=(beta_klett-P_1.loc['beta_t'])/P_1.loc['beta_t']
#    delta_klett2=(beta_klett2-P_1.loc['beta_t'])/P_1.loc['beta_t']
#    betakeys=['Data','Fernald','Klett','Klett2']
#    deltakeys=['Fernald','Klett','Klett2']
#    betadict=dict(zip(betakeys,[P_1.loc['beta_t'],beta_fern,beta_klett,beta_klett2]))
#    deltadict=dict(zip(deltakeys,[delta_fern,delta_klett,delta_klett2]))
#    sigmadict=dict(zip(betakeys,[P_1.loc['sigma_t'],sigma_fern,sigma_klett,sigma_klett2]))
#    df_beta=pan.DataFrame.from_dict(betadict)
#    df_delta=pan.DataFrame.from_dict(deltadict)
#    df_sigma=pan.DataFrame(sigmadict)
    
#    prat = p_rangecor1/p_rangecor0
#    
#    p_slope0 = calc_slope(p_rangecor0)    
#    p_slope1 = calc_slope(p_rangecor1)    
#    prat_slope = calc_slope(prat)

#    fig1 = plt.figure()
#    ax1 = fig1.add_subplot(1,3,1)
#    ax1.plot(P_mol.loc['Temp'],z)
#    ax1.set_xlabel('Temperature [K]')
#    ax1.set_ylabel('Height [m]')
#
#    ax2 = fig1.add_subplot(1,3,2)
#    ax2.plot(P_mol.loc['Press'],z)
#    ax2.set_xlabel('Pressure [Pa]')
#
#    ax3 = fig1.add_subplot(1,3,3)
#    ax3.plot(P_mol.loc['Density'],z)
#    ax3.set_xlabel('Density [kg/m^3]')
#
#    fig2 = plt.figure()
#    ax1 = fig2.add_subplot(1,3,1)
#    ax1.plot(p_norm1,z,p_norm0,z)
#    ax1.set_xscale('log')
#    ax1.set_xlabel('Range corrected signal multiplier')
#    ax1.set_ylabel('Height [m]')
#
#    ax2 = fig2.add_subplot(1,3,2)
#    ax2.plot((P_1.loc['beta_R']+P_1.loc['beta_p']),z)
#    ax2.set_xlabel('Total backscatter coefficient [1/m/sr')
#
#    ax3 = fig2.add_subplot(1,3,3)
#    ax3.plot((P_1.loc['sigma_R']+P_1.loc['sigma_p']),z)
#    ax3.set_xlabel('Total extinction coefficient [1/m')
#    
#    fig3 = plt.figure()
#    ax1 = fig3.add_subplot(1,3,1)
#    ax1.plot(p_norm1,z,p_noisy,z)
#    ax1.set_xscale('log')
#    ax1.set_xlabel('Range corrected signal multiplier')
#    ax1.set_ylabel('Height [m]')
#
#    ax2 = fig3.add_subplot(1,3,2)
#    ax2.plot((P_1.loc['beta_p']/P_1.loc['beta_R']),z)
#    ax2.set_xlabel('Attenuated Backscatter Ratio [beta_P/beta_R]')
#
#    ax3 = fig3.add_subplot(1,3,3)
#    ax3.plot((P_1.loc['sigma_p']),z)
#    ax3.set_xlabel('Total extinction coefficient [1/m')
#    
#    fig4 = plt.figure()
#    ax1 = fig4.add_subplot(1,3,1)
#    ax1.plot(prat,z)
#    ax1.set_xlabel('Signal ratio')
#    ax1.set_ylabel('Height [m]')
#
#    ax2 = fig4.add_subplot(1,3,2)
#    ax2.plot(p_slope0,z)
#    ax2.set_xlabel('Backscatter Slope')
#    
#    ax3 = fig4.add_subplot(1,3,3)
#    ax3.plot(p_slope1,z)
#    ax3.set_xlabel('Backscat Rat Slope')    
    
#    fig5 = plt.figure()
#    ax1 = fig5.add_subplot(1,2,1)
#    ax1.plot(((P_1.loc['beta_t'].values-beta_klett.values)/P_1.loc['beta_t'].values),z)
#    ax1.set_xlabel('Delta Beta [%]')
#    ax1.set_ylabel('Altitude')
#    
#    ax2 = fig5.add_subplot(1,2,2)
#    ax2.plot(beta_klett.values,z,P_1.loc['beta_t'].values,z)
#    ax2.set_xlabel('Original (green) and Klett (blue) Backscatter coeffs')

#    fig6 = plt.figure()
#    ax1 = fig6.add_subplot(1,2,1)
#    ax1.plot((P_1.loc['sigma_t'].values),z)
#    ax1.set_xlabel('Original backscatter coeffs')
#    ax1.set_ylabel('Altitude')
#    
#    ax2 = fig6.add_subplot(1,2,2)
#    ax2.plot(sigma_klett.values,z)
#    ax2.set_xlabel('Fernald Backscatter coeffs')
#    plt.show()
