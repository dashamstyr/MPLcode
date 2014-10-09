import numpy as np
import pandas as pan
from copy import deepcopy
import random
import matplotlib.pyplot as plt

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


def molecular(z,wave):
    """
    Function for generating molecular scattering and extinction coefficients based
    on an altitude and a laser wavelength.  Ozone absorption is ignored.

    Inputs:

    z = altitude [m]
    wave = wavelength [nm]

    Outputs:

    beta = backscatter coefficients [1/m*sr]
    alpha = extinction coefficients [1/m]
    """
    
    #note: make sure z is in the form of a numpy matrix
    
    beta = np.empty_like(z)
    alpha = np.empty_like(z)

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

    alpha = (8*np.pi**3*(m**2-1)**2*N/(3*N_s**2*(wave*1e-9)**4))*((6+3*gamma)/(6-7*gamma))* \
            (P/P_s)*(T_s/T)


    #For Rayleigh scattering the extinction to backscatter ratio is 8*pi/3

    beta = 3*alpha/(8*np.pi)

    return T,P,d,beta,alpha

def molprof(z,wave, T_0 = 1.0):
    """
    Function for generating a theoretical profile of normalized attenuated
    backscatter.  In other words, this provides
    an array that can be multiplied by lidar output power and system constant
    to provide a lidar response profile

    Inputs:

    z = an array of altitudes [m]
    wave = lidar wavelength [nm]
    T_0 = transmissivity to the first z value [unitless], defaults to 1.0 (perfect transmition)

    Outputs:
    P_out = a pandas dataframe with columns representing altitudes and the following index values
    vals = vector of normalized attenuated backscatter [unitless]
    beta_R = vector of molecular backscatter coefficients [1/m*sr]
    alpha_R = vector of molecular extinction coefficients [1/m]
    beta_p = vector of particulate backscatter coefficients (assumed to be zero here)
    alpha_p = vector of particulate extinction coefficients (assumed to be zero here)
    beta_t = vector of total backscatter coefficients [1/m*sr]
    alpha_t = vector of total extinction coefficients [1/m]
    """
    
    T = np.empty_like(z,dtype=float)
    P = np.empty_like(z,dtype=float)
    d = np.empty_like(z,dtype=float)
    beta_R = np.empty_like(z,dtype=float)
    alpha_R = np.empty_like(z,dtype=float)
    P_z = np.empty_like(z,dtype=float)

    for n in range(len(z)):
        [T[n],P[n],d[n],beta_R[n],alpha_R[n]] = molecular(z[n], wave)
    
    T_total = T_0
    P_z[0] = z[0]**-2*beta_R[0]*T_total**2
    
    for n in range(1,len(z)):
        T_step = np.exp(-alpha_R[n]*(z[n]-z[n-1]))
        T_total = T_total*T_step
        P_z[n] = z[n]**-2.0*beta_R[n]*T_total**2
        
    beta_p = np.zeros_like(beta_R,dtype=float)
    alpha_p = np.zeros_like(alpha_R,dtype=float)
    beta_t = beta_R+beta_p
    alpha_t = alpha_R+alpha_p
    keys = ['vals','Temp','Press','Density','beta_R','alpha_R','beta_p','alpha_p','beta_t','alpha_t']
    vals = [P_z,T,P,d,beta_R,alpha_R,beta_p,alpha_p,beta_t,alpha_t]
    P_out = pan.DataFrame(data=vals,index=keys,columns=z)
    
    return P_out


def addlayer(P_in, layer, lrat, inplace=True):
    """
    Function that adds a layer of known backscatter coefficient
    and lidar ratio onto an existing lidar response profile

    Inputs:

    P_in = a pandas dataframe object containing the input profile with indexed values defined in mplprof

    layer = a pandas series object showing backscatter coefficients [1/m*sr] with altitude   
    lrat = lidar (extinction to backscater) ratio of the added layer [1/sr]

    Outputs:
    P_out = dataframe object like P_in with layer added

    """    
    if inplace:
        P_out=P_in
    else:
        P_out = deepcopy(P_in)
    
    beta_in = layer.values
    z_in = np.array(layer.index.values,dtype='float64')

    z_min = min(z_in)
    z_max = max(z_in)
    
    z_old = z_min
    for z in P_out.columns:
        if z < z_min:
            T_total = (P_out.loc['vals',z]/P_out.loc['vals'].iloc[0])* \
            (P_out.loc['beta_t'].iloc[0]/P_out.loc['beta_t',z])* \
            (z/P_out.columns[0])**2
        elif z <= z_max:
            P_out.loc['beta_p',z] += np.interp(z,z_in,beta_in)
            P_out.loc['alpha_p',z] += P_out.loc['beta_p',z]*lrat
            P_out.loc['beta_t',z] = P_out.loc['beta_p',z] + P_out.loc['beta_R',z]
            P_out.loc['alpha_t',z] = P_out.loc['alpha_p',z] + P_out.loc['alpha_R',z]
            T_step = np.exp(-P_out.loc['alpha_t',z]*(z-z_old))
            T_total = T_total*T_step
            P_out.loc['vals',z] = z**-2*P_out.loc['beta_t',z]*T_total**2
            z_old = z
        else:
            T_step = np.exp(-P_out.loc['alpha_t',z]*(z-z_old))
            T_total = T_total*T_step
            P_out.loc['vals',z] = z**-2*P_out.loc['beta_t',z]*T_total**2
            z_old=z

    return P_out

def backandnoise(P_in,background = 0,stdev = [],inplace=True):
    """
    Adds gaussian random noise and background signal to any profile
    
    Inputs:
    P_in = a pandas series with altitude index
    background = a float defining the backgorund signal to be applied(defaults to 0)
    stdev = the standard deviation of the gaussian noise component
    
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
    
    if stdev:
        P_out.values[:] = [v+random.gauss(background,stdev) for v in P_out.values]
    else:
        P_out.values[:] = [v+random.gauss(background,np.sqrt(v+background)) for v in P_out.values]
    return P_out

def background_subtract(P_in,back_avg=[],z_min=9000):
    #subtracts background from signal, without background
    #takes advantage of inverse square law to calculate background signal for
    #a profile

    #start by selecting a region of sufficient altitude that the r squared law
    #will make the ratio of signal to background sufficiently low

    P_out = deepcopy(P_in)
    
    if back_avg:
        P_out.values=P_out.values-back_avg
    else:    
        #select data from altitudes higher than z_min and muliply full signal by
        #range squared, if this gives less than 500 values, take uppermost 500
    
        z = P_out.index[z_min:]
    
        if len(z) <= 500:
            z = P_out.z[-500:]
    
        r_sq = P_out.values[-len(z):]*z**2
    
        #since background is constant while signal is reduced as z^2, first
        #coefficient is equal to background
        
        coeffs = np.polyfit(z,r_sq,2,full=False)
        
    
        background = coeffs[0]
    
        P_out.backsub_vals = P_out.vals-background
        P_out.back = background

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

def fernald(P_in, lrat, wave = 532, E = 1.0, calrange = []):
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
        alpha_out: a pandas series containing extinction coefficients in units of [1/m]
    
    In the langauge of Fernald's original paper, the lidar equation is reduced to the form:
    P(Z) = E*C*Z^-2*[beta_R(Z)+beta_P(Z)]*T_R(Z)^2*T_P(Z)^2
    Where:
    P(Z)=Power at the receiver
    E=laser output power
    C=system constant
    beta_R & beta_P=backscatter ceofficients for Rayleigh and particulate scattering,respectively
    T_R & T_P = integrated extinction (optical depth) from Rayleigh and particulate scattering
    and  
    T_R = exp[-int_0-Z(alpha_R(z)dz)]
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
    
        
    if not calrange:
        for alt in reversed(altitudes):
            if alt == altitudes[-1]:
                #at calibration altitude assume no aerosols
                beta_total.loc[alt] = Rayleigh_coeffs.loc['beta_R',alt]
                oldalt = alt
            
            else:
                X1 = P_in.loc[oldalt]*E
                X0 = P_in.loc[alt]*E
                beta_R1 = Rayleigh_coeffs.loc['beta_R',oldalt]
                beta_R0 = Rayleigh_coeffs.loc['beta_R',alt]
                delta_Z = oldalt-alt
                A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                
                beta_total.loc[alt] = X0*np.exp(A)/((X1/beta_total[oldalt])+lrat* \
                np.abs(X1+X0*np.exp(A))*delta_Z)
                
                oldalt = alt
    else:
        #first assign Rayleigh coefficients to calibration range
        
        minalt = calrange[0]
        maxalt = calrange[1]
        calalts = [z for z in altitudes if minalt < z < maxalt]
        
        for alt in calalts:
            beta_total.loc[alt] = Rayleigh_coeffs.loc['beta_R',alt]
        
        #next calculate coefficients for range below minalt
        if calalts[0] > altitudes[0]:
            oldalt = calalts[0]
            below = [z for z in altitudes if z <= minalt]
            
            for alt in reversed(below):
                X1 = P_in.loc[oldalt]*E
                X0 = P_in.loc[alt]*E
                beta_R1 = Rayleigh_coeffs.loc['beta_R',oldalt]
                beta_R0 = Rayleigh_coeffs.loc['beta_R',alt]
                delta_Z = oldalt-alt
                A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                
                beta_total.loc[alt] = X0*np.exp(A)/((X1/beta_total[oldalt])+lrat* \
                np.abs(X1+X0*np.exp(A))*delta_Z)
                
                oldalt = alt
        
        #then  calculate for altitudes above maxalt
        elif calalts[-1] < altitudes[-1]:
            oldalt = calalts[-1]
            above = [z for z in altitudes if z >= maxalt]
            
            for alt in above:
                X1 = P_in.loc[oldalt]*E
                X0 = P_in.loc[alt]*E
                beta_R1 = Rayleigh_coeffs.loc['beta_R',oldalt]
                beta_R0 = Rayleigh_coeffs.loc['beta_R',alt]
                delta_Z = oldalt-alt
                A = (lrat-Rayleigh_lrat)*(beta_R1+beta_R0)*delta_Z
                
                beta_total.loc[alt] = X0*np.exp(-A)/((X1/beta_total[oldalt])-lrat* \
                np.abs(X1+X0*np.exp(-A))*delta_Z)
                
                oldalt = alt
    
    return beta_total


def klett(P_in,r_m,lrat,sigma_m,k=1):
    """
    Function that calculates backscatter and extinction coefficients based on
    variable lidar ratios using the Klett algorithm 
    
    Inputs:
    P_in = a pandas series with values of signal strength and altitude index
    r_m = the reference altitude - maximum altitude for which calculations are done
            and the point at which the extinction coefficeint is assumed ot be known
    lrat = a Pandas series with values of lidar ratio and altitude index
    sigma_m = the extinction coefficient at altitude r_m
    k = the power coefficient in the power law relationship bwetween backscatter
        and extinction (defaults to 1)
        
    Outputs:
    beta = pandas series of backscatter coefficients
    sigma = pandas series of extinction coefficients
    
    """

    altitudes = P_in.index.values
    
    sigma = pan.Series(index=altitudes)
    beta = pan.Series(index=altitudes)
    
    for alt in reversed(altitudes):
       if alt > r_m:
           sigma.loc[alt] = 'na'
           beta.loc[alt] = 'na'      
       else:
           P_new = P_in.loc[:alt]
           break
    
    newalts = P_new.index.values
    
    for alt in reversed(newalts):
        if alt == newalts[-1]:            
            sigma.loc[alt] = sigma_m
            beta.loc[alt] = lrat.loc[alt]*sigma_m**k
            oldalt = alt
        
        else:
            X1 = (lrat.loc[oldalt]/lrat.loc[alt])**(1/k)   
            sigma.loc[alt] = X1*np.exp((P_new.loc[alt]-P_new.loc[oldalt])/k)/ \
            (sigma.loc[oldalt]**-1+(1.0/k)*(1+X1*np.exp((P_new.loc[alt]-P_new.loc[oldalt])/k))*(oldalt-alt))
            
            beta.loc[alt] = lrat.loc[alt]*sigma.loc[alt]**k
            oldalt = alt
    
    
    return beta, sigma  
    
if __name__ == '__main__':

    z = np.arange(100,15000,3)

    wave = 532.0  #nm
    lrat_p = 30.0
    lrat_m = 8.0*np.pi/3.0

    P_mol = molprof(z,wave)

    z_layer = np.arange(1000,2000,5,dtype=np.float)

    beta_layer = np.ones_like(z_layer)*5e-6

    layer = pan.Series(data=beta_layer, index=z_layer)
    
    P_1 = addlayer(P_mol,layer,lrat_p)
    
    p_rangecor0 = pan.Series(P_mol.loc['vals'].values*(z**2),index=z)
    p_norm0 = p_rangecor0/p_rangecor0.iloc[0]
    p_rangecor1 = pan.Series(P_1.loc['vals'].values*(z**2) ,index=z)  
    p_norm1 = p_rangecor1/p_rangecor1.iloc[0]
    
    p_noisy = backandnoise(p_norm1,inplace=False)
    
    beta_fern = fernald(p_norm1,lrat_p,wave,1.0)
    
    #Klett's lrat is the inverse of Fernald's
    
    lrat_klett = P_1.loc['beta_t']/P_1.loc['alpha_t']
    
    lrat_klett=lrat_klett*0.99
    
    r_m = z[-1]
    
    sigma_m = P_mol.loc['alpha_R'].loc[r_m]
    
    S_klett = np.log(p_norm1)
    
    beta_klett, sigma_klett = klett(S_klett,r_m,lrat_klett,sigma_m)
    
    prat = p_rangecor1/p_rangecor0
    
    p_slope0 = calc_slope(p_rangecor0)    
    p_slope1 = calc_slope(p_rangecor1)    
    prat_slope = calc_slope(prat)

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
#    ax3.plot((P_1.loc['alpha_R']+P_1.loc['alpha_p']),z)
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
#    ax3.plot((P_1.loc['alpha_p']),z)
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
    
    fig5 = plt.figure()
    ax1 = fig5.add_subplot(1,2,1)
    ax1.plot(((P_1.loc['beta_t'].values-beta_klett.values)/P_1.loc['beta_t'].values),z)
    ax1.set_xlabel('Original backscatter coeffs')
    ax1.set_ylabel('Altitude')
    
    ax2 = fig5.add_subplot(1,2,2)
    ax2.plot(beta_klett.values,z,P_1.loc['beta_t'].values,z)
    ax2.set_xlabel('Fernald Backscatter coeffs')

#    fig6 = plt.figure()
#    ax1 = fig6.add_subplot(1,2,1)
#    ax1.plot((P_1.loc['alpha_t'].values),z)
#    ax1.set_xlabel('Original backscatter coeffs')
#    ax1.set_ylabel('Altitude')
#    
#    ax2 = fig6.add_subplot(1,2,2)
#    ax2.plot(sigma_klett.values,z)
#    ax2.set_xlabel('Fernald Backscatter coeffs')
#    plt.show()
