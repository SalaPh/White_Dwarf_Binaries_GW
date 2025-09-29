############################################################################
# this script calculates information for White Dwarf Binaries (WDB)
# such as the cutoff frequency due to Roche lobe overflow
# and checks if the masses inputted are valid for a WD binary system
############################################################################

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import scipy.constants as phconst

from . import constants as const
from . import waveform_LO

###################################################
# WD CONSTANTS
###################################################

minWDMass=0.1 # [Msol]
maxWDMass=1.44 # [Msol], ChandrasekharMass

###################################################
# FUNCTIONS
###################################################

def WD_RadiusMass_Relation(Mass):
    '''
    white dwarf mass radius relation, 
    Eq. (24) of astro-ph/0312577
    
    input:
        Mass - mass of the WD in [Msol]
    output:
        radius of the WD in [m]
    '''
    #Formula from https://arxiv.org/pdf/astro-ph/0312577.pdf
    Mp=0.00057

    # Make the input an array if it is not already
    if not isinstance(Mass, np.ndarray):
        Mass = np.array([Mass])
    Mass = np.where(Mass == maxWDMass, Mass * (1 - 1e-10), Mass)
    Mass = np.where(Mass == Mp, Mass * (1 + 1e-10), Mass)
    invalid_mass_indices = np.where((Mass < Mp) | (Mass > maxWDMass))
    # Warn if any masses are out of range
    if invalid_mass_indices[0].size > 0:
        print(f'WARNING:\n The mass {Mass[invalid_mass_indices]} ' + r'$M_\odot$' + ' is not in the range of the WD mass. A radius cannot be computed')

    radius = 0.0114*const.SolarRadius*((Mass/maxWDMass)**(-2/3)-(Mass/maxWDMass)**(2/3))**(1/2)*(1+3.5*(Mass/Mp)**(-2/3)+(Mass/Mp)**(-1))**(-2/3)
    return radius if radius.size > 1 else radius[0]

#------------------------------------------------#

def WD_Density(Mass,Radius):
    '''
    Density of a white dwarf given its mass and radius
    input:
        Mass - mass of the WD in [Msol]
        Radius - radius of the WD in [m]
    output:
        density of the WD in [kg/m^3]
    '''
    return Mass*const.SolarMass/((4*np.pi/3)*Radius**3)

#------------------------------------------------#

def RLOF1(q):
    '''
    Roche lobe overflow radius of the more massive body in a binary system
    input:
        q - mass ratio m1/m2, with m1>m2
    output:
        Roche lobe overflow radius in units of the orbital separation
    '''
    #Formula from https://arxiv.org/pdf/1403.4754.pdf
    return 0.49*q**(2/3)/(0.6*q**(2/3) + np.log(1+q**(1/3)))

#------------------------------------------------#

def RLOF2(q):
    '''
    Roche lobe overflow radius  of the less massive body in a binary system
    input:
        q - mass ratio m1/m2, with m1>m2
    output:
        Roche lobe overflow radius in units of the orbital separation
    '''
    return RLOF1(1/q) 

#------------------------------------------------#

def WD_RLOFCutoffFrequency(Mass1,Mass2):
    '''
    Frequency cutoff due to Roche lobe overflow of the less massive body in a binary system
    input:
        Mass1 - mass of the more massive body in [Msol]
        Mass2 - mass of the less massive body in [Msol]
    output:
        Frequency cutoff in [Hz]
    '''
    return waveform_LO.Frequency_Distance(WD_RadiusMass_Relation(Mass2)/RLOF2(Mass1/Mass2),Mass1,Mass2)

#------------------------------------------------#

def WD_TauatRLOF(Mass1,Mass2):
    '''
    Time to merger at Roche lobe overflow of the less massive body in a binary system
    input:
        Mass1 - mass of the more massive body in [Msol]
        Mass2 - mass of the less massive body in [Msol]
    output:
        Time to merger at Roche lobe overflow in [s]
    '''
    return waveform_LO.tau_fGW(waveform_LO.Mc(Mass1,Mass2),WD_RLOFCutoffFrequency(Mass1,Mass2)) 

#------------------------------------------------#

def Delta(M1, M2, a):
    '''
    Roche lobe overflow of the donor,
    Eq. (10) of astro-ph/0312577
    input:
        M1 - mass of the accretor in [Msol]
        M2 - mass of the donor in [Msol]
        a - orbital separation in [m]
    output:
        Delta in [m]
    '''
    return WD_RadiusMass_Relation(M2) - a*RLOF2(M1/M2)

#------------------------------------------------#

def aRLOF(M1, M2):
    '''
    returns the orbital separation for which the donor fills its Roche Lobe
    input:
        M1 - mass of the accretor in [Msol]
        M2 - mass of the donor in [Msol]
    output:
        orbital separation at Roche Lobe overflow [m]
    '''
    return WD_RadiusMass_Relation(M2)/RLOF2(M1/M2)

#------------------------------------------------#
# Analytical derivatives required for computations

def DWD_RadiusMass_Relation(Mass):
    '''
    Derivative over mass of the white dwarf mass radius relation, 
    Eq. (24) of astro-ph/0312577
    
    input:
        Mass - mass of the WD in [Msol]
    output:
        radius of the WD in [m/Msol]
    '''

    Mp=0.00057
    return 0.0114*const.SolarRadius*((1/2)*((Mass/maxWDMass)**(-2/3)-(Mass/maxWDMass)**(2/3))**(-1/2)
                                        *((-2/3)*(Mass/maxWDMass)**(-2/3)-(2/3)*(Mass/maxWDMass)**(2/3))
                                        *(1+3.5*(Mass/Mp)**(-2/3)+(Mass/Mp)**(-1))**(-2/3)
                                        +(-2/3)*((Mass/maxWDMass)**(-2/3)-(Mass/maxWDMass)**(2/3))**(1/2)
                                        *(1+3.5*(Mass/Mp)**(-2/3)+(Mass/Mp)**(-1))**(-5/3)
                                        *((-2/3)*3.5*(Mass/Mp)**(-2/3)+(-1)*(Mass/Mp)**(-1))
                                        )/(Mass)

#------------------------------------------------#

def DMcWD_RadiusMass_RelationR2(Mc,q):
    '''
    Derivative of the mass-radius relation over q for Mc
    input:
        Mc - chirp mass in [Msol]
        q - mass ratio m1/m2, with m1>m2
    output:
        derivative of the radius of the WD in [m/Msol]
    '''
    return waveform_LO.DMcm2(q)*DWD_RadiusMass_Relation(waveform_LO.m2(Mc,q))

#------------------------------------------------#

def DqWD_RadiusMass_RelationR2(Mc,q):
    '''
    Derivative of the mass-radius relation over q for m2
    input:
        Mc - chirp mass in [Msol]
        q - mass ratio m1/m2, with m1>m2
    output:
        derivative of the radius of the WD in [m]    
    '''
    return waveform_LO.Dqm2(Mc,q)*DWD_RadiusMass_Relation(waveform_LO.m2(Mc,q))

#------------------------------------------------#

def DRLOF1(q):
    '''
    Derivative of the Roche lobe overflow radius of the more massive body in a binary system over q
    input:
        q - mass ratio m1/m2, with m1>m2
    output:
        derivative of the Roche lobe overflow radius in units of the orbital separation
    '''
    return (
        0.49*2/3*q**(-1/3)/(0.6*q**(2/3) + np.log(1+q**(1/3)))
        -(0.49*q**(2/3)/(0.6*q**(2/3) + np.log(1+q**(1/3)))**2)
        *(0.6*2/3*q**(-1/3) + 1/3*q**(-2/3)/(1+q**(1/3)))
            )

#------------------------------------------------#

def DRLOF2(q):
    '''
    Derivative of the Roche lobe overflow radius of the less massive body in a binary system over q
    input:
        q - mass ratio m1/m2, with m1>m2
    output:
        derivative of the Roche lobe overflow radius in units of the orbital separation
    '''
    return (
        0.49*(-2/3)*q**(-5/3)/(0.6*q**(-2/3) + np.log(1+q**(-1/3)))
        -(0.49*q**(-2/3)/(0.6*q**(-2/3) + np.log(1+q**(-1/3)))**2)
        *(0.6*(-2/3)*q**(-5/3) - 1/3*q**(-4/3)/(1+q**(-1/3)))
            )

#------------------------------------------------#

def qunstabilitylimit(Mass1,Mass2):
    '''
    Limit in mass ratio q for unstable mass transfer in a binary system with two WDs
    input:
        Mass1 - mass of the more massive body in [Msol]
        Mass2 - mass of the less massive body in [Msol]
    output:
        limit in mass ratio q=m1/m2, with m1>m2, for unstable mass transfer
        if q>qunstabilitylimit(Mass1,Mass2) the mass transfer is unstable
        if q<qunstabilitylimit(Mass1,Mass2) the mass transfer is stable
    '''
    
    # the formulas required to compute the unstable regime are from this article: https://arxiv.org/pdf/astro-ph/0312577.pdf
    Q=Mass1/Mass2

    csi2=Mass2*DWD_RadiusMass_Relation(Mass2)/WD_RadiusMass_Relation(Mass2)    
    csiRL=(1+1/Q)/3*(2*np.log(1+Q**(-1/3))-Q**(-1/3)/(1+Q**(-1/3)))/(0.6*Q**(-2/3)+np.log(1+Q**(-1/3)))
    k=0.1939*(1.44885-Mass1)**(0.1917)
    DkM1=-0.1917*0.1939*(1.44885-Mass1)**(0.1917-1)
    Lambda=1+2*Mass1*DWD_RadiusMass_Relation(Mass1)/WD_RadiusMass_Relation(Mass1)+Mass1*DkM1/k

    r1=WD_RadiusMass_Relation(Mass1)/(WD_RadiusMass_Relation(Mass2)/RLOF2(Q))

    # q delimiting unstable scenario is given by equation (30)
    return 1+(csi2-csiRL)/2-k*r1**2*(1+1/Q)*Lambda

#------------------------------------------------#

def qstabilitylimit(Mass1,Mass2):
    '''
    Limit in mass ratio q for stable mass transfer in a binary system with two WDs
    input:
        Mass1 - mass of the more massive body in [Msol]
        Mass2 - mass of the less massive body in [Msol]
    output:
        limit in mass ratio q=m1/m2, with m1>m2, for stable mass transfer
        if q>qstabilitylimit(Mass1,Mass2) the mass transfer is unstable
        if q<qstabilitylimit(Mass1,Mass2) the mass transfer is stable
    '''
    
    # the formulas required to compute the unstable regime are from this article: https://arxiv.org/pdf/astro-ph/0312577.pdf
    Q=Mass1/Mass2
    csi2=Mass2*DWD_RadiusMass_Relation(Mass2)/WD_RadiusMass_Relation(Mass2)
    csiRL=(1+1/Q)/3*(2*np.log(1+Q**(-1/3))-Q**(-1/3)/(1+Q**(-1/3)))/(0.6*Q**(-2/3)+np.log(1+Q**(-1/3)))
    rh=0.0883+ 0.04858*np.log10(Q)+0.11489*np.log10(Q)**2-0.020475*np.log10(Q)**3
    
    # q delimiting unstable scenario is given by equation (31)
    return 1+(csi2-csiRL)/2-np.sqrt((1+1/Q)*rh)

######################################################
# CLASS FOR CHECKING ADDITIONAL INFORMATION OF WDB
######################################################

class WDclass:
    # class to compute properties of a WDB
    def __init__(self,
                 source_Mc,
                 source_q,
                 detector_fGWmin,
                 detector_fGWmax,
                 max_measurement_time,
                 detector_orbit_R_cb=8.44e6,
                 detector_orbit_R_sat=2.e7,
                 detector_orbit_period=2.8e4,
                 SolarMass=const.SolarMass,
                 SolarRadius=const.SolarRadius,
                 year=const.year,
                 m_per_Mpc=const.m_per_Mpc,
                 Mpc_per_m=const.Mpc_per_m,
                 minWDMass=minWDMass,
                 maxWDMass=maxWDMass, #ChandrasekharMass
                 MaxNSMass=2.1, #TOVMass 
                 EarthOrbitRadius = const.EarthOrbitRadius,
                 detector="MAGIS",
                 check = True,
                 WDBcheck = True
                 ):
        '''
        GW source inputs:
            source_Mc -- detector-frame chirp mass [Msol]
            source_q -- mass ratio q=m_1/m_2
        Measurement inputs:
            detector_fGWmin -- smallest frequency considered [Hz]
            detector_fGWmax -- largest frequency considered [Hz]
            max_measurement_time -- maximal time of measurement considered [seconds] 
        Detector inputs:
            detector_orbit_R_cb -- radius of the orbit of the center of the baseline around Earth [m]
                                    over which the computation of SNR and Fisher matrix computed
            detector_orbit_R_sat -- radius of the orbit of the satellite around Earth [m], for period 
            detector_orbit_period -- period of the orbit of the satellite around Earth [seconds]
        Parameters and constants for computation:
            SolarMass -- Solar mass [kg]
            SolarRadius -- Solar radius [m]
            year -- Year [s]
            m_per_Mpc -- Mega parsec [m]
            Mpc_per_m -- Mega parsec [m]
            minWDMass -- minimum mass for a WD [Msol]
            maxWDMass -- maximum mass for a WD [Msol], ChandrasekharMass
            MaxNSMass -- maximum mass for a NS [Msol], TOVMass
            EarthOrbitRadius -- radius of the Earth's orbit around the Sun [m]
            detector -- string to choose the detector, currently only "MAGIS" implemented
            check -- boolean Flag. Set to True to check that the RLOF frequency is measurable (with the parameters inputted and the detector)
            WDBcheck -- boolean Flag. Set to True to check that the binary is in fact a WDB
        '''    
        
        self.source_Mc = source_Mc
        self.source_q = source_q
        # compute individual masses
        self.source_M1 = waveform_LO.m1(source_Mc,source_q)
        self.source_M2 = waveform_LO.m2(source_Mc,source_q)

        self.detector_fGWmin = detector_fGWmin
        self.detector_fGWmax = detector_fGWmax
        self.detector_orbit_R_cb = detector_orbit_R_cb
        self.detector_orbit_R_sat = detector_orbit_R_sat
        self.detector_orbit_period = detector_orbit_period
        self.SolarMass = SolarMass 
        self.SolarRadius = SolarRadius 
        self.year = year  
        self.m_per_Mpc = m_per_Mpc
        self.Mpc_per_m = Mpc_per_m  
        self.minWDMass = minWDMass
        self.maxWDMass = maxWDMass
        self.maxNSMass = MaxNSMass
        self.max_measurement_time=max_measurement_time
        self.EarthOrbitRadius=EarthOrbitRadius
        self.check = check
        self.WDBcheck = WDBcheck

        # initialize sensitivity curve of the detector
        self.NoiseSpace = np.loadtxt('Analysis_Scripts/Instruments_Sensitivity/NoiseCurve_Space.dat', usecols=(0,1))
        self.minFrequencySpaceNoise=min(self.NoiseSpace[:,0])
        if detector=="MAGIS":
            self.Noise=self.NoiseSpace
            self.minFrequencyDetector=self.minFrequencySpaceNoise
        else:  
            print("WARNING:\n The detector chosen is not valid. The MAGIS detector will be used.")
            self.Noise=self.NoiseSpace
            self.minFrequencyDetector=self.minFrequencySpaceNoise

        # compute cutoff frequency and check if the binary is a WDB
        SingularCriticalPointX=((1+2*np.pi*self.EarthOrbitRadius/(phconst.c*self.year))/(1-2*np.pi*self.EarthOrbitRadius/(phconst.c*self.year)))
        self.SunSingularCriticalPointTime=(self.year/2+self.EarthOrbitRadius/phconst.c*((SingularCriticalPointX)**(8/3)+1))/((SingularCriticalPointX)**(8/3)-1)
        self.DetectorSingularCriticalPointTime=3/8*((phconst.c-2*np.pi*self.EarthOrbitRadius/self.year)/(4*np.pi**2*self.detector_orbit_R_sat/(self.detector_orbit_period**2)))

    ###################################################
    # Class Functions
    ###################################################

    def CutoffFrequency_func(self):
        '''
        Function to compute the cutoff frequency due to Roche lobe overflow of the less massive body in a binary system
        and check if the masses inputted are valid for a WD binary system
        input:
            self.source_M1 -- mass of the more massive body [Msol]
            self.source_M2 -- mass of the less massive body [Msol]
            self.detector_fGWmin -- smallest frequency considered [Hz]
            self.detector_fGWmax -- largest frequency considered [Hz]
            self.minWDMass -- minimum mass for a WD [Msol]
            self.maxWDMass -- maximum mass for a WD [Msol], ChandrasekharMass
            self.maxNSMass -- maximum mass for a NS [Msol], TOVMass
            self.minFrequencyDetector -- minimal frequency the detector is sensitive to [Hz]
        output:
            Note: the function does not return any value, but sets the attributes of the class
            1) self.CutoffFrequency
            2) self.Binarytype
            3) self.WD_detector_fGWmax
            4) self.check
            5) self.WDBcheck

            Warning messages are printed if the masses are not valid for a WD binary system
            1) if at least one mass is not a WD mass, the possible object types are printed
            2) if both masses are not WD masses, no cutoff frequency is computed
            3) if the chosen minimal frequency is higher than maximal one, self.check is set to False
            4) if the maximal frequency is lower than the sensibility of the interferometer, self.check is set to False
        '''

        # Order the masses, so that M1>M2
        if self.source_M2>self.source_M1:
            smallM=self.source_M1
            self.source_M1=self.source_M2
            self.source_M2=smallM
        
        # Check if M2 is a WD mass
        if self.source_M2>=self.minWDMass and self.source_M2<=self.maxWDMass:
            # If M2 is a valid WD mass, then the frequency cutoff will be given by RLOF: compute it
            self.CutoffFrequency=WD_RLOFCutoffFrequency(self.source_M1,self.source_M2)
            # Check if the binary is a WD binary
            if self.source_M1>=self.minWDMass and self.source_M1<=self.maxWDMass:
                self.Binarytype = "WD-WD"
                self.WDBcheck = True
            # If the second mass is not a WD mass, print warnings and possible object types
            # The second mass can only be a WD, NS or BH (or other exotic very compact objects), otherwise the binary would have merged at lower frequencies 
            elif self.source_M1>self.maxWDMass and self.source_M1<=self.maxNSMass:
                print("The heavier body (", "%.2f"% self.source_M1,") is too massive to be a WD, but is in the NS mass range.")
                print("This scenario can still have the distruption of the WD, the less massive body (", "%.2f"% self.source_M2,"), as cutoff-frequency.")
                self.Binarytype = "WD-NS"
                self.WDBcheck = False
            elif self.source_M1>self.maxNSMass:
                print("The heavier body (", "%.2f"% self.source_M1,") is too massive to be a WD, but can be a BH.")
                print("This scenario can still have the distruption of the WD, the less massive body (", "%.2f"% self.source_M2,"), as cutoff-frequency.")
                self.Binarytype = "WD-BH"
                self.WDBcheck = False
        # If M2 is lighter than a WD mass, print warnings and possible object types
        # This M2 can only be a PBH, because have mass too small to be a NS or BH, and if it is not a WD, nor NS, nor BH, it must be a PBH (or other exotic very compact objects beyond the standard astrophysical objects)
        elif self.source_M2<self.minWDMass:
            print("WARNING:\n  The lighter body (", "%.2f"% self.source_M2,") does not have enough mass to be a WD. For this merger to not have happened at lower frequencies this body must be a PBH.")
            # Identify the possible object types for the more massive body
            if self.source_M1<self.minWDMass:
                print("The heavier body (", "%.2f"% self.source_M1,") also does not have enough mass to be a WD, hence must be a PBH for this binary to not have merged yet.")
                print("This scenario has no cutoff-frequency.")
                self.CutoffFrequency=self.detector_fGWmax
                self.Binarytype = "PBH-PBH"
                self.WDBcheck = False
            elif self.source_M1<=self.maxWDMass and self.source_M2<self.minWDMass:
                print("This scenario has the distruption of the WD, the more massive body (", "%.2g"% self.source_M1,"), as cutoff-frequency.")
                self.CutoffFrequency=WD_RLOFCutoffFrequency(self.source_M2,self.source_M1)
                self.Binarytype = "PBH-WD"
                self.WDBcheck = False
            elif self.source_M1>self.maxWDMass and self.source_M1<=self.maxNSMass:
                print("The heavier body (", "%.2f"% self.source_M1,") is too massive to be a WD, but is in the NS mass range.")
                print("This scenario has no cutoff-frequency.")
                self.CutoffFrequency=self.detector_fGWmax
                self.Binarytype = "PBH-NS"
                self.WDBcheck = False
            elif self.source_M1>self.maxNSMass:
                print("The heavier body (", "%.2f"% self.source_M1,") is too massive to be a WD, but can be a BH.")
                print("This scenario has no cutoff-frequency.")
                self.CutoffFrequency=self.detector_fGWmax
                self.Binarytype = "PBH-BH"
                self.WDBcheck = False
        # If M2 is heavier than a WD mass, print warnings and possible object types
        # This M2 can only be a NS or BH
        elif self.source_M2>self.maxWDMass:
            print("WARNING:\n  The lighter body (", "%.2f"% self.source_M2,") has too much mass to be a WD. For this merger to not have happened at lower frequencies this body must either a NS or a BH, and so must be the more massive body (", "%.2g"% self.source_M1,").")
            print("This scenario has no cutoff-frequency.")
            self.CutoffFrequency=self.detector_fGWmax
            if self.source_M2<self.maxNSMass:
                if self.source_M1>self.maxWDMass and self.source_M1<=self.maxNSMass:
                    self.Binarytype = "NS-NS"
                    self.WDBcheck = False
                elif self.source_M1>self.maxNSMass:
                    self.Binarytype = "NS-BH"
                    self.WDBcheck = False
            elif self.source_M2>self.maxNSMass and self.source_M1>self.maxNSMass:
                    self.Binarytype = "BH-BH"
                    self.WDBcheck = False
        
        # Compute the maximal frequency to be considered for the WD binary system
        self.WD_detector_fGWmax=min(self.CutoffFrequency,self.detector_fGWmax)
        
        # Check if the frequencies chosen allow a correct evaluation
        # If not, set self.check to False and print warnings
        if self.detector_fGWmin >= self.WD_detector_fGWmax:
            self.check = False
            print("WARNING:\n  The chosen minimal frequency is higher than maximal one.")
        if self.WD_detector_fGWmax <= self.minFrequencyDetector:
            self.check = False
            print("WARNING:\n  The maximal frequency is lower than the sensibility of the interferometer.")
        return

    #------------------------------------------------#

    def INFO(self):
        '''
        Function to print additional information about the WD binary system and the measurement
        First, it is checked that the function CutoffFrequency_func has been run, if not it is run

        input:
            self.source_M1 -- mass of the more massive body [Msol]
            self.source_M2 -- mass of the less massive body [Msol]
            self.source_Mc -- detector-frame chirp mass [Msol]
            self.detector_fGWmin -- smallest frequency considered [Hz]
            self.detector_fGWmax -- largest frequency considered [Hz]
            self.WD_detector_fGWmax -- largest frequency considered [Hz], min(self.CutoffFrequency,self.detector_fGWmax)
            self.max_measurement_time -- maximal time of measurement considered [seconds]
        output:
            prints information about the WD binary system and the measurement
            1) Type of binary
            2) Masses of the two bodies
            3) Cutoff Frequency either due to RLOF or the chosen maximal frequency
            4) Frequencies between which the evaluation is done
            5) If the frequencies allow a correct evaluation, the measurement duration is printed, else a warning is printed
            6) Frequencies at which multiple critical points could appear 

            Warning messages are printed
            1) if the choosen maximal frequency is lower than the cutoff
            2) if the time between cutoff-frequency and minimal given frequency is less than the measurement time
            3) if multiple critical points could appear in the frequency range considered
        '''
        
        print("\n-------------------------------------------")
        # Check if CutoffFrequency_func has been run, if not run it
        try:
            self.WD_detector_fGWmax
        except:
            self.CutoffFrequency_func()

        # Print the information about the WD binary system and the measurement time and frequency range
        print("Type of binary: ", self.Binarytype)
        minfGW=max(self.detector_fGWmin,waveform_LO.fGW(self.source_Mc,
            waveform_LO.tau_fGW(self.source_Mc, self.WD_detector_fGWmax) + self.max_measurement_time))
        print("The masses of the two bodies is ", "%.3g"% self.source_M2, " and ", "%.3g"% self.source_M1, " Solar Masses")
        # Print the cutoff frequency due to RLOF if it is lower than the chosen maximal frequency
        if self.CutoffFrequency < self.detector_fGWmax:
            print("Cutoff Frequency due to RLOF: ",  "%.3g"% self.CutoffFrequency, " Hz")
        elif self.CutoffFrequency == self.detector_fGWmax:
            print("No cutoff for WD has been computed hence the chosen one is used: ",  "%.3g"% self.detector_fGWmax, " Hz")
        # Print a warning if the chosen maximal frequency is lower than the cutoff
        else:
            print("WARNING:\n  The choosen maximal frequency is lower than the cutoff ", "%.3g"% self.CutoffFrequency)
        # Print the frequencies between which the evaluation is done
        print("Evaluation for frequencies between ", "%.3g"% minfGW, " and ", "%.3g"% self.WD_detector_fGWmax,   " Hz")

        # Print the measurement duration if the frequencies allow a correct evaluation, else print a warning
        if self.check:
            if minfGW == self.detector_fGWmin:
                print("WARNING:\n  Time between cutoff-frequency and minimal given frequency is less than ", "%.3g"% (self.max_measurement_time/self.year), " yr")
                print("In fact, the binary covers the above frequency range in only ", "%.3g"% ((waveform_LO.tau_fGW(self.source_Mc, minfGW)-waveform_LO.tau_fGW(self.source_Mc, self.CutoffFrequency))/self.year), " yr")
            else:
                print("Measurement duration is ", "%.3g"% (self.max_measurement_time/self.year), " yr")
        else:
            print("WARNING:\n  The chosen parameters do not allow a correct evaluation (see previous warnings).")
        
        # Print the frequencies at which multiple critical points could appear
        # If these frequencies are in the range which lead to multiple critical points, print a warning
        print("Multiple critical point could appear at:")
        SingularCriticalPointFrequency=(waveform_LO.fGW(self.source_Mc,self.SunSingularCriticalPointTime))
        if minfGW < SingularCriticalPointFrequency:
            print("- WARNING:\n  From Earth rotation around the Sun from frequency ", "%.3g"% SingularCriticalPointFrequency, " Hz and above")
        else:
            print("- From Earth rotation around the Sun from frequency ", "%.3g"% SingularCriticalPointFrequency, " Hz  and above")

        SingularCriticalPointFrequency2=(waveform_LO.fGW(self.source_Mc,self.DetectorSingularCriticalPointTime))
        if minfGW < SingularCriticalPointFrequency2:
            print("- WARNING:\n  From Satellite rotation around Earth from frequency ", "%.3g"% SingularCriticalPointFrequency2, " Hz and above")
        else:
            print("- From Satellite rotation around Earth from frequency ", "%.3g"% SingularCriticalPointFrequency2, " Hz  and above")

        print("-------------------------------------------\n")
        return

    #------------------------------------------------#  