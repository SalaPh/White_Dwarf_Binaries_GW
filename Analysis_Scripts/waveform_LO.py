###################################################
# leading order waveform computation following Maggiore 2008, chapter 4.1
###################################################

import numpy as np
import scipy.constants as phconst

from . import constants as const

###################################################
# FUNCTIONS
###################################################

def Mc(m1, m2):
   '''
   returns the chirp mass knowing both binary masses, Eq. 4.2 of Maggiore
   input:
      m1 -- binary larger mass [Msol]
      m2 -- binary smaller mass [Msol]
   output:
      Mc -- detector frame chirp mass [Msol]
   '''
   return (m1*m2)**(3/5)/(m1+m2)**(1/5)

#------------------------------------------------#

def m1(Mc, q):
   '''
   returns the larger mass of the binary, could be found from Eq. 4.2 of Maggiore
   input:
      Mc -- detector frame chirp mass [Msol]
      q -- binary masses relation (m1/m2)
   output:
      m1 -- binary larger mass [Msol]
   '''
   return Mc*q**0.4*(1+q)**0.2

#------------------------------------------------#

def m2(Mc, q):
   '''
   returns the larger mass of the binary, could be found from Eq. 4.2 of Maggiore
   input:
      Mc -- detector frame chirp mass [Msol]
      q -- binary masses relation (m1/m2)
   output:
      m2 -- binary smaller mass [Msol]
   '''
   return Mc*q**-0.6*(1+q)**0.2

#------------------------------------------------#

def mu(m1, m2):
    '''
    mass parameter defined in Eq. (11) of astro-ph/0312577
     input:
        m1 - mass of the accretor in [Msol]
        m2 - mass of the donor in [Msol]
    output:
        mass parameter mu
    '''
    return m2/(m1+m2)

#------------------------------------------------#

def GWPhase(Mc, Theta0, tau):
   '''
   Returns the phase of the GW signal in radians, using Eq. 4.193 of Maggiore
   input:
      Mc -- detector frame chirp mass [Msol]
      Theta0 - GW phase at coalesence
      tau -- time to coalesence [seconds] 
   output:
      phase -- phase of the GW signal [rad]
   '''
   return -2.*(5.*(phconst.G*const.SolarMass)*Mc/phconst.c**3)**(-5./8.)*tau**(5./8.)+Theta0

#------------------------------------------------#

def fGW(Mc, tau):
   '''
   Returns the frequency of the GW signal in Hz, using Eq. 4.195 of Maggiore
   input:
      Mc -- detector frame chirp mass [Msol]
      tau -- time to coalesence [seconds] 
   output:
      fGW -- frequency of the GW signal [Hz]
   '''
   return 1./np.pi*(5./256./tau)**(3./8.)*((phconst.G*const.SolarMass)*Mc/phconst.c**3)**(-5./8.)

#------------------------------------------------#

def tau_fGW(Mc, fGW):
   '''
   Returns the time before merger in [sec] at which the GW signal has the 
   frequency fGW in Hz, using Eq. 4.195 of Maggiore
   input:
      Mc -- detector frame chirp mass [Msol]
      fGW -- frequency of the GW signal [Hz] 
   output:
      tau -- time to coalesence [seconds]
   '''
   return 5./256. * (np.pi*fGW)**(-8./3.)*(phconst.c**3/(phconst.G*const.SolarMass)/Mc)**(5./3.)

#------------------------------------------------#

def Porb(m1, m2, a):
    '''
    orbital period of the binary, from Kepler's third law
    input:
        m1 - mass of the accretor in [Msol]
        m2 - mass of the donor in [Msol]
        a - orbital separation in [m]
    output:
        orbital period of the binary in [s]
    '''

    return 2.*np.pi*np.sqrt(a**3 / (phconst.G*const.SolarMass) / (m1+m2))

#------------------------------------------------#

def forb(M1, M2, a):
    '''
    Simple computation of the ORBITAL (not GW) frequency given the masses and distances at a certain time, and considering the orbit to still be circular
    input:
        M1 - mass of the accretor in [Msol]
        M2 - mass of the donor in [Msol]
        a - orbital separation in [m]
    output:
        f in [Hz]
    '''
    
    return np.sqrt((phconst.G*const.SolarMass)*(M1+M2)/(a**3))/(2*np.pi)

#------------------------------------------------#

def Bodyhf(fGW,Mc,dL):
   '''
   returns the common body of the strain for the GW signal
   input:
      fGW -- GW frequency [Hz]
      Mc -- detector frame chirp mass [Msol]
      dL -- luminosity distance in Mpc
   output:
      hfac -- common body of the strain for the GW signal
   '''
   return (2*Mc**(5/3)*np.pi**(2/3)*fGW**(2/3)/dL*const.Mpc_per_m)*((phconst.G*const.SolarMass)**(5/3)/phconst.c**4)

#------------------------------------------------#

def get_hp_hc(Mc, iota, dL, Theta0, tau):
   '''
   return a tuple of plus-polarization, cross polarization strain signals. 
   input:
      Mc -- detector frame chirp mass [Msol]
      iota -- angle between the sky-localization and the orbital momentum of the binary [rad]
      dL -- luminosity distance in Mpc
      Theta0 - GW phase at coalesence
      tau -- time to coalesence [seconds] 
   output:
      (hp, hc) -- tuple of plus-polarization, cross polarization strain signals
   '''
   fGW_vec = fGW(Mc, tau)
   phase_vec = GWPhase(Mc, Theta0, tau)
   hfac_vec = Bodyhf(fGW_vec, Mc, dL)
   return hfac_vec*(1.+np.cos(iota)**2)*np.cos(phase_vec), hfac_vec*np.cos(iota)*np.sin(phase_vec)

#------------------------------------------------#

def Distance_Frequency(Frequency,Mass1,Mass2):
   '''
   returns the Distance of a binary system given its GW frequency and masses
   input:
      Frequency -- GW frequency [Hz]
      Mass1 -- mass of the first body [Msol]
      Mass2 -- mass of the second body [Msol]
   output:
      Distance -- Distance of the binary system [m]
   '''
   return (phconst.G*((Mass1+Mass2)*const.SolarMass)/(np.pi**2*(Frequency)**2))**(1/3)

#------------------------------------------------#

def Frequency_Distance(Distance,Mass1,Mass2):#GW frequency =  2*orbital frequency
   '''
   returns the GW frequency of a binary system given its Distance and masses
   input:
      Distance -- Distance of the binary system [m]
      Mass1 -- mass of the first body [Msol]
      Mass2 -- mass of the second body [Msol]
   output:
      Frequency -- GW frequency [Hz]
   '''
   return np.sqrt(phconst.G*(Mass1+Mass2)*const.SolarMass)/(np.pi*Distance**(3/2))

#------------------------------------------------#
#Derivatives of a selection of elements, which are then needed in some computations

def dfGWdt(fGW,Mc):
   '''
   returns the derivative of the GW frequency with respect to time
   input:
      fGW -- GW frequency [Hz]
      Mc -- detector frame chirp mass [Msol]
   output:
      derivative of the GW frequency with respect to time [Hz/s]
   '''
   return (96/5)*np.pi**(8/3)*Mc**(5/3)*abs(fGW)**(11/3)*((phconst.G*const.SolarMass)**(5/3)/phconst.c**5)

#------------------------------------------------#
def DMcm2(q):
   '''
   derivative over Mc of M2
   input:
      q -- binary masses relation (m1/m2)
   output:
      derivative of m2 over Mc
   '''
   return q**-0.6*(1+q)**0.2

#------------------------------------------------#

def Dqm2(Mc, q):
   '''
   derivative over q of M2
   input:
      Mc -- detector frame chirp mass [Msol]
      q -- binary masses relation (m1/m2)
   output:
      derivative of m2 over q [Msol]
   '''
   return Mc*(-0.6*q**-1.6*(1+q)**0.2+0.2*q**-0.6*q**-0.8)

#------------------------------------------------#
