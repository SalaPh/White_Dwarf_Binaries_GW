############################################################################
# Main module for computing SNR and Fisher matrix of GW signals from inspiraling binaries
# as seen by a space-based detector consisting of satellites orbiting Earth.
#
# The code uses the PyCBC package to generate the frequency-domain waveforms
# and the SciPy package for interpolation and root finding.
############################################################################

import numpy as np
import time
import sys

from scipy.optimize import root_scalar
from scipy import interpolate

# import pycbc
import pycbc.waveform
import pycbc.types

# import my modules
from Analysis_Scripts import helper_funs
from Analysis_Scripts import waveform_LO
from Analysis_Scripts import waveform_PN
from Analysis_Scripts import antennaFuns_satellites

###################################################
# import the noise curve data

noise_curve_data = np.loadtxt('Analysis_Scripts/Instruments_Sensitivity/NoiseCurve_Space.dat').T
noise_fun_PSD_Hz = interpolate.interp1d(
   noise_curve_data[0], 
   noise_curve_data[1]**2, 
   bounds_error=False, 
   fill_value = np.inf
   )

###################################################
# FUNCTIONS
###################################################

def fGW_tau_PN(tau, nu, m, PN_order):
   '''
   GW frequency (in [Hz]) of the signal
   at time tau before the merger, using the PN expressions
   from the waveform_PN module. 
   input:
      tau -- time before merger [sec]
      nu -- symmetric mass ratio of the binary
      m -- total mass of the binary [Msol]
      PN_order -- order of the PN expression
   output:
      fGW -- GW frequency at time tau before merger [Hz]
   '''
   theta = waveform_PN.theta(-tau, nu, m)
   x = waveform_PN.x_PN(theta, nu, PN_order)
   return waveform_PN.c_m_s**3/(waveform_PN.G_m3_Msol_s2*m*np.pi)*x**1.5

#------------------------------------------------#


def fGW_tau_PN_root(tau, fGWtarget, nu, m, PN_order):
   '''
   Same as fGW_tau_PN, except that fGWtarget is subtracted
   such that this function can be used in a root finder
   to solve for tau. 
   input:
      tau -- time before merger [sec]
      fGWtarget -- target frequency to subtract [Hz]
      nu -- symmetric mass ratio of the binary
      m -- total mass of the binary [Msol]
      PN_order -- order of the PN expression
   output:
      fGW - fGWtarget [Hz]
   '''
   return fGW_tau_PN(tau, nu, m, PN_order) - fGWtarget

#------------------------------------------------#

def get_htilde(
   source_Mc,           
   source_q,            
   source_iota,            
   source_phi0,            
   source_tc,           
   source_dL,           
   source_RA,           
   source_DEC,          
   source_psi,          
   detector_orbit_R,          
   detector_orbit_period,           
   detector_orbit_t0,            
   detector_orbit_RA0,           
   detector_orbit_DEC0,
   eval_fvec,
   clock_tD,          
   PN_order_phase,            
   PN_order_amplitude,            
   fast_tau_evaluation   
   ):
   '''
   Returns array of frequency domain waveform [1/Hz])
   GW source input:
      source_Mc -- detector-frame chirp mass [Msol]
      source_q -- mass ratio q=m_1/m_2
      source_iota -- angle between sky localization and orbital momentum 
                     of the binary [rad]
      source_phi0 -- reference phase of the GW signal
      source_tc -- time of merger, measured from the solar equinox [seconds]
      source_dL -- luminosity distance [Mpc]
      source_RA --  right ascension [rad]
      source_DEC -- declination [rad]
      source_psi -- polarization angle [rad]
   Detector input:
      detector_orbit_R -- radius of the orbit of the satellites around 
                          Earth [m]
      detector_orbit_period -- period of the detector around Earth [seconds]
      detector_orbit_t0 -- reference time for fixing orbit of satellite 
                           around Earth [seconds]
      detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
      detector_orbit_DEC0 -- declination of the satellite at t0 [rad]
   Parameters for computation:
      eval_fvec -- array of frequencies for which htilde is evaluated [Hz]
      clock_tD -- reference time for clock location
      PN_order_phase -- order of PN expansion of phase
      PN_order_amplitude -- order of PN expansion of waveform
      fast_tau_evaluation -- Boolean flag. Set to True to use leading-order
                             expression for calculation of tau(f_GW), 
                             if False the PN expressions will be used
   output:
      htilde -- array of frequency-domain waveform [1/Hz]
   '''
   # compute mass params
   source_m1 = source_Mc*source_q**0.4*(1+source_q)**0.2
   source_m2 = source_Mc*source_q**-0.6*(1+source_q)**0.2
   source_m = source_Mc*(1.+source_q)**1.2/source_q**0.6
   source_nu = source_q/(1.+source_q)**2
   #  polarization basis waveforms
   hptilde, hctilde = pycbc.waveform.get_fd_waveform_sequence(
      approximant = 'TaylorF2',
      phase_order = 2*PN_order_phase,
      amplitude_order = 2*PN_order_amplitude,
      mass1 = source_m1,
      mass2 = source_m2,
      inclination = source_iota,
      coa_phase = source_phi0,
      distance = source_dL,
      sample_points = pycbc.types.array.Array(eval_fvec)
      )
   # compute the time-before-merger tau for the frequencies
   # from the leading order time-frequency relation
   tau_vec = waveform_LO.tau_fGW(source_Mc, eval_fvec)
   # Improve tau(fGW) computation with PN equations if asked for
   if fast_tau_evaluation==False:
      for i in range(len(tau_vec)):
         tau_low = tau_vec[i]*0.9
         tau_high = tau_vec[i]*1.1
         while (
            fGW_tau_PN_root(
               tau_low, 
               eval_fvec[i], 
               source_nu, 
               source_m, 
               PN_order_phase
               )
            * fGW_tau_PN_root(
               tau_high, 
               eval_fvec[i], 
               source_nu, 
               source_m, 
               PN_order_phase
               )
            ) > 0:
            tau_low *= 0.9
            tau_high *=1.1
         tau_vec[i] = root_scalar(
            fGW_tau_PN_root,
            args=(eval_fvec[i], source_nu, source_m, PN_order_phase),
            bracket=[tau_low, tau_high]
            ).root
   # get antenna functions and time delays as functions of tau/time
   timeDelay, Fp, Fc = antennaFuns_satellites.timeDelay_antennaFuns(
      source_tc - tau_vec,
      source_RA,
      source_DEC,
      source_psi,
      detector_orbit_R,
      detector_orbit_t0,
      detector_orbit_RA0,
      detector_orbit_DEC0,
      clock_tD,
      detector_period = detector_orbit_period
      )
   # and return the frequency-domain waveform function
   return (
      np.exp(2.*np.pi*1j*np.mod(eval_fvec*source_tc, 1.)) # t_c phase
      * np.exp(2.*np.pi*1j*np.mod(timeDelay*eval_fvec, 1.)) # Doppler Phase
      * (Fp*np.array(hptilde) + Fc*np.array(hctilde))
      )

#------------------------------------------------#

def SNR_FisherMatrix(
   source_Mc,
   source_q,
   source_iota,
   source_phi0,
   source_tc,
   source_dL,
   source_RA,
   source_DEC,
   source_psi,
   detector_orbit_R,
   detector_orbit_period,
   detector_orbit_t0,
   detector_orbit_RA0,
   detector_orbit_DEC0,
   eval_fvec,
   clock_tD,
   PN_order_phase,
   PN_order_amplitude,
   Delta_Mc_rel,
   Delta_q_rel,
   Delta_iota_abs,
   Delta_phi0_abs,
   Delta_tc_abs,
   Delta_dL_rel,
   Delta_angleSL_abs,
   Delta_psi_abs,
   fast_tau_evaluation
   ):
   '''
   Big wrapper function which returns the 
   SNR^2 and (upper triangle of the) FisherMatrix of the GW signal
   GW source input:
      source_Mc -- detector-frame chirp mass [Msol]
      source_q -- mass ratio q=m_1/m_2
      source_iota -- angle between sky localization and orbital momentum 
                     of the binary [rad]
      source_phi0 -- reference phase of the GW signal at fGWmin
      source_tc -- time of merger, measured from the solar equinox [seconds]
      source_dL -- luminosity distance [Mpc]
      source_RA -- right ascension [rad]
      source_DEC -- declination [rad]
      source_psi -- polarization angle [rad]
   Detector input:
      detector_orbit_R -- radius of the orbit of the satellites around 
                          Earth [m]
      detector_orbit_period -- period of the detector around Earth [seconds]
      detector_orbit_t0 -- reference time for fixing orbit of satellite 
                           around Earth [seconds]
      detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
      detector_orbit_DEC0 -- declination of the satellite at t0 [rad]
   Parameters for computation:
      eval_fvec -- array of frequencies for which htilde is evaluated [Hz]
      clock_tD -- reference time for clock location
      PN_order_phase -- order of PN expansion of phase
      PN_order_amplitude -- order of PN expansion of amplitude
      Delta_Mc_rel -- relative change of chirp mass parameter for finite
                      difference computation of the Fisher matrix
      Delta_q_rel -- relative change of q parameter for the finite
                     difference computation of the Fisher Matrix
      Delta_iota_abs -- absolute change of iota parameter for finite
                        difference computation of the Fisher matrix [rad]
      Delta_phi0_abs -- ansolute change of phi0 parameter for finite
                        difference computation of the Fisher matrix [rad]
      Delta_tc_abs -- absolute change of t_c parameter [seconds] for finite
                      difference computation of the Fisher matrix
      Delta_dL_rel -- relative change of dL parameter for finite
                      difference computation of the Fisher matrix
      Delta_angleSL_abs -- absolute change in sky localization angle 
                           parameters for finite difference computation 
                           of Fisher matrix [rad]
      Delta_psi_abs -- absolute change of psi parameter for finite
                       difference computation of the Fisher matrix [rad]
      fast_tau_evaluation -- boolean Flag. Set to True to use leading-order
                             expression for calculation of tau(f_GW), 
                             if False the PN expressions will be used
   output:
      SNR2 -- signal-to-noise ratio squared of the GW signal  
      FisherMat -- upper triangle of the Fisher matrix of the GW signal
   '''
   # make array of source parameters
   source_params = np.array([
      source_Mc,
      source_q,
      source_iota,
      source_phi0,
      source_tc,
      source_dL,
      source_RA,
      source_DEC,
      source_psi
      ])
   # make an array of the stops for each finite difference computation
   source_params_steps = np.array([
      source_Mc*Delta_Mc_rel,
      source_q*Delta_q_rel,
      Delta_iota_abs,
      Delta_phi0_abs,
      Delta_tc_abs,
      source_dL*Delta_dL_rel,
      Delta_angleSL_abs,
      Delta_angleSL_abs,
      Delta_psi_abs
      ])
   # make a new get_htilde with fixed detector and computational params
   func_htilde = lambda params: get_htilde(
      params[0],           
      params[1],            
      params[2],            
      params[3],            
      params[4],           
      params[5],           
      params[6],           
      params[7],          
      params[8],          
      detector_orbit_R,          
      detector_orbit_period,           
      detector_orbit_t0,            
      detector_orbit_RA0,           
      detector_orbit_DEC0,
      eval_fvec,
      clock_tD,          
      PN_order_phase,            
      PN_order_amplitude,            
      fast_tau_evaluation = fast_tau_evaluation          
      )
   # get the reference signal
   htilde_ref = func_htilde(source_params)
   # get the derivatives
   Dhtilde_list = []
   for i in range(len(source_params)):
      temp_steps = np.zeros(source_params.size)
      temp_steps[i] = source_params_steps[i]
      Dhtilde_list.append(
         helper_funs.numDerivative_5s(
            func_htilde, 
            source_params, 
            temp_steps
            )
         )
   # get noise curve
   noise_curve = noise_fun_PSD_Hz(eval_fvec)
   # get SNR
   SNR2 = helper_funs.inner_product(
      htilde_ref,
      htilde_ref, 
      noise_curve, 
      eval_fvec
      )
   # get entries of Fisher matrix
   # create output array for Fisher matrix
   FisherMat = np.zeros((
      len(source_params),
      len(source_params)
      ))
   # fill upper triangle of Fisher matrix
   for i in range(FisherMat.shape[0]):
      for j in range(i,FisherMat.shape[1]):
         FisherMat[i,j] = helper_funs.inner_product(
            Dhtilde_list[i],
            Dhtilde_list[j],
            noise_curve,
            eval_fvec
            )
   return SNR2, FisherMat

#------------------------------------------------#

def SNR2(
   source_Mc,
   source_q,
   source_iota,
   source_phi0,
   source_tc,
   source_dL,
   source_RA,
   source_DEC,
   source_psi,
   detector_orbit_R,
   detector_orbit_period,
   detector_orbit_t0,
   detector_orbit_RA0,
   detector_orbit_DEC0,
   eval_fvec,
   clock_tD,
   PN_order_phase,
   PN_order_amplitude,
   fast_tau_evaluation
   ):
   '''
   Big wrapper function which returns the 
   SNR^2 only of the GW signal
   GW source input:
      source_Mc -- detector-frame chirp mass [Msol]
      source_q -- mass ratio q=m_1/m_2
      source_iota -- angle between sky localization and orbital momentum 
                     of the binary [rad]
      source_phi0 -- reference phase of the GW signal at fGWmin
      source_tc -- time of merger, measured from the solar equinox [seconds]
      source_dL -- luminosity distance [Mpc]
      source_RA -- right ascension [rad]
      source_DEC -- declination [rad]
      source_psi -- polarization angle [rad]
   Detector input:
      detector_orbit_R -- radius of the orbit of the satellites around 
                          Earth [m]
      detector_orbit_period -- period of the detector around Earth [seconds]
      detector_orbit_t0 -- reference time for fixing orbit of satellite 
                           around Earth [seconds]
      detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
      detector_orbit_DEC0 -- declination of the satellite at t0 [rad]
   Parameters for computation:
      eval_fvec -- array of frequencies for which htilde is evaluated [Hz]
      clock_tD -- reference time for clock location
      PN_order_phase -- order of PN expansion of phase
      PN_order_amplitude -- order of PN expansion of amplitude
      fast_tau_evaluation -- boolean Flag. Set to True to use leading-order
                             expression for calculation of tau(f_GW), 
                             if False the PN expressions will be used
   output:
      SNR2 -- signal-to-noise ratio squared of the GW signal
   '''
   # make array of source parameters
   source_params = np.array([
      source_Mc,
      source_q,
      source_iota,
      source_phi0,
      source_tc,
      source_dL,
      source_RA,
      source_DEC,
      source_psi
      ])

   # make a new get_htilde with fixed detector and computational params
   func_htilde = lambda params: get_htilde(
      params[0],           
      params[1],            
      params[2],            
      params[3],            
      params[4],           
      params[5],           
      params[6],           
      params[7],          
      params[8],          
      detector_orbit_R,          
      detector_orbit_period,           
      detector_orbit_t0,            
      detector_orbit_RA0,           
      detector_orbit_DEC0,
      eval_fvec,
      clock_tD,          
      PN_order_phase,            
      PN_order_amplitude,            
      fast_tau_evaluation = fast_tau_evaluation          
      )
   # get the reference signal
   htilde_ref = func_htilde(source_params)
   # get noise curve
   noise_curve = noise_fun_PSD_Hz(eval_fvec)
   # get SNR
   SNR2 = helper_funs.inner_product(
      htilde_ref,
      htilde_ref, 
      noise_curve, 
      eval_fvec
      )
   return SNR2

######################################################
# CLASS FOR GW EVENT
######################################################

class GW_event:
   # GW_event class which wraps all the previous functions in one big analysis
   def __init__(
      self,
      source_Mc,
      source_q,
      source_iota,
      source_phi0,
      source_tc,
      source_dL,
      source_RA,
      source_DEC,
      source_psi,
      detector_fGWmin,
      detector_fGWmax,
      detector_orbit_R,
      detector_orbit_period,
      detector_orbit_t0,
      detector_orbit_RA0,
      detector_orbit_DEC0,
      PN_order_phase,
      PN_order_amplitude,
      min_waveform_length = 1e4,
      max_waveform_length = 1e7,
      max_time_spacing = np.nan,
      max_measurement_time = 3.155814954e7,
      Delta_Mc_rel = 1e-12,
      Delta_q_rel = 1e-5,
      Delta_iota_abs = 1e-6,
      Delta_phi0_abs = 1e-5,
      Delta_tc_abs = 1e-6,
      Delta_dL_rel = 1e-4,
      Delta_angleSL_abs = 1e-8,
      Delta_psi_abs = 1e-5,
      Delta_autoadjust = True,
      fast_tau_evaluation = True,
      verbose = False,
     
      ):
      '''
      GW source input:
         source_Mc -- detector-frame chirp mass [Msol]
         source_q -- mass ratio q=m_1/m_2
         source_iota -- angle between sky localization and orbital momentum 
                        of the binary [rad]
         source_phi0 -- reference phase of the GW signal
         source_tc -- time of merger, measured from the solar equinox [seconds]
         source_dL -- luminosity distance [Mpc]
         source_RA --  right ascension [rad]
         source_DEC -- declination [rad]
         source_psi -- polarization angle [rad]
      Detector input:
         detector_fGWmin -- smallest frequency considered [Hz]
         detector_fGWmax -- smallest frequency considered [Hz]
         detector_orbit_R -- radius of the orbit of the satellites around 
                             Earth [m]
         detector_orbit_period -- period of the detector around Earth [seconds]
         detector_orbit_t0 -- reference time for fixing orbit of satellite 
                              around Earth [seconds]
         detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
         detector_orbit_DEC0 -- declination of the satellite at t0 [rad]
      Parameters for computation:
         PN_order_phase -- order of PN expansion of phase
         PN_order_amplitude -- order of PN expansion of waveform
         min_waveform_length -- min number of points in frequency-domain waveforms
         max_waveform_length -- max number of points in frequency-domain waveforms
         max_time_spacing -- smallest time spacing of waveform evaluation in [s]. 
                             for the default nan, the code will set this automatically
                             to 1e-2 of the detector_orbit_period
         max_measurement_time -- if the lifetime of the signal in the band
                                 is larger than max_measurement_time, only the
                                 part of the signal within max_measurement_time 
                                 at frequencies below detector_fGWmax will be 
                                 considered. 
                                 Units of max_measurement_time are [seconds]
         Delta_Mc_rel -- relative change of chirp mass parameter for finite
                           difference computation of the Fisher matrix
         Delta_q_rel -- relative change of q parameter for the finite
                        difference computation of the Fisher Matrix
         Delta_iota_abs -- absolute change of iota parameter for finite
                           difference computation of the Fisher matrix [rad]
         Delta_phi0_abs -- ansolute change of phi0 parameter for finite
                           difference computation of the Fisher matrix [rad]
         Delta_tc_abs -- absolute change of t_c parameter [seconds] for finite
                         difference computation of the Fisher matrix
         Delta_dL_rel -- relative change of dL parameter for finite
                         difference computation of the Fisher matrix
         Delta_angleSL_abs -- absolute change in sky localization angle 
                              parameters for finite difference computation 
                              of Fisher matrix [rad]
         Delta_psi_abs -- absolute change of psi parameter for finite
                          difference computation of the Fisher matrix [rad]
         Delta_autoadjust -- boolean to turn on/off the automatic adjustment
                             of Delta_Mc_rel with the chirp mass
         fast_tau_evaluation -- boolean Flag. Set to True to use leading-order
                                expression for calculation of tau(f_GW), \
                                if False the PN expressions will be used
         verbose -- boolean Flag. Set to True to get the code to write the 
                    progress of the computation to std_out
      '''
      self.source_Mc = source_Mc
      self.source_q = source_q
      self.source_iota = source_iota
      self.source_phi0 = source_phi0
      self.source_tc = source_tc
      self.source_dL = source_dL
      self.source_RA = source_RA
      self.source_DEC = source_DEC
      self.source_psi = source_psi
      self.detector_fGWmin = detector_fGWmin
      self.detector_fGWmax = detector_fGWmax
      self.detector_orbit_R = detector_orbit_R
      self.detector_orbit_period = detector_orbit_period
      self.detector_orbit_t0 = detector_orbit_t0
      self.detector_orbit_RA0 = detector_orbit_RA0
      self.detector_orbit_DEC0 = detector_orbit_DEC0
      self.PN_order_phase = PN_order_phase
      self.PN_order_amplitude = PN_order_amplitude
      self.min_waveform_length = min_waveform_length
      self.max_waveform_length = max_waveform_length
      self.max_time_spacing = max_time_spacing
      self.max_measurement_time = max_measurement_time
      self.Delta_Mc_rel = Delta_Mc_rel
      self.Delta_q_rel = Delta_q_rel
      self.Delta_iota_abs = Delta_iota_abs
      self.Delta_phi0_abs = Delta_phi0_abs
      self.Delta_tc_abs = Delta_tc_abs
      self.Delta_dL_rel = Delta_dL_rel
      self.Delta_angleSL_abs = Delta_angleSL_abs
      self.Delta_psi_abs = Delta_psi_abs
      self.Delta_autoadjust = Delta_autoadjust
      self.fast_tau_evaluation = fast_tau_evaluation
      self.verbose = verbose
      if np.isnan(self.max_time_spacing):
         self.max_time_spacing = 1e-2*self.detector_orbit_period
      if self.Delta_autoadjust:
         self.Delta_Mc_rel *= self.source_Mc**(5/3)

   ###################################################
   # Class Functions
   ###################################################

   def get_SNR2(self):
      '''
      Appends the class with SNR^2 of the GW signal
      '''
      # compute the frequency array on which to compute the GW signal
      # This computation is based on the LO time-frequency relation
      #
      # start by calculating the time-before-merger when the signal
      # is at detector_fGWmax
      eval_taumin = waveform_LO.tau_fGW(
         self.source_Mc, 
         self.detector_fGWmax
         )
      # calculate time-before-merger when the signal is at 
      # detector_fGWmin or one year before eval_taumin
      eval_taumax = np.min([
         waveform_LO.tau_fGW(self.source_Mc, self.detector_fGWmin),
         eval_taumin + self.max_measurement_time
         ])
      # check that fvec would not be longer that max_waveform_length
      if (eval_taumax-eval_taumin)/self.max_time_spacing > self.max_waveform_length:
         print('''
            your choices for 'max_time_spaceing'
            and 'max_measurement_time' lead to a 
            signal longer than 'max_waveform_length'
            ''')
         print("stopping calculation")
         sys.exit()
      # make linearly spaced time array
      eval_tauvec = np.linspace(
         eval_taumax, 
         eval_taumin, 
         int(np.ceil(
            (eval_taumax-eval_taumin)/self.max_time_spacing
            ))
         )
      # check if the time_array is longer than min_waveform_length
      if len(eval_tauvec) < self.min_waveform_length:
         eval_tauvec = np.linspace(
         eval_taumax, 
         eval_taumin, 
         int(self.min_waveform_length)
         )
      # convert to frequency array
      eval_fvec = waveform_LO.fGW(self.source_Mc, eval_tauvec)

      # get the time and location of the clock
      self.clock_tD = np.double(
         self.source_tc 
         - 0.5*(eval_tauvec[0]+eval_tauvec[-1])
         )
      self.clock_r = antennaFuns_satellites.orbit_earth_sun(np.array(self.clock_tD))
      # make tuple of run_params
      run_params = (
         self.source_Mc,
         self.source_q,
         self.source_iota,
         self.source_phi0,
         self.source_tc,
         self.source_dL,
         self.source_RA,
         self.source_DEC,
         self.source_psi,
         self.detector_orbit_R,
         self.detector_orbit_period,
         self.detector_orbit_t0,
         self.detector_orbit_RA0,
         self.detector_orbit_DEC0,
         eval_fvec,
         self.clock_tD,
         self.PN_order_phase,
         self.PN_order_amplitude,
         self.fast_tau_evaluation
         )
      self.SNR2 = SNR2(*run_params)
      return

#-----------------------------------------------------------------------------------------------------------------------------#   
 
   def get_SNR_FisherMatrix(self):
      '''
      Appends the class with SNR^2 and FisherMatrix of the GW signal
      '''
      # compute the frequency array on which to compute the GW signal
      # This computation is based on the LO time-frequency relation
      #
      # start by calculating the time-before-merger when the signal
      # is at detector_fGWmax
      eval_taumin = waveform_LO.tau_fGW(
         self.source_Mc, 
         self.detector_fGWmax
         )
      # calculate time-before-merger when the signal is at 
      # detector_fGWmin or one year before eval_taumin
      eval_taumax = np.min([
         waveform_LO.tau_fGW(self.source_Mc, self.detector_fGWmin),
         eval_taumin + self.max_measurement_time
         ])
      # check that fvec would not be longer that max_waveform_length
      if (eval_taumax-eval_taumin)/self.max_time_spacing > self.max_waveform_length:
         print('''
            your choices for 'max_time_spaceing'
            and 'max_measurement_time' lead to a 
            signal longer than 'max_waveform_length'
            ''')
         print("stopping calculation")
         sys.exit()
      # make linearly spaced time array
      eval_tauvec = np.linspace(
         eval_taumax, 
         eval_taumin, 
         int(np.ceil(
            (eval_taumax-eval_taumin)/self.max_time_spacing
            ))
         )
      # check if the time_array is longer than min_waveform_length
      if len(eval_tauvec) < self.min_waveform_length:
         eval_tauvec = np.linspace(
         eval_taumax, 
         eval_taumin, 
         int(self.min_waveform_length)
         )
      # convert to frequency array
      eval_fvec = waveform_LO.fGW(self.source_Mc, eval_tauvec)

      # get the time and location of the clock
      self.clock_tD = np.double(
         self.source_tc 
         - 0.5*(eval_tauvec[0]+eval_tauvec[-1])
         )
      self.clock_r = antennaFuns_satellites.orbit_earth_sun(np.array(self.clock_tD))
      # make tuple of run_params
      run_params = (
         self.source_Mc,
         self.source_q,
         self.source_iota,
         self.source_phi0,
         self.source_tc,
         self.source_dL,
         self.source_RA,
         self.source_DEC,
         self.source_psi,
         self.detector_orbit_R,
         self.detector_orbit_period,
         self.detector_orbit_t0,
         self.detector_orbit_RA0,
         self.detector_orbit_DEC0,
         eval_fvec,
         self.clock_tD,
         self.PN_order_phase,
         self.PN_order_amplitude,
         self.Delta_Mc_rel,
         self.Delta_q_rel,
         self.Delta_iota_abs,
         self.Delta_phi0_abs,
         self.Delta_tc_abs,
         self.Delta_dL_rel,
         self.Delta_angleSL_abs,
         self.Delta_psi_abs,
         self.fast_tau_evaluation
         )
      # run the Fisher matrix calculation
      starttime = time.time() 
      if self.verbose:
         print("Starting calculation of FisherMat.")
      self.SNR2, self.FisherMat = SNR_FisherMatrix(*run_params)
      if self.verbose:
         print("Finished calculation of SNR and FisherMat. Time elapsed: {:.1f} sec".format(time.time()-starttime))
      # make Fisher matrix symmetric
      for i in range(self.FisherMat.shape[0]):
         for j in range(i+1, self.FisherMat.shape[1]):
            self.FisherMat[j,i] = self.FisherMat[i,j]
      # and add dimensionless version of Fisher matrix
      self.FisherMat_dimless = np.copy(self.FisherMat)
      source_params_dimful = np.ones(self.FisherMat.shape[0])
      source_params_dimful[0] = self.source_Mc
      source_params_dimful[5] = self.source_dL
      for i in range(self.FisherMat.shape[0]):
         self.FisherMat_dimless[:,i] *= source_params_dimful[i]
         self.FisherMat_dimless[i,:] *= source_params_dimful[i]
      # key for Fisher/Covariance matrix entries
      self.FisherMat_key = {
         0 : 'source_Mc [Msol]',
         1 : 'source_q',
         2 : 'source_iota [rad]',
         3 : 'source_phi0 [rad]',
         4 : 'source_tc [s]',
         5 : 'source_dL [Mpc]',
         6 : 'source_RA [rad]',
         7 : 'source_DEC [rad]',
         8 : 'source_psi [rad]',
         }
      # matrix to be added to Fisher matrix to include
      # Gaussian priors with std (2)pi for angles.
      self.PriorMat = np.diag([
         0., 
         0., 
         1./np.pi**2, 
         1./(4.*np.pi**2),
         0.,
         0.,
         1./(4.*np.pi**2),
         1./np.pi**2, 
         1./(4.*np.pi**2)
         ])
      return

   #------------------------------------------------#   

   def get_CoVaMat(self):
      '''
      Appends class by covariance matrix.
      If Fisher Matrix not computed, executing this
      function will run the FisherMatrix computation
      '''
      try:
         self.FisherMat
      except:
         self.get_SNR_FisherMatrix()
      self.CoVaMat = np.linalg.inv(np.double(self.FisherMat))
      return
      
   #------------------------------------------------#   

   def get_CoVaMat_dimless(self):
      '''
      Appends class by dimensionless covariance matrix.
      If Fisher Matrix not computed, executing this
      function will run the FisherMatrix computation
      '''
      try:
         self.FisherMat_dimless
      except:
         self.get_SNR_FisherMatrix()
      self.CoVaMat_dimless = np.linalg.inv(
         np.double(
            self.FisherMat_dimless
            )
         )
      return   
      
   #------------------------------------------------#   

   def get_CoVaMat_priors(self):
      '''
      Appends class by covariance matrix
      including Gaussian priors with std (2)pi for angles.
      If Fisher Matrix not computed, executing this
      function will run the FisherMatrix computation
      '''
      try:
         self.FisherMat
      except:
         self.get_SNR_FisherMatrix()
      self.CoVaMat_priors = np.linalg.inv(
         np.double(
            self.FisherMat + self.PriorMat
            )
         )
      return
         
   #------------------------------------------------#   

   def get_CoVaMat_priors_dimless(self):
      '''
      Appends class by dimensionless covariance matrix
      including Gaussian priors with std (2)pi for angles.
      If Fisher Matrix not computed, executing this
      function will run the FisherMatrix computation
      '''
      try:
         self.FisherMat_dimless
      except:
         self.get_SNR_FisherMatrix()
      self.CoVaMat_priors_dimless = np.linalg.inv(
         np.double(
            self.FisherMat_dimless + self.PriorMat
            )
         )
      return
         
   #------------------------------------------------#   

   def get_angular_resolution(self):
      '''
      Appends class by angular resolution in [sr]
      If Fisher CoVaMat (and FisherMatrix) is not computed,
      executing this function will run the CoVaMat 
      (and FisherMatrix) computation
      '''
      try:
         self.CoVaMat
      except:
         self.get_CoVaMat()
      self.angular_resolution = (
         2.*np.pi
         * np.sin(np.pi/2. - self.source_DEC)
         * np.sqrt(
            self.CoVaMat[6,6]*self.CoVaMat[7,7]
            - self.CoVaMat[6,7]**2
            )
         )
      return
         
   #------------------------------------------------#   

   def get_angular_resolution_priors(self):
      '''
      Appends class by angular resolution in [sr]
      including Gaussian priors with std (2)pi for angles.
      If Fisher CoVaMat (and FisherMatrix) is not computed,
      executing this function will run the CoVaMat 
      (and FisherMatrix) computation
      '''
      try:
         self.CoVaMat_priors
      except:
         self.get_CoVaMat_priors()
      self.angular_resolution_priors = (
         2.*np.pi
         * np.sin(np.pi/2. - self.source_DEC)
         * np.sqrt(
            self.CoVaMat_priors[6,6]*self.CoVaMat_priors[7,7]
            - self.CoVaMat_priors[6,7]**2
            )
         )
      return
