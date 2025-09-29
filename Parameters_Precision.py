############################################################################
# Code to compute a grid of parameter estimation precisions with code ParamEstimator_SpaceAI.py
# The SNR and fisher matrix results (which can be used to compute the parameter precision) are stored the results in the folder Precision_Output/
############################################################################

import pickle
import numpy as np
import time
import os
import sys

import multiprocessing as mp

from Analysis_Scripts import constants as const
from Analysis_Scripts import waveform_LO
from Analysis_Scripts import WD

######################################################
# load parameter estimation code
import ParamEstimator_SpaceAI as AI_PE

###################################################
# Multiprocessing settings
###################################################

Ncores = 8

###################################################
# CLASS to compute parameter reconstruction precision
###################################################

class Precisions:
   # Class to compute a grid of parameter estimation precisions with code ParamEstimator_SpaceAI.py It stores the results in the folder Precision_Output/
   ''' 
   This class can be called for computing the precision of several GW parameters for a grid of values 
   of either the two WD masses (M1,M2) or the chirp mass and luminosity distance (Mc,dL).
   The other parameters of the source are kept fixed.
   '''
   def __init__(self,
                 run_type,
                 rangelength, 
                 source_Mc, 
                 source_q,
                 source_dL,
                 source_tau0,
                 source_tc,
                 source_phi0, 
                 source_iota,
                 source_psi,
                 source_RA,
                 source_DEC_title,
                 Mcmin=0.35,
                 Mcmax=1.24,
                 minMass=0.4,
                 maxMass=WD.maxWDMass*0.95,
                 dLmin=0.1,
                 dLmax=500,
                 minint1=0.001,
                 maxint1=10.,
                 duration=1.,
                 cbar_min=1000,
                 cbar_max=1e9, 
                 detector_orbit_R_cb = 8.44e6,
                 detector_orbit_R_sat = 2e7, 
                 detector_orbit_t0 = 0.,
                 detector_orbit_RA0 = np.pi/2, 
                 detector_orbit_DEC0 = -23.4*np.pi/180, 
                 PN_order_phase = 3.5,
                 PN_order_amplitude = 3.0,
                 fast_tau_evaluation = True 
                ):
      '''
      Grid setting inputs:
         run_type: "M1M2" or "McdL"
            "M1M2": grid in the two WD masses (M1,M2)
            "McdL": grid in the chirp mass and luminosity distance (Mc,dL)
         rangelength: number of bins in the chosen grid axes
      GW source inputs:
         source_Mc: detector-frame chirp mass [Msol]
         source_q: mass ratio (m_1/m_2)
         source_dL: luminosity distance [Mpc]
         source_tau0: decided to be the same as coalescence time
         source_tc: time of merger, measured from the solar equinox [seconds]
         source_phi0: [0,2*np.pi] phase of the GW signal at fGWmin
         source_iota: [0,np.pi] angle between sky localization and orbital momentum of the binary [rad]
         source_psi: [0,2*np.pi] polarization angle [rad]
         source_RA: [0,2*np.pi] source right ascesnsion [rad], max at np.pi/2 for other values maximised
         source_DEC_title: [-np.pi/2,np.pi/2] source declination [rad], max at np.pi/2
      Setting regarding the ranges to explore:
         Mcmin: minimum chirp mass considered [Msol]
         Mcmax: maximum chirp mass considered [Msol]
         minMass: minimum mass considered [Msol]
         maxMass: maximum mass considered [Msol]
         dLmin: minimum luminosity distance considered [Mpc]
         dLmax: maximum luminosity distance considered [Mpc]
         minint1: detector_fGWmin, smallest frequency considered [Hz]
         maxint1: detector_fGWmax, largest frequency considered [Hz]
         duration: detector_max_measurement_time, max length of signal considered [years]
         cbar_min: minimum time considered for detecting the chosen SNR [s]
         cbar_max: maximum time considered for detecting the chosen SNR [s]
      Fixed GW parameters:
         detector_orbit_R_cb: radius of the orbit of the center of the baseline around Earth [m],
                                       over which the DNR and Fisher matrix is computed
         detector_orbit_R_sat: radius of the orbit of the satellite around Earth [m], for period
         detector_orbit_t0: reference time at which to define the location of the detector [seconds after vernal equinox]
         detector_orbit_RA0: right ascension of the satellite at t0 [rad]
         detector_orbit_DEC0: declination of the satellite at t0 [rad]
      Parameters for computation:
         PN_order_phase: post-newtonian order for the phase of the waveform
         PN_order_amplitude: post-newtonian order for the amplitude of the waveform
         fast_tau_evaluation: boolean Flag. Set to True to use leading-order
                           expression for calculation of tau(f_GW), \
                           if False the PN expressions will be used
      '''    
      self.run_type=run_type 
      self.rangelength=rangelength
      self.minMass=minMass
      self.maxMass=maxMass
      self.Mcmin=Mcmin
      self.Mcmax=Mcmax
      # Initialize the grid ranges   
      self.M=np.linspace(minMass,maxMass,rangelength)
      self.Mc = np.linspace(max(Mcmin,waveform_LO.Mc(source_q*WD.minWDMass,WD.minWDMass)*1.01),
                              min(Mcmax,waveform_LO.Mc(WD.maxWDMass,WD.maxWDMass/source_q)/1.01),rangelength)
      self.dL = np.geomspace(dLmin,dLmax,rangelength)

      self.cbar_min=cbar_min
      self.cbar_max=cbar_max
      self.source_Mc=source_Mc 
      self.source_q=source_q 
      self.source_dL=source_dL 
      self.source_tau0=source_tau0 
      self.source_tc=source_tc 
      self.source_phi0=source_phi0
      self.source_iota=source_iota 
      self.source_psi=source_psi 
      self.source_RA=source_RA 
      self.source_DEC_title=source_DEC_title 
      self.source_DEC=source_DEC_title + detector_orbit_DEC0
      self.detector_fGWmin = minint1 
      self.detector_fGWmax = maxint1 
      self.detector_max_measurement_time = duration*const.year 
      self.detector_orbit_R_cb = detector_orbit_R_cb 
      self.detector_orbit_R_sat = detector_orbit_R_sat 
      self.detector_orbit_t0 = detector_orbit_t0 
      self.detector_orbit_RA0 = detector_orbit_RA0
      self.detector_orbit_DEC0 = detector_orbit_DEC0
      self.PN_order_phase = PN_order_phase
      self.PN_order_amplitude = PN_order_amplitude
      self.fast_tau_evaluation = fast_tau_evaluation
      ######################################################
      # get some derived parameters
      self.detector_orbit_period = AI_PE.antennaFuns_satellites.period_satellite_earth(detector_orbit_R_sat)
      self.max_time_spacing=self.detector_orbit_period

      ######################################################
      # set up path to store results
      self.filename=f"a{source_RA/np.pi:.2f}d{source_DEC_title/np.pi:.2f}i{source_iota/np.pi:.2f}p{source_psi/np.pi:.2f}-"+run_type+f'-{rangelength}/'
      self.fpath_out = "Precision_Output/"+self.filename
      if not os.path.exists(self.fpath_out):
         os.makedirs(self.fpath_out)

      ######################################################
      # Initialize parameters to run
      self.run_params_list = []
      
      # In the case of a grid in (M1,M2)
      if run_type=="M1M2":
         for i in range(self.rangelength):
            j=0
            while j<=i:
               source_McM1M2 = waveform_LO.Mc(self.M[i],self.M[j])
               source_qM1M2 = self.M[i]/self.M[j]

               self.run_params_list.append((
                  source_McM1M2,
                  source_qM1M2,
                  self.source_iota,
                  self.source_phi0,
                  self.source_tc,
                  self.source_dL,
                  self.source_RA,
                  self.source_DEC,
                  self.source_psi,
                  self.detector_fGWmin,
                  self.detector_fGWmax,
                  self.detector_orbit_R_cb,
                  self.detector_orbit_period,
                  self.detector_orbit_t0,
                  self.detector_orbit_RA0,
                  self.detector_orbit_DEC0,
                  self.PN_order_phase,
                  self.PN_order_amplitude,
                  f"M1{self.M[i]:.2f}M2{self.M[j]:.2f}"
                  ))
               j+=1

      # In the case of a grid in (Mc,dL)
      elif self.run_type=="McdL":
         self.run_params_list = []
         for i in range(self.rangelength):
            source_McMcdL = self.Mc[i]

            self.run_params_list.append((
               source_McMcdL,
               self.source_q,
               self.source_iota,
               self.source_phi0,
               self.source_tc,
               self.source_dL,
               self.source_RA,
               self.source_DEC,
               self.source_psi,
               self.detector_fGWmin,
               self.detector_fGWmax,
               self.detector_orbit_R_cb,
               self.detector_orbit_period,
               self.detector_orbit_t0,
               self.detector_orbit_RA0,
               self.detector_orbit_DEC0,
               self.PN_order_phase,
               self.PN_order_amplitude,
               f"Mc{self.Mc[i]:.2f}"
               ))


      # If the run_type is not recognised, exit
      else:
         print("run_type not recognised")
         exit()
         return

   ###################################################
   # Functions
   ###################################################

   def run_event(self, i, run_params):
      '''
      Function to run a single event, given the parameters in run_params
      input:
         i: index of the run in the list of runs
         run_params: list of parameters to run the event
      output:
         saves the event in a pickle file in the folder Precision_Output/
      '''

      # Initialize the filename
      datafile = self.fpath_out+f'classdata_in_progress/data'+run_params[-1]+'.pkl'
      # check if the event was already computed, if not initialize the WD class and check if it is a measurable WDB
      if not os.path.exists(datafile):
         copy_run_params=list(run_params[:-1])
         WDparam = WD.WDclass(
                  copy_run_params[0], # source_Mc
                  copy_run_params[1], # source_q
                  copy_run_params[9], # detector_fGWmin
                  copy_run_params[10], # detector_fGWmax
                  self.detector_max_measurement_time, # max_measurement_time
                  detector_orbit_R_cb=self.detector_orbit_R_cb,
                  detector_orbit_R_sat=self.detector_orbit_R_sat,
                  detector_orbit_period=self.detector_orbit_period,
                  )
         WDparam.CutoffFrequency_func()

         # Compute the SNR and Fisher matrix of the GW signal of the initialized WDB 
         # only if the source is in the detector band and that is a WDB
         if WDparam.check and WDparam.WDBcheck: 
            copy_run_params[10]=WDparam.WD_detector_fGWmax

            event = AI_PE.GW_event(
               *copy_run_params, 
               **{
                  "max_measurement_time" : self.detector_max_measurement_time,
                  # "max_time_spacing" : max_time_spacing,
                  "fast_tau_evaluation" : self.fast_tau_evaluation,
                  "verbose" : True,
                  }
               )
            event.get_SNR_FisherMatrix()
            # save the event in a pickle file
            with open(self.fpath_out+f"classdata_in_progress/data"+run_params[-1]+".pkl", 'wb') as outp:
               pickle.dump(event, outp, pickle.HIGHEST_PROTOCOL)
            del event
         del WDparam
      return

   #------------------------------------------------#

   def track_mp_progress(self, job, mp_inputs, starttime, update_interval=1*60):
      '''
      Function to track the progress of a multiprocessing job
      input:
         job: multiprocessing job
         mp_inputs: list of inputs to the multiprocessing job
         starttime: time when the job started
         update_interval: time interval between updates in seconds 
      output:
         prints the progress of the job every update_interval seconds
      '''
      # function to track multiprocessing progress
      sleep_time = 10
      my_t0 = time.time()
      while job._number_left > 0:
         time.sleep(sleep_time)
         if time.time()-my_t0 > update_interval:
            print(
               str(int(100 * job._number_left * job._chunksize / len(mp_inputs))),
               "% of jobs remain to be computed; runtime:",
               str(int((time.time() - starttime)/60)),
               "min"
               )
            my_t0 = time.time()
            sys.stdout.flush()
      return

   def run(self):
      '''
      Function to run the grid of parameter estimation precisions.
      Runs the computation of function run_event in parallel using multiprocessing.
      output:
         saves the results in the folder Precision_Output/ through the function run_event      
      '''
      
      ######################################################
      #Print parameters
      print("\nComputation for ", self.run_type, " case\n")
      print(f"Computation at PN order: {self.PN_order_phase} and {self.PN_order_amplitude}\n")
      if self.run_type=="McdL":
         print("%.3g"% self.source_q, " source_q")        
      elif self.run_type=="M1M2":
         print("%.3g"% self.source_dL, " source_dL [Mpc]")
      print("%.3g"% self.source_iota, " source_iota [rad]")
      print("%.3g"% self.source_phi0, " source_phi0 [rad]")
      print("%.3g"% self.source_tc, " source_tc [s] ")
      print("%.3g"% self.source_RA, " source_RA [rad]")
      print("%.3g"% self.source_DEC_title, " source_DEC [rad] (compared to Earth equator)")
      print("%.3g"% self.source_psi, " source_psi [rad]\n")

      print("%.3g"% self.detector_orbit_RA0, " detector_RA [rad]")
      print("%.3g"% self.detector_orbit_DEC0, " detector_DEC [rad]\n")

      # Start the timer of the computation
      time_start = time.time()

      # Check if the computation was already done
      if not os.path.exists(self.fpath_out+'classdata'):
         # create a folder to store the results in progress
         if not os.path.exists(self.fpath_out+'classdata_in_progress'):
            os.mkdir(self.fpath_out+'classdata_in_progress')

         print("###########################################")
         print("starting to run")
         print("###########################################")
         
         # run the computation in parallel using multiprocessing (mp) with Ncores cores
         pool = mp.Pool(Ncores)
         p = pool.starmap_async(self.run_event, enumerate(self.run_params_list))
         self.track_mp_progress(p, self.run_params_list, time_start)
         p.get()
         pool.close()
         pool.join()

         #######################
         # write info file for future reference
         fo = open(self.fpath_out+'/info_Parameters_Precision.txt', 'w')
         fo.write(f"#WDB precision on variables {self.filename} \n")

         fo.write('# inputs for the disappearance computation\n')
         fo.write(self.run_type+'# Measurement Type\n')
         fo.write(str(self.rangelength)+' # number of scan points\n')

         fo.write('# inputs for GW source\n')
         if self.run_type=="M1M2":
               fo.write('Variable WD masses, hence chirp mass and source_q\n')
         elif self.run_type=="McdL":
               fo.write(f'Variable between {self.Mc[0]} and {self.Mc[-1]} # detector-frame chirp mass [Msol]\n')
               fo.write(str(self.source_q)+' # source_q, mass ratio (m_1/m_2)\n')
         fo.write(str(self.source_phi0)+' # source_phi0, phase of the GW signal at fGWmin\n')
         fo.write(str(self.source_tc)+' # source_tc, time of merger, measured from the solar equinox [seconds]\n')
         fo.write(str(self.source_dL)+' # source_dL, luminosity distance [Mpc]\n')
         fo.write(str(self.source_iota)+' # source_iota, angle between sky localization and orbital momentum of the binary [rad]\n')
         fo.write(str(self.source_psi)+' # source_psi, polarization angle [rad]\n')
         fo.write(str(self.source_RA)+' # source right ascension [rad]\n')
         fo.write(str(self.source_DEC)+' # source declination [rad]\n')

         fo.write('# inputs for detector\n')
         fo.write(str(self.detector_fGWmin)+' # detector_fGWmin, smallest frequency considered [Hz]\n')
         fo.write(str(self.detector_fGWmax)+' # detector_fGWmax, largest frequency considered [Hz]\n')
         fo.write(str(self.detector_max_measurement_time)+' # detector_max_measurement_time, max length of signal considered [seconds]\n')
         fo.write(str(self.detector_orbit_R_cb)+' # detector_orbit_R_cb, radius of the orbit of the center of the baseline around Earth [m]\n')
         fo.write(str(self.detector_orbit_R_sat)+' # detector_orbit_R_sat, radius of the orbit of the satellite around Earth [m], for period\n')
         fo.write(str(self.detector_orbit_t0)+' # detector_orbit_t0, reference time for fixing orbit of satellite around Earth [seconds]\n')
         fo.write(str(self.detector_orbit_RA0)+' # detector_orbit_RA0, right ascension of the satellite at t0 [rad]\n')
         fo.write(str(self.detector_orbit_DEC0)+' # detector_orbit_DEC0, declination of the satellite at t0 [rad]\n')

         fo.write('# inputs computation\n')
         fo.write(str(self.PN_order_phase)+' # PN_order_phase\n')
         fo.write(str(self.PN_order_amplitude)+' # PN_order_amplitude\n')
         fo.write(str(self.fast_tau_evaluation)+ ' # Fast tau(fGW) evaluation flag\n')

         fo.write('# inputs on the ranges to explore\n')
         if self.run_type=="M1M2":
               fo.write(str(self.minMass)+' # minMass, minimum mass considered [Msol]\n')
               fo.write(str(self.maxMass)+' # maxMass, maximum mass considered [Msol]\n')
         elif self.run_type=="McdL":
               fo.write(str(self.Mcmin)+' # Mcmin, minimum chirp mass considered [Msol]\n')
               fo.write(str(self.Mcmax)+' # Mcmax, maximum chirp mass considered [Msol]\n')
               fo.write(str(self.dLmin)+' # dLmin, minimum luminosity distance considered [Mpc]\n')
               fo.write(str(self.dLmax)+' # dLmax, maximum luminosity distance considered [Mpc]\n')
         fo.write(str(self.cbar_min)+' # cbar_min, minimum time considered for detecting the chosen SNR [s]\n')
         fo.write(str(self.cbar_max)+' # cbar_max, maximum time considered for detecting the chosen SNR [s]\n')

         fo.write(f"# total computation time: {time.time() - time_start} seconds\n")
         fo.close()

         # the computation is finished, rename the folder for storing the results in progress to the final name
         os.rename(self.fpath_out+'classdata_in_progress', self.fpath_out+'classdata')

         print("###########################################")
         print("Finished!")
         print("###########################################")
        
      else:
         print('Computation was previously already done.')
      return

   #------------------------------------------------#
