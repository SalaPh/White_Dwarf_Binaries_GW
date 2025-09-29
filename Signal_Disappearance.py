############################################################################################################
# Code to compute the time after which a certain SNR is reached after RLOF frequency (signaling the SNR expected if a binary has not merged)
# The computation considers no mass transfer
# The time to reach a certain SNR is computed with a multiprocessor code and stored in the folder Precision_Output/
############################################################################################################

import pickle
import numpy as np
import time
import os
import sys

import multiprocessing as mp

from scipy.optimize import root_scalar

from Analysis_Scripts import constants as const
from Analysis_Scripts import waveform_LO
from Analysis_Scripts import WD

######################################################
# load parameter estimation code
import ParamEstimator_SpaceAI as AI_PE

###################################################
# Multiprocessing settings
###################################################

Ncores = 4

###################################################
# CLASS to compute the disappearance time of the signal
###################################################

class Disappearance:
    # Class to compute a grid of times required to reach a certain SNR after RLOF frequency with the code ParamEstimator_SpaceAI.py
    # It stores the results in the folder Precision_Output/ 
    '''
    This class can be called to estimate the time required to recognise the disappearance of a WD binary signal
    by computing the time after which a certain SNR (signal expected) is reached after RLOF frequency.
    The computation considers no mass transfer
    This computation is performed for a for a grid of values of either the two WD masses (M1,M2) 
    or the chirp mass and luminosity distance (Mc,dL). The other parameters are kept fixed.
    '''
    def __init__(self,
                 run_type,
                 selected_SNR, 
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
            selected_SNR: minimum SNR to reach to consider the signal as disappeared
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
        Range setting inputs:
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
            PN_order_phase: post-Newtonian order for the phase of the waveform
            PN_order_amplitude: post-Newtonian order for the amplitude of the waveform
            fast_tau_evaluation: boolean Flag. Set to True to use leading-order
                        expression for calculation of tau(f_GW), \
                        if False the PN expressions will be used
        '''
        self.run_type=run_type 
        self.selected_SNR=selected_SNR
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
        if self.run_type=="M1M2":
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
                        self.detector_orbit_R_sat,
                        self.detector_orbit_period,
                        self.detector_orbit_t0,
                        self.detector_orbit_RA0,
                        self.detector_orbit_DEC0,
                        self.PN_order_phase,
                        self.PN_order_amplitude,
                        self.detector_max_measurement_time,
                        self.fast_tau_evaluation
                        ))
                    j+=1

        # In the case of a grid in (Mc,dL)
        elif self.run_type=="McdL":
            for i in range(self.rangelength):
                source_McMcdL = self.Mc[i]
                for j in range(self.rangelength):
                    source_dLMcdL = self.dL[j]

                    self.run_params_list.append((
                        source_McMcdL,
                        self.source_q,
                        self.source_iota,
                        self.source_phi0,
                        self.source_tc,
                        source_dLMcdL,
                        self.source_RA,
                        self.source_DEC,
                        self.source_psi,
                        self.detector_fGWmin,
                        self.detector_fGWmax,
                        self.detector_orbit_R_cb,
                        self.detector_orbit_R_sat,
                        self.detector_orbit_period,
                        self.detector_orbit_t0,
                        self.detector_orbit_RA0,
                        self.detector_orbit_DEC0,
                        self.PN_order_phase,
                        self.PN_order_amplitude,
                        self.detector_max_measurement_time,
                        self.fast_tau_evaluation
                        ))

        # If the run_type is not recognised, exit
        else:
            print("run_type not recognised")
            exit()
        
    ###################################################
    # Functions
    ###################################################

    def SNR_after_RLOF(self,
        source_Mc,
        source_q,
        source_iota,
        source_phi0,
        source_tc,
        source_dL,
        source_RA,
        source_DEC,
        source_psi,
        starting_frequency,
        duration_noSignal,
        detector_orbit_R_cb,
        detector_orbit_period,
        detector_orbit_t0,
        detector_orbit_RA0,
        detector_orbit_DEC0,
        PN_order_phase,
        PN_order_amplitude,
        detector_max_measurement_time,
        fast_tau_evaluation
        ):
        '''
        Redefinition of the ending frequencies of the computation of the SNR, and computation of it with
        ParamEstimator_SpaceAI.SNR2
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
            starting_frequency -- frequency at which the SNR is computed [Hz]
            duration_noSignal -- time after which the SNR is computed (and the signal could have disappeard) [seconds]
            detector_orbit_R -- radius of the orbit of the satellites around 
                                Earth [m]
            detector_orbit_period -- period of the detector around Earth [seconds]
            detector_orbit_t0 -- reference time for fixing orbit of satellite 
                                around Earth [seconds]
            detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
            detector_orbit_DEC0 -- declination of the satellite at t0 [rad]
        Parameters for computation:
            PN_order_phase -- order of PN expansion of phase
            PN_order_amplitude -- order of PN expansion of amplitude
            fast_tau_evaluation -- boolean Flag. Set to True to use leading-order
                                    expression for calculation of tau(f_GW), 
                                    if False the PN expressions will be used
        output:
            sol_SNR -- SNR computed between starting_frequency and the frequency
                        after a time duration_noSignal
        '''
        # Given a starting frequency and a duration, compute ending frequency at leading order
        fGW_start=starting_frequency
        Tau_starting=np.double(waveform_LO.tau_fGW(source_Mc,starting_frequency))
        fGW_after_noSignal=np.double(waveform_LO.fGW(source_Mc,Tau_starting-duration_noSignal))
            
        # Compute the SNR between starting_frequency and fGW_after_noSignal
        event = AI_PE.GW_event(
            source_Mc,
            source_q,
            source_iota,
            source_phi0,
            source_tc,
            source_dL,
            source_RA,
            source_DEC,
            source_psi,
            fGW_start,
            fGW_after_noSignal,
            detector_orbit_R_cb,
            detector_orbit_period,
            detector_orbit_t0,
            detector_orbit_RA0,
            detector_orbit_DEC0,
            PN_order_phase,
            PN_order_amplitude,
            max_measurement_time = detector_max_measurement_time,
            fast_tau_evaluation = fast_tau_evaluation,
            verbose = False
            )
        event.get_SNR2()
        sol_SNR=np.sqrt(event.SNR2)
        del event
        return sol_SNR

    #------------------------------------------------#

    def Time_until_SNR(self,
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
        detector_orbit_R_cb,
        detector_orbit_R_sat,
        detector_orbit_period,
        detector_orbit_t0,
        detector_orbit_RA0,
        detector_orbit_DEC0,
        PN_order_phase,
        PN_order_amplitude,
        detector_max_measurement_time,
        fast_tau_evaluation,
        SNR_choice = 2,
        starting_frequency = "RLOF"
        ):
        '''
        Function to compute the time after which a certain SNR is reached after RLOF frequency
        This is achieved by a root finding algorithm
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
            detector_fGWmin -- smallest frequency considered [Hz]
            detector_fGWmax -- largest frequency considered [Hz]
            detector_orbit_R -- radius of the orbit of the satellites around 
                                Earth [m]
            detector_orbit_period -- period of the detector around Earth [seconds]
            detector_orbit_t0 -- reference time for fixing orbit of satellite 
                                around Earth [seconds]
            detector_orbit_RA0 -- right ascension of the satellite at t0 [rad]
        Parameters for computation:
            PN_order_phase -- order of PN expansion of phase
            PN_order_amplitude -- order of PN expansion of amplitude
            detector_max_measurement_time -- max length of signal considered [seconds]
            fast_tau_evaluation -- boolean Flag. Set to True to use leading-order
                                    expression for calculation of tau(f_GW), 
                                    if False the PN expressions will be used
        Disappearance criteria input:                            
            SNR_choice -- minimum SNR to reach to consider the signal as disappeared
            starting_frequency -- frequency at which the SNR is computed, 
                                either "RLOF" or a number in Hz
        output:
            sol_time -- time after RLOF after which the SNR_choice is reached [seconds]
        '''

        # initialize the WD class and check if it is a measurable WDB
        WDprop=WD.WDclass(source_Mc,
                        source_q,
                        detector_fGWmin,
                        detector_fGWmax,
                        detector_max_measurement_time,
                        detector_orbit_R_cb=detector_orbit_R_cb,
                        detector_orbit_R_sat=detector_orbit_R_sat,
                        detector_orbit_period=detector_orbit_period
                        )
        WDprop.CutoffFrequency_func()

        # If the binary is not measurable, return NaN (we can compute this quantity even if the binary is not a WDB)
        if not WDprop.check:
            return np.nan

        # Pick the initial frequency
        if starting_frequency == "RLOF":
            fGW_start=WDprop.WD_detector_fGWmax
        else:
            fGW_start=detector_fGWmin
            print('No RLOF frequency was used')

        # Initialize the function for root finding
        fun = lambda x: SNR_choice-self.SNR_after_RLOF(
                                        source_Mc,
                                        source_q,
                                        source_iota,
                                        source_phi0,
                                        source_tc,
                                        source_dL,
                                        source_RA,
                                        source_DEC,
                                        source_psi,
                                        fGW_start,
                                        x,
                                        detector_orbit_R_cb,
                                        detector_orbit_period,
                                        detector_orbit_t0,
                                        detector_orbit_RA0,
                                        detector_orbit_DEC0,
                                        PN_order_phase,
                                        PN_order_amplitude,
                                        detector_max_measurement_time,
                                        fast_tau_evaluation
                                        )
        
        #Choose a possible range of solution time, minimum is cbar_min, maximum is either cbar_max or the coalescence time
        coalescence_time=waveform_LO.tau_fGW(source_Mc,fGW_start) 
        possible_max=min(self.cbar_max,coalescence_time/1.01)

        # Check if the solution time is in the range, otherwise return cbar_min or cbar_max
        if fun(possible_max) > 0 :
            sol_time=self.cbar_max*1.01
            if coalescence_time < self.cbar_max :
                sol_time=np.nan
            elif fun(coalescence_time/1.01) > 0 :
                sol_time=np.nan
        elif fun(self.cbar_min) < 0 :
            sol_time=self.cbar_min/1.01
        # If the solution time is in the range, find it with root_scalar
        else:
            sol_time = root_scalar(fun, bracket=[self.cbar_min, possible_max], method="ridder").root
        del WDprop
        return sol_time

    #------------------------------------------------#

    def run_event(self,run_params,Solution):
        '''
        Function to run a single event of the grid and store the result in a multiprocessor list
        input:
            run_params -- tuple of parameters to run the event
            Solution -- multiprocessor list to store the result
        output:
            the result is stored in the multiprocessor list Solution
        '''
        # Depending on the run_type, store in the multiprocessor list the result of the run_event

        # In the case of M1M2, store (m1,m2,time)
        if self.run_type=="M1M2":
            sol_time = self.Time_until_SNR(*run_params,SNR_choice = self.selected_SNR)
            Solution.append([waveform_LO.m1(run_params[0],run_params[1]),
                            waveform_LO.m2(run_params[0],run_params[1]),
                            sol_time])
            return 

        # In the case of McdL, store (Mc,dL,time)
        elif self.run_type=="McdL":
            sol_time = self.Time_until_SNR(*run_params,SNR_choice = self.selected_SNR)
            Solution.append([run_params[0],
                            run_params[5],
                            sol_time])
            return 
        
        # If the run_type is not recognised, exit
        else:
            print("run_type not recognised")
            exit()
            return

    #------------------------------------------------#

    def track_mp_progress(self,job, mp_inputs, starttime, update_interval=1*60):
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
        if not os.path.exists(self.fpath_out + 'Signal_to_SNR_after_RLOF.pkl'):
            
            print("###########################################")
            print("starting to run")
            print("###########################################")
            
            # run the computation in parallel using multiprocessing (mp) with Ncores cores
            with mp.Manager() as manager:
                Solution = manager.list()  # This is the managed list
                pool = mp.Pool(Ncores)
                p = pool.starmap_async(self.run_event, [(run_params, Solution) for run_params in self.run_params_list])
                self.track_mp_progress(p, self.run_params_list, time_start)
                p.get()
                pool.close()
                pool.join()
                # Convert the managed list to a regular numpy array
                solution=np.array(Solution)

            #######################
            # write info file for future reference
            fo = open(self.fpath_out+'/info_Signal_Disappearance.txt', 'w')
            fo.write(f"#WDB precision on variables {self.filename} \n")

            fo.write('# inputs for the disappearance computation\n')
            fo.write(self.run_type+'# Measurement Type\n')
            fo.write(str(self.selected_SNR)+' # selected_SNR, minimum SNR to consider the signal as disappeared\n')
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

            print("###########################################")
            print("Finished!")
            print("###########################################")

            # Convert lists to numpy arrays for storage

            # Create the grid for the parameters
            if self.run_type=="M1M2":
                X,Y = np.meshgrid(self.M,self.M)
            elif self.run_type=="McdL":
                X,Y = np.meshgrid(self.Mc,self.dL)

            # Filter the parameters to avoid numerical issues in finding the right place in the grid
            parameter_precision =1e-5

            # Create the Z array with the same shape as X and Y and fill it with filtered results
            Z = np.full(np.shape(X),np.nan)
            for i in range(self.rangelength):
                for j in range(self.rangelength):
                    #Place the value for a certain combination of m1 and m2 (or Mc and dL) in the right position
                    parameters_check=((abs(solution[:,0]-X[i,j])<=parameter_precision)*(abs(solution[:,1]-Y[i,j])<=parameter_precision))
                    if sum(parameters_check)>=1:
                        Z[i,j] = solution[parameters_check,2][0]

            # Save the results as a pickle file
            with open(os.path.join(self.fpath_out, 'Signal_to_SNR_after_RLOF.pkl'), 'wb') as f:
                            pickle.dump((X, Y, Z), f, pickle.HIGHEST_PROTOCOL)        
        else:
            print('Computation was previously already done.')
        return