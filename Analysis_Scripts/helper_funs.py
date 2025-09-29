############################################################################
# Helper functions for Fisher matrix calculations
############################################################################

import numpy as np

###################################################
# FUNCTIONS
###################################################

def inner_product(
   htilde1, 
   htilde2, 
   noise, 
   freqs
   ):
   '''
   Inner product of two waveforms htilde1 and htilde2
   input:  
      htilde1, htilde2 : waveforms in frequency domain
      noise : noise power spectral density
      freqs : frequency array
   output:
      inner product <h1|h2>
   '''
   integrand = np.real(htilde1*np.conj(htilde2))/noise
   return 4.*np.trapz(integrand, x=freqs)

#------------------------------------------------#

def numDerivative_3s(fun, params, steps):
   '''
   Returns a 3-stencil central finite difference derivative
   of the function fun evaluated at the parameters given in 
   params shifted by the steps input.
   Params and steps must be of the same size, and organized 
   such that fun(params) returns the function at the parameter 
   values, and fun(params+steps) the function at the shifted
   parameter values
   input:
      fun : function to differentiate
      params : parameters at which to evaluate the derivative
      steps : step sizes for each parameter
   output:
      derivative of fun at params
   '''
   return (
      - fun(params-steps)
      + fun(params+steps)
      )/(2*np.linalg.norm(steps))

#------------------------------------------------#

def numDerivative_5s(fun, params, steps):
   '''
   Returns a 5-stencil central finite difference derivative
   of the function fun evaluated at the parameters given in 
   params shifted by the steps input.
   Params and steps must be of the same size, and organized 
   such that fun(params) returns the function at the parameter 
   values, and fun(params+steps) the function at the shifted
   parameter values
   input:
      fun : function to differentiate
      params : parameters at which to evaluate the derivative
      steps : step sizes for each parameter
   output:
      derivative of fun at params
   '''
   return (
      fun(params-2*steps)
      - 8*fun(params-steps)
      + 8*fun(params+steps)
      - fun(params+2*steps)
      )/(12*np.linalg.norm(steps))

#------------------------------------------------#
