############################################################################
# custom format functions for matplotlib plots
############################################################################

import numpy as np

#------------------------------------------------#
# normal format with 1, 2 or 3 decimal places
fmt1 = lambda x,pos: '{:.1f}'.format(x)
fmt2 = lambda x,pos: '{:.2f}'.format(x)
fmt3 = lambda x,pos: '{:.3f}'.format(x)


#------------------------------------------------#
# scientific format with 1, 2 or 3 decimal places
emt1 = lambda x,pos: '{:.1e}'.format(x)
emt2 = lambda x,pos: '{:.2e}'.format(x)
emt3 = lambda x,pos: '{:.3e}'.format(x)

###################################################
# FUNCTIONS
###################################################
 
def scientific_formatter(x, pos):
    '''
    Format that shows only 10^exponent or 1 if x=1, does not show e.g. 2*10^3 or any other number
    '''
    absx=abs(x)
    if x >=0:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            if coeff==1:
                string = '${:.0f}$'.format(coeff) 
            else:
                string = ''
        else:
            if coeff==1:
                string = r'$10^{{{}}}$'.format(exponent)
            else:
                string = '' 
    else:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '$-{:.0f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$-10^{{{}}}$'.format(exponent)
            else:
                string = '' 
    return string

#------------------------------------------------#

def scientific_formatter0(x, pos):
    '''
    Format that shows only 10^exponent or 1 if x=1, does show any coefficient with 0 decimal places
    '''
    absx=abs(x)
    if x >=0:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '${:.0f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$10^{{{}}}$'.format(exponent)
            else:
                string = '{:.0f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    else:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '$-{:.0f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$-10^{{{}}}$'.format(exponent)
            else:
                string = '-{:.0f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    return string

#------------------------------------------------#

def scientific_formatter1(x, pos):
    '''
    Format that shows only 10^exponent or 1.0 if x=1, does show any coefficient with 1 decimal places
    '''
    absx=abs(x)
    if x >=0:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '${:.1f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$10^{{{}}}$'.format(exponent)
            else:
                string = '{:.1f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    else:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '$-{:.1f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$-10^{{{}}}$'.format(exponent)
            else:
                string = '-{:.1f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    return string

#------------------------------------------------#

def scientific_formatter2(x, pos):
    '''
    Format that shows only 10^exponent or 1.00 if x=1, does show any coefficient with 2 decimal places
    '''
    absx=abs(x)
    if x >=0:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '${:.2f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$10^{{{}}}$'.format(exponent)
            else:
                string = '{:.2f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    else:
        exponent = int(np.floor(np.log10(absx)))
        coeff = absx / 10**exponent
        if exponent == 0:
            string = '$-{:.2f}$'.format(coeff) 
        else:
            if coeff==1:
                string = r'$-10^{{{}}}$'.format(exponent)
            else:
                string = '-{:.2f}'.format(coeff) + r'$\cdot10^{{{}}}$'.format(exponent)
    return string

#------------------------------------------------#
