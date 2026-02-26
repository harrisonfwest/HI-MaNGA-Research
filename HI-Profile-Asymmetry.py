# Python rewriting of profile_asymmetry.pro

import numpy as np

# IDL has a robust_sigma() function taking in an array-like, which we will recreate:
def robust_sigma(vals) -> float:
    x = np.array(x)
    median = np.median(x)
    mad = np.median(np.abs(x - median))
    return mad * 1.4826

def profile_asymmetry(velocity,flux,v0,vlow,vhigh,rms=None,plot=None,inspect=None):
    if rms == None:
        sel = np.sort(np.concatenate(np.argwhere(), np.argwhere()))