# Python rewriting of profile_asymmetry.pro
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from astropy.io import fits

# IDL has a robust_sigma() function taking in an array-like, which we will recreate:
def robust_sigma(vals) -> float:
    x = np.array(x)
    median = np.median(x)
    mad = np.median(np.abs(x - median))
    return mad * 1.4826

def profile_asymmetry(velocity: ArrayLike, flux: ArrayLike, v0: float, vlow: float, vhigh: float, rms: float = None, plot:bool = False): # IDL function also had an 'inspect' keyword which halts the function before returning?
    velocity = np.array(velocity)
    flux = np.array(flux)
    
    if rms == None:
        sel = np.sort(np.concatenate(np.argwhere(velocity < vlow), np.argwhere(velocity > vhigh)))
        rms = robust_sigma(flux[sel])
        
    sel1 = np.sort(np.concatenate(np.argwhere(velocity >= vlow), np.argwhere(velocity <= v0)))
    sel2 = np.sort(np.concatenate(np.argwhere(velocity >= v0), np.argwhere(velocity <= vhigh)))
    count1 = len(sel1)
    count2 = len(sel2)
    
    if count1 >= 2 and count2 >= 2:
        flow = np.sum(flux[sel1])
        eflow = rms * np.sqrt(count1)
        fhigh = np.sum(flux[sel2])
        efhigh = rms*np.sqrt(count2)
        
        ratio = flow/fhigh
        eratio = np.sqrt((eflow/fhigh)**2 + (flow*efhigh/fhigh**2)**2)
        
        asym = np.log10(ratio)
        easym = 1/np.log(10)*eratio/ratio
    
        if plot == True:
            plt.plot(velocity, flux, c='black')
            plt.axvline(v0, c = 'tab:green', lw = 3)
            plt.plot(velocity[sel1], flux[sel1], c='tab:red')
            plt.plot(velocity[sel2], flux[sel2], c='tab:blue')
            plt.axhline(1.5*rms, linestyle = 2)
            plt.show()
    else:
        asym = -999
        easym = -999
    
    return [asym, easym]

def asymmetry_himanga(id, asym, easym, v0 = None, path = None, plot = False):
    # id should be SDSS ID
    if path == None:
        path = './'
    par = path + id + '_par.sav' # TODO What is this file? IDL uses 'restore' which is not built into python the same way
    spec = path + 'mangaHI-' + id + '.fits'
    
    s = fits.open(spec)
    hdr = s[0].header
    if v0 == None:
        v0 = hdr['OBJ_VEL']
    rms = np.sqrt(np.mean(par**2)) # TODO this probably doesn't work, see earlier todo note
    
    if len(s[1]) >= 1000:
        print('Error in spectrum size. Skipping')
        return
    
    ### Remaining IDL in this procedure (not function, how to deal with?):
    # awvinds = [par.awvinfo.xmin,par.awvinfo.xmax] - round(150/4)
    # if awvinds[0] lt 0 or awvinds[1] lt 0 then return

    # vmin = min([s.vhi[awvinds[0]],s.vhi[awvinds[1]]])
    # vmax = max([s.vhi[awvinds[0]],s.vhi[awvinds[1]]])

    # out = profile_asymmetry(s.vhi,s.fhi,v0,vmin,vmax,rms=rms,plot=1)
    # asym = out[0]
    # easym = out[1]

def find_edges(velocity: ArrayLike, flux: ArrayLike, v0, rms = None, threshold = None, maxcount = None, reverse = False, plot = False, center_offset_start = None):
    if maxcount == None:
        maxcount = 3
    if threshold == None:
        threshold = 1.5
    if center_offset_start == None:
        center_offset_start = 50
    
    # Find index corresponding to v0
    ind0 = np.argwhere(velocity == v0)
