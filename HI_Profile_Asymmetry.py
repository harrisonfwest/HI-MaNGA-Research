# Python script to calculate global HI profile asymmetry, loosely based on that in /HI Profile/profile_asymmetry.pro
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.integrate import trapezoid
from astropy.io import fits

def profile_asymmetry(plateIFU: str, velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, width: float, plot: bool = False) -> ArrayLike:
    """
    Function to calculate asymmetry parameters of global gas profile

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        VHI (float): Gas-derived central velocty of the galaxy
        VOPT (float): Optically-derived central velocty of the galaxy
        width (ArrayLike): Effective widths of the galaxy along the velocity axis. Include WM50, WP50, WP20, W2P50, and WF50
        plot (bool, optional): Argument to determine if the chosen global gas profile will be plotted. Defaults to False.

    Returns:
        ArrayLike: List of parameters passed into function call, as well as integrated fluxes for each of left (low velocity) and right (high velocity) sides of the galaxy for each of the provided central velocities. Asymmetry itself is not directly calculated to allow for flexibility in subsequent analysis
    """
    # instead of float parameter width, use ArrayLike widths corresponding to (in order) WM50, WP50, WP20, W2P50, and WF50
    
    # for width, widthname in zip(widths, ['WM50', 'WP50', 'WP20', 'W2P50', 'WF50']):
        # and show/calculate values for all of these options? waiting for KLM feedback
    
    # Calculate rms to subtract from flux
    mean_v0 = (VHI + VOPT)/2
    rms_sel = np.transpose(np.argwhere(abs(velocity - mean_v0) >= width))[0][20:-20] # Trim first and last 20 velocities of baseline to remove noise
    rms = np.sqrt(np.mean(flux[rms_sel]**2))
    subtracted_flux = flux - rms
    
    lo_sel_HI = np.transpose(np.argwhere((velocity > (VHI - width)) & (velocity < VHI)))[0]
    hi_sel_HI = np.transpose(np.argwhere((velocity < (VHI + width)) & (velocity > VHI)))[0]

    lo_sel_OPT = np.transpose(np.argwhere((velocity > (VOPT - width)) & (velocity < VOPT)))[0]    
    hi_sel_OPT = np.transpose(np.argwhere((velocity < (VOPT + width)) & (velocity > VOPT)))[0]
    
    
    # We will integrate using the Trapezoidal rule, found in scipy.integrate.trapezoid
    
    # Low velocity side integrated flux centered at VHI
    lo_flux_HI = trapezoid(subtracted_flux[lo_sel_HI], velocity[lo_sel_HI])
    # High velocity side integrated flux centered at VHI
    hi_flux_HI = trapezoid(subtracted_flux[hi_sel_HI], velocity[hi_sel_HI])
    
    # Low velocity side integrated flux centered at VOPT
    lo_flux_OPT = trapezoid(subtracted_flux[lo_sel_OPT], velocity[lo_sel_OPT])
    # High velocity side integrated flux centered at VOPT 
    hi_flux_OPT = trapezoid(subtracted_flux[hi_sel_OPT], velocity[hi_sel_OPT])
    
    
    if plot == True:
        plt.figure(figsize= (10, 8))
        plt.plot(velocity, subtracted_flux, color = 'black')
        plt.axhline(0, c='gray', lw = 1)
        
        plt.axvline(VHI - width, lw = 2, color = 'tab:red', label = 'VHI-based width')
        plt.axvline(VHI + width, lw = 2, color = 'tab:red')
        plt.axvline(VHI, lw = 2, color = 'tab:red', linestyle = '--', label = 'VHI central velocity')
        
        plt.axvline(VOPT - width, lw = 2, color = 'tab:blue', label = 'VOPT-based width')
        plt.axvline(VOPT + width, lw = 2, color = 'tab:blue')
        plt.axvline(VOPT, lw = 2, color = 'tab:blue', linestyle = '--', label = 'VOPT central velocity')
        
        plt.legend()
        plt.title('Global HI profile of Plate-IFU: ' + plateIFU)
        plt.xlabel('Velocity')
        plt.ylabel('Flux')
        
        plt.show()

    return np.array([plateIFU, VHI, VOPT, width, lo_flux_HI, hi_flux_HI, lo_flux_OPT, hi_flux_OPT])