# Python script to calculate global HI profile asymmetry, loosely based on that in /HI Profile/profile_asymmetry.pro
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.integrate import trapezoid
from astropy.io import fits

def find_edges(velocity: ArrayLike, flux: ArrayLike, v_central: float, rms: float = None, threshold: float = 3, max_bins: int = 7) -> ArrayLike:
    """
    Algorithm to find edges of 1D spectrum (i.e. flux over a range of velocities) by detecting where the spectrum is reduced to noise

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        v_central (float): Central velocity of the galaxy
        rms (float, optional): Optional parameter of pre-calculated RMS noise. If no value is provided, it will be calculated from the spectrum. Defaults to None.
        threshold (float, optional): Sigma level below which spectrum is considered noisy (i.e. above how many multiples of the RMS is the data a detection?). Defaults to 3.
        max_bins (int, optional): Number of consecutive points required to be within the noise threshold for an edge to be reported. Defaults to 7.

    Returns:
        ArrayLike: Returns a two-element array of low and high velocity side edges corresponding indices along the velocity axis.
    """
    
    low_edge = None
    high_edge = None
    
    if rms == None:
        rms_sel = np.transpose(np.argwhere(abs(velocity - v_central) >= width))[0][20:-20] # Trim first and last 20 velocities of baseline to remove noise>
        rms = np.sqrt(np.mean(flux[rms_sel]**2))
        
    cen_vel = np.argwhere(velocity == v_central)
    low_vel = cen_vel - 200 # Assuming 200 indices is a reasonable range within which to find the edge of the galaxy
    high_vel = cen_vel + 200
    
    for i in range(cen_vel-6, low_vel):
        if flux[i:i+max_bins - 1].all <= rms * threshold:
            low_edge = i+max_bins-1
            break
    
    for i in range(cen_vel+6, high_vel):
        if flux[i:i-max_bins + 1].all <= rms * threshold:
            low_edge = i-max_bins+1
            break
    
    # If no suitable edges are found, choose 200 indices (~1000 km/s) from central velocity
    if low_edge == None:
        low_edge = low_vel
    if high_edge == None:
        high_edge = high_vel
    
    return [low_edge, high_edge]


def profile_asymmetry(plateIFU: str, velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, widths: ArrayLike, plot: bool = False) -> ArrayLike:
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
        ArrayLike: List of parameters passed into function call, as well as integrated fluxes for each of left (low velocity) and right (high velocity) sides of the galaxy
        for each of the provided central velocities. Asymmetry itself is not directly calculated to allow for flexibility in subsequent analysis
    """
    
    # instead of float parameter width, use ArrayLike widths corresponding to (in order) WM50, WP50, WP20, W2P50, and WF50
    
    color_count = 1 # For automatic coloring of same width-v0 paired lines (i.e. v0, v0-width, v0+width). See plotting block below
    integrated_fluxes = np.empty([0]) # List of pairs of integrated fluxes for low and high velocity sides about central vel (not grouped within list)

    if plot == True:
        plt.figure(figsize= (10, 8))
        plt.plot(velocity, flux, color = 'black')
        plt.axhline(0, c='gray', lw = 2)
        
        plt.title('Global HI profile of Plate-IFU: ' + plateIFU)
        plt.xlabel('Velocity')
        plt.ylabel('Flux')
        
        plt.axvline(VHI, lw = 2, color = 'lime', linestyle = '--', label = 'VHI central velocity')
        plt.axvline(VOPT, lw = 2, color = 'magenta', linestyle = '--', label = 'VOPT central velocity')
        
    for width, widthname in zip(widths, ['WM50', 'WP50', 'WP20', 'W2P50', 'WF50']):    
        # Calculate rms to subtract from flux
        mean_v0 = (VHI + VOPT)/2
        rms_sel = np.transpose(np.argwhere(abs(velocity - mean_v0) >= width))[0][20:-20] # Trim first and last 20 velocities of baseline to remove noise
        rms = np.sqrt(np.mean(flux[rms_sel]**2))
        subtracted_flux = flux - rms
        
        # low_edge_index, high_edge_index = find_edges(velocity, flux, mean_v0, rms)
        
        lo_sel_HI = np.transpose(np.argwhere((velocity > (VHI - width)) & (velocity < VHI)))[0]
        hi_sel_HI = np.transpose(np.argwhere((velocity < (VHI + width)) & (velocity > VHI)))[0]

        lo_sel_OPT = np.transpose(np.argwhere((velocity > (VOPT - width)) & (velocity < VOPT)))[0]    
        hi_sel_OPT = np.transpose(np.argwhere((velocity < (VOPT + width)) & (velocity > VOPT)))[0]
        
        
        # We will integrate using the Trapezoidal rule, implemented in scipy.integrate.trapezoid
        
        # Low velocity side integrated flux centered at VHI
        lo_flux_HI = trapezoid(subtracted_flux[lo_sel_HI], velocity[lo_sel_HI])
        integrated_fluxes = np.append(integrated_fluxes, lo_flux_HI)
        # High velocity side integrated flux centered at VHI
        hi_flux_HI = trapezoid(subtracted_flux[hi_sel_HI], velocity[hi_sel_HI])
        integrated_fluxes = np.append(integrated_fluxes, hi_flux_HI)
        
        # Low velocity side integrated flux centered at VOPT
        lo_flux_OPT = trapezoid(subtracted_flux[lo_sel_OPT], velocity[lo_sel_OPT])
        integrated_fluxes = np.append(integrated_fluxes, lo_flux_OPT)
        # High velocity side integrated flux centered at VOPT 
        hi_flux_OPT = trapezoid(subtracted_flux[hi_sel_OPT], velocity[hi_sel_OPT])
        integrated_fluxes = np.append(integrated_fluxes, hi_flux_OPT)
        
        
        if plot == True: # Put this inside loop of widths
            plt.axvline(VHI - width, lw = 2, color = 'C' + str(color_count), label = 'VHI-centered ' + widthname + ' width')
            plt.axvline(VHI + width, lw = 2, color = 'C' + str(color_count))
            color_count += 1
            
            plt.axvline(VOPT - width, lw = 2, color = 'C' + str(color_count), label = 'VOPT-centered ' + widthname + ' width')
            plt.axvline(VOPT + width, lw = 2, color = 'C' + str(color_count))
            color_count += 1
    
    if plot == True: # Leave this outside loop of widths
        plt.legend()
        plt.show()
        
    res = np.append(np.array([plateIFU, VHI, VOPT]), widths)
    res = np.append(res, integrated_fluxes)
    return res