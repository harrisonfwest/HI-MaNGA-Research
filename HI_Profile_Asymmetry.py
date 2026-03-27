# Python script to calculate global HI profile asymmetry, loosely based on that in /HI Profile/profile_asymmetry.pro
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.integrate import trapezoid

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_edges(velocity: ArrayLike, flux: ArrayLike, v0: float, max_bins: int = 12) -> ArrayLike:
    """
    Algorithm to find edges of 1D spectrum (i.e. flux over a range of velocities) by detecting where the spectrum is reduced to noise
    for a given number of consecutive points

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        v0 (float): Central velocity of the galaxy
        max_bins (int, optional): Number of consecutive points required to be within the noise threshold to define an edge. Defaults to 12.

    Returns:
        ArrayLike: Returns a two-element array of low and high velocity side edge indices corresponding indices along the velocity axis.
    """
    
    low_edge = None
    high_edge = None
    
    cen_vel = np.argwhere(velocity == find_nearest(velocity, v0))[0, 0] # Index of velocity closest to v0
    
    low_rms_sel = np.arange(20, cen_vel - 100)
    hi_rms_sel = np.arange(cen_vel + 100, len(velocity) - 20)
    rms_sel = np.append(low_rms_sel, hi_rms_sel)
    rms = np.sqrt(np.mean(flux[rms_sel]**2))
            
    low_vel = cen_vel - 100 # Assumes that the edges of the galaxy are each within 100 indices of v0 such that beyond 100 in either direction is all baseline
    high_vel = cen_vel + 100

    for i in np.arange(low_vel, cen_vel - max_bins)[::-1]:
        if all(x <= (rms * 2) for x in flux[i: i + max_bins - 1]):
            low_edge = i + max_bins - 1
            break
    
    for i in np.arange(cen_vel + max_bins, high_vel):
        if all(x <= (rms * 2) for x in flux[i - max_bins + 1:i]):
            high_edge = i - max_bins + 1
            break
    
    res = [low_edge, high_edge]
    
    return res

def spectrum_analysis(plateIFU: str, velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, widths: ArrayLike, MA: bool, row_ind: int, LOGMSTARS: float, LOGMHI: float, HI_is_lim: bool, plot: bool = False) -> ArrayLike:
    """ Function to calculate relative integrated fluxes about central velocity of galaxy from global HI profile. Preserves and returns input parameters and their derived values.

    Args:
        plateIFU (str): Galaxy observation plate IFU number for identification in HI-MaNGA database
        velocity (ArrayLike): Array of velocities; x-axis of global HI profile
        flux (ArrayLike): Array of fluxes; y-axis of global HI profile
        VHI (float): HI-derived central velocity of galaxy along velocity axis
        VOPT (float): Optically derived central velocity of galaxy along velocity axis
        widths (ArrayLike): List of different measured galaxy widths along velocity axis. Found in HI-MaNGA data table
        plot (bool, optional): Option to display a plot of the global profile, along with vertical lines indicating where galaxy edges are reported or have been found. Defaults to False.

    Returns:
        ArrayLike: Returns an array of Plate IFU number; HI and Optically derived central velocities; 
        pairs of velocities indices corresponding to edges both from given widths and calculated edges, for both central velocites;
        and pairs of relative integrated fluxes for the provided edges.
        Edge index and flux pairs are ordered as follows: VHI-WM50, VHI-WP50, VHI-WP20, VHI-W2P50, VHI-WF50, VHI-calculated,
        VOPT-WM50, VOPT-WP50, VOPT-WP20, VOPT-W2P50, VOPT-WF50, VOPT-calculated
    """

    # Apparently, the direction of the velocity array varies by file. Must make sure it's uniform to correctly analyze
    if velocity[0] > velocity[1]:
        velocity = velocity[::-1]
        flux = flux[::-1]
    
    v0_arr = [VHI, VOPT]
    
    flux_pairs = []
    velocity_index_pairs = []
    
    for v0 in v0_arr:
        for w in widths:
            if w == -999.0:
                arr = [None, None]
            else:
                lo_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 - w)))[0, 0]
                hi_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 + w)))[0, 0]
                arr = [lo_vel_index, hi_vel_index]
            velocity_index_pairs.append(arr)
        velocity_index_pairs.append(find_edges(velocity, flux, v0))
    
    '''
    velocity_index_pairs are ordered as follows:
    [VHI based WM50 indices], [VHI based WP50 indices], [VHI based WP20 indices], [VHI based W2P50 indices], [VHI based WF50 indices], [VHI based analytical indices]
    ... [VOPT based WM50 indices], [VOPT based WP50 indices], [VOPT based WP20 indices], [VOPT based W2P50 indices], [VOPT based WF50 indices], [VOPT based analytical indices]
    '''

    c = 0
    for vel_pair in velocity_index_pairs:
        if any(x == None for x in vel_pair):
            flux_pairs.append([None, None])
            c += 1
            continue
        
        if c <= 5:
            cen_vel = np.argwhere(velocity == find_nearest(velocity, VHI))[0, 0]
            c += 1
        else:
            cen_vel = np.argwhere(velocity == find_nearest(velocity, VOPT))[0, 0]
        
        lo_sel = np.arange(vel_pair[0], cen_vel)
        hi_sel = np.arange(cen_vel + 1, vel_pair[1])
        
        lo_flux = np.trapz(flux[lo_sel], velocity[lo_sel])
        hi_flux = np.trapz(flux[hi_sel], velocity[hi_sel])
        
        flux_pairs.append([lo_flux, hi_flux])
    
        
    res = [row_ind, plateIFU, v0_arr]
    res += velocity_index_pairs
    res += flux_pairs
    res += [MA, LOGMSTARS, LOGMHI, HI_is_lim]
        
    width_names = ['WM50', 'WP50', 'WP20', 'W2P50', 'WF50']
    
    if plot:
        color_count = 1
        plt.figure(figsize= (10, 8), dpi= 500)
        plt.plot(velocity, flux, color = 'black')
        
        plt.title('Global HI profile of Plate-IFU: ' + plateIFU)
        plt.xlabel(r'Velocity [km s$^{-1}$]')
        plt.ylabel('Flux [Jy]')
        
        plt.axvline(VHI, lw = 2, color = 'lime', linestyle = 'dashdot', label = 'VHI central velocity')
        plt.axvline(VOPT, lw = 2, color = 'magenta', linestyle = 'dotted', label = 'VOPT central velocity')
                
        for vels, width_name in zip(velocity_index_pairs[0:5], width_names):
            if not any(x == None for x in vels):
                plt.axvline(velocity[vels[0]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count), label = 'VHI-based ' + width_name + ' edges')
                plt.axvline(velocity[vels[1]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count))
            color_count += 1
        for vels, width_name in zip(velocity_index_pairs[6:11], width_names):
            if not any(x == None for x in vels):
                plt.axvline(velocity[vels[0]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count), label = 'VOPT-based ' + width_name + ' edges')
                plt.axvline(velocity[vels[1]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count))
            color_count += 1
        
        if not any(x == None for x in velocity_index_pairs[5]):
            vels = velocity_index_pairs[5]
            plt.axvline(velocity[vels[0]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count), label = 'VHI-based analytically calculated edges')
            plt.axvline(velocity[vels[1]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count))
        color_count += 1
        
        if not any(x == None for x in velocity_index_pairs[11]):
            vels = velocity_index_pairs[11]
            plt.axvline(velocity[vels[0]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count), label = 'VOPT-based analytically calculated edges')
            plt.axvline(velocity[vels[1]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count))
        
        plt.legend()
        plt.show()
    
    return res