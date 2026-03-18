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

def find_edges(velocity: ArrayLike, flux: ArrayLike, v0: float, threshold: float = 2, max_bins: int = 12) -> ArrayLike:
    """
    Algorithm to find edges of 1D spectrum (i.e. flux over a range of velocities) by detecting where the spectrum is reduced to noise
    for a given number of consecutive points

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        v0 (float): Central velocity of the galaxy
        threshold (float, optional): Sigma level below which spectrum is considered noisy (i.e. above how many multiples of the RMS is the data a detection?). Defaults to 2.
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
            
    cen_vel = np.argwhere(velocity == find_nearest(velocity, v0))[0, 0] # Index of velocity closest to v0
    low_vel = cen_vel - 100 # Assumes that the edges of the galaxy are each within 100 indices of v0 such that beyond 100 in either direction is all baseline
    high_vel = cen_vel + 100
        
    for i in np.arange(low_vel, cen_vel - max_bins)[::-1]:
        if all(x <= (rms * threshold) for x in flux[i : i + max_bins - 1]) and ((i + max_bins - 1) < cen_vel):
            low_edge = i# + max_bins - 1
            break
    
    for i in np.arange(cen_vel + max_bins, high_vel):
        if all(x <= (rms * threshold) for x in flux[i - max_bins + 1:i]):
            high_edge = i# - max_bins + 1
            break
    
    # If either a low or high edges fails to be found, choose 100 indices from central velocity as the edge
    if low_edge == None:
        low_edge = low_vel
    if high_edge == None:
        high_edge = high_vel
    
    res = [low_edge, high_edge]
    return res

def spectrum_analysis(plateIFU: str, velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, widths: ArrayLike, plot: bool = False) -> ArrayLike:
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
    v0_arr = [VHI, VOPT]
    
    flux_pairs = []
    velocity_index_pairs = []
    
    for v0 in v0_arr:
        for w in widths:
            lo_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 - w)))[0, 0]
            hi_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 + w)))[0, 0]
            arr = [lo_vel_index, hi_vel_index]
            velocity_index_pairs.append(arr)
        velocity_index_pairs.append(find_edges(velocity, flux, v0))
    
    # velocity_index_pairs now ordered as follows:
    # [VHI based WM50 indices], [VHI based WP50 indices], [VHI based WP20 indices], [VHI based W2P50 indices], [VHI based WF50 indices],
    # ... [VOPT based WM50 indices], [VOPT based WP50 indices], [VOPT based WP20 indices], [VOPT based W2P50 indices], [VOPT based WF50 indices]
    
    for i in range(5):
        lo_sel = np.transpose(np.argwhere((velocity > velocity[velocity_index_pairs[i][0]]) & (velocity < VHI)))[0]
        hi_sel = np.transpose(np.argwhere((velocity < velocity[velocity_index_pairs[i][1]]) & (velocity > VHI)))[0]
        lo_flux = trapezoid(flux[lo_sel], velocity[lo_sel])
        hi_flux = trapezoid(flux[hi_sel], velocity[hi_sel])
        flux_pairs.append([lo_flux, hi_flux])
        
    for i in range(5):
        lo_sel = np.transpose(np.argwhere((velocity > velocity[velocity_index_pairs[i+5][0]]) & (velocity < VOPT)))[0]
        hi_sel = np.transpose(np.argwhere((velocity < velocity[velocity_index_pairs[i+5][1]]) & (velocity > VOPT)))[0]
        lo_flux = trapezoid(flux[lo_sel], velocity[lo_sel])
        hi_flux = trapezoid(flux[hi_sel], velocity[hi_sel])
        flux_pairs.append([lo_flux, hi_flux])
    
    # flux_pairs now contains left and right integrated fluxes for ordering described above
    
    res = [plateIFU, v0_arr]
    res.append(velocity_index_pairs)
    res.append(flux_pairs)
    
    if plot:
        color_count = 1
        plt.figure(figsize= (10, 8), dpi= 500)
        plt.plot(velocity, flux, color = 'black')
        plt.axhline(0, c='gray', lw = 2)
        
        plt.title('Global HI profile of Plate-IFU: ' + plateIFU)
        plt.xlabel('Velocity')
        plt.ylabel('Flux')
        
        plt.axvline(VHI, lw = 2, color = 'lime', linestyle = '--', label = 'VHI central velocity')
        plt.axvline(VOPT, lw = 2, color = 'magenta', linestyle = '--', label = 'VOPT central velocity')
        
        width_names = ['WM50', 'WP50', 'WP20', 'W2P50', 'WF50', 'analytically calculated']
        
        for vels, width_name in zip(velocity_index_pairs[0:6], width_names):
            plt.axvline(velocity[vels[0]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count), label = 'VHI-based ' + width_name + ' edges')
            plt.axvline(velocity[vels[1]], linestyle = 'dashdot', lw = 1, color = 'C' + str(color_count))
            color_count += 1
        for vels, width_name in zip(velocity_index_pairs[6::], width_names):
            plt.axvline(velocity[vels[0]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count), label = 'VOPT-based ' + width_name + ' edges')
            plt.axvline(velocity[vels[1]], linestyle = 'dotted', lw = 1, color = 'C' + str(color_count))
            color_count += 1
        
        plt.legend()
        plt.show()
    
    return res