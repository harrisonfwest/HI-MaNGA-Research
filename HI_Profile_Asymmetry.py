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

def find_edges(velocity: ArrayLike, flux: ArrayLike, v_central: float, threshold: float = 5, max_bins: int = 12) -> ArrayLike:
    """
    Algorithm to find edges of 1D spectrum (i.e. flux over a range of velocities) by detecting where the spectrum is reduced to noise

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        v_central (float): Central velocity of the galaxy
        rms (float, optional): Optional parameter of pre-calculated RMS noise. If no value is provided, it will be calculated from the spectrum. Defaults to None.
        threshold (float, optional): Sigma level below which spectrum is considered noisy (i.e. above how many multiples of the RMS is the data a detection?). Defaults to 5.
        max_bins (int, optional): Number of consecutive points required to be within the noise threshold to define an edge. Defaults to 12.

    Returns:
        ArrayLike: Returns a two-element array of low and high velocity side edges corresponding indices along the velocity axis.
    """
    
    low_edge = None
    high_edge = None
    
    rms_sel = np.transpose(np.argwhere(abs(velocity - v_central) >= 750))[0][20:-20] # Trim first and last 20 velocities of baseline to remove artifacts, and assume 750 km/s away from center is baseline noise
    rms = np.sqrt(np.mean(flux[rms_sel]**2))
        
    cen_vel = np.argwhere(velocity == find_nearest(velocity, v_central))[0, 0]
    low_vel = cen_vel - 200 # Assuming 200 indices is a reasonable range within which to find the edge of the galaxy
    high_vel = cen_vel + 200
    
    for i in range(cen_vel - max_bins + 1, low_vel):
        if flux[i:i+max_bins - 1].all() <= rms * threshold:
            low_edge = i+max_bins-1
            break
    
    for i in range(cen_vel + max_bins - 1, high_vel):
        if flux[i:i-max_bins + 1].all() <= rms * threshold:
            low_edge = i-max_bins+1
            break
    
    # If no suitable edges are found, choose 200 indices (~1000 km/s) from central velocity
    if low_edge == None:
        low_edge = low_vel
    if high_edge == None:
        high_edge = high_vel
    
    return [low_edge, high_edge]

def spectrum_analysis(plateIFU: str, velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, widths = ArrayLike, plot: bool = False) -> ArrayLike:
    v0_arr = [VHI, VOPT]
    flux_pairs = []
    
    velocity_index_pairs = []
    for v0 in v0_arr:
        for w in widths:
            lo_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 - w)))[0, 0]
            hi_vel_index = np.argwhere(velocity == find_nearest(velocity, round(v0 + w)))[0, 0]
            arr = [lo_vel_index, hi_vel_index]
            velocity_index_pairs.append(arr)
        arr = find_edges(velocity, flux, v0)
        velocity_index_pairs.append(arr)
    
    # velocity_index_pairs now ordered as follows:
    # [VHI based WM50 indices], [VHI based WP50 indices], [VHI based WP20 indices], [VHI based W2P50 indices], [VHI based WF50 indices],
    # [VOPT based WM50 indices], [VOPT based WP50 indices], [VOPT based WP20 indices], [VOPT based W2P50 indices], [VOPT based WF50 indices]
    
    # TODO: resume debugging here
    for i in range(5):
        lo_sel = np.transpose(np.argwhere((velocity > velocity_index_pairs[i][0]) & (velocity < VHI)))[0]
        hi_sel = np.transpose(np.argwhere((velocity < velocity_index_pairs[i][1]) & (velocity > VHI)))[0]
        lo_flux = trapezoid(flux[lo_sel], velocity[lo_sel])
        hi_flux = trapezoid(flux[hi_sel], velocity[hi_sel])
        flux_pairs.append([lo_flux, hi_flux])
        
    for i in range(5):
        lo_sel = np.transpose(np.argwhere((velocity > velocity_index_pairs[i+5][0]) & (velocity < VOPT)))[0]
        hi_sel = np.transpose(np.argwhere((velocity < velocity_index_pairs[i+5][1]) & (velocity > VOPT)))[0]
        lo_flux = trapezoid(flux[lo_sel], velocity[lo_sel])
        hi_flux = trapezoid(flux[hi_sel], velocity[hi_sel])
        flux_pairs.append([lo_flux, hi_flux])
    
    # flux_pairs now contains left and right integrated fluxes for ordering described above
    
    res = [plateIFU, v0_arr]
    res.append(velocity_index_pairs)
    res.append(flux_pairs)
    
    if plot:
        color_count = 1
        plt.figure(figsize= (10, 8))
        plt.plot(velocity, flux, color = 'black')
        plt.axhline(0, c='gray', lw = 2)
        
        plt.title('Global HI profile of Plate-IFU: ' + plateIFU)
        plt.xlabel('Velocity')
        plt.ylabel('Flux')
        
        plt.axvline(VHI, lw = 2, color = 'lime', linestyle = '--', label = 'VHI central velocity')
        plt.axvline(VOPT, lw = 2, color = 'magenta', linestyle = '--', label = 'VOPT central velocity')
        
        width_names = ['WM50', 'WP50', 'WP20', 'W2P50', 'WF50', 'analytically calculated']
        
        for vels, width_name in zip(velocity_index_pairs, width_names):
            if color_count < 6:
                cen = 'VHI'
            elif color_count >= 6:
                cen = 'VOPT'
            plt.axvline(velocity[vels[0]], color = 'C' + str(color_count), label = cen + '-based ' + width_name + ' edges')
            plt.axvline(velocity[vels[1]], color = 'C' + str(color_count))
            color_count += 1
        
        plt.legend()
        plt.show()
    
    return res