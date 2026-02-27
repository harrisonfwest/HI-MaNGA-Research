# Python script to calculate global HI profile asymmetry, loosely based on that in /HI Profile/profile_asymmetry.pro
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from astropy.io import fits

def profile_asymmetry(velocity: ArrayLike, flux: ArrayLike, VHI: float, VOPT: float, width: float | ArrayLike, plot: bool = False) -> ArrayLike:
    """
    Function to calculate asymmetry parameters of global gas profile

    Args:
        velocity (ArrayLike): Array of velocities (i.e. x-axis of gas profile spectrum)
        flux (ArrayLike): Array of fluxes (i.e. y-axis of gas profile spectrum)
        VHI (float): Gas-derived central velocty of the galaxy
        VOPT (float): Optically-derived central velocty of the galaxy
        width (float | ArrayLike): Effective width of the galaxy along the velocity axis. Taken as the largest of the galaxy widths given in the MaNGA data table (i.e. WM50, WP50, etc.)
        plot (bool, optional): Argument to determine if the chosen global gas profile will be plotted. Defaults to False.

    Returns:
        ArrayLike: List of parameters passed into function call, as well as integrated fluxes for each of left (low velocity) and right (high velocity) sides of the galaxy for each of the provided central velocities. Asymmetry itself is not directly calculated to allow for flexibility in subsequent analysis
    """
    v0 = np.array([VHI, VOPT])
    low_int_flux = -999 # Low velocity side integrated flux
    hi_int_flux = -999 # High velocity side integrated flux
    