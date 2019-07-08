# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:17:06 2019

@author: Louis
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def snr(
        time_data,
        signal_data,
        noise_data,
        time_range,
        line_shape='gaussian',
        sub_zero=True,
        plot=False,
        plot_path=None,
        title=""
        ):
    """
    Function to automatically calculate the signal-to-noise ratio according to
    European Pharmacopoeia guidelines. A Lorentzian or Gaussian lineshape is
    fitted to the signal in order to determine the half-height peak width and
    height of the peak. The amplitude of the noise is calculated over twenty
    peak widths around the fitted signal maximum, and the ratio of fitted peak
    height to noise amplitude is returned.

    Parameters
    ----------
    time_data : array_like
        Time points.
    signal_data : array_like
        Values representing the signal data
    noise_data: array-like
        Values representing noise data in the region of the signal, taken from
        a 'blank'.
    time_range: tuple or list of two values
        Time window in which the signal can be found. Should not contain any
        other signals.
    line_shape: string, optional
        If this is set to 'gaussian' (default), a Gaussian curve is fitted. If
        set to 'lorentzian', a lorentzian curve is fitted.
    sub_zero: bool, optional
        If the noise data can take values below zero, this should be set to
        true. If the noise can only take positive values, this should be false.
        True by default.
    plot: bool, optional
        If set to true, plots the fitted Gaussian and Lorenztian curves, along
        with the experimental data.
    plot_path: string, optional.
        If a string containing a path is provided, and plot=True, the plot is
        saved to this path.
    title: string, optional
        If a plot is generated, it will be given this title.

    Returns
    -------
    snr : float
        The signal-to-noise ratio.
    """
    window = (np.argmin(np.abs(time_data-time_range[0])),   # Calculate window indices
              np.argmin(np.abs(time_data-time_range[1])))

    x_data = time_data[window[0]:window[1]]                 # Time points within the window
    y_data = signal_data[window[0]:window[1]]               # Signal points within the window
    
    # Define lineshape functions
    
    def lorentzian(x_data, max_pos, width, height):
        """
        Returns Lorentzian lineshape as an array of y-values
        """
        x_data = np.array(x_data)
        x = 2*(max_pos-x_data)/width
        return height/(1+np.square(x))
    
    def gaussian(x_data, max_pos, width, height):
        """
        Returns a Gaussian lineshape as an array of y-values
        """
        x_data = np.array(x_data)
        x = 2*(max_pos-x_data)/width
        return height * np.exp(-0.693*np.square(x))
    
    initial_params = (np.average(time_range), 1, np.max(y_data))                            # Estimate initial parameters for lineshape fitting
    params = {'lorentzian' : curve_fit(lorentzian, x_data, y_data, p0=initial_params)[0],   # Fit Lorentzian curve
              'gaussian' : curve_fit(gaussian, x_data, y_data, p0=initial_params)[0]}       # Fit Gaussian curve
       
    [max_pos, width, height] = params[line_shape]                                           # Choose appropriate parameters for SNR calculation
    
    width = abs(width)                                                                      # Sometimes negatives widths can be fitted
    begin = max_pos - width*10                                                              # Calculate window for noise amplitude measurement
    end = max_pos + width*10
    
    noise_window = (noise_data[np.argmin(np.abs(time_data-begin)):np.argmin(np.abs(time_data-end))])    # Convert timepoint window into indices window
    
    signal = height                                 # Use fitted peak height for SNR calculation; this appears to be implied but not stated by European Pharmacopoeia guidelines
    noise = max(noise_window) - min(noise_window)   # Peak to peak height of the noise within the 20 half-height peak widths
    snr = (sub_zero+1)*signal/noise                  # Calculate SNR
    
    if plot:                                        # Plot fitted curves and data, if required
        fig, axis = plt.subplots(dpi=72)
        axis.plot(x_data, y_data)
        [max_pos, width, height] = params['lorentzian']
        axis.plot(x_data, lorentzian(x_data, max_pos, width, height))
        [max_pos, width, height] = params['gaussian']
        axis.plot(x_data, gaussian(x_data, max_pos, width, height))
        axis.set_title(title + ' S/N: %8.2f' % snr)
        plt.legend(['Experimental', 'Lorentzian fit', 'Gaussian fit'],  frameon=False)
        axis.set_xlabel("Retention time / min")
        axis.set_ylabel('Intensity')
        axis.set_xlim([min(x_data), max(x_data)])
        axis.set_ylim([min(y_data), max(y_data)*1.2])
        if plot_path is not None:                   # Save plot, if required
            fig.savefig(plot_path, dpi=300)

    return snr