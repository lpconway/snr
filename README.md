# snr

    Calculate the signal-to-noise ratio according to European Pharmacopoeia
    guidelines. A Lorentzian or Gaussian lineshape is fitted to the signal in
    order to determine the half-height peak width and height of the peak. The
    amplitude of the noise is calculated over twenty peak widths around the
    fitted signal maximum, and the ratio of fitted peak height to noise amplitude
    is returned.

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
