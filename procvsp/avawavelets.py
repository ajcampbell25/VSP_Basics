def ricker( length,sample_rate, c_freq):
    """
    Generate time and amplitude values for a zero-phase wavelet. The second derivative of the Gaussian function
    or the third derivative of the normal-probability density function. A Ricker wavelet is often used as a zero-phase
    embedded wavelet in modeling and synthetic seismogram manufacture.

    Reference:
        Norman H. Ricker (1896–1980), American geophysicist.

    :param sample_rate: sample rate in seconds (float, int)
    :param length: length of time (dt) in seconds (float, int)
    :param c_freq: central frequency of wavelet (cycles/seconds or Hz). (float, int)
    :return: ndarray

    Usage:
        time, wavelet = (sample_rate,duration,c_freq)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    t_min = -length / 2
    t_max = (length - sample_rate) / 2
    wv_time = np.linspace(t_min, t_max, int(length / sample_rate))
    wv_amp = (1. - 2. * (np.pi ** 2) * (c_freq ** 2) * (wv_time ** 2)) * np.exp(
        -(np.pi ** 2) * (c_freq ** 2) * (wv_time ** 2))
    ''' 
    print (' \nricker :')
    print ( ' wv_amp.shape :', wv_amp.shape)
    print ( ' argmax wv_amp :', wv_amp.argmax())
    plt.figure(figsize=(7, 4))
    plt.plot(wv_time, wv_amp, lw=2, color='black', alpha=0.5)
    plt.fill_between(wv_time, wv_amp, 0, wv_amp > 0.0, interpolate=False, color='blue', alpha=0.5)
    plt.fill_between(wv_time, wv_amp, 0, wv_amp < 0.0, interpolate=False, color='red', alpha=0.5)
    plt.title('%d Hz Ricker wavelet' % c_freq, fontsize=14, pad=12)
    plt.xlabel('TWT (s)', fontsize=12)
    plt.ylabel('Amplitude', fontsize=12)
    plt.ylim((-1.1, 1.1))
    plt.xlim((min(wv_time), max(wv_time)))
    plt.grid()
    plt.show()
    '''
    return wv_time, wv_amp