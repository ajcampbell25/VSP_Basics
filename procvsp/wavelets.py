def ormsby(duration, dt, f):
    """
    https://github.com/agile-geoscience/bruges/blob/master/bruges/filters/wavelets.py
    
    The Ormsby wavelet requires four frequencies which together define a
    trapezoid shape in the spectrum. The Ormsby wavelet has several sidelobes,
    unlike Ricker wavelets.
    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (usually 0.001, 0.002,
            or 0.004).
        f (ndarray): Sequence of form (f1, f2, f3, f4), or list of lists of
            frequencies, which will return a 2D wavelet bank.
    Returns:
        ndarray: A vector containing the Ormsby wavelet, or a bank of them.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    f = np.asanyarray(f).reshape(-1, 1)
    
    print ('\nOrmsby Info')
    print (' duration :',duration,' dt :',dt, ' f :',f)

    try:
        f1, f2, f3, f4 = f
    except ValueError:
        raise ValueError("The last dimension must be 4")

    def numerator(f, t):
        return (np.sinc(f * t)**2) * ((np.pi * f) ** 2)

    pf43 = (np.pi * f4) - (np.pi * f3)
    pf21 = (np.pi * f2) - (np.pi * f1)

    t = np.arange(-duration/2, duration/2, dt)

    w = ((numerator(f4, t)/pf43) - (numerator(f3, t)/pf43) -
         (numerator(f2, t)/pf21) + (numerator(f1, t)/pf21))
    
    print (' presqueeze w shape :', w.shape)

    w = np.squeeze(w) / np.amax(w)
    
    # make filter plots
    
    plt.figure(figsize=(15, 5))    

    ax1 = plt.subplot(111)    
    ax1.plot(t, w, c = 'red')  # using fftfreq to get x axis    
    ax1.set_title('Ormsby wave low cut %s low pass %s high pass %s high cut %s'
                  %(f[0], f[1],f[2],f[3]))    
    ax1.set_xlabel('Time')
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    plt.show() 

    return w

def ricker(duration, dt, f):
    """
    https://github.com/agile-geoscience/bruges/blob/master/bruges/filters/wavelets.py
    
    Also known as the mexican hat wavelet, models the function:
    
    .. math::
        A =  (1 - 2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}
    If you pass a 1D array of frequencies, you get a wavelet bank in return.
    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (ndarray): Centre frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.

    Returns:
        ndarray. Ricker wavelet(s) with centre frequency f sampled on t.
    .. plot::
        plt.plot(bruges.filters.ricker(.5, 0.002, 40))
     
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    print ('\nRicker Info')
    print (' duration :',duration,' dt :',dt, ' f :',f)
    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2, duration/2, dt)
    pft2 = (np.pi * f * t)**2
    w = np.squeeze((1 - (2 * pft2)) * np.exp(-pft2))
    
    ############   make single trace plots   ################
    
    plt.figure(figsize=(15,5))
    
    ax1 = plt.subplot(111)    
    ax1.plot(t, w, c = 'red')  # using fftfreq to get x axis    
    ax1.set_title('Ricker wave %s Hertz'%(f))    
    ax1.set_xlabel('Time')    
#    ax1.set_xlim(0,200)    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    plt.show() 

    return w
    
def butterworth(lowcut, highcut, fs, N, order):

    '''Generate a butterworth wavelet.
     For description of using Second Order Section (sos) instead of b,a
     see https://stackoverflow.com/questions/12093594/
     how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter
    '''     
    print("\u0332".join('\nButterworth Wavelet Parameters :'))    
    print ('\n fs :',fs, ' N :',N,' order :',order)

    '''
    # The standard way 
    b, a = butter(order, [low, high], btype='band')    
    #get a discrete time impulse respnse    
    x = np.zeros(N)     
    x[center] = 1     
    coeff = filtfilt(b,a, x)
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import sosfiltfilt, butter

    dt = 1/fs
    t = np.arange((-N*dt)/2, (N*dt)/2, dt)
    nyq = 0.5 * fs        
    low = lowcut / nyq    
    high = highcut / nyq        

    sos_out = butter(order, [low, high], analog=False, btype='band', output='sos')
    
    center = (N//2)+1  # if N odd, keeps spike at middle of window    
    x = np.zeros(N)     
    x[center] = 1    
    print (' center :', center,' a few values of x :','\n', x[center -5:center+5])    
    coeff = sosfiltfilt(sos_out, x)
        
    # make single trace plots    
    plt.figure(figsize=(15,5))
    
    ax1 = plt.subplot(111)    
    ax1.plot(t,coeff, c = 'red')  # using fftfreq to get x axis
    ax1.set_title('Butterworth Wave %sHz lowcut, %sHz highcut order= %s'
                  %(lowcut, highcut,order))    
    ax1.set_xlabel('Time')    
    ax1.xaxis.grid()   
    ax1.yaxis.grid()

    plt.show()
    
    return coeff    