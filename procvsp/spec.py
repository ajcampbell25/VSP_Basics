
def spec_1d(data, timerange, frange, thead, trace, fs,twin, title_spec):    
    ''' Frequency Analysis of whole trace
        Uses scipy fft and avoids using powers of 2 for number of taps
        
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n'
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    
    import math 
    from math import ceil
    import procvsp.utils as Utils
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds
    
    print("\u0332".join('\nFrAn Parameters :'))    
    print ('fs :', fs,)

    # extract analysis trace
    data_single, thead_single = Utils.chosetrace(data, thead, trace)    
    data_single = data_single.T # samples in a column

    # useful headers
    TTobs_single = thead_single[:,8]
    zrcv_select = thead_single[:,2]
    trnum_single = thead_single[:,0]

    ####### optimize number of samples for transform  #############

    # for a segment of trace
    if (twin =='y')or(twin=='Y'):
        start = int(timerange[0]*(fs/1000))
        stop = int(start+(timerange[1]*(fs/1000)))
        data_trimd = data_single[start:stop,]
        
        # Apply window to segment
        w=sig.tukey(data_trimd.shape[0], alpha=.1)    
        data_win = data_trimd[:,0]*w        # Multiply trace by window
        N = scipy.fft.next_fast_len(int(timerange[1]-timerange[0]*(fs/1000)))
    
    # for a full trace
    else:
        # Apply window to whole trace
        w=sig.tukey(data_single.shape[0], alpha=.1)    
        data_win = data_single[:,0]*w        # Multiply trace by window    
        N = scipy.fft.next_fast_len(data_win.shape[0]) # in samples, best for scipy fft
    
    # pad with zeros if optimal N greater than trace or segment
    if (N > data_win.shape[0]):                    
        pad = N - data_win.shape[0]
        data_win = np.pad(data_win, (0,int(pad)), 'constant')

    ####### do the fft  #############
    
    X = scipy.fft.fft(data_win[:])      # from 0, to TT plus window 
                                                 # [:] if a window is used    
    X_db = 20*np.log10(np.abs(X)/np.max(np.abs(X)))# db=20*np.log10(S/np.max(S))
    
    freq = scipy.fft.fftfreq(X.shape[0], d=samprate)    # Generate plot frequency axis     
#    f = np.arange(0, N)*fs/N                # alternate method
                                                           
    ####### Only keep positive frequencies #########
    
    keep = freq>=0    
    X = X[keep]    
    X_db = X_db[keep]    
    freq = freq[keep]

    ############   make Spectral plots   ################
    
    plt.figure(figsize=(15,5))    
    ax1 = plt.subplot(121)    

    ax1.plot(freq, np.absolute(X), c = 'red')  # using fftfreq to get x axis    

    ax1.set_title('Amplitude Spectrum of %s at Depth %s'
                  %(title_spec, zrcv_select))    
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(frange[0], frange[1]) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    ax2 = plt.subplot(122)
    
    ax2.plot(freq,X_db, c='blue') #using number of samples and sample rate to get x axis
    
    ax2.set_title('Power Spectrum in db of %s at Depth %s'
                  %(title_spec, zrcv_select))    
    ax2.set_xlabel('Frequency hz')
    ax2.set_xlim(frange[0], frange[1]) # extents must be set   
    ax2.set_ylabel('Power(db)')    
    ax2.xaxis.grid()    
    ax2.yaxis.grid()

    plt.show() 

def bandpass_filter(data, lowcut, highcut, fs, order, N, QCP):
    '''
    for description of using Second Order Section (sos) instead of b,a
    see https://stackoverflow.com/questions/12093594/
    how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter
    recently updated
    '''
    
    import scipy.signal as sig
    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    #import scipy.fft
    #from scipy.signal import butter, lfilter, freqz
    
    # frequencies need to be normalized by nyquist frequency
    nyq = 0.5 * fs
    low = lowcut / nyq    
    high = highcut / nyq        
    sos_out = sig.butter(order, [low, high], analog=False, btype='band', output='sos')

    # create a time axis symmetric around zero
    dt = 1/fs
    t = np.arange(-N*dt/2,N*dt/2,dt)

    buttfilt = sig.sosfiltfilt(sos_out, data)
    buttfilt = np.float32(buttfilt) #temporary fix to type getting changed bug
    ordertest=order
    
    if (QCP == 'y') or (QCP =='Y'):
        # get a discrete time impulse response 
        center = N//2  #seems important to keep spike at middle of window    
        x = np.zeros(N)     
        x[center] = 1     
        coeff = sig.sosfiltfilt(sos_out, x)    

        fig = plt.figure(figsize=(15,6))    
        gs = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace = .15)
       
        ax1 = plt.subplot(gs[0])    
        for ordertest in [ordertest, ordertest*2, ordertest *3]:        
            sostest = sig.butter(ordertest, [low, high], analog=False, btype='band', output='sos')        
            w, h = sig.sosfreqz(sostest, worN=512)        
            ax1.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % ordertest)

        ax1.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')        
        ax1.set_xlim(0, highcut*4)       # half of nyquist should be good enough    
        ax1.set_xlabel('Frequency (Hz)')    
        ax1.set_ylabel('Gain')    
        ax1.grid(True)    
        ax1.legend(loc='best')        
        ax1.set_title('Butterworth Filter %sHz low cut, %sHz high cut, order =%s'
                      %(lowcut, highcut, order),fontsize=14)
    
        ax2 = plt.subplot(gs[1])    
        ax2.plot(t,coeff, c='red', label = ' %s Order Zero-phase Butterworth'
                 %(order) )                                              
        ax2.set_title('Butterworth Impulse Response %sHz lowcut, %sHz highcut, order=%s'
                      %(lowcut, highcut, order),fontsize=14)    
        ax2.set_xlabel('Time (s)')        
        ax2.legend(loc='lower left',borderaxespad=0, fontsize = 10)        
        ax2.grid(True)
        
        plt.show()

    return buttfilt    