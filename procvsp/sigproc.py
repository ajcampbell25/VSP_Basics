def tar(VSP, thead, fs, exp):
    # Apply geometric spreading correction in the form of
    # T**exp
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    print("\u0332".join('\nTAR Parameters :'))

    numsamp = VSP.shape[1]
#    time =  np.arange(0,numsamp*(1000/fs),(1000/fs))    
    time =  np.arange(0,numsamp*(1/fs),(1/fs))    

    gainfunc = time**exp
    gainfunc = gainfunc.reshape(1,-1)
    
    tarred = VSP*gainfunc
    print (' VSP.shape :', VSP.shape,' gainfunc.shape :', gainfunc.shape)
    
    plt.figure(figsize=(15,7))    
    ax1 = plt.subplot(111)
   
    ax1.plot(time, gainfunc.reshape(-1,1), c = 'red')  # using fftfreq to get x axis    
    ax1.set_title('T**%s Gain Function'%(exp))   
    ax1.set_xlabel('Time(ms)')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()
    plt.show()
    
    return tarred
    
def medfilt (x, k):
    
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating zeros. Seems to be better edge results   
    see https://stackoverflow.com/questions/24585706/scipy-medfilt-wrong-result 

    """
    import scipy.signal as sig
    import numpy as np
   
    mfilt  = np.zeros(shape = (x.shape[0], x.shape[1]),dtype=np.float32)
    x = x.T    
    assert k % 2 == 1, "Median filter length must be odd."
    
    for b, trace in enumerate(x[:,:]):
        mfilt[:,b] = sig.medfilt(x[b], k) #scipy, in future test need for transpose
        
    return mfilt
    
def medfilt_across (x, k):
    
    """ from https://gist.github.com/bhawkins/3535131
        
    Apply a length-k median filter to a VSP data array.
    Boundaries are extended by repeating endpoints. at edge, input and output     
    traces are the same.     
    This can lead all zeros at edge traces after differencing.    
    Traces with all zeros lead to division errors when trying to normalize

    """
    import numpy as np

    out  = np.zeros(shape = (x.shape[0], x.shape[1]))
    x = x.T
    
    assert k % 2 == 1, "Median filter length must be odd."

    for b, trace in enumerate(x[:, :]):
        k2 = (k - 1) // 2                
        y = np.zeros ((len (x[b,]), k), dtype=x.dtype)
        y[:,k2] = x[b,]        
        row=b
        
        for i in range (k2):
            j = k2 - i
            y[j:,i] = x[row,:-j]
            y[:j,i] = x[row, 0]
            y[:-j,-(i+1)] = x[row,j:]
            y[-j:,-(i+1)] = x[row,-1]
            
        out[:,row] = np.median(y, axis = 1)
        
    return out   
    
def medfilt_along (VSP, k):
    
    """ from https://gist.github.com/bhawkins/3535131        
    Apply a length-k median filter to a VSP data array x.
    Boundaries are extended by repeating endpoints.

    """
 
    import numpy as np
    
    out = VSP
    
    assert k % 2 == 1, "Median filter length must be odd."

    for b, trace in enumerate(VSP[::1, :]):
        k2 = (k - 1) // 2                
        y = np.zeros ((len (VSP[b,]), k), dtype=VSP.dtype)        
        y[:,k2] = VSP[b,]        
        row=b

        for i in range (k2):
            j = k2 - i
            y[j:,i] = VSP[row,:-j]
            y[:j,i] = VSP[row, 0]
            y[:-j,-(i+1)] = VSP[row,j:]
            y[-j:,-(i+1)] = VSP[row,-1]
            
        out[row,] = np.median(y, axis = 1)

    return out
    
def normalize (Seismic, norm, thead, scal):
    """
    data normalization
    Frobenius norm is the default, L1 norm optional using data.sum
    row_sums = data.sum(axis=1)

    This method takes the whole trace, zeros before TT included so result 
    can be unexpected
    
    I need to create a method for normalizing in a window then apply 
    that norm factor to the whole trace
    stackoverflow.com/questions/8904694/how-to-normalize-a-2-dimensional-numpy-array-in-python-less-verbose
    """
    import numpy as np
    import matplotlib.pyplot as plt
     
#    print("\u0332".join('\nNormalization Stats :'))
    
    rdepth = thead[:,2]
    
    data2 = np.zeros(shape = (Seismic.shape[0], Seismic.shape[1]), dtype=np.float32)    
    #data1 = Seismic
    
    if (norm == 'Y') or (norm =='y'):        
        #row_sums = np.linalg.norm(data1, axis=1)        
        #print (' row_sums shape', row_sums.shape)        
        #data2 = (data1 / row_sums[:, np.newaxis])        
        #datascaled = data2 * scal

        ampmax = np.nanmax(np.abs(Seismic), axis=1) 
        data2 = (Seismic / ampmax[:, np.newaxis])        
        datascaled = data2 * scal
        
        amp_pre=np.amax(Seismic, axis=1) # for before-after plot
        amp_post=np.amax(datascaled, axis = 1) # for before-after plot
            
        plt.figure(figsize=(14,5))    
        ax1 = plt.subplot(121)
        ax1.plot(rdepth, ampmax, c = 'red')  # using fftfreq to get x axis    
        ax1.set_title('Normalization Factors')   
        ax1.set_xlabel('TVD Depth')    
        ax1.xaxis.grid()    
        ax1.yaxis.grid()
   
        ax2 = plt.subplot(122)
        ax3 = ax2.twinx()
    
        ax2.plot(rdepth, amp_pre, c = 'red', label='Pre Norm')         
        ax2.legend(loc='lower left',borderaxespad=0)#, fontsize = 7
        ax2.set_title('Maximum Amplitudes')   
        ax2.set_xlabel('TVD Depth')    
        ax2.set_ylabel('Pre Normalization (red)')    
        ax2.xaxis.grid()    
        ax2.yaxis.grid()
    
        ax3.plot(rdepth, amp_post, c = 'blue', label='Post Norm')  
        ax3.legend(loc='best',borderaxespad=0)#, fontsize = 7
        ax3.set_ylabel('Post Normalization (blue)') 
        ax3.yaxis.set_label_position('right')
        
        plt.show()
        
    else:        
        datascaled = Seismic * scal

    return datascaled
    
def normalize2(a, *, datarange = None):
    # not tested
    # from https://community.opengroup.org/osdu/platform/domain-data-mgmt-services
    # /seismic/open-zgy/-/blob/fff61f457d7d3ed90f239cafc5dc74a2d64400f8/python/openzgy/tools/viewzgy.py
    
    a = a.astype(np.float32)
    dead = np.isnan(a)
    amin, amax = datarange or (np.nanmin(a), np.nanmax(a))
    # Zero should be at the center
    if amin * amax < 0:
        x = max(abs(amin), abs(amax))
        amin, amax = (-x, x)
    # NaN and Inf show as smallest number
    a[dead] = amin
    if amin == amax:
        a *= 0
    else:
        # Avoid underflow, because app might have np.seterr(all='raise')
        a = a.astype(np.float64)
        a = (a - amin) / (amax - amin)
    a = (a * 255).astype(np.uint8)
    return a, dead

def cstack(upvsp, thead, fs, cwin, ctrnum, repeat):
    # Generate the a windowed data set and it's stack (corridor) 
    # A  tail mute is first applied at TWT plus cwin to get the corridor of data
    # The bottom ctrnum traces are left untouched
    # Headers must include twt in column 9
    
    import numpy as np
    import warnings

    cwinsamps= int(cwin *1000/fs)
    premtime = thead[:,8]*2*1000/fs
    mtime = (thead[:,8]*2*1000/fs)+cwinsamps
    
    # convert TWT and TWT + window to integer sample numbers 
    premtime = premtime.astype(int)
    mtime = mtime.astype(int)
    
#    corrmute = np.zeros(shape = (upvsp.shape[0], upvsp.shape[1]))#,dtype=np.float32) 
    corrmute = np.empty(shape = (upvsp.shape[0], upvsp.shape[1]))#,dtype=np.float32)
    corrmute = corrmute*np.nan
    
    emuted = corrmute.shape[0]-ctrnum
    snomute = corrmute.shape[0]-ctrnum
    enomute = corrmute.shape[0]

    # for each trace (k), copy the data between the time headers 
    # premtime and mtime  to an array of nans    
    for k in range(0,emuted):
        # copy data array up to break time+window into array of zeros   
        corrmute[k,premtime[k]:mtime[k]] = upvsp[k,premtime[k]:mtime[k]]

    for k in range(snomute,enomute):
        # copy data array from bottom traces   
        corrmute[k,premtime[k]:-1] = upvsp[k,premtime[k]:-1]

    # there will always be nans above the shallowest receiver, calculating the mean or median
    # of all nan samples will generate a warning which we can suppress
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
        corrstk= np.nanmedian(corrmute, axis=0).reshape(-1,1)
#        corrstk= np.nanmean(corrmute, axis=0).reshape(-1,1)
#        corrstk= np.sum(corrmute, axis=0).reshape(-1,1)
    
    print (' np.nanmax(corrmute) :', np.nanmax(corrmute), ' np.nanmax(corrstk) :', np.nanmax(corrstk))

    corrstk = np.repeat(corrstk,repeat,axis=1)

    return corrmute, corrstk.T

def attributes(data, fs):
    ''' Instantaneous attributes of seismic trace
    
    The analytic signal x_a(t) of signal x(t) is:

    where F is the Fourier transform, U the unit step function, 
    and y the Hilbert transform of x. [1]
    In other words, the negative half of the frequency spectrum
    is zeroed out, turning the real-valued signal into a complex signal. 
    The Hilbert transformed signal can be obtained from 
    np.imag(hilbert(x)), and the original signal from np.real(hilbert(x)).
    
    Use the scipy hilbert transform to get the complex valued analytic 
    trace. 
    The magnitude of the analytic signal gives the amplitude envelope.
    The instantaneous phase corresponds to the phase angle of the 
    analytic signal.
    The instantaneous frequency can be obtained by differentiating the 
    instantaneous phase in respect to time.
    
    Data - the VSP traces
    fs - sampling rate in hertz
        
    '''    
    import numpy as np
    
    from scipy.signal import hilbert
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    data_nrm = data//np.max(data)
    
    # generate a time axis vector    
    t = np.arange(0,data.shape[1], fs/1000).reshape(-1)
    tim1, tim2 = t.min(),t.max()
    
    analytic_signal = hilbert(data)
    amp_env = np.abs(analytic_signal)
    inst_phase = np.rad2deg(np.angle(analytic_signal))
    inst_phase_unwrap = np.unwrap(np.angle(analytic_signal))
    inst_frequency = (np.diff(inst_phase_unwrap) / (2.0*np.pi) * fs)
    
    print ("\u0332".join('\nAttribute Parameters :'))
    print (' fs :', fs, '\n',
           ' max trace amplitude : ', np.max(data),'\n',
           ' max amplitude envelope : ', np.max(amp_env),'\n' ,
           ' max inst frequency : ', np.max(inst_frequency),'\n', 
           ' max inst phase : ', np.max(inst_phase))
    print (' inst freq shape :', inst_frequency.shape)    
    print (' inst amp shape :', amp_env.shape)
    
    pad_freq = np.pad(inst_frequency, [(0, 0), (0, 1)], mode='constant', constant_values=0) # fix lost sample    
    
    mid_row = int((data.shape[0]/3))
    
    data_mid = data[mid_row:mid_row+1,].reshape(-1)
    amp_mid = amp_env[mid_row:mid_row+1,].reshape(-1)
    freq_mid = pad_freq[mid_row:mid_row+1,].reshape(-1)
    phase_mid = inst_phase[mid_row:mid_row+1,].reshape(-1)
    print (' freq_mid shape: ', freq_mid.shape,' amp_mid shape: ', amp_mid.shape)
    
    fig = plt.figure(figsize=(14,12))    
    gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,1], hspace = .25)    
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.plot(t,data_mid)
    ax1.plot(t,amp_mid,'k', linewidth=1)
    ax1.plot(t,-amp_mid,'k', linewidth=1)
    ax1.set_xlim(tim1,tim2)
    ax1.set_xlabel('Time(ms)')    
    ax1.set_title("raw amplitude and envelope for trace %s of %s"%(mid_row,data.shape[0]))
    ax1.grid('on')
    ax1.legend(labels=['seismic','envelope'])

    ax2.plot(t,phase_mid,'k', linewidth=1)
    ax2.set_xlim(tim1,tim2)
    ax2.set_xlabel('Time(ms)')
    ax2.set_title("unwrapped instantaneous phase for trace %sof %s"%(mid_row,data.shape[0]))    
    ax2.grid('on')

    ax3.plot(t,freq_mid,'k', linewidth=1)
    ax3.set_xlim(tim1,tim2)
    ax3.set_xlabel('Time(ms)')
    ax3.set_title("instantaneous frequency for trace %s of %s"%(mid_row,data.shape[0]))
    ax3.grid('on')

    ''' A test section based on 
    https://curvenote.com/@stevejpurves/geoscience/phaseandhilbert
    
    '''
    '''
    hilbert_trace = np.zeros_like(data)
    hilbert_trace = np.imag(hilbert(data))

    z = data + 1j*hilbert_trace
    env = np.abs(z)
    phase2 = np.angle(z)
    
    data_mid = data[mid_row:mid_row+1,].reshape(-1)
    hil_tr_mid = hilbert_trace[mid_row:mid_row+1,].reshape(-1)
    env_mid = env[mid_row:mid_row+1,].reshape(-1)
    phase2_mid = phase2[mid_row:mid_row+1,].reshape(-1)
    print (' hil_tr_mid shape: ', hil_tr_mid.shape)

    plt.subplot(5,1,4)
    plt.plot(t,data_mid)
    plt.plot(t, hil_tr_mid, ':')
    plt.plot(t, env_mid, 'k', linewidth=1)
    plt.plot(t, -env_mid, 'k', linewidth=1)
    plt.title('Trace Amplitudes')
    plt.grid('on')
    plt.xlim(tim1,tim2)
    plt.legend(labels=['seismic','hilbert','envelope'])

    plt.subplot(5,1,5)
    plt.plot(t,phase2_mid, 'k', linewidth=1)
    plt.title('Instantaneous Phase')
    plt.grid('on')
    plt.xlim(tim1,tim2)
    plt.xlabel('Time (ms)' )
    
    plt.tight_layout()
    plt.show()
    '''
 
    return amp_env, pad_freq, inst_phase 

