def wavtrace(vsp, trheaders, wave, aligntime, numsamp, fs):
    
    """
    vsp: is the data
    ttime: is the travel time of down waves to be aligned with
    wave: is the wavelet - generated by a previous run
    aligntime: gives the time in ms to align wavelet on
    numsamp: number of samples in wavelet
    fs: is sample rate in hertz
    
    remember data arrays are in sample numbers - different sample rates have to be compensated
    """
    import numpy as np
    import matplotlib as plt
#    import procvsp.utils as utilvsp
    from procvsp.utils import shift  
    
    thead_copy = np.copy(trheaders)    
    ttime = thead_copy[:,8]    
    dt =1/fs *1000   # sample rate in ms    
    alignsamp = int(aligntime /dt)     # convert time to sample number    
    shifttime = int((numsamp/2)*dt)   #  array of times to align wavelets upon

    thead_copy[:,8] = shifttime

    wave = wave.T
    wave = wave.reshape(-1,wave.shape[0])
    
    # add zeros after the wavelet to make the trace the same length as the seismic trace
    zer  = np.zeros(vsp.shape[1] - wave.shape[1])    
    zer = zer.reshape(-1,zer.shape[0])    
    wavelets = np.hstack((wave, zer))         # zero pad the kernel to same length    
    wavelets_plot = wavelets.T.reshape(-1)
    plot_time = np.arange(0,wavelets_plot.shape[0],dt)

    # align wavelets with downgoing arrival time
    updown = 'down'
    waveshift,_ = shift(wavelets, thead_copy, updown, aligntime, fs) 
    waveshift_plot = waveshift.T.reshape(-1)
    waveshift_plot = waveshift_plot[:wavelets_plot.shape[0],]
    
    print("\u0332".join('\nWavetrace Parameters :'))
    print(' fs :',fs,' numsamp :',numsamp)
    print(' plot_time shape :',plot_time.shape)       
    print(' wave shape :',wave.shape, ' wavelets shape :', wavelets.shape,
          ' wavelets_plot shape :',wavelets_plot.shape,
          ' waveshift shape :' , waveshift.shape,
          ' waveshift_plot shape :', waveshift_plot.shape
          )    
    '''        
    plt.figure(figsize=(12,5))    

    ax1=plt.subplot(2,1,1)
    ax1.set_title("wavelets unshifted")
    ax1.plot(plot_time,wavelets_plot, c='black')
    ax1.set_xlabel('Time (ms)')
    ax1.xaxis.grid()
    ax1.yaxis.grid()
    ax1.set_xlim(0,plot_time.max())
    
    ax2=plt.subplot(2,1,2)
    ax2.set_title("wave shifted")
    ax2.plot(plot_time,waveshift_plot, c='black')
    ax2.set_xlabel('Time (ms)')
    ax2.xaxis.grid()
    ax2.yaxis.grid()
    ax2.set_xlim(0,plot_time.max())

    plt.tight_layout()

    plt.show()
    '''    
    return waveshift.T
    
def Weiner_waveshape_decon(trin_up,trdsign_dwn,theaders_all,**kwargs):    
    '''
    def deconw(trin=None,trdsign=None,n=None,stab=None,wndw=None,*args,**kwargs):
        varargin = deconw.varargin
        nargin = deconw.nargin
    DECONW: Wiener (Robinson) deconvolution    
    [trout,d]=deconw(trin,trdsign,n,stab,wndw)    
    Wiener (Robinson) deconvolution of the input trace. The operator is designed from the second
    input trace.    
    trin= input trace to be deconvolved
    trdsign= input trace to be used for operator design
    n= number of autocorrelogram lags to use (and length of inverse operator)
    stab= stabilization factor expressed as a fraction of the
         zero lag of the autocorrelation.
        ********* default= .0001 **********
    wndw= the type of window for the autocorrelation. 1 for boxcar, 2 for triangle, 
          3 for Gaussian
    ********** default =1 ************    
    trout= output trace which is the deconvolution of trin
    d= output deconvolution operator used to deconvolve trin
    The estimated wavelet is w=ifft(1./fft(d));    
    by: G.F. Margrave, May 1991
    
    My adaptations:
    Input Data:
    1. Align the downgoing waves at an arbitrary time (another function)
    2. Shift the upgoing waves to two-way time
    Decon: 
    1. Use a zero-phase Butterworth wavelet as the desired output.
    2. Design the filters on the downgoing waves, one operator per trace.
    3. Apply the operators to the corresponding upgoing traces.
    
    aligntime = downgoing wave alignment time
    backup = backup time before direct arrival
    winlngth = time widow starting at backup time for filter design
    stab = pre-whitening stabilization factor expressed as a fraction of the
           zero lag of the autocorrelation.
    lowcut = low freq. for Butterworth wavelet and post decon bandpass filter
    highcut = high freq. for Butterworth wavelet and post decon bandpass filter
    order = order for Butterworth wavelet and post decon bandpass filter
    wndw = apply a window to autocorr 1= do nothing 
           2 and 3 to be developed
           
    Useage comments:
    1. Stab should be in range .0001 to .001. 
    2. Stab affects relative amplitudes, stab = 0 gets deconned down similar 
       to wavelet
    3. Butterworth order should be 2 or 3, ringing gets bad for high orders
    4. The wavelet is normalized to have a max amplitude of 1, regardless of 
       bandpass parameters
    
    '''
    import math 
    from math import ceil
    import scipy.signal as sig
    from scipy.linalg import solve_toeplitz   
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    import procvsp.spec as specvsp
    import procvsp.decon as dec
    import procvsp.sigproc as sgp
    
    aligntime = kwargs['aligntime']
    wndw = kwargs['wndw']
    stab = kwargs['stab']
    fs = kwargs['fs']
    N =  kwargs['N']
    backup = kwargs['backup']
    winlngth = kwargs['winlngth']
    lowcut = kwargs['lowcut'] 
    highcut = kwargs['highcut'] 
    order = kwargs['order']
    norm_wav = kwargs['norm_wav']
    qc_plot = kwargs['qc_plot']
    qc_trc = kwargs['qc_trc'] - 1
    endplot = kwargs['endplot']

    print("\u0332".join('\nWeiner_waveshape Parameters :'))
    print (' fs  :', fs, '\n', 'N :', N,
           '\n','order :', order,'\n', 'aligntime :', aligntime)
            
    theaders_all_copy = np.copy(theaders_all)
    TTobs_all = theaders_all_copy[:,8]

    dt =1/fs *1000   # sample rate in ms
    nyq=0.5*fs       # nyquist frequency
    
    # transpose data matrices to 1 trace per column
    trdsign = trdsign_dwn.T    
    trin_up = trin_up.T    

    ##################### generate wavelet  ##########################
    #get Butterworth coefficients
    low = lowcut / nyq    
    high = highcut / nyq        
    sos_wave = sig.butter(order, [low, high], analog=False, 
                      btype='band', output='sos')
    # convolve coefficients with spike
    center = N//2  #seems important to keep spike at middle of window    
    x = np.zeros(N)     
    x[center] = 1     
    coeff = sig.sosfiltfilt(sos_wave, x)
    
    # apply a normalization to wavelet, so max amplitude  = 1 (optional)      
    nrm_factor = 1
    if (norm_wav =='y')or(norm_wav == 'Y'):
        wavemax = np.abs(coeff).max()
        nrm_factor = 1/wavemax   
    coeff = coeff * nrm_factor
    
    ################## window the downgoing data   #############
    
    backup = int(backup/dt)    
    winlngth = int(winlngth/dt)
    start = (int(ceil(aligntime) / dt))-backup    # direct arrival plus a window 
    stop = int(start + winlngth) *2               # Do for VSProwess comparison
    
    if (stop > trdsign.shape[0]):       # if window is longer than trace, pad        
        pad_win = stop-trdsign.shape[0]        
        trin_up = np.pad(trin_up, [(0, pad_win),(0,0)], 'constant', constant_values = 0)        
        trdsign = np.pad(trdsign, [(0, pad_win),(0,0)], 'constant', constant_values = 0)        
        wavetrace = dec.wavtrace(trdsign.T,theaders_all,coeff,aligntime, N, fs)        
        wave = wavetrace
    else:
        wavetrace = dec.wavtrace(trdsign_dwn,theaders_all_copy,
                         coeff,aligntime, N, fs)
        wave = wavetrace
        
    # Apply window to downgoing trace without chopping direct arrival
    trdwin = trdsign[start:stop,]
    # Window wavelet to same length as downgoing waves
    wave = wave[start:stop,]             
    wave = wave.reshape(-1)
    # find the max amplitude of the wavelet for qc plot annotation
    wavind = np.argmax(wave)
    wavmax = wave[wavind]
    wavelabelxy = (wavind, np.max(wavmax))
    wavelabel = 'Wavelet Amplitude %s'%(wavmax)

    ################## window and reshape downgoing   #############

    shape = (trdsign.shape[0], trdsign.shape[1])
    decon_dwn_all = np.zeros(shape,dtype=np.float32)
    decon_up_all = np.zeros(shape,dtype=np.float32)
    cc_all = np.zeros(shape,dtype=np.float32)
    ac_all = np.zeros(shape,dtype=np.float32)
    
    for k in range(0,(trdsign.shape[1])):
        ac = np.correlate(trdwin[:,k], trdwin[:,k],mode ='full')[len(trdwin)-1:]
        cc = np.correlate(wave,trdwin[:,k], mode = 'same')            
        #   stabilize the autocorrelation
        ac[0]=ac[0]*(1.0 + stab)
        # do the levinson recursion
        x = solve_toeplitz(ac,cc) # from scipy linalg       
        # apply the filter 
        decon_dwn = np.convolve(trdsign[:,k], x, mode = 'same')    
        decon_up = np.convolve(trin_up[:,k], x, mode = 'same')            
        decon_dwn_all[:,k] = decon_dwn[0:trdsign.shape[0]]        
        decon_up_all[:,k] = decon_up[0:trdsign.shape[0]]
        ac_all[:int(ac.shape[0]),k] = ac
        cc_all[:int(cc.shape[0]),k] = cc
        
    decon_dwn_filt=np.zeros(shape)    
    decon_up_filt=np.zeros(shape)
    for k in range(0,(decon_dwn_all.shape[1])):        
        decon_dwn_filt[:,k] = specvsp.bandpass_filter(decon_dwn_all[:,k], lowcut, 
                                                     highcut, fs, order, N, QCP='n')        
        decon_up_filt[:,k] = specvsp.bandpass_filter(decon_up_all[:,k], lowcut, 
                                                    highcut, fs, order, N, QCP='n')

    # find the max amplitude of the deconvolved downgoing for qc plot annotation
    decdwn_ind = np.argmax(decon_dwn_all[:,qc_trc])
    decdwn_max = decon_dwn_all[decdwn_ind,qc_trc]
    declabelxy = (decdwn_ind, np.max(wavmax))    
    declabel = 'Dec. Down Max. Amp %s'%(decdwn_max)

    # plot a single trace of data before and after decon as well as acor and xcor
    if (qc_plot == 'y')or(qc_plot == 'Y'):
        
        x_limit = int(endplot/dt) # plotting limit for test trace
        
        # create a dictionary so that names and data sets are tied together
        plot_data = {'Aligned Downgoing':trdsign[:x_limit,qc_trc],
                   'Upgoing ':trin_up[:x_limit,qc_trc],
                   'Windowed Downgoing':trdwin[:x_limit,qc_trc],
                   'Windowed Wavelet': wave[:x_limit],
                   'Auto correlation':ac_all[:x_limit,qc_trc],
                   'Cross correlation':cc_all[:int(winlngth/dt)*2,qc_trc],
                   'Deconvolved Up withBPF':decon_up_filt[:x_limit,qc_trc],                
                   'Deconvolved down':decon_dwn_all[:x_limit,qc_trc] }

        scal = np.max((decon_dwn_all[:x_limit,qc_trc]))/np.max(trdsign[:x_limit,qc_trc])
        
        # set up variables for reading dictionary
        titles = list(plot_data.keys())
        values = list(plot_data.values())
        
        #fix the centering of cc for plotting with symmetric time
        print (' np.argmax(cc_all[:,qc_trc]) :',np.argmax(cc_all[:,qc_trc]))
        values[5] = values[5][:int(winlngth/dt)*2]
        cc_start = values[5].shape[0]//2 - int(winlngth/dt)//2
        cc_end = values[5].shape[0]//2 + int(winlngth/dt)//2
        print (' values[5].shape[0] :',values[5].shape[0],' cc_start, cc_end : ', cc_start, cc_end )
        
        gs = gridspec.GridSpec(8, 1, height_ratios=[1,1,1,1,1,1,1,3], hspace = .5)
        fig = plt.figure(figsize=(14,21))
        
        for n in range(7):
            ax = fig.add_subplot(gs[n])
            ax.set_xlabel('Time(ms)')
            # make xtime longer for xcorr plot only
            if n==5:
                values[n] =values[n][cc_start:cc_end,]
                xtime = np.arange(-values[n].shape[0]//2,values[n].shape[0]//2,dt)
#                values[n] =values[n]
#                xtime = np.arange(0,values[n].shape[0],dt)
            else:
                xtime = np.arange(0,x_limit,dt)
                values[n] =values[n]

            ax.plot(xtime,values[n], c='black')
            ax.set_title(titles[n] + " For Trace %s"%(qc_trc))
            ax.xaxis.grid()
            ax.yaxis.grid()
            j=n
            if j==3:
                ax.annotate(wavelabel, xy=wavelabelxy, ha='center', va='center')
                                         
        ax = fig.add_subplot(gs[j+1])        

        ax.set_xlabel('Time(ms)') 
        ax.plot(xtime,values[j+1], c='black')
        ax.plot(xtime,values[j-j]*scal, c='brown')
        ax.annotate(declabel, xy=declabelxy, ha='center', va='center')
        ax.set_title(titles[j+1] + ' For Trace %s'
                     ', Pre-decon trace scaled for display'%(qc_trc))
        ax.set_xlim(0,x_limit//2)
        ax.xaxis.grid()
        ax.yaxis.grid()        
        plt.show()
            
    return decon_dwn_filt.T,decon_up_filt.T
    
def spike_decon_1trace(trin,trdsign,thead_single,fs,**kwargs):    

#def deconw(trin=None,trdsign=None,n=None,stab=None,wndw=None,*args,**kwargs):
#    varargin = deconw.varargin
#    nargin = deconw.nargin
# DECONW: Wiener (Robinson) deconvolution    
# [trout,d]=deconw(trin,trdsign,n,stab,wndw)    
# Wiener (Robinson) deconvolution of the input trace. The operator is designed from the second
# input trace.    
# trin= input trace to be deconvolved
# trdsign= input trace to be used for operator design
# n= number of autocorrelogram lags to use (and length of inverse operator)
# stab= stabilization factor expressed as a fraction of the
#       zero lag of the autocorrelation.
#      ********* default= .0001 **********
# wndw= the type of window for the autocorrelation. 1 for boxcar, 2 for triangle, 3 for Gaussian
# ********** default =1 ************    
# trout= output trace which is the deconvolution of trin
# d= output deconvolution operator used to deconvolve trin
# The estimated wavelet is w=ifft(1./fft(d));    
# by: G.F. Margrave, May 1991 

    import math 
    from math import ceil
    import scipy.signal as sig
    from scipy.linalg import solve_toeplitz   
    import numpy as np
    import matplotlib.pyplot as plt
    
    import procvsp.sigproc as sgp
    import procvsp.spec as specvsp
    
    Tmin = kwargs['stime'] 
    Tmax = kwargs['etime'] 
    lowcut = kwargs['lowcut']          # Butterworth filter parameters
    highcut = kwargs['highcut'] 
    order = kwargs['order'] 
    N = kwargs['numfsamp']  #= 1024
    wintap = kwargs['window']       # apply a window to prevent spectral leackage in fft
    wndw = kwargs['wndw']           # window type    
    stab = kwargs['stab'] 
    
    
    thead_one = np.copy(thead_single)
    TTobs_single = thead_one[:,8]
    
    Data_norm = 'n'
    DScalar = 1

    trdsign_nrm = sgp.normalize(trdsign, Data_norm,thead_one, DScalar)    
    trdsign = trdsign_nrm.T    
    trdwin  = trdsign[0:2000]      # could add a section to window the downgoing as done in waveshape decon
    trdwin = trdwin.reshape(-1)
    trdsign = trdsign.reshape(-1)

    ac = np.correlate(trdwin, trdwin,mode='full')[len(trdwin)-1:]     

    # window the autocorrelation
    if (wndw == 2):
        w=np.linspace(1,0,n)
        ac=multiply(ac,w)

    else:
        if (wndw == 3):
            #make g 2 sigma down at n
            sigma=n / 2
            g=exp(- (arange(0,n - 1)) ** 2 / sigma ** 2)
            a=multiply(ac,g)
            
    #   stabilize the autocorrelation
    ac[0]=ac[0]*(1.0 + stab)    
    # generate the right hand side of the normal equations
    b = np.hstack((np.ones(1), np.zeros(len (ac)-1)))    
    # do the levinson recursion
    x = solve_toeplitz(ac,b)
    # apply the filter
    decon = np.convolve(trdsign, x, mode = 'full')    
    decon = decon[0:len(trdsign)]
    
    buttfilt_decon = specvsp.bandpass_filter(decon, lowcut , highcut, fs, order, N, QCP='y')
    
    plt.figure(figsize=(14, 12))    

    ax1=plt.subplot(3,1,1)
    ax1.plot(trdsign)
    ax1.set_title("down")
    ax2=plt.subplot(3,1,2)
    ax2.plot(ac)
    ax2.set_title("auto correlation")
    ax3=plt.subplot(3,1,3)
    ax3.set_title("BPF spiking deconvolved Down and Input Down")
    ax3.plot(trdsign)
    ax3.plot(buttfilt_decon)
    plt.tight_layout()
    plt.show()

    return decon,buttfilt_decon
      