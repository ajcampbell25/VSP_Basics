def q_simulation(VSP, timerange, frange, thead, trace, fs,twin, title_spec,**kwargs):
    ''' Run a Frequency Analysis on a selected trace and then apply a 
    'forward Q filter' to the spectrum to simulate the effect of Q
        
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n'
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter
        
        The frequency spectrum is 

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    
    import math 
    from math import ceil
    import procvsp.utils as Utils
    
    Q=kwargs["Q"]
    vel=kwargs["Vel"]
    Zint=kwargs["Zint"]
    save=kwargs["savepng"]
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds
    
    print("\u0332".join('\nFrAn Parameters :'))    
    print ('fs :', fs,)

    # extract analysis trace
    data_single, thead_single = Utils.chosetrace(VSP, thead, trace)    
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
        print (' start:',start,' stop:',stop)        
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

    freq = scipy.fft.fftfreq(X.shape[0], d=samprate)    # Generate plot frequency axis     
#    f = np.arange(0, N)*fs/N                # alternate method
                                                           
    ####### Only keep positive frequencies #########
    
    keep = freq>=0    
    X = X[keep]    

    freq = freq[keep]
    
    ############ apply Q to spectrum  ##########################
    
    Z0 = zrcv_select
    Z=Z0+Zint
    
    Qfwd= X*np.exp(-1*(np.pi*freq)/(Q*vel)*(Z-Z0))
    
    ########### use spectral ratios to get Q from Q filtered spectrum

    specrat=np.log(np.absolute(Qfwd)/np.absolute(X))  
    # get sample number for max frequency to control plots
    fmax_ind=np.where(freq<=frange[1])
    specrat = specrat[fmax_ind]
    freqrat = freq[fmax_ind]
    
    slopew = np.polyfit(freqrat,specrat,  1)# 1 means first order - linear regression
    
    T=Z/vel
    T0=Z0/vel
    
    Qeff=-1*(np.pi*((T-T0)/slopew[0]))
    print (' slopew :', slopew, ' slopew.shape :',slopew.shape, 'Qeff :',Qeff)
    
    ############   make Q demonstration plots   ################
    
    fig=plt.figure(figsize=(15,5))    
    ax1 = plt.subplot(121)    

    ax1.plot(freq, np.absolute(X), c = 'red', label = "Raw Spectrum")  # using fftfreq to get x axis    
    ax1.plot(freq, np.absolute(Qfwd), c = 'blue', label = "Q Filtered Spectrum")  # using fftfreq to get x axis    
    ax1.legend(loc='upper right',borderaxespad=0, fontsize = 8)        

    ax1.set_title('Spectrum Before and after Q Filtering - Reference depth %s'
                  %(Z0[0]))
    ax1.text(.5,.95,"Q input = %s,\nSimulated target depth %s"%(Q,Z[0]),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(frange[0], frange[1]) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')   
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    ax2 = plt.subplot(122)
    
    ax2.scatter(freqrat,specrat, c='blue', label = "Spectral Ratio") #using number of samples and sample rate to get x axis
    
    ax2.set_title('Spectral Ratio from Reference Depth %s to Simulated Depth %s'
                  %(Z0[0],Z[0]))
    ax2.text(.7,.95,"Slope = %.4f, Q calculated = %.0f"%(slopew[0],Qeff[0]),
                  horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    ax2.set_xlabel('Frequency hz')
    ax2.set_xlim(frange[0], frange[1]) # extents must be set   
    ax2.set_ylabel('Log Spectral Ratio')    
    ax2.xaxis.grid()    
    ax2.yaxis.grid()

    plt.show()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\Q _%s_modelling.png' 
        %(Q), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)

def qest(VSP, thead,fs, **kwargs):
    ''' Run a Frequency Analysis on two selected traces.
    Estimate the spectral ratios in a frequency band. From the
    slope of the psectral ratios, calculate Q
    
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n'
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter
        
        The frequency spectrum is 

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    
    import math 
    from math import ceil
    import procvsp.utils as Utils
    
    twin = kwargs["apply_window"]
    tmin = kwargs["time_min"]
    tmax = kwargs["time_max"]
    fmin = kwargs["freq_min"]
    fmax= kwargs["freq_max"]
    srfmin = kwargs["srfreq_min"]
    srfmax= kwargs["srfreq_max"]    
    trace_1=kwargs["trace_1"]
    trace_2=kwargs["trace_2"]
    save=kwargs["savepng"]
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds
    
    print("\u0332".join('\nFrAn Parameters :'))    
    print ('fs :', fs,)

    # extract analysis traces
    data_ref, thead_ref = Utils.chosetrace(VSP, thead, trace_1)    
    data_ref = data_ref.T # samples in a column

    data_two, thead_two = Utils.chosetrace(VSP, thead, trace_2)    
    data_two = data_two.T # samples in a column

    # useful headers
    TTobs_ref = thead_ref[:,8]
    zrcv_ref = thead_ref[:,2]
    trnum_ref = thead_ref[:,0]
    
    TTobs_two = thead_two[:,8]
    zrcv_two = thead_two[:,2]
    trnum_two = thead_two[:,0]


    ####### optimize number of samples for transform  #############

    # for a segment of trace
    if (twin =='y')or(twin=='Y'):
        start = int(tmin*(fs/1000))
        stop = int(start+(tmax*(fs/1000)))
        data_trimd_ref = data_ref[start:stop,]
        data_trimd_two = data_two[start:stop,]
        
        # Apply window to segment
        w=sig.tukey(data_trimd_ref.shape[0], alpha=.1)    
        data_win_r = data_trimd_ref[:,0]*w        # Multiply trace by window
        data_win_2 = data_trimd_two[:,0]*w        # Multiply trace by window
        N = scipy.fft.next_fast_len(int(tmax-tmin*(fs/1000)))
    
    # for a full trace
    else:
        # Apply window to whole trace
        w=sig.tukey(data_ref.shape[0], alpha=.1)    
        data_win_r = data_ref[:,0]*w        # Multiply trace by window    
        data_win_2 = data_two[:,0]*w        # Multiply trace by window    
        N = scipy.fft.next_fast_len(data_win_r.shape[0]) # in samples, best for scipy fft
    
    # pad with zeros if optimal N greater than trace or segment
    if (N > data_win_r.shape[0]):                    
        pad = N - data_win_r.shape[0]
        data_win_r = np.pad(data_win_r, (0,int(pad)), 'constant')
        data_win_2 = np.pad(data_win_2, (0,int(pad)), 'constant')

    ####### do the fft  #############
    
    X_ref = scipy.fft.fft(data_win_r[:])      # from 0, to TT plus window 
                                                 # [:] if a window is used    
    X_two = scipy.fft.fft(data_win_2[:])      # from 0, to TT plus window 
                                                 # [:] if a window is used    
    freq = scipy.fft.fftfreq(X_ref.shape[0], d=samprate)    # Generate plot frequency axis     
#    f = np.arange(0, N)*fs/N                # alternate method
                                                           
    ####### Only keep positive frequencies #########
    
    keep = freq>=0    
    X_ref = X_ref[keep]       
    X_two = X_two[keep]    
    
    freq = freq[keep]
    
    ############ get time and depth between receivers  ##########################

    Z0 = zrcv_ref
    Z = zrcv_two

    T=TTobs_two/1000
    T0=TTobs_ref/1000    

    ########### use spectral ratios to get Q from Q filtered spectrum

    specrat=np.log(np.absolute(X_two)/np.absolute(X_ref))  
    # get sample number for max frequency to control plots
    fmax_ind=np.where((freq>=fmin)&(freq<=fmax))
    srfmax_ind=np.where((freq>=srfmin)&(freq<=srfmax))

    slopew = np.polyfit(freq[srfmax_ind],specrat[srfmax_ind],  1)# 1 means first order - linear regression
    plot_slope = slopew[0]*freq[srfmax_ind]+slopew[1] 

    Qeff=-1*(np.pi*((T-T0)/slopew[0]))
    
    # trendpoly = np.poly1d(slopew) # an alternative for plotting instead of y=mx+b 
    
    print (' slopew :', slopew, ' slopew.shape :',slopew.shape, 'Qeff :',Qeff)
    
    ############   make Q demonstration plots   ################
    
    fig=plt.figure(figsize=(15,5))    
    ax1 = plt.subplot(121)    

    ax1.plot(freq, np.absolute(X_ref), c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis    
    ax1.plot(freq, np.absolute(X_two), c = 'blue', label = "Target Depth Spectrum")  # using fftfreq to get x axis    
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        

    ax1.set_title('Spectra at Reference Depth %s and at Target Depth %s'
                  %(zrcv_ref[0],zrcv_two[0]))
#    ax1.text(.5,.95,"Q = %s"%(Q),
#                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(fmin, fmax) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    ax2 = plt.subplot(122)
    
    ax2.plot(freq,specrat, c='blue', label = "Spectral Ratio") #using number of samples and sample rate to get x axis
    ax2.plot(freq[srfmax_ind],plot_slope, color='red', linestyle='--', label='fitted slope') #using number of samples and sample rate to get x axis
#    ax2.plot(freq[srfmax_ind],trendpoly(freq[srfmax_ind]),color='red', linestyle='--', label='fitted slope')    
    ax2.set_title('Spectral Ratio from Reference Depth %s to %s Depth'
                  %(zrcv_ref[0],zrcv_two[0]))
    ax2.text(.1,.05,"Slope = %.4f, Q calculated = %.0f"%(slopew[0],Qeff[0]),
                  horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes)
    ax2.legend(loc='upper right',borderaxespad=0, fontsize = 8)        

    ax2.set_xlabel('Frequency hz')
    ax2.set_xlim(fmin,fmax) # extents must be set   
    ax2.set_ylabel('Log Spectral Ratio')    
    ax2.xaxis.grid()    
    ax2.yaxis.grid()

    plt.show()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\Q _%s_estimation.png' 
        %(Qeff), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)   
        
