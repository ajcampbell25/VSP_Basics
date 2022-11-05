def q_wiggles(axq,thead, VSPdata, fs, **kwargs):
    
    """Make a wiggle plot of seismic traces.
    
    Crossplot x (amplitude) and y (time). Add amplitude to receiver depth to 
    get trace deflection. Alternatively add amplitude to receiver number to get 
    trace deflection. Scaling in X direction (amplitude) is different in each 
    case
    
    Trace deflection is based on sample value. Plots are spaced by receiver 
    number or trace number. 
    
    A scalar may need to be applied to make reasonable deflections, dependent 
    on data amplitudes and plot spacing
    
    Plot parameter definitions:
    
    pol = polarity 'n' for normal or tape polarity, 'r' to flip polarity
    Tmax, Tmin = start and end time of plot    
    first_rcv = first receiver in plot  
    spacing =  'Z' for traces spread by receiver depth
    skiplabel =  plot every nth header
    norm = plot trace normalization 'n' or 'y'         
    plot_polarity = 'n'     # n for normal or tape polarity, r to flip polarity 
    scal = scale plot amplitudes by this value
    info_wig = print diagnostic information to terminal
    Title_plot = 'plot title '

    """
      
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import matplotlib.cm as cm
    
#    pol = kwargs['pol']
    Tmax = kwargs["time_max"]
    Tmin = kwargs["time_min"]
#    first_rcv = kwargs['first_rcv']
#    spacing = kwargs['spacing']    
    skiplabel = kwargs['skiplabel']
    norm = kwargs['norm']
    scal = kwargs['scal']
    title_top = kwargs['title_top']
    info_wig = kwargs['info_wig']
#    tframe = kwargs['timframe']

    # trace header info for main (decon up) plot          
    TT = thead[:,8]
    rcv_depth = thead[:,2]
    trace_num = thead[:,0]

    # get the y axis as a time array    
    numsamp = VSPdata.shape[1]
    y = np.arange(0,numsamp*(1000/fs),(1000/fs)  )

    #create an empty array to put normalized data into     
    data2 = np.zeros(shape = (VSPdata.shape[0], VSPdata.shape[1]), dtype = np.float32)    
    data1 = VSPdata
    
    # apply a trace normalization to main plot also a scale factor         
    if (norm == 'Y') or (norm =='y'):        
         #row_sums = np.linalg.norm(VSPdata, axis=1)
        #data2 = (VSPdata / row_sums[:, np.newaxis]) # problem for traces of all zeros,ie. after median and subtraction
        #datascaled = data2
        amax = np.nanmax(np.abs(VSPdata), axis=1) 
        data2 = (VSPdata / amax[:, np.newaxis])        
        datascaled = data2 * scal        
        
    else:        
        datascaled = data1 * scal

    ## Flip polarity if requested        
#    if (pol == 'r') or (pol =='R'):        
#        datascaled = datascaled * -1
        
    ##### Set up the baselines for each trace #####
    ##### Either receiver depth or trace number ###
#    if (spacing == 'Z') or (spacing == 'z'):        
#        dscaler, pad = (rcv_depth, 10)        
#        dlabel = 'Receiver Depth'
        
#    else:    
    dscaler, pad = (trace_num, 1)        
    dlabel = 'Receiver Number'
    
    # for labeling trace number on top of main track
    dscaler_tracenum, pad = (trace_num, 1)        
    dlabel_tracenum = 'Receiver Number'
   
    for i, trace in enumerate(datascaled[::1, :]):

        #add sample values to either receiver number or trace number     
        x = trace + dscaler[i]    
#        ax2.plot(x, y, 'k-',  linewidth = .5)
        axq.plot(x, y, 'k-',  linewidth = .5)
        axq.fill_betweenx(y, dscaler[i], x, where=(x > dscaler[i]), color='k')
        axq.set_xlim(dscaler[0]-pad, dscaler[-1]+pad )
        axq.set_ylim(Tmax, Tmin)
        axq.set_xticks(dscaler[:-1:1])        
        axq.set_xlabel(dlabel)
        
    axq.plot(dscaler,TT[0:],c='red',linewidth=2, label='Travel Time' )    
    for n, label in enumerate(axq.xaxis.get_ticklabels()):
        label.set_rotation(90)
        if n % skiplabel != 0:
            label.set_visible(False)            
        
    axq.yaxis.grid()

 
    if(info_wig=='y')or(info_wig=='Y'):    
        print("\u0332".join('\nWiggle Plot Global Information :')) 
        print (' VSPdata.shape :', VSPdata.shape,' (traces,samples)')
        print(' VSPdata type :', VSPdata.dtype)
        print (' Max an Min Amplitude VSPdata :',np.nanmax(VSPdata),np.nanmin(VSPdata))
        print (' datascaled.shape ',datascaled.shape,' (traces,samples)')    
        print (' thead shape :', thead.shape,' (traces,header columns)')
        
def q_datafran(VSP, thead,fs,**kwargs):
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
    from matplotlib  import gridspec

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
    
    # get the magnitude of the complex spectral signal 
    Xref_mag=np.abs(X_ref)
    Xtwo_mag=np.abs(X_two)
    # max and min for plotting relative spectra
    Xmax=np.max(Xref_mag) 
    Xmin=np.min(Xref_mag)

    ####### Only keep positive frequencies #########
    
    keep = freq>=0    
    Xref_mag = Xref_mag[keep]       
    Xtwo_mag = Xtwo_mag[keep]    
    
    freq = freq[keep]
    
    ############ get time and depth between receivers  #####################

    Z0 = zrcv_ref
    Z = zrcv_two

    T=TTobs_two/1000
    T0=TTobs_ref/1000    

    ########### use spectral ratios to get Q from Q filtered spectrum

    specrat=np.log(Xtwo_mag/Xref_mag)  
    # get sample number for max frequency to control plots
    fmax_ind=np.where((freq>=fmin)&(freq<=fmax))
    srfmax_ind=np.where((freq>=srfmin)&(freq<=srfmax))

    slopew = np.polyfit(freq[srfmax_ind],specrat[srfmax_ind],  1)# 1 means first order - linear regression
    plot_slope = slopew[0]*freq[srfmax_ind]+slopew[1] 

    Qeff=-1*(np.pi*((T-T0)/slopew[0]))
    
    # trendpoly = np.poly1d(slopew) # an alternative for plotting instead of y=mx+b 
    
    print (' slopew :', slopew, ' slopew.shape :',slopew.shape, 'Qeff :',Qeff)
    ############   make spectra and wiggle plots   ################
    
    fig = plt.figure(figsize=(18,8)) # usually needs adjusting    
    gs = gridspec.GridSpec(1, 4,width_ratios=[.3,1,.3,.5], hspace = .05)

    ax1 = plt.subplot(gs[0])    

    ax1.plot(Xref_mag,freq,  c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis
    ax1.fill_betweenx(freq, Xref_mag, 0, where=(Xref_mag > 0), color='red')
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Spectra at Reference Trace %s'
                  %(trnum_ref[0]))
#    ax1.text(.5,.95,"Q = %s"%(Q),
#                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_ylabel('Frequency hz')    
    ax1.set_ylim(fmax, fmin) # extents must be set
    ax1.set_xlim(Xmin, Xmax) # extents must be set   

    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_xlabel('Amplitude')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()
    for label in ax1.xaxis.get_ticklabels(): #rotate the x axis labels by 45 deg
	    label.set_rotation(90)


    ax2 = plt.subplot(gs[1])    
    q_wiggles(ax2,thead,VSP,fs,**kwargs)

    ax3 = plt.subplot(gs[2])
    ax3.plot(Xref_mag,freq,  c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis    
    ax3.plot(Xtwo_mag, freq, c = 'blue', label = "Target Depth Spectrum")  # using fftfreq to get x axis
    ax3.fill_betweenx(freq, Xtwo_mag,0, where=(Xtwo_mag > 0), color='blue')
    ax3.fill_betweenx(freq, Xref_mag, Xtwo_mag, where=(Xref_mag > Xtwo_mag), color='red')
    ax3.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax3.set_title('Reference Spectra and \nSpectra at Target Trace %s'
                  %(trnum_two[0]))
    ax3.set_ylabel('Frequency hz')    
    ax3.set_ylim(fmax, fmin) # extents must be set
    ax3.set_xlim(Xmin, Xmax) # extents must be set   

    ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax3.set_xlabel('Amplitude')    
    ax3.xaxis.grid()    
    ax3.yaxis.grid()
    for label in ax3.xaxis.get_ticklabels(): #rotate the x axis labels by 45 deg
	    label.set_rotation(90)

    ax4 = plt.subplot(gs[3])

 
    ax4.plot(specrat,freq, c='blue', label = "Spectral Ratio") #using number of samples and sample rate to get x axis
    ax4.plot(plot_slope,freq[srfmax_ind], color='red', linestyle='--', label='fitted slope') #using number of samples and sample rate to get x axis
    ax4.set_title('Spectral Ratios \nBetween Trace %s to Trace %s'
                  %(trnum_ref[0],trnum_two[0]))
    ax4.text(.5,.55,"Slope = %.4f, \nQ calculated = %.0f"%(slopew[0],Qeff[0]),
                  horizontalalignment='left', verticalalignment='center', transform=ax4.transAxes)
    ax4.legend(loc='upper right',borderaxespad=0, fontsize = 8)        

    ax4.set_ylabel('Frequency hz')
    ax4.set_ylim(fmax,fmin) # extents must be set   
    ax4.set_xlabel('Log Spectral Ratio')    
    ax4.xaxis.grid()    
    ax4.yaxis.grid()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\Q_wiggle_spec.png',
        dpi=DPI, bbox_inches = 'tight', pad_inches = .1)
            

    plt.show()    
