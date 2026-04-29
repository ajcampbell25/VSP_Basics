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
    
    Tmax = kwargs["time_max"]
    Tmin = kwargs["time_min"]
    skiplabel = kwargs['skiplabel']
    norm = kwargs['norm']
    scal = kwargs['scal']
    title_top = kwargs['title_top']
    info_wig = kwargs['info_wig']

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

        amax = np.nanmax(np.abs(VSPdata), axis=1) 
        data2 = (VSPdata / amax[:, np.newaxis])        
        datascaled = data2 * scal        
        
    else:        
        datascaled = data1 * scal
 
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
        
    for n, label in enumerate(axq.xaxis.get_ticklabels()):
        label.set_rotation(90)
        if n % skiplabel != 0:
            label.set_visible(False)
    axq.set_ylabel('Time (ms)')            
    axq.yaxis.set_label_position("right")
            
    axq.yaxis.grid()

    if(info_wig=='y')or(info_wig=='Y'):    
        print("\u0332".join('\nWiggle Plot Global Information :')) 
        print (' VSPdata.shape :', VSPdata.shape,' (traces,samples)')
        print(' VSPdata type :', VSPdata.dtype)
        print (' Max an Min Amplitude VSPdata :',np.nanmax(VSPdata),np.nanmin(VSPdata))
        print (' datascaled.shape ',datascaled.shape,' (traces,samples)')    
        print (' thead shape :', thead.shape,' (traces,header columns)')
        
def q_datafran(VSP,thead,theadspec,aspec, freq,fs,**kwargs):
    ''' Run a Frequency Analysis on two selected traces.
    Estimate the spectral ratios in a frequency band. From the
    slope of the spectral ratios, calculate Q
    
    VSP : for wiggle trace plot - all traces in data set
    thead : headers corresponding to wiggle trace data
    theadspec : headers corresponding to spectral traces
    aspec : spectral traces F-Z
    tr_ref : reference or shallow trace
    tr_targ : target or deep trace
    specave : average n spectra around the reference and target
    srfreq_min : min frequency for spectral ratio calc
    srfreq_max : max frequency for spectral ratio calc

    fs - sample rate in hertz
    scale - scales amplitude of trace plot, >1 makes plot hotter
        
    The frequency spectrum is magnitude 

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib  import gridspec

    import numpy as np
    import procvsp.qtools as qvsp

    fmin = kwargs["freq_min"]
    fmax= kwargs["freq_max"]
    srfmin = kwargs["srfreq_min"]
    srfmax= kwargs["srfreq_max"]    
    trace_1=kwargs["tr_ref"]-1 # reference trace
    trace_2=kwargs["tr_targ"]-1 # target trace
    save=kwargs["savepng"]

    print("\u0332".join('\nQ_datafran Parameters :'))    
    print ('fs :', fs,)

    # extract analysis traces
    thead_ref = theadspec[trace_1,:]
    thead_two = theadspec[trace_2,:]

    # useful headers
    zrcv_ref = thead_ref[2,]
    trnum_ref = thead_ref[0,]
    
    zrcv_two = thead_two[2,]
    trnum_two = thead_two[0,]
    print (' zrcv_ref :',zrcv_ref,' zrcv_two :',zrcv_two )
    ############create reference and target spectra #############################

    # apply an averaging to spectra (specave=1 for no averaging)
    Xref_mag,Xtwo_mag=qvsp.specave(aspec,**kwargs)
    # find indices on the frequency vector for the reference frequency limits
    srfmax_ind=np.where((freq>=srfmin)&(freq<=srfmax))

    #############  Calculate spectral ratios, Q, and slope of ratios ###########
    
    Qeff,specrat,slopew,plot_slope=qvsp.specratio(Xref_mag,Xtwo_mag,freq,thead_ref,thead_two,**kwargs)

    ############   make spectra and wiggle plots   ################
    
    Xmax=np.max(Xref_mag) 
    Xmin=np.min(Xref_mag)
    
    ############ define some font sizes
    #fparams = {'legend.fontsize': 6,
    #     'axes.labelsize': 16,
    #     'axes.titlesize':10,
    #     'xtick.labelsize':6,
    #    'ytick.labelsize':6}
    #plt.rcParams.update(fparams)
    ftitle=10
    flabel = 7
    
    fig = plt.figure(figsize=(13,6)) # usually needs adjusting    
    gs = gridspec.GridSpec(1, 4,width_ratios=[.3,1,.3,.5], hspace = .35)

    ax1 = plt.subplot(gs[0])       
    ax2 = plt.subplot(gs[1])
    ax4 = plt.subplot(gs[3])
    ax3 = plt.subplot(gs[2])    

    ax1.plot(Xref_mag,freq,  c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis
    ax1.fill_betweenx(freq, Xref_mag, 0, where=(Xref_mag > 0), color='red')
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Spectra at Reference Trace %s'
                  %(trnum_ref),fontsize = ftitle)
#    ax1.text(.5,.95,"Q = %s"%(Q),
#                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_ylabel('Frequency hz',fontsize = flabel)    
    ax1.set_ylim(fmax, fmin) # extents must be set
    ax1.set_xlim(Xmin, Xmax) # extents must be set   

    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_xlabel('Amplitude',fontsize = flabel)    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()
    for label1 in ax1.xaxis.get_ticklabels(): #rotate the x axis labels by 45 deg
        label1.set_rotation(90)

    q_wiggles(ax2,thead,VSP,fs,**kwargs)

    ax3.plot(Xref_mag,freq,  c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis    
    ax3.plot(Xtwo_mag, freq, c = 'blue', label = "Target Depth Spectrum")  # using fftfreq to get x axis
    ax3.fill_betweenx(freq, Xtwo_mag,0, where=(Xtwo_mag > 0), color='blue')
    ax3.fill_betweenx(freq, Xref_mag, Xtwo_mag, where=(Xref_mag > Xtwo_mag), color='red')
    ax3.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax3.set_title('Reference Spectra and \nSpectra at Target Trace %s'
                  %(trnum_two),fontsize = ftitle)

    #ax3.set_ylabel('Frequency hz')
    ax3.yaxis.tick_right()    
    ax3.set_ylim(fmax, fmin) # extents must be set
    ax3.set_xlim(Xmin, Xmax) # extents must be set   

    ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax3.set_xlabel('Amplitude',fontsize = flabel)    
    ax3.xaxis.grid()    
    ax3.yaxis.grid()
    for label2 in ax3.xaxis.get_ticklabels(): #rotate the x axis labels by 45 deg
        label2.set_rotation(90)

    ax4.plot(specrat,freq, c='blue', label = "Spectral Ratio") #using number of samples and sample rate to get x axis
    ax4.plot(plot_slope,freq[srfmax_ind], color='red', linestyle='--', label='fitted slope') #using number of samples and sample rate to get x axis
    ax4.set_title('Spectral Ratios \nBetween Trace %s to Trace %s'
                  %(trnum_ref,trnum_two),fontsize = ftitle)
    ax4.text(.5,.55,"Slope = %.4f, \nQ calculated = %.0f"%(slopew[0],Qeff),
                  horizontalalignment='left', verticalalignment='center', transform=ax4.transAxes)
    ax4.legend(loc='upper right',borderaxespad=0, fontsize = 8)        
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.set_ylabel('Frequency hz',fontsize = flabel)
    ax4.set_ylim(fmax,fmin) # extents must be set   
    ax4.set_xlabel('Log Spectral Ratio',fontsize = flabel)    
    ax4.xaxis.grid()    
    ax4.yaxis.grid()


    DPI=200

    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\Q_wiggle_spec.png',
        dpi=DPI, bbox_inches = 'tight', pad_inches = .1)
            
    plt.show()    

def specfz_qvals(thead,aspec,fdomref,fdom,Q,trnum,freq,fs,Ptitle, **kwargs):

    ''' Plot the dominant frequencies for reference and target traces on top of a 
    2D spectrum image.
    
    Plot the estimated Q values in a separate track
        
        aspec - corresponding 2d amplitude spectrum
        thead - aspec's corresponding header array 
        fmin,fmac - plot frequency range
        fdomref - dominant reference frequencies vector
        fdom - dominant target traces frequencies vector
        Q - calculate Q values vector
        trnum - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scalar applied to amplitude of spectral plot
    '''
    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    fmin = kwargs["freq_min"]
    fmax = kwargs["freq_max"]    
    scale= kwargs['scaler'] # to 'heat up' plots of spectra
    meth = kwargs["method"]
    
    TTobs = thead[:,8]
    rcvdepth = thead[:,2]
    trnum_all = thead[:,0]
 
    print("\u0332".join('\nQ centroid overlay Parameters :'))    
    print (' aspec shape :',aspec.shape, 'thead shape :',thead.shape, 
           'TTobs shape :', TTobs.shape)
    print (' fdom.shape :',fdom.shape)
    print (' Q.shape :',Q.shape)
    print (' trnum_all.shape :',trnum_all.shape,' trnum.shape :',trnum.shape)   
    print (' fs :', fs,)
    print (' scale :', scale)
    print (' trnum_all :',trnum_all)
    print (' meth :',meth)
    
    ############   make amplitude-magnitude plots   ################
    
    fig = plt.figure(figsize=(12,10))
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[.2,1],wspace = .25, hspace=.05)
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    #ax1.plot(trnum,Q, c='red',linewidth = .5, 
    #         label = 'Q', drawstyle = 'steps-pre')
    ax1.scatter(trnum,Q, c='red',s=50,marker='x', 
             label = 'Q')
    #ax1.scatter(trnum,Qint, c='black',s=50,marker='x', 
    #         label = 'Q')
    ax1.set_xlim(np.min(trnum_all), np.max(trnum_all) )
    ax1.set_ylim(0, 500)        
    ax1.grid(which='both', axis='y')
    ax1.set_title(Ptitle,fontsize=14,y=1.3) 
    ax1.set_ylabel('Quality Factor (Q)')
    ax1.set_xlabel('Trace Number')    
    ax1.xaxis.set_label_position('top')    
    ax1.xaxis.set_ticks_position('top') 
    
    plot2 =  ax2.imshow(np.abs(aspec), cmap="gist_ncar_r", interpolation='none',aspect = 'auto',
               vmin = np.min(np.abs(aspec.T))/scale,
               vmax = np.max(np.abs(aspec.T))/scale,
               extent = [trnum_all.min(), trnum_all.max(), freq.max(), freq.min()],zorder=0)
    if meth=="cent":
        ax2.scatter(trnum,fdom,s=50,marker='$\\bigoplus$', color='yellow',linestyle='-', zorder=2)
        ax2.plot(trnum,fdomref, linestyle='-',marker="+", markersize=10,color='white', zorder=1)    
    ax2.set_ylim(fmax, fmin) # extents must be set    
    ax2.set_xlim(np.min(trnum_all), np.max(trnum_all) )    
    ax2.yaxis.grid()    
    ax2.set_xlabel('Trace Number')    
    ax2.set_ylabel('Frequency (hz)')    
    #ax2.set_title(Ptitle)    
    textstr =  'Spectra max : %s'%(aspec.max())    
    ax2.text(0.05, 0.05, textstr, transform=ax2.transAxes, fontsize=12,
             verticalalignment='top')
    
    #plot a colorbar    
    pad = 0.03    
    width = 0.02    
    pos = ax2.get_position()
    axcol = fig.add_axes([pos.xmax + pad, pos.ymin, width, \
                          0.9*(pos.ymax-pos.ymin) ])
    cb2 = fig.colorbar(plot2, label = 'Amplitude', cax = axcol, aspect = 40,format='%.0e')
           
    plt.show()
    
def gauss_demo(aspec,freq,fc1,fc2,fs,ptitle, **kwargs):

    ''' Plot 2 spectra and overlay them with their centroid frequencies
    We want the spectra to look gaussian for the centroid method to work
    1. Reference trace spectra
    2. Reference trace plus trace window spectrum
    
    aspec : 2D F-Z spectra
    freq : frequency vector for aspec
    fc1: centroid frequencies for reference traces
    fc2: centroid frequencies for target traces

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib.patches import Rectangle
    
    import numpy as np

    #tmin = kwargs["time_min"]
    #tmax = kwargs["time_max"]
    fmin = kwargs["freq_min"]
    fmax = kwargs["freq_max"]
    qcfmin = kwargs["qcfreq_min"]
    qcfmax = kwargs["qcfreq_max"]
    trace_1 = kwargs["tr_ref"] # reference trace
    trace_2 = kwargs["tr_win"] # target trace
    save = kwargs["savepng"]
    
    print("\u0332".join('\nGauss Demo Parameters :'))    
    print ('fs :', fs,)
    print (' aspec.shape :',aspec.shape)
    print (' trace_1 :',trace_1,' trace_2 :',trace_2)

    ############create reference and target spectra #############################

    #X_ref,X_two,freq=two_spectra(data_ref,data_two,fs,**kwargs)
    
    X_ref=aspec[:,trace_1]
    X_two=aspec[:,trace_2]
    # get dominant frequencies for the 2 traces
    
    fc_ref=fc1[trace_1]
    fc_two=fc2[trace_2]
    
    # get maximum amplitude in the 2 spectrum arrays
    spmax=np.max((X_ref,X_two)) #for plotting a freq range rectangle
    
    ############  make spectra plots   ################
    
    fig=plt.figure(figsize=(12,5))    
    ax1 = plt.subplot(111)    

    ax1.plot(freq, X_ref, c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis    
    ax1.plot(freq, X_two, c = 'blue', label = "Target Depth Spectrum")  # using fftfreq to get x axis 
    ax1.scatter(fc_ref,0,c = 'red',label = "Reference Centroid")
    ax1.axvline(x=fc_ref, color='red',ls='--')
    ax1.scatter(fc_two,0,c = 'blue',label = "Target Centroid")
    ax1.axvline(x=fc_two, color='blue',ls='--')

    ax1 = plt.gca()
    cbox=ax1.add_patch(Rectangle((qcfmin, 0), qcfmax-qcfmin, spmax,
                      alpha=.4, facecolor='pink',label='Computation Limits'))    
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        

    ax1.set_title('Spectra at Reference Trace %s and at Target Trace %s'
                  %(trace_1,trace_2))
#    ax1.text(.5,.95,"Q = %s"%(Q),
#                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(fmin, fmax) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    plt.show()
    
def average_demo(aspec,freq,fs, **kwargs):

    ''' Plot 2 spectra and overlay them with their centroid frequencies
    1. Reference trace spectra
    2. Reference trace plus trace window spectrum
    
    
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n'
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter
        
        The frequency spectrum is 

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib.patches import Rectangle
    
    import numpy as np
    import procvsp.qtools as qvsp

    fmin = kwargs["freq_min"]
    fmax = kwargs["freq_max"]
    trace_1 = kwargs["tr_ref"]-1 # reference trace
    trace_2 = kwargs["tr_targ"]-1 # target trace
    avenum=kwargs["specave"]
    
    print("\u0332".join('\Spec averaging Demo Parameters :'))    
    print ('fs :', fs,)
    print (' aspec.shape :',aspec.shape)
    print (' trace_1 :',trace_1,' trace_2 :',trace_2)

    ############average the reference and target spectra #############################

    # calculate the averaged spectral traces
    in_ref,in_two,ave_ref,ave_two=qvsp.specave(aspec,**kwargs)
    print (' X_ref.shape :',ave_ref.shape,' X_two.shape :',ave_two.shape)

    # get maximum amplitude in the 2 spectrum arrays
    spmax=np.max((ave_ref,ave_two)) #for plotting a freq range rectangle
    
    ############  make spectra plots   ################
    
    fig=plt.figure(figsize=(12,8))    
    ax1 = plt.subplot(211)    

    ax1.plot(freq, in_ref, c = 'blue')#, label = "Input_Reference Spectrum")  # using fftfreq to get x axis 
    ax1.plot(freq, ave_ref, c = 'red', label = "Mean Reference Spectrum")  # using fftfreq to get x axis    
    #ax1.plot(freq, X_two, c = 'blue', label = "Target Spectrum")  # using fftfreq to get x axis 
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Spectra at Reference Trace %s,  %s level averaging'
                  %(trace_1,avenum))
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(fmin, fmax) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')    
    ax1.xaxis.grid()    

    
    ax2 = plt.subplot(212)
    ax2.plot(freq, in_two, c = 'blue')#, label = "Input_Target Spectrum")  # using fftfreq to get x axis 
    ax2.plot(freq, ave_two, c = 'red', label = "MeanTarget Spectrum")  # using fftfreq to get x axis    
    ax2.legend(loc='best',borderaxespad=0, fontsize = 8)
    ax2.set_title('Spectra at Target Trace %s, %s level averaging'
                  %(trace_2,avenum))
    ax2.set_xlabel('Frequency hz')    
    ax2.set_xlim(fmin, fmax) # extents must be set   
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax2.set_ylabel('Amplitude')     
    ax2.xaxis.grid()

    plt.tight_layout()
    plt.show()

def plotqint(qave,q_int,theadq,**kwargs):
    '''plot the qave and qint versus depth
    inputs:
    qave : q from a single reference trace, cumlaative or average q
    q_int : calculated interval q from qave, not sliding window q
    '''
    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    print("\u0332".join('\nInterval Q test Parameters :'))    
    print (' qave.shape :',qave.shape,' qint.shape :',q_int.shape)
    print (' theadq.shape :',theadq.shape)
    
    qmin = kwargs["q_min"]
    qmax = kwargs["q_max"]    
    Ptitle = kwargs["plot_title"]
    
    TTobs = theadq[:,8]
    rcvdepth = theadq[:,2]
    trnum = theadq[:,0]
 
    ############   make amplitude-magnitude plots   ################
    
    fig = plt.figure(figsize=(8,3))
    
    gs = gridspec.GridSpec(1, 1)#,wspace = .15, hspace=.05)
    
    ax1 = plt.subplot(gs[0])
    
    ax1.scatter(trnum,qave, c='red',s=50,marker='x', 
             label = 'Q ave')
    #ax1.scatter(trnum,qint, c='blue',s=50,marker='x', 
    #         label = 'Q int sliding win.')
    ax1.scatter(trnum,q_int, c='black',s=50,marker='o', 
             label = 'Q int calc')
    ax1.set_xlim(np.min(trnum), np.max(trnum) )    
    ax1.set_ylim(qmin, qmax)    
    ax1.grid(which='both', axis='y')
    ax1.set_title(Ptitle,fontsize=14,y=1.3) 
    ax1.set_ylabel('Quality Factor (Q)')
    ax1.set_xlabel('Trace Number')    
    ax1.xaxis.set_label_position('top')    
    ax1.xaxis.set_ticks_position('top')
    
    plt.show()