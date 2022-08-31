
def wiggle_plot(thead, VSPdata, **kwargs):
    
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
    
    pol = kwargs['pol']
    Tmax = kwargs['Tmax']
    Tmin = kwargs['Tmin']
    first_rcv = kwargs['first_rcv']
    spacing = kwargs['spacing']    
    skiplabel = kwargs['skiplabel']
    fs = kwargs['fs']
    norm = kwargs['norm']
    scal = kwargs['scal']
    title_top = kwargs['title_top']
    info_wig = kwargs['info_wig']
    tframe = kwargs['timframe']

    # trace header info for main (decon up) plot          
    TVDSRD = thead[:,9]
    TT = thead[:,8]
    rcv_depth = thead[:,2]
    trace_num = thead[:,0]
    intVel = thead[:, -1]    

    # change time header to relate to data orientation
    if tframe == "flat":
        TT = TT*0
    if tframe == "twt":
        TT = thead[:,-2]

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
    if (pol == 'r') or (pol =='R'):        
        datascaled = datascaled * -1
        
    ##### Set up the baselines for each trace #####
    ##### Either receiver depth or trace number ###
    if (spacing == 'Z') or (spacing == 'z'):        
        dscaler, pad = (rcv_depth, 10)        
        dlabel = 'Receiver Depth'
        
    else:    
        dscaler, pad = (trace_num, 1)        
        dlabel = 'Receiver Number'
    
    # for labeling trace number on top of main track
    dscaler_tracenum, pad = (trace_num, 1)        
    dlabel_tracenum = 'Receiver Number'                  

    fig = plt.figure(figsize=(15,12))    
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.2, 2], hspace = .05)
    
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])

    ax1.plot(TVDSRD, intVel, c='red',linewidth = .5, 
             label = 'Interval Velocity', drawstyle = 'steps-pre')
    ax1.set_xlim(np.min(TVDSRD)-pad, np.max(TVDSRD) + pad )    
    ax1.set_title('Interval Velocity and %s'%(title_top),fontsize=14)            
    for i, trace in enumerate(datascaled[::1, :]):
        #add sample values to either receiver number or trace number     
        x = trace + dscaler[i]    
        ax2.plot(x, y, 'k-', linewidth = .5)
        ax2.fill_betweenx(y, dscaler[i], x, where=(x > dscaler[i]), color='k')
        ax2.set_xlim(dscaler[0]-pad, dscaler[-1]+pad )
        ax2.set_ylim(Tmax, Tmin)
        ax2.set_xticks(dscaler[:-1:1])        
        ax2.set_xlabel(dlabel)
        
    ax2.plot(dscaler,TT[0:],c='red',linewidth=2, label='Travel Time' )    
    for n, label in enumerate(ax2.xaxis.get_ticklabels()):
        label.set_rotation(90)
        if n % skiplabel != 0:
            label.set_visible(False)            
        
    ax2.yaxis.grid()
    
    plt.show()

    if(info_wig=='y')or(info_wig=='Y'):    
        print("\u0332".join('\nWiggle Plot Global Information :')) 
        print (' Number of traces in plot :', VSPdata.shape[0], 
           ' Number of samples per trace :', VSPdata.shape[1])
        print(' VSPdata type :', VSPdata.dtype)
        print (' sample rate (hz)',fs, 
           ' min max plot time',np.min(y), np.max(y))    
        print (' thead shape :', thead.shape)
        print (' Min TVDSRD - pad', np.min(TVDSRD)-pad, ' Pad :', pad)    
        print (' Max TVDSRD + pad', np.max(TVDSRD)+pad, ' Pad :', pad)    
        print (' min max intvel :', np.min(intVel), np.max(intVel))

def composite_plot(thead, VSPdata, Cstack, **kwargs):
    
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
    
    pol = kwargs['pol']
    Tmax = kwargs['Tmax']
    Tmin = kwargs['Tmin']
    first_rcv = kwargs['first_rcv']
    spacing = kwargs['spacing']    
    skiplabel = kwargs['skiplabel']
    fs = kwargs['fs']
    norm = kwargs['norm']
    scal = kwargs['scal']
    title_top = kwargs['title_top']
    info_wig = kwargs['info_wig']
    tframe = kwargs['timframe']
    
    # trace header info for main (decon up) plot   
    TVDSRD = thead[:,9]
    TT = thead[:,8]
    rcv_depth = thead[:,2]
    trace_num = thead[:,0]
    intVel = thead[:, -1]
    # get trace number array for corridor stack
    trace_num_cstk = np.arange(0, Cstack.shape[0])    

    # change time header to relate to data orientation
    if tframe == "flat":
        TT = TT*0
    if tframe == "twt":
        TT = thead[:,-2]
    
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
       
    # flip polarity if requested
    if (pol == 'r') or (pol =='R'):        
        datascaled = datascaled * -1
        datascaled_cstk = datascaled_cstk*-1

    ##### Set up the baselines for each trace #####
    ##### Either receiver depth or trace number ###
    if (spacing == 'Z') or (spacing == 'z'):        
        dscaler, pad = (rcv_depth, 10)        
        dlabel = 'Receiver Depth'
        
    else:    
        dscaler, pad = (trace_num, 1)        
        dlabel = 'Receiver Number'

    # Baselines for corridor stack will always be trace number
    dscaler_cstk, padc = (trace_num_cstk, 5)       
    dlabel_cstk = 'Trace Number'

    # scale the corridor stack to somewhat similar amplitude range as input to cstack
    # ie. each trace amplitude is added to trace number    
    datascaled_cstk = Cstack * abs(scal/((trace_num.max()+pad)/(trace_num_cstk.max()-padc*2)))
    print (' pol :',pol, ' comp scalar :', scal/((trace_num.max()+pad)/(trace_num_cstk.max()-padc*2)) )
    
    fig = plt.figure(figsize=(17,12))

    gs1 = gridspec.GridSpec(10, 10, hspace = .3)#, wspace=0.01) # make a row by col grid
    ax1 = fig.add_subplot(gs1[0:1, 0:9])            # combine rows or columns
    ax2 = fig.add_subplot(gs1[1:, 0:9])
    ax3 = fig.add_subplot(gs1[1:, 9:10])
     
    ax1.plot(TVDSRD, intVel, c='red',linewidth = .5, 
             label = 'Interval Velocity', drawstyle = 'steps-pre')
    ax1.set_xlim(np.min(TVDSRD)-pad, np.max(TVDSRD) + pad )    
    ax1.set_title('Interval Velocity and %s'%(title_top),fontsize=14) 
    for i, trace in enumerate(datascaled[::1, :]):
        #add sample values to either receiver number or trace number     
        x = trace + dscaler[i]    
        ax2.plot(x, y, 'k-', linewidth = .5)
        ax2.fill_betweenx(y, dscaler[i], x, where=(x > dscaler[i]), color='k')
        ax2.set_xlim(dscaler[0]-pad, dscaler[-1]+pad )
        ax2.set_ylim(Tmax, Tmin)
        ax2.set_xticks(dscaler[:-1:1])        
        ax2.set_xlabel(dlabel)   
    ax2.plot(dscaler,TT[0:],c='red',linewidth=2, label='Travel Time' )    
    for n, label in enumerate(ax2.xaxis.get_ticklabels()):
        label.set_rotation(90)
        if n % skiplabel != 0:
            label.set_visible(False)
    ax2.yaxis.grid()

    for i, tracecstk in enumerate(datascaled_cstk[::1, :]):
        x = tracecstk + dscaler_cstk[i]    
        ax3.plot(x, y, 'k-', linewidth = .5)
        ax3.fill_betweenx(y, dscaler_cstk[i], x, where=(x > dscaler_cstk[i]), color='k')
        ax3.set_xlim(dscaler_cstk[0]-padc, dscaler_cstk[-1]+padc )
        ax3.set_ylim(Tmax, Tmin)
        ax3.set_xticks(dscaler_cstk[:-1:1])        
        ax3.set_xlabel(dlabel_cstk)               
        
    ax3.yaxis.grid()
    
    plt.show()

    if(info_wig=='y')or(info_wig=='Y'):    
        print("\u0332".join('\nWiggle Plot Global Information :')) 
        print (' Number of traces in plot :', VSPdata.shape[0], 
           ' Number of samples per trace :', VSPdata.shape[1])
        print(' VSPdata type :', VSPdata.dtype)
        print (' datascaled shape [0]',datascaled.shape[0], 
           ' datascaled shape [1]',datascaled.shape[1])    
        print (' thead shape :', thead.shape)
        print (' Min TVDSRD - pad', np.min(TVDSRD)-pad, ' Pad :', pad)    
        print (' Max TVDSRD + pad', np.max(TVDSRD)+pad, ' Pad :', pad)    
        print (' min max intvel :', np.min(intVel), np.max(intVel))

        
def plotsingletrace( VSP1, Tmin, Tmax, thead, spacing, scal1, title):

#    import numpy as np
    import matplotlib.pyplot as plt
    
    rcv_depth = thead[0:1,2]
    trace_num = thead[0:1,0]
    
    data1scaled = VSP1 * scal1
    y = np.arange(0, data1scaled.shape[0] )
        
    print("\u0332".join('\nSingle Trace Plot Global Information :'))    
    print (' VSP1 shape :', VSP1.shape)
    print (' VSP1 type :', VSP1.dtype)
    print (' data1scaled shape :', data1scaled.shape)
     
    if (spacing == 'Z') or (spacing == 'z'):        
        dscaler, pad = (rcv_depth, rcv_depth/10)        
        dlabel = 'Receiver Depth'
        
    else:    
        dscaler, pad = (trace_num, 2)        
        dlabel = 'Receiver Number'
        
    x = data1scaled + dscaler    
    xflat = x.ravel()
    
    print (' data1scaled shape [0] :', data1scaled.shape[0], 
           ' Number of samples per trace [1] :', data1scaled.shape[1])    
    print (' x shape :', x.shape, ' x flat shape :', xflat.shape)
            
    fig = plt.figure(figsize=(15,10))    
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,1], wspace = .25)
    
    ax1 = plt.subplot(gs[0])
    
    ax1.plot(y, xflat, 'k-', linewidth = .5)
    ax1.fill_between(y, dscaler, xflat, where=(xflat> dscaler), color='k')    
    ax1.set_xlim(Tmin, Tmax)
    ax1.set_yticks(dscaler[-1:1:]) #careful with the first and last    
    ax1.set_ylabel(dlabel)    
    ax1.set_title('%s'%(title), fontsize = 14)
   
    for label in ax1.yaxis.get_ticklabels():  
        label.set_rotation(90)
        
    ax1.xaxis.grid()
    
    plt.show()
        
def four_plots(VSP1, VSP2, VSP3, VSP4,fs, thead, **kwargs):    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec

  
    txt1 = kwargs['ss_title1']
    txt2 = kwargs['ss_title2']    
    txt3 = kwargs['ss_title3']
    txt4 = kwargs['ss_title4']
    save = kwargs['savePNG']  
    png_txt = kwargs['png_name']       
    scale = kwargs['scalar']              # scaling to plot 1,2,3,4 amplitudes    
    trange = kwargs['time_range']   
    drange = kwargs['depth_range']
    info_gray4 = kwargs['info_gray4']
    
    rcvdepth = thead[:, 2]
    numsamp = VSP1.shape[1]
    
    cmap = 'seismic'
    
    ################ create plot tick labeling and indexing ##############
    
    tindex1 = np.arange(0, numsamp*(1000/fs),(1000/fs) )  # convert fs to msec.
    tindex2 = np.arange(0, VSP2.shape[1]*(1000/fs),(1000/fs))
    tindex3 = np.arange(0, VSP3.shape[1]*(1000/fs),(1000/fs))  
    tindex4 = np.arange(0, VSP4.shape[1]*(1000/fs),(1000/fs))
    
    rindex  = np.stack([rcvdepth for _ in range(VSP1.shape[0])], axis=1)
    
    ############### make plots of input and SVD filtered output ################
    
    fig = plt.figure(figsize=(16,6) )    
    fig.subplots_adjust(wspace=0.125, hspace=0.5)

    from matplotlib.gridspec import GridSpec
    
    gs1 = GridSpec(1, 25, hspace = .25, wspace=0.01) # make a row by col grid
    ax1 = fig.add_subplot(gs1[0:1, 0:6])            # combine rows or columns
    ax2 = fig.add_subplot(gs1[0:1, 6:12])
    ax3 = fig.add_subplot(gs1[0:1, 12:18])
    ax4 = fig.add_subplot(gs1[0:1, 18:24])

    ax1.imshow(VSP1.T, cmap=cmap, interpolation='none', 
               vmin = -VSP1.max()/scale[0],vmax = VSP1.max()/scale[0],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
                tindex1.min()], aspect = 'auto')

    ax1.yaxis.grid(c = 'black', lw = .1)    
    ax1.set_ylim(trange[1], trange[0]) # extents must be set    
    ax1.set_xlim(drange[0], drange[1])    
    ax1.set_xlabel('Receiver Depth')    
    ax1.set_title('%s'%(txt1))

    ax2.imshow(VSP2.T, cmap=cmap, interpolation='none',
               vmin = -VSP2.max()/scale[1],vmax = VSP2.max()/scale[1],
               extent = [rindex.min(), rindex.max(), tindex2.max(), 
                tindex2.min()], aspect = 'auto')    
    ax2.yaxis.grid(c = 'black', lw = .1)    
    ax2.set_ylim(trange[1], trange[0]) # extents must be set    
    ax2.set_xlim(drange[0], drange[1])    
    ax2.set_xlabel('Receiver Depth')    
    ax2.set_title('%s '%(txt2))    
    ax2.set_yticklabels([])
 
    ax3.imshow(VSP3.T, cmap=cmap, interpolation='none',
               vmin = -VSP3.max()/scale[2],vmax = VSP3.max()/scale[2],
               extent = [rindex.min(), rindex.max(), tindex3.max(), 
               tindex3.min()], aspect = 'auto')
    
    ax3.yaxis.grid(c = 'black', lw = .1)    
    ax3.set_ylim(trange[1], trange[0]) # extents must be set    
    ax3.set_xlim(drange[0], drange[1])    
    ax3.set_xlabel('Receiver Depth')    
    ax3.set_title('%s '%(txt3))    
    ax3.set_yticklabels([])  
 
    ax4.imshow(VSP4.T, cmap=cmap, interpolation='none',
               vmin = -VSP4.max()/scale[3],vmax = VSP4.max()/scale[3],
               extent = [rindex.min(), rindex.max(), tindex4.max(), 
               tindex4.min()], aspect = 'auto')    
    ax4.yaxis.grid(c = 'black', lw = .1)    
    ax4.set_ylim(trange[1], trange[0]) # extents must be set    
    ax4.set_xlim(drange[0], drange[1])    
    ax4.set_xlabel('Receiver Depth')    
    ax4.set_title('%s '%(txt4))
    ax4.yaxis.set_label_position('right')    
    ax4.yaxis.set_ticks_position('right')

    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\procflow_gray_%s.png' 
        %(png_txt), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)
            
    plt.show()
    
    if(info_gray4=='y')or(info_gray4=='Y'):        
        print("\u0332".join('\nFour Box Global Information :'))
        print (' Number of traces in plot :', VSP1.shape[0], 
           ' Number of samples per trace :', ' tindex1 shape :,',tindex1.shape,
           ' tindex1 min max :', tindex1.min(),tindex1.max(),            
           ' tindex2 min max :', tindex2.min(),tindex1.max(), 
           ' tindex3 min max :', tindex3.min(),tindex1.max(), 
           ' tindex4 min max :', tindex4.min(),tindex1.max())
 
def plotcolor(thead, VSPdata,**kwargs):

    """Make an image plot of seismic traces.
    
    Plot samples as an image. 
    
    Amplitude is controlled by the max_amp, min_amp parameters
    
    Plot parameter definitions:
    
    pol = polarity 'n' for normal or tape polarity, 'r' to flip polarity
    Tmax, Tmin = start and end time of plot    
    first_rcv = first receiver in plot  
    spacing =  'Z' for traces spread by receiver depth
    skiplabel =  plot every nth header
    norm = plot trace normalization 'n' or 'y'         
    plot_polarity = 'n'     # n for normal or tape polarity, r to flip polarity 
    samp_decimate(ds) = decimation factor to speed plotting - ie. 0.25ms sample rate 
                        is a lot of data to plot
    min_amp, max_amp = user defined to heat up or cool down plot
                       if either value is set to 'data' these values are global,
                       calculated from the complete data set
    color(colr) = chose a useful color map for data being plotted                  
    Title_plot = 'plot title '

    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import matplotlib.cm as cm    

    # generate 2d arrays for the time index and the receiver depths
    
    pol = kwargs['pol']
    Tmin, Tmax = kwargs['Tmin'],kwargs['Tmax']
    shotime = kwargs['show_time']   
    skiplabel = kwargs['skiplabel']
    ds = kwargs['samp_decimate']
    fs = kwargs['fs']
    norm = kwargs['norm']
    plotmin, plotmax = kwargs["min_amp"], kwargs["max_amp"]
    colr = kwargs['color']
    title_top = kwargs['title_top']

    TT = thead[:,8]
    rcv_depth = thead[:,1] # raw measured depth better for deviated well plots
    trace_num = thead[:,0]
    intVel = thead[:, -1]

    # chose to display or not the arrival time line overlay
    if (shotime=='n')or(shotime=='N'):
        TT = thead[:,8]*0

    # make a time vector and receiver depth vector for displays
    numsamp = VSPdata.shape[1]
    tindex1 = np.arange(0, numsamp*(1000/fs),(1000/fs) )  # convert fs in hertz to milliseconds
    tindex  = np.stack([tindex1 for _ in range(VSPdata.shape[0])], axis=0)
    rindex  = np.stack([rcv_depth for _ in range(VSPdata.shape[1])], axis=1)
        
    print('VSPdata shape:', VSPdata.shape, 'tindex shape :', tindex.shape,'rindex shape :', rindex.shape, )

    if (norm == 'Y') or (norm =='y'):
        #row_sums = np.linalg.norm(VSPdata, axis=1)
        #data2 = (VSPdata / row_sums[:, np.newaxis]) # problem for traces of all zeros,ie. after median and subtraction
        #datascaled = data2
        amax = np.nanmax(np.abs(VSPdata), axis=1) 
        data2 = (VSPdata / amax[:, np.newaxis])        
        datascaled = data2 * scal        
    else:
        datascaled = VSPdata 
        
    if (pol == 'r') or (pol =='R'):        
        datascaled = datascaled * -1
    
    # scale the plot amplitude to max in data or by user supplied values    
    vmin = plotmin
    vmax = plotmax
    
    if (plotmin == "data") or ( plotmax == "data"):
        vmin = np.min(datascaled)
        vmax = np.max(datascaled)

    z1min = np.min(rcv_depth)
    z1max = np.max(rcv_depth)
    
    ####### add a trace number track ##########################
    yaxis = rcv_depth[::skiplabel]*0+.2
    xaxis = rcv_depth[::skiplabel]
    value = trace_num[::skiplabel].astype(int) # may need to leave as float

    fig = plt.figure(figsize=(14,12))    
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.15, 2], hspace = .05)    
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])
    
    for i, txt in enumerate(value):
        ax1.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)

    ax1.set_ylabel('Trace Number' )
    ax1.set_yticklabels([])    
    ax1.set_xticks(rcv_depth[:-1:skiplabel])
    ax1.tick_params(axis = 'x', direction = 'in')
    ax1.set_xticklabels([])
    ax1.set_title('%s'%(title_top),fontsize=14)
    
    ax1.set_xlim(z1min,z1max) # comment out for a plot that fits data limits    

    image2 = ax2.pcolormesh(rindex[::ds,::ds], tindex[::ds,::ds], datascaled[::ds,::ds], \
                           vmin = vmin,
                           vmax = vmax,                       
                           cmap = colr, label= 'VSP', shading = 'gouraud')
    ax2.plot(rcv_depth[0::ds],TT[0::ds],c='red',linewidth=1, label='Arrival Time' )                               
    ax2.set_ylim(Tmax,Tmin) # comment out for a plot that fits data limits    
    ax2.grid()    
    ax2.set_xlabel('Depth (ft)')#,fontsize=18)    
    ax2.set_ylabel('Time (ms)')#,fontsize=18)
    ax2.legend(loc='best')    

    pad = 0.03    
    width = 0.02    
    pos2 = ax2.get_position()    
    axcol2 = fig.add_axes([pos2.xmax + pad, pos2.ymin, width, 0.7*(pos2.ymax-pos2.ymin) ])
    fig.colorbar(image2, label = '%s'%(title_top), cax = axcol2, aspect = 40)

    
    plt.show()
'''    
def plot_param_top(self):    
        # create a fancy frame for plot parameter entry in column 0, row 1 of master window
        self.plot_paramframe = ttk.LabelFrame(self.input_frame, text = 'Trace Display Options', width = 400)
        self.plot_paramframe.grid(row = 1, column = 0, padx = 10, pady=5)#, sticky='ns')
        self.plot_paramframe.grid_columnconfigure(0, minsize = 100)#, weight=1)
        self.plot_paramframe.grid_columnconfigure(1, minsize = 50)#, weight=1)
        self.plot_paramframe.grid_columnconfigure(2, minsize = 100)#, weight=1)
        self.plot_paramframe.grid_columnconfigure(3, minsize = 50)#, weight=1)    
'''
def plot_param_top_grid(self):
    # create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk


    # if plot_param frame already exists, delete it and create a new one in it's place
    # using pack it seems to be hard to just overwrite  as I can with grid
    if self.plot_menu_exists ==1:
        self.plot_paramframe.destroy() 

    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.plot_paramframe = ttk.LabelFrame(self.mframe, text = 'Trace Display Options')#, height = 5, width = 800)
    self.plot_paramframe.grid(row = 0, column = 0, padx = 2, pady=2, sticky='nw')

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.plot_paramframe,text="Amp. Multiplier" )
    labelfsplt = tk.Label(master = self.plot_paramframe,text="Samp. Rate (Hz)" )
    labelsplt= tk.Label(master = self.plot_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.plot_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.plot_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.plot_paramframe,text="End Trace" )
    labelips = tk.Label(master = self.plot_paramframe,text="Time Scale." )
    labeltpi = tk.Label(master = self.plot_paramframe,text="Trace Scale" )
    labelskip = tk.Label(master = self.plot_paramframe,text="Header Deci." ) 

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labelips.grid(row=0, column=6,ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=1, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=0, column=8, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryips = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.plot_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=0, column=5, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=1, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryips.grid(row=0, column=7, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=1, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=0, column=9, ipadx=0, padx=5, pady=5, sticky='w')

    # data plotting frame - create action  buttons
    if self.load_type == 1:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.raw_plot)
    elif self.load_type == 2:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.xcorr_plot)
    elif self.load_type == 3:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.stack_plot)        
    buttontrplt.grid(row=0, column=12, columnspan=1, sticky='w',padx=10, pady=5 )
    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')

    cknorm = ttk.Checkbutton(self.plot_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.plot_paramframe, text='Flip Pol.', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.plot_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.plot_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    cknorm.grid(row=0, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=1, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=0, column=11, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=1, column=11, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    # Set some default values which can be manually updated
    self.entryscal.insert(0, '1')
    self.entryips.insert(0,'5')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')    
    # Insert default values that are read from SEG-D data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entryeplt.delete(0, tk.END)
    self.entryeplt.insert(0, '%s'%(self.plottindex.max()))
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))

    self.plot_menu_exists = self.plot_paramframe.winfo_exists() 
def plot_param_top(self):
    # create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk

    # if plot_param frame already exists, delete it and create a new one in it's place
    # using pack it seems to be hard to just overwrite  as I can with grid
    if self.plot_menu_exists ==1:
        self.plot_paramframe.destroy()
    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.plot_paramframe = ttk.LabelFrame(self.mframe, text = 'Trace Display Options', height = 5, width = 800)
    self.plot_paramframe.pack( anchor = "w",side = "top",fill="x", expand=False, padx = 5, pady = 0, )      

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.plot_paramframe,text="Amp. Multiplier" )
    labelfsplt = tk.Label(master = self.plot_paramframe,text="Samp. Rate (Hz)" )
    labelsplt= tk.Label(master = self.plot_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.plot_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.plot_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.plot_paramframe,text="End Trace" )
    labelips = tk.Label(master = self.plot_paramframe,text="Time Scale." )
    labeltpi = tk.Label(master = self.plot_paramframe,text="Trace Scale" )
    labelskip = tk.Label(master = self.plot_paramframe,text="Header Deci." ) 

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=0, column=4,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=1, column=4, ipadx=0,padx=5, pady=5, sticky='w')
    labelips.grid(row=0, column=6,ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=1, column=6, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=0, column=8, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryips = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.plot_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=0, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=0, column=5, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=1, column=5, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryips.grid(row=0, column=7, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=1, column=7, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=0, column=9, ipadx=0, padx=5, pady=5, sticky='w')

    # data plotting frame - create action  buttons
    if self.load_type == 1:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.raw_plot)
    elif self.load_type == 2:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.xcorr_plot)
    elif self.load_type == 3:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.stack_plot)        
    buttontrplt.grid(row=0, column=12, columnspan=1, sticky='w',padx=10, pady=5 )
    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')

    cknorm = ttk.Checkbutton(self.plot_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.plot_paramframe, text='Flip Pol.', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.plot_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.plot_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    cknorm.grid(row=0, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=1, column=10, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=0, column=11, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=1, column=11, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    # Set some default values which can be manually updated
    self.entryscal.insert(0, '1')
    self.entryips.insert(0,'5')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')    
    # Insert default values that are read from SEG-D data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entryeplt.delete(0, tk.END)
    self.entryeplt.insert(0, '%s'%(self.plottindex.max()))
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))

    self.plot_menu_exists = self.plot_paramframe.winfo_exists() 
  
def plot_param_side(self):
    # create a sidebar to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk

#    self.plwindow = tk.Toplevel(self.master)
#    self.plwindow = tk.Toplevel(self.mframe)
#
#    self.plotwin = tk.PanedWindow(self.plwindow, orient=tk.HORIZONTAL, width = 1200,
#    relief='groove',borderwidth=2)
#    self.plotwin.grid(row=0, column=0,columnspan = 3 , padx=10, pady=0, sticky ="nsew")

    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.plot_paramframe = ttk.LabelFrame(self.mframe, text = 'Trace Display Options', width = 200)
    self.plot_paramframe.grid(row = 0, column = 0, padx = 10, pady=5, sticky='ns')
    self.plot_paramframe.grid_columnconfigure(0, minsize = 75)#, weight=1)
    self.plot_paramframe.grid_columnconfigure(1, minsize = 40)#, weight=1)
#    self.plot_paramframe.grid_columnconfigure(2, minsize = 75)#, weight=1)
#    self.plot_paramframe.grid_columnconfigure(3, minsize = 40)#, weight=1)        

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.plot_paramframe,text="Amp. Multiplier" )
    labelfsplt = tk.Label(master = self.plot_paramframe,text="Samp. Rate (Hz)" )
    labelsplt= tk.Label(master = self.plot_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.plot_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.plot_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.plot_paramframe,text="End Trace" )
    labelips = tk.Label(master = self.plot_paramframe,text="Time Scale." )
    labeltpi = tk.Label(master = self.plot_paramframe,text="Trace Scale" )
    labelskip = tk.Label(master = self.plot_paramframe,text="Header Deci." ) 

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=2, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=3, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=4, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=5, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelips.grid(row=6, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=7, column=0, ipadx=0,padx=5, pady=5, sticky='w')
    labelskip.grid(row=8, column=0, ipadx=0,padx=5, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryips = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.plot_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=1, column=1, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=2, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=3, column=1, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=4, column=1, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=5, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryips.grid(row=6, column=1, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=7, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=8, column=1, ipadx=0, padx=5, pady=5, sticky='w')

    # data plotting frame - create action  buttons
    if self.load_type == 1:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.raw_plot)
    elif self.load_type == 2:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.xcorr_plot)
    elif self.load_type == 3:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.stack_plot)        
    buttontrplt.grid(row=13, column=0, columnspan=1, sticky='w',padx=10, pady=5 )
    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')

    cknorm = ttk.Checkbutton(self.plot_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.plot_paramframe, text='Flip Pol.', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.plot_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.plot_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    cknorm.grid(row=9, column=0, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=9, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=10, column=0, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=10, column=1, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    # Set some default values which can be manually updated
    self.entryscal.insert(0, '1')
    self.entryips.insert(0,'5')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')    
    # Insert default values that are read from SEG-D data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entryeplt.delete(0, tk.END)
    self.entryeplt.insert(0, '%s'%(self.plottindex.max()))
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))

def plot_param(self):
    # create a toplevel window (pop-out) to get plot parameters
    # self.master is the root window in a form which can be passed to functions
    
    import tkinter as tk
    from tkinter import ttk

#    self.plwindow = tk.Toplevel(self.master)
    self.plwindow = tk.Toplevel(self.mframe)

    self.plotwin = tk.PanedWindow(self.plwindow, orient=tk.HORIZONTAL, width = 1200,
    relief='groove',borderwidth=2)
    self.plotwin.grid(row=0, column=0,columnspan = 3 , padx=10, pady=0, sticky ="nsew")

    # create a fancy frame for plot parameter entry in column 0, row 1 of master window
    self.plot_paramframe = ttk.LabelFrame(self.plotwin, text = 'Trace Display Options', width = 400)
    self.plot_paramframe.grid(row = 0, column = 0, padx = 10, pady=5)#, sticky='ns')
    self.plot_paramframe.grid_columnconfigure(0, minsize = 100)#, weight=1)
    self.plot_paramframe.grid_columnconfigure(1, minsize = 50)#, weight=1)
    self.plot_paramframe.grid_columnconfigure(2, minsize = 100)#, weight=1)
    self.plot_paramframe.grid_columnconfigure(3, minsize = 50)#, weight=1)        

    # data plotting frame - insert text for user input boxes
    labelscal = tk.Label(master = self.plot_paramframe,text="Plot Scalar" )
    labelfsplt = tk.Label(master = self.plot_paramframe,text="Sample Rate (Hz)" )
    labelsplt= tk.Label(master = self.plot_paramframe,text="Start Time (ms)" )
    labeleplt = tk.Label(master = self.plot_paramframe,text="End Time (ms)" )
    labelstr = tk.Label(master = self.plot_paramframe,text="Start Trace" )
    labeletr = tk.Label(master = self.plot_paramframe,text="End Trace" )
    labelips = tk.Label(master = self.plot_paramframe,text="In. Per Sec." )
    labeltpi = tk.Label(master = self.plot_paramframe,text="Trace Per In." )
    labelskip = tk.Label(master = self.plot_paramframe,text="Every nth Header" ) 

    labelscal.grid(row=0, column=0,padx=5, pady=5, sticky='w')
    labelfsplt.grid(row=0, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelsplt.grid(row=1, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeleplt.grid(row=1, column=2,ipadx=0,padx=5, pady=5, sticky='w')
    labelstr.grid(row=2, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeletr.grid(row=2, column=2, ipadx=0,padx=0, pady=5, sticky='w')
    labelips.grid(row=3, column=0,ipadx=0,padx=5, pady=5, sticky='w')
    labeltpi.grid(row=3, column=2, ipadx=0,padx=0, pady=5, sticky='w')
    labelskip.grid(row=4, column=0, ipadx=0,padx=0, pady=5, sticky='w')

   # data plotting frame- get user inputs, insert defaults
    self.entryscal = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryfsplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrysplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryeplt = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrystr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryetr = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryips = tk.Entry(master=self.plot_paramframe, width=8)
    self.entrytpi = tk.Entry(master=self.plot_paramframe, width=8)
    self.entryskip = tk.Entry(master=self.plot_paramframe, width=8)
    
    self.entryscal.grid(row=0, column=1, padx=5, pady=5, sticky='w')
    self.entryfsplt.grid(row=0, column=3, padx=5, pady=5, sticky='w')
    self.entrysplt.grid(row=1, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryeplt.grid(row=1, column=3, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrystr.grid(row=2, column=1, ipadx=0,padx=5, pady=5, sticky='w')
    self.entryetr.grid(row=2, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryips.grid(row=3, column=1, ipadx=0,padx=5, pady=5, sticky='w')
    self.entrytpi.grid(row=3, column=3, ipadx=0, padx=5, pady=5, sticky='w')
    self.entryskip.grid(row=4, column=1, ipadx=0, padx=5, pady=5, sticky='w')

    # data plotting frame - create action  buttons
    if self.load_type == 1:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.raw_plot)
    elif self.load_type == 2:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.xcorr_plot)
    elif self.load_type == 3:
        buttontrplt = ttk.Button(self.plot_paramframe, text = "Make Plot", command = self.stack_plot)        
    buttontrplt.grid(row=6, column=0, columnspan=1, sticky='w',padx=10, pady=5 )
    # data plotting frame - create checkbuttons for plotting controls
    self.normchk = tk.StringVar()
    self.polchk = tk.StringVar()
    self.vachk = tk.StringVar()
    self.decichk = tk.StringVar()
    
    self.normchk.set('y')
    self.polchk.set('n')
    self.vachk.set('n')
    self.decichk.set('n')

    cknorm = ttk.Checkbutton(self.plot_paramframe, text='Normalize', variable=self.normchk, onvalue='y', offvalue = 'n')
    ckpol = ttk.Checkbutton(self.plot_paramframe, text='Flip Polarity', variable=self.polchk, onvalue='r', offvalue = 'n')
    ckva = ttk.Checkbutton(self.plot_paramframe, text='VA plot', variable=self.vachk, onvalue='y', offvalue = 'n')
    ckdeci = ttk.Checkbutton(self.plot_paramframe, text='Decimate', variable=self.decichk, onvalue='y', offvalue = 'n')
    cknorm.grid(row=5, column=0, ipadx=0, padx=5, pady=5, sticky='w')
    ckpol.grid(row=5, column=1, ipadx=0, padx=5, pady=5, sticky='w')
    ckva.grid(row=5, column=2, ipadx=0, padx=5, pady=5, sticky='w')
    ckdeci.grid(row=5, column=3, ipadx=0, padx=5, pady=5, sticky='w')

   # Set some default values which are updated by SEG-D header values
    self.entryfsplt.insert(0, '1000')
    self.entrysplt.insert(0, '0')
    self.entryeplt.insert(0, '2000')
    self.entrystr.insert(0, '1')
    self.entryetr.insert(0, '10')
    # Set some default values which can be manually updated
    self.entryscal.insert(0, '1')
    self.entryips.insert(0,'5')
    self.entrytpi.insert(0,'10')
    self.entryskip.insert(0,'1')    
    # Insert default values that are read from SEG-D data file
    self.entryfsplt.delete(0, tk.END)
    self.entryfsplt.insert(0, '%s'%(self.fs))
    self.entrysplt.delete(0, tk.END)
    self.entrysplt.insert(0, '0')
    self.entryeplt.delete(0, tk.END)
    self.entryeplt.insert(0, '%s'%(self.plottindex.max()))
    self.entrystr.delete(0, tk.END)
    self.entrystr.insert(0, '1')
    self.entryetr.delete(0, tk.END)
    self.entryetr.insert(0, '%s'%(self.plotdata.shape[0]))    

def scroll_plot(self):
    # Used in SEG-D reader apps, requires tkinter
    # Parameters, arguments, variables, data are passed as object (self)
    
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib  import gridspec
    
    import tkinter as tk    

    import numpy as np
   
    # increase size of right_frame to increase scrolled area
    print("\u0332".join('\nWiggle Plot Global Information :'))
    
    self.first = int(self.entrystr.get())-1
    self.last = int(self.entryetr.get())

    print (' self.plotdata.shape :',self.plotdata.shape)
    print (' self.plotdata.dtype :',self.plotdata.dtype)
    print(' self.first :',self.first,' self.last :',self.last)

    # general settings - font size, rasterizing
    plt.rcParams.update({'font.size': 7})
    rasterized = 'True' # True can speed up large data sets at expense of resolution

    # get user input plot parameters, .get returns strings!
    scal = float(self.entryscal.get())
    fs = int(float(self.entryfsplt.get())) # default is read from file
    stime = int(self.entrysplt.get())
    etime = int(float(self.entryeplt.get()))
    pol = self.polchk.get()
    norm = self.normchk.get()
    va = self.vachk.get()
    dec = self.decichk.get()
    skiplabel = int(self.entryskip.get())
    
    print(' fs :',fs)

    ips = int(self.entryips.get())
    tpi = int(self.entrytpi.get())
    
    # label and trace spacing
    skip =1 # trace spacing
    spacing=1 # to avoid spacing = z
#    skiplabel=skip
    
    # scale the plot width by traces/inch, 12 being baseline plot width
    xscal = tpi/10 # 10 tpi standard
    yscal = 5/ips # 5 ips standard

    # decimating samples to speed plotting        
    ssamp_predeci = int(stime*fs/1000)
    esamp_predeci = int(etime*fs/1000)
    
    fs_deci=fs
    
    # in an effort to speed plotting, I tried using 16 bit floating point
    # does not work well after correlation, get a lot of nans 
    # in row_max when normalizing
    
    VSPtrim = self.plotdata[self.first:self.last,
                        ssamp_predeci:esamp_predeci]

    decimate_factor = 1 
    
    if (dec=='y') or (dec =='Y'):
        if fs == 500:
            # decimate by 2 and repeat 4 times
            # convert to float16 to speed plotting
            decimate_factor = 1
            fs_deci = 250
            for n in range(0,decimate_factor):
                VSPtrim=np.float16(VSPtrim[:,::2])
           
        if fs == 1000:
            decimate_factor = 2
            fs_deci = 250
            for n in range(0,decimate_factor):
                VSPtrim=np.float16(VSPtrim[:,::2])
        if fs == 2000:
            decimate_factor = 3
            fs_deci = 250
            for n in range(0,decimate_factor):
                VSPtrim=np.float16(VSPtrim[:,::2])

    print (' VSPtrim.dtype :', VSPtrim.dtype)
    print (' VSPtrim.shape :',VSPtrim.shape)

    # recalculate indices after decimating
    ssamp_deci = int(stime*fs_deci/1000)
    esamp_deci = int(etime*fs_deci/1000)
    
    # x axis if plotting by trace number
    trnum = np.arange(self.first+1,self.last+1,1)

    # y axis in time
    y = np.arange(ssamp_deci*(1000/fs_deci),esamp_deci*(1000/fs_deci),(1000/fs_deci)  )

    data2 = np.zeros(shape = (VSPtrim.shape[0], VSPtrim.shape[1]))
    data1 = VSPtrim

    # apply normalization if requested
    if (norm == 'Y') or (norm =='y'):
        row_max = np.max(np.abs(data1), axis=1)
        where_0 = np.where(row_max == 0) # find traces of all 0s
        row_max[where_0] = 1 # set 0 to 1
        data2 = (data1 / row_max[:, np.newaxis])
        datascaled = data2 * scal

    else:
        datascaled = data1 * scal

    # flip polarity if requested
    if (pol == 'r') or (pol =='R'):
        datascaled = datascaled * -1

    # choose depth spacing or trace number spacing
    if (spacing == 'Z') or (spacing == 'z'):
        dscaler, pad = (rcv_depth, 10)
        dlabel = 'Receiver Depth'

    else:
        dscaler, pad = (trnum, 1)
        dlabel = 'Receiver Number'

    # prepare header track data
    channel_set = self.plotthead[self.first:self.last,0]
    channel_number = self.plotthead[self.first:self.last,1]
    file_number = self.plotthead[self.first:self.last,4]
    
    yaxis = trnum[::skiplabel]*0+.3 # arbitrary y axis height
    xaxis = trnum[::skiplabel]
    
    trace_num = trnum[::skiplabel]
    ch_set = channel_set[::skiplabel]
    ch_num = channel_number[::skiplabel]
    file_num = file_number[::skiplabel]
    
    track_labels = ['Seq_Trace_number','Chan_set', 'Trace_number', 'File Number']
    
    ##############################################################################
   # create a tkinter canvas and scroll bars
    canvas = tk.Canvas(self.plot_frame,borderwidth=1,relief=tk.RIDGE)#highlightthickness=1, highlightbackground="black")#plot_frame)
    canvas.grid(row=0, column=0, sticky=tk.NSEW)#, padx = 0, pady =0,

    #Trace Plot Scroll bars
    xScrollbar = tk.Scrollbar(self.plot_frame, orient=tk.HORIZONTAL)
    yScrollbar = tk.Scrollbar(self.plot_frame)        
#    xScrollbar = tk.Scrollbar(self.master, orient=tk.HORIZONTAL)
#    yScrollbar = tk.Scrollbar(self.master)        
        
    xScrollbar.grid(row=1, column=0, sticky=tk.EW)
    yScrollbar.grid(row=0, column=1, sticky=tk.NS)

    canvas.config(xscrollcommand=xScrollbar.set)#,xscrollincrement=0)
    canvas.config(yscrollcommand=yScrollbar.set)#,yscrollincrement=0)
#    speed = 2
#    canvas.xview_scroll(speed , 'units')
#    canvas.yview_scroll(speed , 'units')

    xScrollbar.config(command=canvas.xview)
    yScrollbar.config(command=canvas.yview)
    
    ##############################################################################

    # set up plot axes
    fig = Figure(figsize=(12,9)) # use Figure, not pyplot for tkinter to work
#    fig.suptitle('%s'%(self.ptitle), fontsize=14)
    gs = fig.add_gridspec(5, 1, height_ratios=[.15,.15,.15,.15, 2], wspace = 0, hspace = .02)

    # add 4 tracks for headers
    for n in range(4):
        ax = fig.add_subplot(gs[n])

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel(track_labels[n], rotation = 0 )
        ax.set_yticklabels([])
        ax.set_xticks(trnum[:-1:1])
        ax.tick_params(axis = 'x', direction = 'in')
        ax.set_xticklabels([])
        if (n==0):
            for i, txt in enumerate(trace_num):
                ax.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)
                ax.set_xlim((np.min(trnum)-pad), (np.max(trnum) + pad)*xscal )
                
        if (n==1):
            for i, txt in enumerate(ch_set):
                ax.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)
                ax.set_xlim((np.min(trnum)-pad), (np.max(trnum) + pad)*xscal )

        if (n==2):
            for i, txt in enumerate(ch_num):
                ax.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)
                ax.set_xlim((np.min(trnum)-pad), (np.max(trnum) + pad)*xscal )                    
        if (n==3):
            for i, txt in enumerate(file_num):
                ax.annotate(txt, (xaxis[i], yaxis[i]), rotation = 90)
                ax.set_xlim((np.min(trnum)-pad), (np.max(trnum) + pad)*xscal )
                
    # add the main track with traces 
    ax4 = fig.add_subplot(gs[4])

    for i, trace in enumerate(datascaled[::skip, :]):
        x = trace + dscaler[i]

        ax4.plot(x, y, 'k-', linewidth = .5, rasterized = rasterized)
        if (va=='y')or(va=='Y'):
            ax4.fill_betweenx(y, dscaler[i], x, where=(x > dscaler[i]), color='k',
            rasterized = rasterized)
        ax4.set_xlim((dscaler[0]-pad), (dscaler[-1]+pad)*xscal)
        ax4.set_ylim(etime*yscal, stime)
        ax4.set_ylabel(' Time (ms)')
        ax4.set_xticks([])

    #  remove whitespace at top and bottom of plot, adjust if a title is required
    fig.subplots_adjust( top=0.99, bottom=0.01)
    
    plt.show()

    figAgg = FigureCanvasTkAgg(fig, canvas)
    mplCanvas = figAgg.get_tk_widget()
    canvas.create_window(0, 0, window=mplCanvas, anchor=tk.NW)
    canvas.config(scrollregion=canvas.bbox(tk.constants.ALL))

