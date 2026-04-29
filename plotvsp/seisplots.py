
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
#        ax2.plot(x, y, 'k-',  linewidth = .5)
        ax2.plot(x, y, 'k-',  linewidth = .5)
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
        print (' VSPdata.shape :', VSPdata.shape,' (traces,samples)')
        print(' VSPdata type :', VSPdata.dtype)
        print (' Max an Min Amplitude VSPdata :',np.nanmax(VSPdata),np.nanmin(VSPdata))
        print (' datascaled.shape ',datascaled.shape,' (traces,samples)')    
        print (' thead shape :', thead.shape,' (traces,header columns)')
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
    timframe = set the time header to 0 if data is aligned, or chose the TWT header,
               defaults to OWT
    Title_plot = 'plot title '
    
    Example Usage:
    
    plot_params = {"pol":'n', 
                    "Tmax":3500, "Tmin":0, 
                    "first_rcv":first_rcv, 
                    "spacing":'z', 
                    "skiplabel":4, 
                    "fs":fs, 
                    "norm":'n',
                    "scal":1000, 
                    "title_top":'Corridor Mute ',
                    "info_wig":'y',
                    "timframe":'twt'} 
    composite_plot(thead_dec_edit,corr_in,corr_stk, **plot_params)

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
    scalv = kwargs['scal_main']
    scalc = kwargs['scal_stk']
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
    trace_num_cstk = np.arange(1, Cstack.shape[0]+1)    

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
        datascaled = data2 * scalv        
    else:        
        datascaled = data1 * scalv
    
    # we don't trace normalize the cstack, but we do scale it up
    datascaled_cstk = Cstack * scalc
       
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
    '''
    print(' np.nanmax(np.abs(VSPdata)) :',np.nanmax(np.abs(VSPdata)))
    
    # alternative normalization using amplitude ratios
    cstak_rat= np.nanmax(np.abs(datascaled_cstk))/np.nanmax(np.abs(datascaled))
    datascaled_cstk=datascaled_cstk/cstak_rat
    print (' first np.nanmax(np.abs(datascaled_cstk)) :',np.nanmax(np.abs(datascaled_cstk)))
    print(' np.nanmax(np.abs(datascaled)) :',np.nanmax(np.abs(datascaled)),' cstak_rat :',cstak_rat )
    
    # scale the corridor stack to somewhat similar amplitude range as input to cstack
    # ie. each trace amplitude is added to trace number    
    datascaled_cstk = Cstack * abs(scalv/((trace_num.max()+pad)/(trace_num_cstk.max()-padc*2)))
    print (' cstack scalar :', scalv/((trace_num.max()+pad)/(trace_num_cstk.max()-padc*2)) )
    print (' second np.nanmax(np.abs(datascaled_cstk)) :',np.nanmax(np.abs(datascaled_cstk)) ,
    ''' 
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
        print("\u0332".join('Composite Plot Global Information :')) 
        print (' VSPdata.shape :', VSPdata.shape,' (traces,samples)')
        print(' VSPdata type :', VSPdata.dtype)
        print (' Max an Min Amplitude VSPdata :',np.nanmax(VSPdata),np.nanmin(VSPdata))
        print (' Max an Min Amplitude Cstack :',np.nanmax(Cstack),np.nanmin(Cstack))
        print (' datascaled.shape ',datascaled.shape,' (traces,samples)')    
        print (' thead shape :', thead.shape,' (traces,header columns)')
        print (' Min TVDSRD - pad', np.min(TVDSRD)-pad, ' Pad :', pad)    
        print (' Max TVDSRD + pad', np.max(TVDSRD)+pad, ' Pad :', pad)    
        print (' min max intvel :', np.min(intVel), np.max(intVel))

        
def plotsingletrace( VSP1, Tmin, Tmax, thead, spacing, scal1, title):

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
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
    scale = kwargs['scal']              # scaling to plot 1,2,3,4 amplitudes    
    trange = kwargs['time_range']   
    drange = kwargs['depth_range']
    info_gray4 = kwargs['info_plot']
    
    rcvdepth = thead[:, 2]
    numsamp = VSP1.shape[1]
    
    cmap = 'seismic'
    
    ################ create plot tick labeling and indexing ##############
    
    tindex1 = np.arange(0, numsamp*(1000/fs),(1000/fs) )  # convert fs to msec.

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

    ax1.imshow(VSP1[:,:numsamp-1].T, cmap=cmap, interpolation='none', 
               vmin = -np.nanmax(VSP1)/scale[0],vmax = np.nanmax(VSP1)/scale[0],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
                tindex1.min()], aspect = 'auto')

    ax1.yaxis.grid(c = 'black', lw = .1)    
    ax1.set_ylim(trange[1], trange[0]) # extents must be set    
    ax1.set_xlim(drange[0], drange[1])    
    ax1.set_xlabel('Receiver Depth')    
    ax1.set_title('%s'%(txt1))

    ax2.imshow(VSP2[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP2)/scale[1],vmax = np.nanmax(VSP2)/scale[1],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
                tindex1.min()], aspect = 'auto')    
    ax2.yaxis.grid(c = 'black', lw = .1)    
    ax2.set_ylim(trange[1], trange[0]) # extents must be set    
    ax2.set_xlim(drange[0], drange[1])    
    ax2.set_xlabel('Receiver Depth')    
    ax2.set_title('%s '%(txt2))    
    ax2.set_yticklabels([])
 
    ax3.imshow(VSP3[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP3)/scale[2],vmax = np.nanmax(VSP3)/scale[2],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')
    
    ax3.yaxis.grid(c = 'black', lw = .1)    
    ax3.set_ylim(trange[1], trange[0]) # extents must be set    
    ax3.set_xlim(drange[0], drange[1])    
    ax3.set_xlabel('Receiver Depth')    
    ax3.set_title('%s '%(txt3))    
    ax3.set_yticklabels([])  
 
    ax4.imshow(VSP4[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP4)/scale[3],vmax = np.nanmax(VSP4)/scale[3],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')

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
        print("\u0332".join('\nColor Multi Plot Global Information :'))
        print (' Number of traces in plot :', VSP1.shape[0], 
           ' Number of samples per trace :', ' tindex1 shape :,',tindex1.shape,
           ' tindex1 min max :', tindex1.min(),tindex1.max(),)            
 
        print (' Max an Min Amplitude Box 1. :',np.nanmax(VSP1),np.nanmin(VSP1))
        print (' Max an Min Amplitude Box 2 :',np.nanmax(VSP2),np.nanmin(VSP2))
        print (' Max an Min Amplitude Box 3 :',np.nanmax(VSP3),np.nanmin(VSP3))
        print (' Max an Min Amplitude Box 4 :',np.nanmax(VSP4),np.nanmin(VSP4))
        
def decon_plots(VSP1, VSP2, VSP3, VSP4, VSP5,VSP6,fs, thead, **kwargs):    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec

  
    txt1 = kwargs['ss_title1']
    txt2 = kwargs['ss_title2']    
    txt3 = kwargs['ss_title3']
    txt4 = kwargs['ss_title4']
    txt5 = kwargs['ss_title5']
    txt6 = kwargs['ss_title6']

    save = kwargs['savePNG']  
    png_txt = kwargs['png_name']       
    scale = kwargs['scal']              # scaling to plot 1,2,3,4 amplitudes    
    trange = kwargs['time_range']   
    drange = kwargs['depth_range']
    info_gray4 = kwargs['info_plot']
    
    rcvdepth = thead[:, 2]
    
    # Set the number of samples to that of the aligned downgoing
    # This avoids creating a time index for every track as the TWT data 
    # sets are longer than necessary due to aligning/unaligning process
    # This requires all data sets to have the same number of samples
    numsamp = VSP1.shape[1]
    
    # Choose a red-blue color map for the seismic plots
    cmap = 'seismic'
    
    ################ create plot tick labeling and indexing ##############
    
    tindex1 = np.arange(0, numsamp*(1000/fs),(1000/fs) )  # convert fs to msec.    
    rindex  = np.stack([rcvdepth for _ in range(VSP1.shape[0])], axis=1)
    
    ############### make plots of input and SVD filtered output ################
    
    fig = plt.figure(figsize=(16,6) )    
#    fig.subplots_adjust(wspace=0.125, hspace=0.5)

    from matplotlib.gridspec import GridSpec
    
    gs1 = GridSpec(1, 32, hspace = .25, wspace=0.55) # make a row by col grid
    ax1 = fig.add_subplot(gs1[0:1, 0:6])            # combine rows or columns
    ax2 = fig.add_subplot(gs1[0:1, 6:12])
    ax3 = fig.add_subplot(gs1[0:1, 12:18])
    ax4 = fig.add_subplot(gs1[0:1, 18:24])
    ax5 = fig.add_subplot(gs1[0:1, 24:30])
    ax6 = fig.add_subplot(gs1[0:1, 31:32])

    ax1.imshow(VSP1[:,:numsamp-1].T, cmap=cmap, interpolation='none', 
               vmin = -np.nanmax(VSP1)/scale[0],vmax = np.nanmax(VSP1)/scale[0],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
                tindex1.min()], aspect = 'auto')

    ax1.yaxis.grid(c = 'black', lw = .1)    
    ax1.set_ylim(trange[1], trange[0]) # extents must be set    
    ax1.set_xlim(drange[0], drange[1])    
    ax1.set_xlabel('Receiver Depth')    
    ax1.set_title('%s'%(txt1))

    ax2.imshow(VSP2[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP2)/scale[1],vmax = np.nanmax(VSP2)/scale[1],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
                tindex1.min()], aspect = 'auto')    
    ax2.yaxis.grid(c = 'black', lw = .1)    
    ax2.set_ylim(trange[1], trange[0]) # extents must be set    
    ax2.set_xlim(drange[0], drange[1])    
    ax2.set_xlabel('Receiver Depth')    
    ax2.set_title('%s '%(txt2))    
    ax2.set_yticklabels([])
 
    ax3.imshow(VSP3[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP3)/scale[2],vmax = np.nanmax(VSP3)/scale[2],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')
    
    ax3.yaxis.grid(c = 'black', lw = .1)    
    ax3.set_ylim(trange[1], trange[0]) # extents must be set    
    ax3.set_xlim(drange[0], drange[1])    
    ax3.set_xlabel('Receiver Depth')    
    ax3.set_title('%s '%(txt3))    
    ax3.set_yticklabels([])  
 
    ax4.imshow(VSP4[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP4)/scale[3],vmax = np.nanmax(VSP4)/scale[3],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')

    ax4.yaxis.grid(c = 'black', lw = .1)    
    ax4.set_ylim(trange[1], trange[0]) # extents must be set    
    ax4.set_xlim(drange[0], drange[1])    
    ax4.set_xlabel('Receiver Depth')    
    ax4.set_title('%s '%(txt4))
    ax4.set_yticklabels([])    
#    ax4.yaxis.set_label_position('right')    
#    ax4.yaxis.set_ticks_position('right')
    
    ax5.imshow(VSP5[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP5)/scale[4],vmax = np.nanmax(VSP5)/scale[4],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')

    ax5.yaxis.grid(c = 'black', lw = .1)    
    ax5.set_ylim(trange[1], trange[0]) # extents must be set    
    ax5.set_xlim(drange[0], drange[1])    
    ax5.set_xlabel('Receiver Depth')
    ax5.set_yticklabels([])    

    ax5.set_title('%s '%(txt5))
#    ax5.yaxis.set_label_position('right')    
#    ax5.yaxis.set_ticks_position('right')
    
    ax6.imshow(VSP6[:,:numsamp-1].T, cmap=cmap, interpolation='none',
               vmin = -np.nanmax(VSP6)/scale[5],vmax = np.nanmax(VSP6)/scale[5],
               extent = [rindex.min(), rindex.max(), tindex1.max(), 
               tindex1.min()], aspect = 'auto')

    ax6.yaxis.grid(c = 'black', lw = .1)    
    ax6.set_ylim(trange[1], trange[0]) # extents must be set    
    ax6.set_xlim(drange[0], drange[1])    
    ax6.set_xlabel('Trace')    
    ax6.set_title('%s '%(txt6))
    ax6.yaxis.set_label_position('right')    
    ax6.yaxis.set_ticks_position('right')

    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\procflow_gray_%s.png' 
        %(png_txt), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)
            
    plt.show()
    
    if(info_gray4=='y')or(info_gray4=='Y'):        
        print("\u0332".join('\nDecon Plot Global Information :'))
        print (' Number of traces in plot :', VSP1.shape[0], 
           ' Number of samples per trace :', ' tindex1 shape :,',tindex1.shape,
           ' tindex1 min max :', tindex1.min(),tindex1.max())            

        print (' Max an Min Amplitude Box 1. :',np.nanmax(VSP1),np.nanmin(VSP1))
        print (' Max an Min Amplitude Box 2. :',np.nanmax(VSP2),np.nanmin(VSP2))
        print (' Max an Min Amplitude Box 3. :',np.nanmax(VSP3),np.nanmin(VSP3))
        print (' Max an Min Amplitude Box 4. :',np.nanmax(VSP4),np.nanmin(VSP4))
        print (' Max an Min Amplitude Box 5. :',np.nanmax(VSP5),np.nanmin(VSP5))
        print (' Max an Min Amplitude Box 6. :',np.nanmax(VSP6),np.nanmin(VSP6))
        
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
    
    Example usage:
    
    cplot_params = {"pol":'n', 
                    "Tmax":3500, "Tmin":500, 
                    "show_time":'y', 
                    "skiplabel":4,
                    "samp_decimate":1,
                    "fs":fs, 
                    "norm":'n',
                    "min_amp":-.1, "max_amp":.1,
                    "color":'seismic',
                    "title_top":'Decon Upgoing TWT '}
    plotcolor(thead_dec_edit,corr_in,**cplot_params)
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
    scal = kwargs['scal']
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
                           cmap = colr, shading = 'gouraud')
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
