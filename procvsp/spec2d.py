
    
def fktran_numpy(seis,headers, fs, comment1, linevel):
    '''2D fft to generate fk plot
    No padding of receivers or time samples

    seis - VSP data matrix
    headers - need receiver depths
    fs - sampling rate in hertz
    linevel - velocity line for annotating f/k plot

    This assumes constant receeiver spacing!
    '''

    import numpy as np
    import scipy.signal as sig

    import procvsp.spec as Spec
    import procvsp.spec2d as Spec2d
    
    numsamp, ntrace = seis.T.shape            # data needs to be 1 trace/column    
    t = np.arange(0, numsamp*(1/fs),(1/fs) )  # convert fs in hertz to seconds    
    dt = t[1]-t[0]                            # could have used 1/fs b
    
    rcvz = headers[:,1] #rcvdepth               # needs constant rcv separation
    dx = rcvz[1] - rcvz[0]    
    nf = int(2 ** Spec.nextpow2(len(t)))
#    nk = int(2 ** nextpow2(len(x)))          # use for minimizing array size     
    nk  = nf                                  # use for better resolution in k 
    plotnf = int(nf/2+1)                      # plot only positive frequencies
    
    print("\u0332".join('\nfktran_numpy info and parameters before fft') )            
    print(' seis.T shape :',seis.T.shape,'\n numsamp :', numsamp, \
          '\n ntrace :', ntrace)       
    print(' rcvz shape :',rcvz.shape,'\n t shape :', t.shape, '\n dx :', \
          dx,'\n dt :', dt)    
    print(' nf :',nf,'\n nk :', nk)
    
    fk = np.fft.fftshift(np.fft.fft2(seis,s=(nk,nf))) # fftshift puts 0 freq and 
                                                      # 0 k in middle of output          
    freq = np.fft.fftfreq(fk.shape[1], dt)            # get the freq. bin values 
                                                      #from axis 1 of fk array
    freq_unwrap = np.fft.fftshift(freq)               # re-arrange so 0 freq is 
                                                      # the middle of the array
    
    k = np.fft.fftfreq(fk.shape[0], dx)    
    k_unwrap = np.fft.fftshift(k)
    print(' fk shape :',fk.shape)    
    print ( ' freq min', freq.min(),' freq max', freq.max())    
    print ( ' k min', k.min(),' k max', k.max())
            
    fkplot_params={'fs':fs,
               'Fmin':0,'Fmax':freq_unwrap.max(), # 1 for display of positive freqs.
               'Kmin':k_unwrap.min(),'Kmax':k_unwrap.max(), # 1 for display of positive freqs.
               'scalar':1,
               'DBmax':-160,
               'line_velocity':-6000,
               'title_fk':'FK Transform of BPF Raw Z'}
    Spec2d.fk_plot(fk, freq_unwrap, k_unwrap, rcvz,**fkplot_params)
    
    return fk, freq_unwrap, k_unwrap, numsamp,ntrace
    
def fktrans(seis,headers, fs, comment1, linevel, xpad,tpad):
    '''2D fft to generate fk plot
    seis : VSP data matrix
    headers : need receiver depths
    fs : sampling rate in hertz
    linevel : velocity line for annotating f/k plot
    xpad : size (in x units) of spatial zero pad to be afixed to seis
    tpad : size (in t units) of temporal zero pad to be afixed to seis
    
    Returns : 
    fk: 2D frquency-wavenumber spectrum
    freq_unwrap : 1d frequency vector, symmetric about 0 hz  
    k_unwrap : 1d wavenumber vector, symmetric about 0 
    seis.shape[1] : number of traces after padding of seismic
                    Used to truncate inverse f/k
    seis.shape[0] : number of samples after padding of seismic
                    Used to truncate inverse f/k
    padval : the input padding parameters for QC

    This assumes constant receiver spacing!
    '''

    import numpy as np
    import scipy.signal as sig

    import procvsp.spec as Spec
    import procvsp.spec2d as Spec2d

    numsamp, ntrace = seis.T.shape            # data needs to be 1 trace/column

    nx=ntrace
    t = np.arange(0, numsamp*(1/fs), (1/fs))  # convert fs in hertz to seconds    
    dt = t[1]-t[0]                            # could have used 1/fs b

    rcvz = headers[:, 1]                       # rcvdepth needs constant rcv separation
    dx = rcvz[1] - rcvz[0]

    padval=[0,0] # num padded samples, num padded traces 

    if(xpad>0):
        # pad deep end with traces of zeros
        seis_pad = np.pad(seis, ((xpad,xpad),(0, 0)), 'constant',constant_values=0)# pad axis 0 with nxpad columns of 0
        seis=seis_pad
        # pad the depth headers with extrapolated depths
#        rcvz_pad = np.arange(rcvz[-1]+dx,(rcvz[-1]+ xpad*dx)+dx,dx )  # get a depth array for padded traces
#        rcvz_new=np.hstack((rcvz,rcvz_pad))  # add pad trace depths to original trace depths
#        rcvz=rcvz_new        
        rcvz_pad_deep = np.arange(rcvz[-1]+dx,(rcvz[-1]+ xpad*dx)+dx,dx )  # get a depth array for deep padded traces
        rcvz_pad_shallow = np.arange(rcvz[0]- xpad*dx,rcvz[0]-dx,dx )  # get a depth array for shallow padded traces
        rcvz_new=np.hstack((rcvz_pad_shallow,rcvz,rcvz_pad_deep))  # add pad trace depths to original trace depths
        rcvz=rcvz_new
    
    ntpad = 0
    if(tpad>0):
        # pad zeros at ends of traces
        ntpad=int(tpad/dt)
        # pad axis 1 with ntpad columns of 0
        time_pad = np.pad(seis, ((0,0),(0, ntpad)),'constant',
                constant_values=0)  
        # get a time array after padding traces
        seis=time_pad
        numsamp_tpad = int(time_pad.shape[1])
        t = np.arange(0, numsamp_tpad*(1/fs),(1/fs) )  # convert fs in hertz to seconds    
        print (' ntpad :',ntpad,' time_pad.shape :',time_pad.shape, ' t.shape :',t.shape)
    
    # save number of padded traces and samples for future removal
    padval = [ntpad,xpad]
    
    nf=int(2**Spec.nextpow2(len(t)))
    # nk=int(2** nextpow2(len(x))) use for minimizing array size
    nk=nf                                  # use for better resolution in k
    plotnf = int(nf/2+1)                     # plot only positive frequencies

    print("\u0332".join('\nfktrans with padding info and parameters before fft'))
    print(' seis.T shape :',seis.T.shape,'\n numsamp :', numsamp, \
          '\n ntrace :', ntrace) 
    print (' seis.min() :',seis.min(),' seis.max() :',seis.max())          
    print(' rcvz shape :',rcvz.shape,'\n t shape :', t.shape, '\n dx :', \
          dx,'\n dt :', dt)
    print (' rcvz min :',rcvz.min(),' rcvz max :', rcvz.max())        
    print(' nf :',nf,'\n nk :', nk)
    
    fk = np.fft.fftshift(np.fft.fft2(seis,s=(nk,nf))) # fftshift puts 0 freq and 
                                                      # 0 k in middle of output          
    freq = np.fft.fftfreq(fk.shape[1], dt)            # get the freq. bin values 
                                                      #from axis 1 of fk array
    freq_unwrap = np.fft.fftshift(freq)               # re-arrange so 0 freq is
                                                      # the middle of the array
    
    k = np.fft.fftfreq(fk.shape[0], dx)    
    k_unwrap = np.fft.fftshift(k)
    print(' fk shape :',fk.shape)
    print (' fk.min() :',fk.min(),' fk.max() :',fk.max())     
    print ( ' freq min', freq.min(),' freq max', freq.max())    
    print ( ' k min', k.min(),' k max', k.max())
            
    fkplot_params={'fs':fs,
               'Fmin':0,'Fmax':freq_unwrap.max(), # 1 for display of positive freqs.
               'Kmin':k_unwrap.min(),'Kmax':k_unwrap.max(), # 1 for display of positive freqs.
               'scalar':1,
               'DBmax':-160,
               'line_velocity':linevel,

               'title_fk':'FK Transform of BPF Raw Z'}
    #Spec2d.fk_plot(fk, freq_unwrap, k_unwrap, rcvz,**fkplot_params)
    
    return fk, freq_unwrap, k_unwrap, seis.shape[1],seis.shape[0], padval
                   
def invfk_numpy(fk, numsamp, ntrace, pads):    
    ''' inverse 2d fft to transform to depth-time
    fk: the forward transformed data
    numsamp : number of samples in input seismic data, with padding
    ntrace : number of traces in input seismic data, with padding
    pads: remove the effect of time and depth padding
    '''        
    import numpy as np

    VSP_invtrans = np.fft.ifft2(np.fft.ifftshift(fk))
    
    print("\u0332".join('\ninv fk info and parameters ') )    
    print (' fk shape :', fk.shape)
    print (' numsamp :',numsamp, ' ntrace :',ntrace)
    print(' VSP_invtran shape :',VSP_invtrans.shape)    

    # get real component
    VSP_invtrans = np.real(VSP_invtrans[pads[1]:ntrace-pads[1],:numsamp-pads[0]])
    # remove padding when saving
    #first_trc=0
    #last trc=
    #if (pads[1]!=0):
    #    first_trc= pads[1]   
    
    return VSP_invtrans
    
def fk_plot(FKdata, f, k, rcvdepth,**kwargs):
    ''' try to plot negative and positive frequencies, with controls
    Inputs
    FKdata : frequency/wavenumber array, rows are wavenumber spectra
    f: 1d frequency array from 2d fft
    k: 1d wavenumber array from 2d fft
    rcvdepth : 1d array with receiver depths
    fmax, fmin : plot frequency limits, can be negative as fft is 2-sided
    
    Useage
    fkplot_params={'fs':fs,
               'Fmin':-60,'Fmax':60, # 1 for display of positive freqs.
               'Kmin':-.01,'Kmax':.002, # 1 for display of positive freqs.
               'scalar':1,
               'DBmax':-160,
               'line_velocity':-6000,
               'title_fk':'FK Transform of BPF Raw Z'}
    fk_plot_test(fktrans, F, K, rcvz,**fkplot_params)
    ''' 
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.ticker as mtick

    fs=kwargs['fs']
    fmax=kwargs['Fmax']
    fmin=kwargs['Fmin']
    kmax=kwargs['Kmax']
    kmin=kwargs['Kmin']
    dbmax=kwargs['DBmax']
    vel = kwargs['line_velocity']
    txt1=kwargs['title_fk']
    
    # get correct start and stop stop indices for frequency limits
    df =f[1] - f[0] # sample rate in hz        
    startf = int(f.shape[0]/2 - abs(fmin/df))  
    stopf = int(f.shape[0]/2 + abs(fmax/df)) 
    f=f[startf:stopf]

    # Find the indices of the mute start and end
    # 1. subtract mute value from every element of slowness array
    # 2. find the index where the difference is smallest
    startk = (np.abs(k-kmin)).argmin()
    stopk = (np.abs(k-kmax)).argmin()
    k=k[startk:stopk]    
    
    print("\u0332".join('\nfk plot info and parameters') )    
    print ('fs  :', fs)       
    print ('FKdata shape :', FKdata.shape)    
    print (' startf :, stopf :',startf,stopf)
    print (' startk :, stopk :',startk,stopk)
    print (' k.min(),k.max(), f.min(), f.max()', k.min(),k.max(), f.min(), f.max())
    print(' df :',df,'\n start freq index :', startf,'stop freq index :', stopf,
          '\n startk :', startk, 'stopk :',stopk)         
    
    FKdata = FKdata.T  # rows in frequency now       
    FKdata =  FKdata[startf:stopf,startk:stopk]

    FKdata = np.absolute(FKdata)    # get magnitude of complex values of FKdata 
                                    # np. real could be tried to get real comp.
    FKdata = 20*np.log10((FKdata)/np.max(FKdata))  # Convert to db 20 * np.log10\
                                                   # (S / np.max(S))
       
    print ('Trimmed FKdata shape :', FKdata.shape)       
    print (' Trimmed FK Min, max amplitude :', FKdata.min(),FKdata.max())       

    # get velocity line end points
    kvel=[0,np.max(f)/vel]
    fvel=[0,np.max(f)]
    slope=np.max(f)/(np.max(f)/vel) # line always strts at 0,0
    print(' kvel, fvel:',kvel,fvel,' slope :', slope)    
    # locate velocity annotation
    xylabel = (kmax/2, (slope*(kmax/2)))
    if (slope < 0):
        xylabel = (-1*kmax/2, -1*(slope*(kmax/2)))

    label = 'Velocity = %s'%(vel)    
    print(' xy label location:',xylabel)

    # plot spectra in dB     
    fig = plt.figure(figsize=(12,10))    
    ax1 = fig.add_subplot(111)    

    plot1 = ax1.imshow(FKdata, cmap="gist_rainbow", interpolation='none', \
            vmin = dbmax,vmax = 0,extent = [k.min(), k.max(), f.max(), f.min()],\
            aspect = 'auto')
    ax1.set_xlim(kmin, kmax)
    ax1.set_ylim(fmax, fmin)
    ax1.plot(kvel,fvel,c='k',linewidth=4, linestyle='--', marker='')
    ax1.yaxis.grid()    
    ax1.set_xlabel('Wavenumber 1\\ft')
    ax1.set_ylabel('Frequency (hz)')    
    ax1.set_title('%s'%(txt1))
    # get rotation angle for velocity annotation
    p1 = ax1.transData.transform_point((kvel[0], fvel[0]))
    p2 = ax1.transData.transform_point((kvel[1], fvel[1]))
    dy = (p2[1] - p1[1])
    dx = (p2[0] - p1[0])
    rotn = np.degrees(np.arctan2(dy, dx))
    if (vel<0):
        rotn=rotn-180
    txtbox = ax1.annotate(label, xy=xylabel, ha='center', va='center', rotation=rotn)
    txtbox.set_bbox(dict(facecolor='white', alpha=1, edgecolor='black'))
    # plot a colorbar    
    pad = 0.03    
    width = 0.02    
    pos = ax1.get_position()    
    axcol = fig.add_axes([pos.xmax + pad, pos.ymin, width, \
                          0.9*(pos.ymax-pos.ymin) ])
    fig.colorbar(plot1, label = 'Amplitude in db', cax = axcol, aspect = 40)
        
    plt.show()
    
def radon_linear(d,theader, **kwargs):
    ''' Linear radon transform using pylops chirp2d
    -Referenced to middle trace of input gather
    -Requires constant depth sampling
    
    Inputs:
    d : VSP data
    theader : trace header array
    minvel : minimum expected velocity in gather, used to calculate max slowness 
    
    Parameters
    npx: number of slopes
    pxmax : maximum slowness in gather
    
    Outputs:
    R2Op : radon operator
    dL_chirp: radon forward transform of d
    dinv_chirp: inverse radon transform to qc edge effects etc.
    '''
    import pylops
    import numpy as np
    
    fs=kwargs['fs']
    minvel=kwargs['minvel']
    pclip=kwargs['plot_scaler']
    # number of slopes in operator
    npx = d.shape[0] # cannot control number of slopes, always number of data traces
    
    # slowness sampling
    pxmax = 1/minvel #1e-5  # s/ft
    px = np.linspace(-pxmax, pxmax, npx)

    # depth sampling-assumed to be constant
    xaxis = theader[:,1]
    dx=(xaxis[1]-xaxis[0])/2
    
    # get time axis for data plots
    dt=1/fs # time sample rate in s.
    numsamp, ntrace = d.T.shape
    print (pxmax, dx, dt)
    taxis = np.arange(0, numsamp*(1/fs),(1/fs) ) 

    #do the transform - calculate the operator
    R2Op = pylops.signalprocessing.ChirpRadon2D(taxis, xaxis, pxmax*dx/dt, dtype="float64")
    
    # this is the radon spectrum
    dL_chirp = R2Op * d
    
    # calculate the transform adjoint
    dadjoint_chirp = R2Op.H * dL_chirp

    #calculate the inverse transform with no filtering as a qc
    dinv_chirp = R2Op.inverse(dL_chirp).reshape(R2Op.dimsd)
    
    print (' VSP data shape :',d.shape,
           ' Operator shape :',R2Op.shape,
           ' Spectrum shape: ',dL_chirp.shape,
           ' Adjoint shape :',dadjoint_chirp.shape, )
    
    return R2Op,dL_chirp,dinv_chirp,px

def radon_filter(data_edit,TauPop,TauPspec,P_slopes, **kwargs):
    ''' Mute rectangular part of the radon transform
    
    Inputs:
    data_edit : VSP data
    TauPop : Tau-P (linear radon) operator
    TauPspec : linear forward transform of data_edit
    P_slopes : array of slowness(slopes) used in linear forward transform
    fs : sampling rate of data_edit in hertz

    Parameters    
    left_mute : minimum slowness to start mute zone
    right_mute : maximum slowness to end mute zone
    top_mute : minimum tau (time) to start mute zone
    base_mute : maximum tau (time) to start mute zone
    
    pxmax : maximum slowness in gather
    '''

    import numpy as np
    
    fs=kwargs['fs']
    pmute=kwargs['pmute']
    tmute=kwargs['taumute']

    # get the mute definition in slowness (p) and time (t)
    left_mute=pmute[0]
    right_mute=pmute[1]
    top_mute= tmute[0]
    base_mute=tmute[1]

    dt=1/fs

    # Find the indices of the mute start and end
    # 1. subtract mute value from every element of slowness array
    # 2. find the index where the difference is smallest
    slow1 = (np.abs(P_slopes-left_mute)).argmin()
    slow2 = (np.abs(P_slopes-right_mute)).argmin()
    
    tau1=int(top_mute//dt)
    tau2 = int(base_mute//dt)

    print (' slow1 index :',slow1, ' slow2 index :', slow2)

    xfilt_arb = np.copy(TauPspec)
    xfilt_arb[slow1:slow2,tau1:tau2]= 0

    invTauP_arbfilt = TauPop.inverse(xfilt_arb).reshape(TauPop.dimsd)

    return xfilt_arb, invTauP_arbfilt

def radon_plot(headers,Pslopes,*args,**kwargs):
    ''' Plot either 2 or 4 data sets
    1. input VSP data
    2. forward transform of 1.
    if 4 data sets are input:
    3. inverse transform of 2., usually after muting
    4. forward transform shown in 2. after muting
    
    Inputs:
    headers : VSP header array for receiver depths
    Pslopes : array of slowness (slopes) used in transform
    
    '''

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import numpy as np
    
    #return args, could be 2 or 4 data items
    print( 'elements in args :',len(args))
    
    VSP=args[0]
    TauPspec=args[1]
    
    # must maintain order of inputs to get this part correctly
    if len(args)==4:
        VSPfilt=args[2]
        TauPspec_filt=args[3]
        pmute=kwargs['pmute']
        tmute=kwargs['taumute']
        #create a box for mute limits
        xs = [pmute[0], pmute[1],pmute[1],pmute[0],pmute[0]]
        ys = [tmute[0],tmute[0],tmute[1],tmute[1],tmute[0]]        
    
    # check if filtered data exists so we can add an 2 axes
    try:
        VSPfilt
    except NameError:
        VSPfilt_exists = False
    else:
        VSPfilt_exists = True

    fs=kwargs['fs']
    #pmute=kwargs['pmute']
    #tmute=kwargs['taumute']
    pclip=kwargs['plot_scaler']

    numsamp, ntrace = VSP.T.shape            # data needs to be 1 trace/column    
    taxis = np.arange(0, numsamp*(1/fs),(1/fs) ) 
    xaxis = headers[:,1]
    

    
    fig = plt.figure(figsize=(12,8) )    
    gs = GridSpec(1, 2, hspace = .22, wspace=0.2) # make a row by col grid
   
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1],sharey=ax1)
     
    if VSPfilt_exists==True:
        ax1.remove()
        ax2.remove() 

        fig = plt.figure(figsize=(12,16) )    
        gs = GridSpec(2,2, hspace = .22, wspace=0.2)
    
        ax1 = plt.subplot(gs[0])    
        ax2 = plt.subplot(gs[1])         
        ax3 = plt.subplot(gs[2])    
        ax4 = plt.subplot(gs[3],sharey=ax3)

    ax1.imshow(
        VSP.T,
        cmap="gray",
        vmin=-pclip * np.abs(VSP).max(),
        vmax=pclip * np.abs(VSP).max(),
        extent=(xaxis[0], xaxis[-1], taxis[-1], taxis[0]),
    )
    ax1.set(xlabel="$x$ [m]", ylabel="$t$ [s]", title="Data")
    ax1.axis("tight")

    ax2.imshow(
        TauPspec.T,
        cmap="jet",
        vmin=-pclip * np.abs(TauPspec).max(),
        vmax=pclip * np.abs(TauPspec).max(),
        extent=(Pslopes[0], Pslopes[-1], taxis[-1], taxis[0]),
    )

    ax2.set(xlabel="$p$ [s/m]",  ylabel="$tau$ [s]", title="Radon")
    ax2.axis("tight")
    
    if VSPfilt_exists==True:
        ax3.imshow(
            TauPspec_filt.T,
            cmap="jet",
            vmin=-pclip * np.abs(TauPspec).max(),
            vmax=pclip * np.abs(TauPspec).max(),
            extent=(Pslopes[0], Pslopes[-1], taxis[-1], taxis[0]),
            )
        #ax3.axvline(left_mute, color="r", linestyle="--")
        #ax3.axvline(right_mute, color="r", linestyle="--")
        ax3.plot(xs,ys,color="r", linestyle="--")
        ax3.set(xlabel="$p$ [s/m]",  ylabel="$t$ [s]",title="Radon with Mute")
        ax3.axis("tight")

        ax4.imshow(
            VSPfilt.T,
            cmap="gray",
            vmin=-pclip * np.abs(VSP).max(),
            vmax=pclip * np.abs(VSP).max(),
            extent=(xaxis[0], xaxis[-1], taxis[-1], taxis[0]),
        )
        ax4.set(xlabel="$x$ [m]", ylabel="$tau$ [s]",title="Inverse transform")
        ax4.axis("tight")
    plt.show() 
   