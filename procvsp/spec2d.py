def fktran_numpy(seis,headers, fs, ishift, comment1, linevel):
    '''2D fft to generate fk plot
    seis - VSP data matrix
    headers - need receiver depths
    fs - sampling rate in hertz

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

    scalar = 1    
    DBmax = -60    
    fk_for_plot = fk    
    if ishift == 1:    
        fk_for_plot = fk_for_plot[:,plotnf:]        
        print (' fk_for_plot shape :',fk_for_plot.shape)
    
    Spec2d.fk_plot(fk_for_plot, freq_unwrap, k_unwrap, fs, rcvz, scalar, DBmax,\
            comment1,linevel)
    
    return fk, freq_unwrap, k_unwrap, numsamp,ntrace
    
def invfk_numpy(fk, numsamp, ntrace):    
    ''' inverse 2d fft to transform to depth-time
    '''        
    import numpy as np

    VSP_invtrans = np.fft.ifft2(np.fft.ifftshift(fk))
    
    print("\u0332".join('\ninv fk numpy info and parameters ') )    
    print (' fk shape :', fk.shape)
    print(' VSP_invtran shape :',VSP_invtrans.shape)    

    # get real component
    VSP_invtrans = np.real(VSP_invtrans[:ntrace,:numsamp])
    
    return VSP_invtrans
    
def fk_plot(FKdata, f, k, fs, rcvdepth,scal, dbmax,txt1,vel):

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.ticker as mtick

#    plt.rcParams.update({'font.size': 16}) # may affect other plots
    
    FKdata = FKdata.T                    # to get frequency on y axis    
    kmin, kmax, fmin, fmax = -.0075, .0075, 0, 80
#    kmin, kmax, fmin, fmax = k.min(),k.max(), f.min(), f.max()
    
    FKdata = np.absolute(FKdata)    # get magnitude of complex values of FKdata 
                                    # np. real could be tried to get real comp.
    FKdata = 20*np.log10((FKdata)/np.max(FKdata))  # Convert to db 20 * np.log10\
                                                   # (S / np.max(S))
    
    df =f[1] - f[0]                                # sample rate in hz    
    stopf = int(fmax/df)                           # get correct stop index    
    dk = k[1]-k[0]    
    startk = int(abs((k.min() - kmin)/dk))    
    stopk = startk + int(kmax/dk*2)

    print("\u0332".join('\nfk plot info and parameters') )    
    print ('fk shape :', FKdata.shape)       
    print(' df :',df,'\n stop freq index :', stopf, '\n f shape :', f.shape, \
          '\n dk :', dk, '\n startk :', startk, '\n stopk :',stopk)     
    
    FKdata =  FKdata[0:stopf,startk:stopk]

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
            vmin = dbmax,vmax = 0,extent = [kmin, kmax, fmax, fmin],\
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