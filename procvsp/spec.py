
def nextpow2(x):
    import numpy as np
 
    ''' find a number that is a power of 2
    from https://stackoverflow.com/questions/14267555/find-the-
    smallest-power-of-2-greater-than-n-in-python
    ''' 
    if x == 0:            
        y = 1
       
    else: y = np.ceil(np.log2(x))
        
    print("\u0332".join('\nnextpow2 info and parameters') )
    print(' x :', x)
    print(' nextpow2 :', y)    

    return int(y)
    
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

    import procvsp.utils as Utils
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds
    
    print("\u0332".join('\nFrAn Parameters :'))    
    print ('fs :', fs,)

    # extract analysis trace
    data_single, thead_single = Utils.chosetrace(data, thead, trace)    

    # useful headers
    TTobs_single = thead_single[:,8]
    zrcv_select = thead_single[:,2]
    trnum_single = thead_single[:,0]
    
    # generate the spectra for the chosen trace
    X, X_db,freq = spectra(data_single,timerange, frange,  fs,twin)

    X=X.T
    X_db=X_db.T

    ############   make Spectral plots   ################
    
    plt.figure(figsize=(15,5))    
    ax1 = plt.subplot(121)    

    ax1.plot(freq, np.absolute(X), c = 'red')  # using fftfreq to get x axis    

    ax1.set_title('Amplitude Spectrum of %s Depth %s'
                  %(title_spec, zrcv_select))    
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(frange[0], frange[1]) # extents must be set   
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.set_ylabel('Amplitude')    
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    ax2 = plt.subplot(122)
    
    ax2.plot(freq,X_db, c='blue') #using number of samples and sample rate to get x axis
    
    ax2.set_title('Amplitude Spectrum (db) of %s Depth %s'
                  %(title_spec, zrcv_select))    
    ax2.set_xlabel('Frequency hz')
    ax2.set_xlim(frange[0], frange[1]) # extents must be set   
    ax2.set_ylabel('dB')    
    ax2.xaxis.grid()    
    ax2.yaxis.grid()
    
    plt.show() 
    
def spectra(idata, timerange, frange, fs,twin):

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    

    ####### window the trace to prevent edge effect in fft  #############
    
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds

    ####### optimize number of samples for transform  #############

    # for a segment of trace
    if (twin =='y')or(twin=='Y'):
        start = int(timerange[0]*(fs/1000))
        stop = int(timerange[1]*(fs/1000))
        idata_trimd = idata[:,start:stop]
        
        # Apply window to segment
        w=sig.tukey(idata_trimd.shape[1], alpha=.1)    
        idata_win = idata_trimd[:,:]*w.reshape(-1, w.shape[0])        # Multiply trace by window
        N = scipy.fft.next_fast_len(int((timerange[1]-timerange[0])*(fs/1000)))
    
    # for a full trace
    else:
        # Apply window to whole trace
        w=sig.tukey(idata.shape[1], alpha=.1)    
        idata_win = idata[:,:]*w.reshape(-1, w.shape[0])        # Multiply trace by window    
        N = scipy.fft.next_fast_len(idata_win.shape[1]) # in samples, best for scipy fft
    
    # pad with zeros if optimal N greater than trace or segment
    #if (N > idata_win.shape[1]):                    
    #    pad = N - idata_win.shape[1]
    #    idata_win = np.pad(idata_win, ((0,0),(0,int(pad))), 'constant')        

    ####################### do the fft  #####################################
    X = np.zeros(shape = (idata_win.shape[0], N),dtype=np.complex_)
    X_db = np.zeros(shape = (idata_win.shape[0],N))
    X_pow = np.zeros(shape = (idata_win.shape[0],N))

    for k in range(0,(idata_win.shape[0])):
        
        X[k,:] = scipy.fft.fft(idata_win[k,:], n = N)      # from 0, to TT plus window
        # normalize db by max of first trace, which is usually the shallowest and has the strongest signal, 
        # but could be changed to max of each trace if desired
        X_db[k,:] = 20*np.log10(np.abs(X[k,:])/(np.abs(np.max(X[0,:])))) # db=20*np.log10(S/np.max(S))
        #X_pow[k,:] = np.abs(X[k,:])**2 # power spectrum
    # Only keep positive frequencies #########    
    freq = scipy.fft.fftfreq(N, d=samprate)    # Generate plot frequency axis
    keep = freq>=0
    freq = freq[keep]
    
    X_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    X_db_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    #X_pow_posfreq = np.zeros(shape = (idata_win.shape[0], freq.shape[0]))
    
    X_posfreq = X[:, keep]
    X_db_posfreq = X_db[:, keep]
    #X_pow_posfreq = X_pow[:, keep]
    
    return X_posfreq, X_db_posfreq,freq

def spec_FZ_SLB(idata, thead, fs,**kwargs):
    ''' Frequency Analysis of every trace (F-Z)    
        Uses scipy fft 
        Number of taps defined using next_fast_len, not next pow 2
        
        timerange - time window to extract from traces
        frange - plot frequency range
        dbrange - plot db range
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scalar applied to amplitude of spectral plot
    '''
    import numpy as np
    
    import scipy.signal as sig
    import scipy.fft    
    import scipy.ndimage
    
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib  import gridspec
    
    window = 'y' # put as arg
    
    db=kwargs['spec_type'] # 'y' for dB, 'n' for amplitude
    spacing=kwargs['spacing'] # Z for traces spread by receiver depth    
    scale=kwargs['scale'] # scale up image apmlitude plot
    title_spec=kwargs['title_fran']
    twin=kwargs['time_win'] # window trace prior to fft
    timerange=kwargs['time_range'] 
    frange=kwargs['freq_range']    
    trrange=kwargs['trace_range'] 
    dbrange=kwargs['db_range'] 
    skiplabel=kwargs['header_skip'] # plot every nth header
    rotate=kwargs['plt_rotate'] # rotate the plot 90 degrees if rotate =='y'
    save = kwargs['savepng']
    interp = kwargs['interp']
    TTobs = thead[:,8]
    rcvdepth = thead[:,2]
    trnum = thead[:,0]
    sdb=dbrange[1]
    edb=dbrange[0]
    
    idata=idata[trrange[0]:trrange[1],:]
    
    trindex  = np.stack([trnum for _ in range(idata.shape[0])], axis=1)
    rcvindex = rcvdepth[trrange[0]:trrange[1]]
    
    print("\u0332".join('\nFrAn image2 Parameters :'))    
    print (' data shape :',idata.shape, ' TTobs shape :',
           TTobs.shape)
    print (' trindex.min():',trindex.min(), ' trindex.max():',trindex.max())
    print ('fs :', fs,)
    # get the spectra, one for every trace
    X_posfreq,X_db_posfreq,freq=spectra(idata,timerange, frange, fs,twin)
    
    # set up the different plotting parameters for db or amplitude
    if (db=='y')or(db=='Y'):
        spec=X_db_posfreq.T
        vminp = edb
        vmaxp = sdb
        scale=1
        zlabel = 'dB'
        textstr =  'Spectra min : %s'%(spec.min())  
    else:
        spec=abs(X_posfreq.T)
        vminp = np.min(np.abs(spec))/scale
        vmaxp = np.max(np.abs(spec))/scale
        zlabel = 'Amplitude'
        textstr =  'Spectra max : %s'%(spec.max())  
    
    print (' freq.min():',freq.min(), ' freq.max():',freq.max())
    print (' trnum.shape :, trindex.shape :, Freq.shape :',trnum.shape,trindex.shape,freq.shape)
    print (' X_db_posfreq.shape :',X_db_posfreq.shape)
    print (' vmin :',vminp,' vmax :',vmaxp)

    ############   make amplitude-magnitude plots   ################
    
    fig = plt.figure(figsize=(10,8))
    
    gs = gridspec.GridSpec(1, 1, wspace = .25)
    
    ax1 = plt.subplot(gs[0])
    
    # set up axes and image extent
    extent = [rcvindex.min(), rcvindex.max(), freq.max(), freq.min()]
    x=rcvindex
    y=freq
    
    if (rotate=='y')or(rotate=='Y'):
        spec=spec.T
        extent = [ freq.min(), freq.max(),rcvindex.max(), rcvindex.min(),]
        y=rcvindex
        x=freq 

    # interpolate spectra if requested
    if (interp=='y')or(interp=='Y'):
        spec = scipy.ndimage.gaussian_filter(spec, 1)

    # plot the spectra as an image
    plot1 =  ax1.imshow(spec, cmap="gist_ncar_r", interpolation='none',aspect = 'auto',
               vmin = vminp,
               vmax = vmaxp,
               extent = extent)
    # plot the contours
    cp = ax1.contour(x,y,spec,20,colors='k',linestyles='solid',linewidths=0.4)   
    ax1.clabel(cp, inline=True, fontsize=12)   
    ax1.set_title('2D Amplitude Spectra from %s'%(title_spec))    

    if (rotate=='n')or(rotate=='N'):     
        ax1.set_ylim(frange[1], frange[0]) # extents must be set    
        ax1.set_xlim(rcvindex.min(), rcvindex.max())                           
        ax1.yaxis.grid()    
        ax1.set_xlabel('TV Depth',fontsize = 14)    
        ax1.set_ylabel('Frequency (hz)',fontsize = 14)    
        #plot all the tick marks         
        ax1.set_xticks(x[:-1:1])
        #force tick marks and labels on top and bottom of plot        
        ax1.tick_params(top='True',labeltop=True)
        #rotate the tick labels 90 degrees, and skip some labels    
        for n, label in enumerate(ax1.xaxis.get_ticklabels()):
            label.set_rotation(270)
            if n % skiplabel != 0:
                label.set_visible(False)        

    else:
        ax1.set_xlim(frange[0], frange[1]) # extents must be set    
        ax1.set_ylim(rcvindex.max(), rcvindex.min())              
        ax1.xaxis.grid()    
        ax1.set_ylabel('TV Depth',fontsize = 14)    
        ax1.set_xlabel('Frequency (hz)',fontsize = 14)       
        #plot all the tick marks
        x=x[::-1]        
        ax1.set_yticks(y[:-1:1])
        #force tick marks and labels on top and bottom of plot        
        ax1.tick_params(left='True',labelleft=True)
        #rotate the tick labels 90 degrees, and skip some labels    
        for n, label in enumerate(ax1.yaxis.get_ticklabels()):
            label.set_rotation(0)
            if n % skiplabel != 0:
                label.set_visible(False)                 

    #plot a colorbar    
    pad = 0.03    
    width = 0.02    
    pos = ax1.get_position()
    axcol = fig.add_axes([pos.xmax + pad, pos.ymin, width, \
                          0.9*(pos.ymax-pos.ymin) ])
    cb1 = fig.colorbar(plot1, cax = axcol, aspect = 40)#,format='%.0e')
    cb1.set_label(label=zlabel,size=12, )
 
    plt.show()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('graphics\\spec_fz_%s.png' 
        %(title_spec), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)
            
 
    return np.abs(X_posfreq.T), X_db_posfreq.T,freq,thead[trrange[0]:trrange[1],:]   


def spec_FZ(idata, timerange, thead, fs,spacing, dbrange,frange, trrange,scale, title_spec, twin):
    ''' Frequency Analysis of every trace (F-Z)    
        Uses scipy fft 
        Number of taps defined using next_fast_len, not next pow 2
        
        timerange - time window to extract from traces
        frange - plot frequency range
        dbrange - plot db range
        fs - sample rate in hertz
        scale - scalar applied to amplitude of spectral plot
    '''
    import numpy as np
    
    import scipy.signal as sig
    import scipy.fft    

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib  import gridspec
    
    window = 'y' # put as arg

    TTobs = thead[:,8]
    rcvdepth = thead[:,2]
    trnum = thead[:,0]
    sdb=dbrange[1]
    edb=dbrange[0]
    print (' idata.shape :',idata.shape)
    print (' trrange :', trrange)
    print (' trnum.min(), trnum.max() :',trnum.min(), trnum.max())
    print (' rcvdepth.min(), rcvdepth.max() :',rcvdepth.min(), rcvdepth.max())
    trnum=trnum[trrange[0]:trrange[1]]
    idata=idata[trrange[0]:trrange[1],:]
    
    trindex  = np.stack([trnum for _ in range(idata.shape[0])], axis=1)
    rcvindex = rcvdepth[trrange[0]:trrange[1]]
    
    print("\u0332".join('\nFrAn image2 Parameters :'))    
    print (' data shape :',idata.shape, ' TTobs shape :',
           TTobs.shape)
    print (' trindex.min():',trindex.min(), ' trindex.max():',trindex.max())
    print ('fs :', fs,)
    
    X_posfreq,X_db_posfreq,freq=spectra(idata,timerange, frange, fs,twin)

    ############   make amplitude-magnitude plots   ################
    
    fig = plt.figure(figsize=(12,14))
    
    gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,.1], wspace = .25)
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    plot1 =  ax1.imshow(np.abs(X_posfreq.T), cmap="gist_ncar_r", interpolation='none',aspect = 'auto',
               vmin = np.min(np.abs(X_posfreq))/scale,
               vmax = np.max(np.abs(X_posfreq))/scale,
               extent = [rcvindex.min(), rcvindex.max(), freq.max(), freq.min()])
    
    ax1.set_ylim(frange[1], frange[0]) # extents must be set    
#    ax1.set_xlim(trrange[0], trrange[1])
    ax1.set_xlim(rcvindex.min(), rcvindex.max())                                                      
    ax1.yaxis.grid()    
    ax1.set_xlabel('Receiver Depth')    
    ax1.set_ylabel('Frequency (hz)')    
    ax1.set_title('Combined Amplitude Spectra from %s'%(title_spec))    
    textstr =  'Spectra max : %s'%(X_posfreq.max())    
    ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top')
    
    #plot a colorbar    
    pad = 0.03    
    width = 0.02    
    pos = ax1.get_position()
    axcol = fig.add_axes([pos.xmax + pad, pos.ymin, width, \
                          0.9*(pos.ymax-pos.ymin) ])
    cb1 = fig.colorbar(plot1, label = 'Amplitude', cax = axcol, aspect = 40,format='%.0e')

    ############   make amplitude decibel plot   ################################

    plot2 = ax2.imshow(X_db_posfreq.T, cmap="gist_ncar_r", interpolation='none',aspect = 'auto',
               vmin = edb,
               vmax = sdb,
               extent = [rcvindex.min(), rcvindex.max(), freq.max(), freq.min()])

    ax2.set_ylim(frange[1], frange[0]) # extents must be set    
    ax2.set_xlim(rcvindex.min(), rcvindex.max())                                                      
    ax2.yaxis.grid()    
    ax2.set_xlabel('Receiver Depth')    
    ax2.set_ylabel('Frequency (hz)')    
    ax2.set_title('Combined Amplitude Spectra (dB) from %s'%(title_spec))    
    textstr =  'Spectra min : %s'%(X_db_posfreq.min())    
    ax2.text(0.05, 0.05, textstr, transform=ax2.transAxes, fontsize=12,
             verticalalignment='top')
    
    #plot a colorbar  
    pad = 0.03    
    width = 0.02    
    pos2 = ax2.get_position()
    axcol2 = fig.add_axes([pos2.xmax + pad, pos2.ymin, width, \
                          0.9*(pos2.ymax-pos2.ymin) ])
    cb2 = fig.colorbar(plot2, label = 'dB', cax = axcol2, aspect = 40,format='%.1f')
    
    # add a separate trace number track
    deci= int((trindex.max() - trindex.min())/20)
    trnum=np.int64(trnum)
    ytracnum=np.copy(trnum)*0
    
 
    ax3.plot(trnum[::deci],ytracnum[::deci],linestyle='None')
    for z,y in zip(trnum[::deci],ytracnum[::deci]):
        trlabel=z
        ax3.annotate(trlabel, (z,y))#, textcoords='data')
    ax3.set_xlim(trnum.min(), trnum.max())
    # turn off bounding box
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    #ax1.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.set_xticks([])
    ax3.set_yticks([])

    ax3.set_xlabel('Trace')
    
    plt.show()
 
    return np.abs(X_posfreq.T), X_db_posfreq.T,freq,thead[trrange[0]:trrange[1],:]   

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

def simple_bpf(VSPdata, tf1,tf2,tf3,tf4, fs, qc ):
#def simple_bpf(VSPdata, f1, f2, transL, transH, fs, qc ):
    
    '''from excellent tutorial at:
    https://tomroelandts.com/articles/how-to-create-simple-band-pass-and-band-reject-filters
    
    1. Create low pass and high pass windowed sync filters   
    2. Convolve the 2 filters to get a band pass filter
    3. Convolve the data trace with band pass filter. Numpy has an easy sinc generator
    
    Original method uses transtion zones:
    - transL is the transition band or the low pass filter - the high frequency roll-off
    - transH is the transition band or the high pass filter - the low frequency roll-off
    To call original method:
    decon_dwn_filt[:,k] = simple_bpf(decon_dwn_all[:,k], lc1, hc1,tz1,tz2,fs,qc='y')
    
    Updated method uses corner frequencies, from which the transition zones are calculated
    tf1 - low cut-off corner in Hertz
    tf2 - low pass corner in Hertz
    tf3 - low pass end corner in Hertz
    tf4 - high cut-off corner in Hertz
    

    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    VSPdata=VSPdata.T
    
    #print("\u0332".join('\nSimple BPF Information :'))
    #print (' VSPdata.shape :',VSPdata.shape)
    
    L=1024 # length of frequency response
    samprate = 1/fs            #sample rate in seconds
    
    dt = 1/fs
    '''  
    # This is the original method from T.Roelandts 
    # I adapted it to use corner frequencies instead of transition zones
    
    fL = f1/fs  # Cutoff frequency as a fraction of the sampling rate 
    fH = f2/fs  # Cutoff frequency as a fraction of the sampling rate 
    bL = transL/fs  # Transition band, as a fraction of the sampling rate 
    bH = transH/fs  # Transition band, as a fraction of the sampling rate 
    NL = int(np.ceil((4 / bH))) # samples in low pass
    NH = int(np.ceil((4 / bL))) # samples in high pass
    print (' fL :', fL,' fH :', fH)    
    print (' NL :', NL,' NH :', NH)
    ''' 
    # This is the adapted method using corner frequencies instead 
    # of transition zones
    
    fL=(((tf2-tf1)/2)+tf1)/fs
    fH=(((tf4-tf3)/2)+tf3)/fs
    bL = (tf2-tf1)/fs
    bH = (tf4-tf3)/fs
    NL=int(np.ceil((4/bH)))
    NH=int(np.ceil((4/bL)))    
    
    if not NL % 2: NL += 1  # Make sure that NL is odd.
    nL = np.arange(NL)
    if not NH % 2: NH += 1  # Make sure that NH is odd.
    nH = np.arange(NH)
 
    # Compute a low-pass filter with cutoff frequency fH.
    hlpf = np.sinc(2 * fH * (nL - (NL - 1) / 2))
    hlpf *= np.blackman(NL)
    hlpf = hlpf / np.sum(hlpf)
 
    # Compute a high-pass filter with cutoff frequency fL.
    hhpf = np.sinc(2 * fL * (nH - (NH - 1) / 2))
    hhpf *= np.blackman(NH)
    hhpf = hhpf / np.sum(hhpf)
    hhpf = -hhpf
    hhpf[(NH - 1) // 2] += 1
    
    # Convolve both filters.    
    h = np.convolve(hlpf, hhpf)

    # Pad filter with zeros.
    h_padded = np.zeros(L)
    h_padded[0:h.shape[0]] = h
    
    # do the fft
    H = np.abs(np.fft.fft(h_padded))
    freq = np.fft.fftfreq(H.shape[0], d=dt)    # Generate plot frequency axis
    ######### get rid of negative frequencies
    keep = freq>=0    
    H = H[keep]    
    freq = freq[keep]
    
    # Test if VSP data is a single trace - number of dimensions after squeeze
    # should be 1
    VSPdata=np.squeeze(VSPdata)
    numdims=VSPdata.ndim
    # print (' VSPdata.shape after squeeze :',VSPdata.shape)
    
    if numdims == 1:
        BPFdata = np.zeros(shape = (VSPdata.shape), dtype=np.float32)
        BPFdata = np.convolve(VSPdata,h, mode='same')
    else:  
        # apply filter to data
        BPFdata = np.zeros(shape = (VSPdata.shape[0], VSPdata.shape[1]), dtype=np.float32)          
        for k in range(0,(VSPdata.shape[1])):        
            #BPFdata[k,:-1] = np.convolve(VSPdata[k,:-1],h, mode='same')
            BPFdata[:,k] = np.convolve(VSPdata[:,k],h, mode='same')

    if (qc=='y')or(qc=='Y'):
    # Plot frequency response (in amplitude and dB) and impulse response
        fig = plt.figure(figsize=(15,5))    
        gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1], wspace = .25)
    
        ax1 = plt.subplot(gs[0])    
        ax1.plot(freq,H)      
        ax1.set_xlim(0, tf4*2)       # subjective choice, f2 for original method
        ax1.set_xlabel('Frequency (Hz)')    
        ax1.set_ylabel('Gain')    
        ax1.grid(True)    
        ax1.set_title('Simple Bandpass Frequency Response\n%s/%s - %s/%s hz corner frequencies'
                      %(tf1, tf2, tf3, tf4),fontsize=12)

        ax2 = plt.subplot(gs[1])    
        ax2.plot(freq, 20 * np.log10(H))      
        ax2.set_xlim(0, tf4*4)      # subjective choice, f2 for original method
        ax2.set_ylim(-200,0)       # subjective choice
        ax2.set_xlabel('Frequency (Hz)')    
        ax2.set_ylabel('Gain [dB]')    
        ax2.grid(True)        
        ax2.set_title('Simple Bandpass Frequency Response\n%s/%s - %s/%s hz corner frequencies'
                      %(tf1, tf2, tf3, tf4),fontsize=12)    
    
        ax3 = plt.subplot(gs[2])    
        x = np.arange((-h.shape[0]*dt)/2, (h.shape[0]*dt)/2, dt)
        ax3.plot(x,h, c='red')                                                 
        ax3.set_title('Simple Bandpass Impulse Response\n%s/%s - %s/%s hz corner frequencies'
                      %(tf1, tf2, tf3, tf4),fontsize=12)    
        ax3.set_xlabel('Time (s)')                
        ax3.grid(True)        
        
        plt.show()

    return BPFdata.T    
    
def FrAn_VSProwess(data, timerange, frange, thead, trace, fs, title_spec):
    ''' Frequency Analysis of section of trace 
        Shift the data prior to fft so direct arrival is center of input to fft
        Apply a window to reduce edge effects
        Use next power of 2 to get number of taps ( scipy does not require this)
        Finally run fft
        
        thead - requires a valid travel time     
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n' -not used
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter
    Example:    
    
    trace = 6
    time_win = 'n' # window trace prior to fft
    time_range = [0, 2000]
    freq_range = [0, 150]
    
    title_fran = 'Raw Z'
    FrAn_VSProwess(data_edit, time_range,freq_range, theader_edit, trace, fs, 
         title_fran)
    '''

    import math
    import scipy 
    from math import ceil
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import procvsp.utils as utilvsp
            
    dt =1/fs *1000             # sample rate in ms
    samprate = 1/fs            #sample rate in seconds
    
    # extract analysis trace
    data_single, thead_single = utilvsp.chosetrace(data, thead, trace)
    data_single = data_single.T
    # get useful headers
    TTobs_single = thead_single[:,8]
    zrcv_select = thead_single[:,2]
    trnum_single = thead_single[:,0]
        
    print("\u0332".join('\nFrAn2 Parameters :'))   
    print ('fs :', fs,)
    
    ###################### section the trace  #######################
    
    data_strim = data_single[int(timerange[0]*(fs/1000)):int(timerange[1]*(fs/1000))]
    
    ############# optimize number of samples for transform  ############

    win = nextpow2(int((timerange[1]-timerange[0])*(fs/1000)))  # get next power of 2 for window
    
    start = (int(ceil(TTobs_single) / dt))  # include data from start - 
                                            # direct arrival - plus a window 
    stop = int(start + win/2)               # Divide window by 2 for VSProwess  
    data_trim = data_single[0:stop]       # Apply window to whole trace so to
                                            # avoid tapering direct arrival
    N = data_trim.size
   
    ######### shift to put direct arrival in center of a window ##########
    
    tshift = thead_single[:,8] * (fs/1000) # get index number for travel time
    atime = int(win/2)                     # get index number to alignment at window center            
    xshift = tshift - atime                # shift to move direct arrival to center of window
    xshift = xshift.astype(int) * -1 
    arr=data_trim
    
    # apply shift to sectioned data
    if (xshift>0):
        arr = np.pad(data_trim, ((0,int(xshift)),(0,0)), 'constant')
    data_shft = np.zeros(shape = (arr.shape[0], arr.shape[1]))    
    data_shft = np.roll(arr,xshift) # carefulwith input array shape

    if xshift > 0:
        data_shft[xshift] = 0            # [1, 4000] need row number 0
        
    ####### window the trace to prevent edge effect in fft  #############

#    w = np.blackman(N) # design a window 
#    w = np.kaiser(N,14) # design a window 
#    w = np.hamming(N) # design a window 
    w = np.hanning(data_shft.shape[0]) # design a window
#    w = np.bartlett(N) # design a window
    
    data_wind = data_shft[:,0]*w        # Multiply trace by window
    print (' start :', start, ' stop :', stop, ' data windowed shape :', 
           data_wind.shape,' N :', N, ' w shape :', w.shape)

    ############################ do the fft  #############################
    
    X = scipy.fft.fft(data_wind)                      # number of taps equals input samples
    X_db = 20*np.log10(np.abs(X)/np.max(np.abs(X)))   # db=20*np.log10(S/np.max(S))    
    freq = scipy.fft.fftfreq(X.shape[0], d=samprate)  # Generate plot frequency axis 

    # Only keep positive frequencies    
    keep = freq>=0    
    X = X[keep]    
    X_db = X_db[keep]    
    freq = freq[keep]

    ############   make single trace plots   ################
    
    plt.figure(figsize=(15,15))    
    ax1 = plt.subplot(211)
    
    ax1.plot(freq, np.absolute(X), c = 'red')  # using fftfreq to get x axis    
    ax1.set_title('Amplitude Spectrum of %s at Depth %s'
                  %(title_spec, zrcv_select))    
    ax1.set_xlabel('Frequency hz')    
    ax1.set_xlim(frange[0], frange[1]) # extents must be set   
    #ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.xaxis.grid()    
    ax1.yaxis.grid()

    ax2 = plt.subplot(212)
    
    ax2.plot(freq,X_db, c='blue') #using number of samples and sample rate to get x axis    
    ax2.set_title('Amplitude Spectrum in db of %s at Depth %s'
                  %(title_spec, zrcv_select))    
    ax2.set_xlabel('Frequency hz')    
    ax2.set_xlim(frange[0], frange[1]) # extents must be set   
    ax2.xaxis.grid()    
    ax2.yaxis.grid()

    plt.show()     