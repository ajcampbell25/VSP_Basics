def specave(aspec,**kwargs):
    ''' calculate the mean of several spectra
    aspec : The input 2D spectra
    trace_1 : shallow trace to be averaged with surrounding traces
    trace_2 : deep trace to be averaged with surrounding traces
    
    return : 
    '''
    import numpy as np
    
    trace_1 = kwargs["tr_ref"]-1 # reference trace
    trace_2 = kwargs["tr_targ"]-1 # target trace
    avenum=kwargs["specave"]
    print_qc=kwargs['print_qc']

    if (print_qc=='y') or (print_qc=='Y'):
        print("\u0332".join('\Average Spec. Rat. Parameters :'))    
        
    # get the trace range around the one to be averaged 
    deltnum=avenum//2
    if print_qc=='y':
        print ('avenum :',avenum,'deltnum :',deltnum)
    # for one trace, need different slicing notation

    w1start=trace_1-deltnum
    w1end=trace_1+deltnum
    w2start=trace_2-deltnum
    w2end=trace_2+deltnum
    
    if w1start<0:
        w1start=0

    if avenum==1:
        A1av=aspec[:,trace_1]    
        A2av=aspec[:,trace_2]    
    
    else:
        A1av=np.mean(aspec[:,w1start:w1end],axis=1)    
        A2av=np.mean(aspec[:,w2start:w2end],axis=1)

    if (print_qc=='y') or (print_qc=='Y'):
        print ('w1start,w1end :',w1start,w1end)
        print ('w2start,w2end :', w2start,w2end)
        print ('A1av.shape :',A1av.shape,' A2av.shape :',A2av.shape)
    
    return aspec[:,w1start:w1end],aspec[:,w2start:w2end],A1av,A2av
    
    
def q_simulation(VSP, timerange, frange, thead, trace, fs,twin, title_spec,**kwargs):
    ''' Run a Frequency Analysis on a selected trace and then apply a 
    'forward Q filter' to the spectrum to simulate the effect of Q
        
        timerange - desired analysis window 
        twin - apply analyis window 'y' or use whole trace 'n'
        trace - trace number to extract from 2D data array
        fs - sample rate in hertz
        scale - scales amplitude of trace plot, >1 makes plot hotter
        
        The frequency spectrum is complex, the absolute value (real) is plotted

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    import numpy as np
    import scipy.signal as sig
    import scipy.fft    
    import math 
    from math import ceil
    import procvsp.utils as Utils
    import procvsp.spec as specvsp
    
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

    # generate the spectra for the chosen trace
    X, X_db,freq = specvsp.spectra(data_single.T,timerange, frange, fs,twin)   
    X=X.T    
    X=np.reshape(X,-1) # make array 1d
 
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
    
    fig=plt.figure(figsize=(14,5))    
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

def qest(thead,aspec,freq,fs, **kwargs):
    ''' Run a Frequency Analysis on two selected traces.
    Estimate the spectral ratios in a frequency band. From the
    slope of the spectral ratios, calculate Q
    
    theadspec : headers corresponding to spectral traces
    aspec : spectral traces F-Z
    tr_ref : reference or shallow trace
    tr_targ : target or deep trace
    specave : average n spectra around the reference and target
    srfreq_min : min frequency for spectral ratio calc
    srfreq_max : max frequency for spectral ratio calc

    fs - sample rate in hertz
        
    The frequency spectrum is magnitude

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
    #tmin = kwargs["time_min"]
    #tmax = kwargs["time_max"]
    fmin = kwargs["freq_min"]
    fmax= kwargs["freq_max"]
    srfmin = kwargs["srfreq_min"]
    srfmax= kwargs["srfreq_max"]    
    trace_1=kwargs["tr_ref"]-1 # reference trace
    trace_2=kwargs["tr_targ"]-1 # target trace
    save=kwargs["savepng"]
    
    print("\u0332".join('\nFrAn Parameters :'))    
    print ('fs :', fs,)

    # extract analysis traces
    thead_ref = thead[trace_1,:]
    thead_two = thead[trace_2,:]
    
    zrcv_ref = thead_ref[2,]
    zrcv_two = thead_two[2,]
    
    ############create reference and target spectra #############################
    
    # calculate the averaged spectral traces
    in_ref,in_two,X_ref,X_two=specave(aspec,**kwargs)
    # find indices on the frequency vector for the reference frequency limits
    srfmax_ind=np.where((freq>=srfmin)&(freq<=srfmax))

    #############  Calculate spectral ratios, Q, and slope of ratios ###########
    
    Qeff,specrat,slopew,plot_slope=specratio(X_ref,X_two,freq,thead_ref,thead_two,**kwargs)

    ############   make Q spectra and spectral ratio plots   ################
    
    fig=plt.figure(figsize=(12,5))    
    ax1 = plt.subplot(121)    

    ax1.plot(freq, X_ref, c = 'red', label = "Reference Spectrum")  # using fftfreq to get x axis    
    ax1.plot(freq, X_two, c = 'blue', label = "Target Depth Spectrum")  # using fftfreq to get x axis    
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        

    ax1.set_title('Spectra at Reference Depth %s \nand at Target Depth %s'
                  %(zrcv_ref,zrcv_two))
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
    ax2.set_title('Spectral Ratio from Reference Depth %s \nto %s Depth'
                  %(zrcv_ref,zrcv_two))
    ax2.text(.1,.05,"Slope = %.4f, Q calculated = %.0f"%(slopew[0],Qeff),
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


def specratio(X_ref,X_two,freq,thead_ref,thead_two,**kwargs):
    ''' calculate the spectral ratios for Q estimation
    X_ref : magnitude spectrum at reference depth
    X_two : magnitude spectrum at target depth
    zrcv_ref : Depth to reference receiver
    zrcv_two : Depth to target receiver
    
    '''   
    import numpy as np
    
    fmin = kwargs["freq_min"]
    fmax= kwargs["freq_max"]
    srfmin = kwargs["srfreq_min"]
    srfmax= kwargs["srfreq_max"]    
    print_qc=kwargs['print_qc']

    if (print_qc=='y') or (print_qc=='Y'):
        print("\u0332".join('\Spectral Ratio Parameters :'))    

    # useful headers
    TTobs_ref = thead_ref[8,]
    zrcv_ref = thead_ref[2,]
    trnum_ref = thead_ref[0,]
    
    TTobs_two = thead_two[8,]
    zrcv_two = thead_two[2,]
    trnum_two = thead_two[0,]
    
    ############ get time and depth between receivers  ##########################

    Z0 = zrcv_ref
    Z = zrcv_two

    T=TTobs_two/1000
    T0=TTobs_ref/1000    

    ########### use spectral ratios to get Q from Q filtered spectrum ############

    #pecrat=np.log(np.absolute(X_two)/np.absolute(X_ref))
    specrat=np.log(X_two/X_ref)  
    
    # get sample number for max frequency to control plots
    # fmax_ind=np.where((freq>=fmin)&(freq<=fmax))
    srfmax_ind=np.where((freq>=srfmin)&(freq<=srfmax))

    slopew = np.polyfit(freq[srfmax_ind],specrat[srfmax_ind],  1)# 1 means first order - linear regression
    plot_slope = slopew[0]*freq[srfmax_ind]+slopew[1] 

    Qeff=-1*(np.pi*((T-T0)/slopew[0]))
    
    # trendpoly = np.poly1d(slopew) # an alternative for plotting instead of y=mx+b 
    
    if (print_qc=='y') or (print_qc=='Y'):
        print (' slopew :', slopew, ' slopew.shape :',slopew.shape, 'Qeff :',Qeff)
    
    return Qeff,specrat,slopew, plot_slope

def qave_specrat(thead,aspec, freq,fs,**kwargs):
    ''' Use the spectral ratio method to calculate Q
    Adapted from Crewes matlab routines

    Calculate an average Q log:
        1. Set a reference trace
    2. Calculate Q between reference trace and every deeper trace
        
    '''   
    import numpy as np
    import procvsp.utils as utilvsp
    import procvsp.qtools as qvsp
    import plotvsp.qplots as qplts
    
    usepow=kwargs['use_power']
    print_qc=kwargs['print_qc']

    ########## find average q from a reference shallow trace ######
    # - to all the other traces

    # create a power spectrum
    if (usepow=='y') or (usepow=='Y'):
        aspec = np.abs(aspec)**2
        
    #if (print_qc=='y') or (print_qc=='Y'):
    print (' thead.shape :',thead.shape,' aspec.shape :',aspec.shape) 
    
    # change reference trace to python indexing    
    tr_ref=kwargs["tr_ref"]-1
    trref_plot=kwargs["tr_ref"]
    
    numlevs=aspec.shape[1]-tr_ref
                      
    qsrave = np.zeros(shape=(numlevs))
    fdempty = np.zeros(shape=(numlevs))# for plotting
    trnumsr = np.zeros(shape=(numlevs))    

    print (' tr_ref :',tr_ref,'qsrave.shape :',qsrave.shape) 
    print (' range :',range(0,numlevs))
    for i in range(0,numlevs-1):       
        # target trace does not need to be re-indexed
        tr_target=tr_ref+i
        # reset dictionary value of target trace number
        kwargs['tr_targ'] = tr_target
        # for plotting trace number, needs to start at 1
        trnumsr[i] = i+tr_ref
                
        # extract reference trace headers
        theadref=thead[tr_ref,:]
        theadtwo=thead[tr_target,:]
        
        # calculate the averaged spectral traces
        in_ref,in_two,Xref_sr,Xtwo_sr=qvsp.specave(aspec,**kwargs)
        
        if (print_qc=='y') or (print_qc=='Y'):
            print (' theadref :',theadref)
            print (' theadtarget :',theadtwo)         
            
        # calculate Q from spectral ratios between reference and target traces
        #qsrave[i],specrat[i],slopew[i],plot_slope[i]=qvsp.specratio(Xref_sr,Xtwo_sr,freq,theadref,theadtwo,**kwargs)
        qsrave[i],_,_,_=qvsp.specratio(Xref_sr,Xtwo_sr,freq,theadref,theadtwo,**kwargs)
        if (print_qc=='y') or (print_qc=='Y'):
            print (' tr_ref :',tr_ref,' tr_target :',kwargs['tr_targ'],
                   ' trnumsr[i] :',trnumsr[i],' qsrave :',qsrave[i]) 
        print (' trnumsr[i] :',trnumsr[i],' numlevs :',numlevs,' i :',i,
               ' tr_ref :',tr_ref,' tr_target :',tr_target,' qsrave[i] :',qsrave[i])
#    trnumsr = trnumsr+thead[0,0]# make receiver number equal to the spectra receiver number
    kwargs['method']='sr'# put dominant frequencies on plot

    ptitle = 'Average Q From Spectral Ratios and Spectra - Reference Trace %s'%(tr_ref)
    qplts.specfz_qvals(thead,aspec, fdempty, fdempty,qsrave,trnumsr,freq,fs,ptitle,**kwargs)
    
    
    return  qsrave,trnumsr 

def q_centroid(thead_ref,thead_two,aspec, freq,fs,**kwargs):
    '''
    This version takes a reference and a target trace as  input
    for analysis
    
    Adapted to python by A. Campbell, from original CREWES matlab code.
    
    function [Q,obj,fd1,fd2]=q_centroid(A1,A2,f,t1,t2,f1,f2,Qmax,Qmin,p)
    % Q_centroid: estimate Q from matching centroid frequecies
    %
    % [Q,obj,fd1,fd2]=Q_centroid(A1,A2,f,t1,t2,f1,f2,Qmax,Qmin,p)
    %
    % It is assumed that A2 is related to A1 through the model
    % A2=A1*T*exp(-pi*f*(t2-t1)/Q), where T and Q are real, positive, scalars
    % to be determined. T represents frequency independent loss while Q
    % parameterizes the frequency-dependent attenuation. The centroid
    % frequency, also known as the dominant frequency is defined by
    % fc=sum(f.*A.^p)/sum(A.^p) for spectrum A where p defaults to 2. Let fc2
    % be the centroid frequency for spectrum A2 and consider the spectrum
    % Aq=A1*exp(-pi*f*(t2-t1)/Q) which has centroid frequency fc. By letting Q
    % range over all integer valies from Qmin to Qmax, we can find the
    % minimimum value of Qtest.*(fc-fc2).^2 where Qtest=Qmin:Qmax. The value of
    % Q at the minimum is returned as the estimated Q. A virtue of this method
    % is that is does not requite estimation of T. The method fails if fc1 is
    % less than fc2.
    %
    % A1 ... amplitude spectrum at time t1 (A1 must never vanish)
    % A2 ... amplitude spectrum at time t2 (Note t2 should be greater than t1)
    % f  ... frequency coordinate vector for A1 and A2. 
    %NOTE: A1, A2, and f must be real-valued vectors of identical size
    % t1 ... time at which the spectrum A1 was measured
    % t2 ... time at which the spectrum t2 was measured
    % Qmax ... maximum value of Q in the search
    % *********** default Qmax=250 ************
    % Qmin ... minimum value of Q in the search
    % *********** default Qmin=5 ***********
    % p ... exponent for dominant frequency calclation
    % ************* default p=2 *********
    %
    % Q       ... estimated Q value
    % obj     ... objective function
    % fc1     ... dominant frequency of spectrum A1
    % fc2     ... dominant frequency of spectrum A2
    %NOTE: obj has length Qmin:Qmax  
    %
    % by G.F. Margrave, 2014
    %
    % NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
    % with the exception of re-selling or re-distributing the SOFTWARE.
    % By using this software, you are agreeing to the terms detailed in this software's
    % Matlab source file.

    % BEGIN TERMS OF USE LICENSE
    %
    % This SOFTWARE is maintained by the CREWES Project at the Department
    % of Geology and Geophysics of the University of Calgary, Calgary,
    % Alberta, Canada.  The copyright and ownership is jointly held by
    % its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
    % project may be contacted via email at:  crewesinfo@crewes.org
    %
    % The term 'SOFTWARE' refers to the Matlab source code, translations to
    % any other computer language, or object code
    %
    % Terms of use of this SOFTWARE
    %
    % 1) This SOFTWARE may be used by any individual or corporation for any purpose
    %    with the exception of re-selling or re-distributing the SOFTWARE.
    %
    % 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
    %    presentations
    %
    % 3) This SOFTWARE is provided "as is" with no warranty of any kind
    %    either expressed or implied. CREWES makes no warranties or representation
    %    as to its accuracy, completeness, or fitness for any purpose. CREWES
    %    is under no obligation to provide support of any kind for this SOFTWARE.
    %
    % 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
    %    notice. New versions will be made available at www.crewes.org .
    %
    % 5) Use this SOFTWARE at your own risk.
    %
    % END TERMS OF USE LICENSE
    '''   
    import numpy as np
    from math import exp,pi
    
    #twin = kwargs["apply_window"]
    fmin = kwargs["freq_min"]
    fmax= kwargs["freq_max"]
    f1 = kwargs["qcfreq_min"]
    f2= kwargs["qcfreq_max"]    
    trace_1=kwargs["tr_ref"]-1
    trace_2=kwargs["tr_targ"]

    Qmin=kwargs["Q_min"]
    Qmax=kwargs["Q_max"]
    p=kwargs['p']
    save=kwargs["savepng"]
    print_qc=kwargs['print_qc']

    # get time headers
    t1 = thead_ref[8,]/1000
    t2 = thead_two[8,]/1000

    # time separation between 2 levels
    delt=(t2-t1) # assuming OWT
    
    if (print_qc=='y')or(print_qc=='Y'):
        print("\u0332".join('\nq_centroid Parameters :'))    
        print ('fs :', fs,)
        print (' aspec.shape :',aspec.shape)
        print (' trace_1 :',trace_1,' trace_2 :',trace_2)
        

    A1=aspec[:,trace_1]
    A2=aspec[:,trace_2]    
    #print ('A1.shape :',A1.shape,'A2.shape :',A2.shape)

    if len(A1)!=len(A2):
        print('A1 and A2 must be the same size')
    if len(freq)!=len(A1):
        print('size of f must be the same as A1 and A2')
            
    # create a testing Q-value vector
    Qtest=np.arange(Qmin,Qmax+1)# step default is 1

    # find the indices for the desired frequency range
    ind=np.where((freq>=f1)&(freq<=f2))
    
    # calculate the dominant frequencies and variances from the measured data
    fd1= np.sum(freq[ind]*A1[ind]**p)/np.sum(A1[ind]**p)
    fd2=np.sum(freq[ind]*A2[ind]**p)/np.sum(A2[ind]**p)

    sig1=np.sum(((freq[ind]-fd1)**2)*(A1[ind]**p))/np.sum(A1[ind]**p)
    sig2=np.sum(((freq[ind]-fd2)**2)*(A2[ind]**p))/np.sum(A2[ind]**p)

    fdtest=np.zeros(shape=Qtest.shape[0])
    sigtest=np.zeros(shape=Qtest.shape[0])

    if(fd1>fd2): 
        for k in range(0,len(Qtest)-1):
            #print (' k :',k,)
            #print (' A1.shape :',A1.shape,' freq.shape : ',freq.shape,' Qtest.shape :',Qtest.shape)
            # test spectrum for different q values            
            Atest=A1[ind]*np.exp(-1*np.pi*freq[ind]*(delt/Qtest[k]))
            # get centroid values for different q values from the test spectra 
            fdtest[k]=np.sum(freq[ind]*Atest**p)/np.sum(Atest**p)
            # calculate spectral spread-variance for each of the test spectra            
            sigtest[k]=np.sum(((freq[ind]-fdtest[k])**2)*(Atest**p))/np.sum(Atest**p)
            
        # find Q value where:
        # a)where fdtest is closest to fd2
        # b)where sigtest (spread) is closest to sig2
        #
        # 1. Square the difference of the test centroid array and the measured centroid
        # 2. Multiply by q value testing array
        obj1=Qtest*(fdtest-fd2)**2 # objective function for fdom
        # 1. Square the difference of the test variance array and the measured variance
        # 2. Multiply by q value testing array        
        obj2=Qtest*(sigtest-sig2)**2 # objective function for sigma
        # normalize the 2 objective functions and multiply them together
        # smaller numbers where modeled values match measured values
        obj=obj1/np.max(obj1)+1*obj2/np.max(obj2)
        # find the index where objective function is minimized
        [omin,imin]=obj.min(0),obj.argmin(0)
        # match the index number to the Q test array to get Q
        Q=Qtest[imin]
    else:
        Q=float('inf')
        obj=[]
    return Q, fd1, fd2        

def qave_centroid(thead,aspec, freq,fs,**kwargs):
    ''' Use the centroid method to calculate Q
    Adapted from Crewes matlab routines
    
    
    '''   
    import numpy as np
    import procvsp.utils as utilvsp
    import procvsp.qtools as qtls
    import plotvsp.qplots as qplts
    
    usepow=kwargs['use_power']
    print_qc=kwargs['print_qc']

    qave = np.zeros(shape=aspec.shape[1])
    qint = np.zeros(shape=aspec.shape[1])
    #q_int = np.zeros(shape=aspec.shape[1])
    ########## find average q from a reference shallow trace ######
    # - to all the other traces
   
    fdomrefa = np.zeros(shape=aspec.shape[1])
    fdoma = np.zeros(shape=aspec.shape[1])
    trnuma = np.zeros(shape=aspec.shape[1])    

    # create a power spectrum
    if (usepow=='y') or (usepow=='Y'):
        aspec = np.abs(aspec)**2
        
    if (print_qc=='y') or (print_qc=='Y'):
        print (' thead.shape :',thead.shape,' aspec.shape :',aspec.shape) 
    
    # change reference trace to python indexing    
    tr_ref=kwargs["tr_ref"]-1
    trref_plot=kwargs["tr_ref"]

    for i in range(0,aspec.shape[1]):       
        # target trace does not need to be re-indexed
        tr_target=i
        # reset dictionary value of target trace number
        kwargs['tr_targ'] = tr_target
        # for plotting trace number, needs to start at 1
        trnuma[i] = i+1 
        # extract reference trace
        theadref=thead[tr_ref,:]
        theadtwo=thead[tr_target,:]
        
        if (print_qc=='y') or (print_qc=='Y'):
            print (' theadref :',theadref)
            print (' theadtarget :',theadtwo)         
            
        # calculate Q from dominant frequency of reference and target traces
        qave[i], fdomrefa[i], fdoma[i]=qtls.q_centroid(theadref,theadtwo,aspec,freq,fs,**kwargs)
        
        if (print_qc=='y') or (print_qc=='Y'):
            print (' tr_ref :',tr_ref,' tr_target :',kwargs['tr_targ'],' qave :',qave[i], 
            ' fdomrefa [i] :', fdomrefa[i],' fdoma[i] :', fdoma[i])
    
    # make a plot of FZ with Q as a header graph
    trnuma = trnuma+thead[0,0]# make receiver number equal to the spectra receiver number
    kwargs['method']='cent'# put dominant frequencies on plot   
    ptitle = 'Average Q and Spectra with Dominant Frequencies Overlay - Reference Trace %s'%(tr_ref)
    qplts.specfz_qvals(thead,aspec, fdomrefa, fdoma,qave,trnuma,freq,fs,ptitle,**kwargs)
    
    # make a plot of 1D spectra with dominant frequencies a shaded analysis window 
    trwin=fdoma.shape[0] # make the window for all traces where qave was calculated
    kwargs['tr_ref'] = trref_plot
    kwargs['tr_win']= trwin-1
    
    ptitle = ' Spectra with Dominant Frequencies  - Window Length %s traces'%(trwin) 
    qplts.gauss_demo(aspec, freq, fdomrefa, fdoma,fs,ptitle,**kwargs)
    
    return  qave,trnuma 

def qint_centroid(thead,aspec, freq,fs,**kwargs):
    ''' Use the centroid method to calculate Q. \
    Interval is defined as a number of traces, not geological
    
    Adapted from Crewes matlab routines
    
    This version uses an overlapping sliding window approach
    to get interval Q
    aspec : F-Z spectral matrix
    thead : trace headers synchronized with aspec
    trref : reference trace, calculated at start of each sliding window
    trwin : sliding window in traces
    
    returns : qint is the interval Q over a trace interval
       
    '''   
    import numpy as np
    import procvsp.utils as utilvsp
    import procvsp.qtools as qtls
    import plotvsp.qplots as qplts
    
    usepow=kwargs['use_power']
    print_qc=kwargs['print_qc']
    trwin=kwargs['tr_win']
    
    qint = np.zeros(shape=aspec.shape[1])
    
    # create a power spectrum
    if (usepow=='y') or (usepow=='Y'):
        aspec = np.abs(aspec)**2
        
    if (print_qc=='y') or (print_qc=='Y'):
        print (' thead.shape :',thead.shape,' aspec.shape :',aspec.shape) 
    
    ######### find interval q between sliding windows of traces. ####### 
    # The windows overlap by one receiver

    # figure out the array shapes based on the for loop parameters
    rlist=list(range(0,aspec.shape[1]-trwin,1))
    shape= len(rlist)

    fdomref = np.zeros(shape=shape)
    fdom = np.zeros(shape=shape)
    trnum=np.zeros(shape=shape)
    qint=np.zeros(shape=shape)
    q_int=np.zeros(shape=shape)

    if (print_qc=='y') or (print_qc=='Y'):
        print (' fdomref.shape :',fdomref.shape,' trnum.shape :',trnum.shape,' shape :',shape) 
    level=0
    for i in range(0,aspec.shape[1]-trwin,1):
        tr_ref=i      
        tr_target=tr_ref+trwin-1
        # add 1 to tr_ref as 1 gets subtracted in q_centroid
        kwargs['tr_ref'] = tr_ref+1
        kwargs['tr_targ'] = tr_target
        
        trnum[level]=i+(trwin+1)/2 #location of Q value, needs to be corrected to midpoint? 
        
        theadref=thead[tr_ref,:]
        theadtwo=thead[tr_target,:]
        
        qint[level], fdomref[level], fdom[level]=qtls.q_centroid(theadref,theadtwo,aspec,freq,fs,**kwargs)
        if (print_qc=='y') or (print_qc=='Y'):
            print (' level :',level,' tr_ref :',tr_ref,' tr_target :',kwargs['tr_targ'], ' trnum :',trnum[level],' qint :',qint[level], 
                   ' fdomref [i] :', fdomref[level],' fdom[i] :', fdom[level])
        level=level+1
    
    # make a plot of FZ with Q as a header graph
    trnum=trnum+thead[0,0] # make receiver number relative to spectrum receiver number    
    kwargs['method']='cent'# put dominant frequencies on plot
    ptitle = 'Interval Q and Spectra with Dominant Frequencies Overlay - Window Length %s traces'%(trwin)
    qplts.specfz_qvals(thead,aspec, fdomref, fdom,qint,trnum,freq,fs,ptitle,**kwargs)
    
    return  qint,trnum

def ave2int(qave_xu,theader,**kwargs):
    ''' see Xu, Stewart,2006, Seismic Attenuation (Q) Estimation from VSP Data, 
        CSEG Expanded abstract
        inputs:
        qave : average or cumulative Q
        
        returns:
        q_int_xu : interval Q calculated from qave

        Example Useage:

        ave2int_params={
             "q_min": 0,"q_max":300, # for spectral plot limits            
             "tr_ref": 1,# reference trace
             "plot_title":"Average and Interval Q",
             "savepng":'n',
               }
        q_int=ave2int(qave,thead_spec,**ave2int_params)
    '''        
    import numpy as np
    import plotvsp.qplots as qplts

    q_int_xu=np.zeros((len(qave_xu),))    
    q_int_xu[0]=qave_xu[0]
    
    deltat=np.zeros((len(qave_xu),))
    twt = theader[:,8]#/1000
    
    ####### fix TWT to be measured from the reference trace
    trref=kwargs['tr_ref']-1
    twtr=twt[trref:len(qave_xu)+1]
    print ( ' twtr.shape :',twtr.shape)
    
    #deltat[1:]=twt[1:] - twt[:-1]

    denom = np.zeros((len(qave_xu),))
    qnum = np.zeros((len(qave_xu),))
    # vectorized method
    #q_int_xu[1:] = (twt[1:] - twt[:-1])/((twt[1:]/qave_xu[1:]) - (twt[:-1]/qave_xu[:-1]))
    
    #for loop as demo of method
    for i in range(1, qave_xu.shape[0]):
        #q_int_xu[i] = (twt[i] - twt[i-1])/(twt[i]/qave[i]  - twt[i-1]/qave[i-1] )
        q_int_xu[i] = (twtr[i] - twtr[i-1])*(twtr[i]/qave_xu[i]  - twtr[i-1]/qave_xu[i-1] )**-1
        deltat[i] = twtr[i] - twtr[i-1]
        denom[i] =  twtr[i]/qave_xu[i]  - twtr[i-1]/qave_xu[i-1]
        qnum[i] = i
        print ('qnum[i] : ',qnum[i], 'twtr[i] :',twtr[i],' i :',i,' qave[i] :',qave_xu[i],'qave[i-1] :',qave_xu[i-1],
               ' q_int[i] :',q_int_xu[i])
    #print (qnum.shape, twtr.shape,qave_xu.shape, q_int.shape)

    # replace negative (unphysical) values with pi
    neg = q_int_xu < 0
    if np.any(neg):
        q_int_xu[neg] = np.pi

    from tabulate import tabulate        
    #thead_tabl = np.vstack((qnum,twtr,deltat,denom,qave, q_int )).T        
    thead_tabl = np.vstack((qnum,deltat,denom,qave_xu, q_int_xu )).T        

    cheaders2 = ["Number","twt\nsec","t2-t1", "denominator","q ave","q int calc"]                       
    #numfmt = (".0f",".4f",".4f",".4f",".2f", ".2f")
    numfmt = (".0f",".4f",".4f",".2f", ".2f")                 
                 
    table = tabulate(thead_tabl, headers = cheaders2,  floatfmt = numfmt)#,tablefmt="pretty")
    print(table)
    
    qplts.plotqint(qave_xu, q_int_xu,theader,**kwargs)
    
    return q_int_xu     