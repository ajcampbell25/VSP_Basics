def medfilt (x, k):
    
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating zeros. Seems to be better edge results   
    see https://stackoverflow.com/questions/24585706/scipy-medfilt-wrong-result 

    """
    import scipy.signal as sig
    import numpy as np
   
    mfilt  = np.zeros(shape = (x.shape[0], x.shape[1]),dtype=np.float32)
    x = x.T    
    assert k % 2 == 1, "Median filter length must be odd."
    
    for b, trace in enumerate(x[:,:]):
        mfilt[:,b] = sig.medfilt(x[b], k) #scipy, in future test need for transpose
        
    return mfilt
    
def medfilt_across (x, k):
    
    """ from https://gist.github.com/bhawkins/3535131
        
    Apply a length-k median filter to a VSP data array.
    Boundaries are extended by repeating endpoints. at edge, input and output     
    traces are the same.     
    This can lead all zeros at edge traces after differencing.    
    Traces with all zeros lead to division errors when trying to normalize

    """
    import numpy as np

    out  = np.zeros(shape = (x.shape[0], x.shape[1]))
    x = x.T
    
    assert k % 2 == 1, "Median filter length must be odd."

    for b, trace in enumerate(x[:, :]):
        k2 = (k - 1) // 2                
        y = np.zeros ((len (x[b,]), k), dtype=x.dtype)
        y[:,k2] = x[b,]        
        row=b
        
        for i in range (k2):
            j = k2 - i
            y[j:,i] = x[row,:-j]
            y[:j,i] = x[row, 0]
            y[:-j,-(i+1)] = x[row,j:]
            y[-j:,-(i+1)] = x[row,-1]
            
        out[:,row] = np.median(y, axis = 1)
        
    return out   
    
def medfilt_along (VSP, k):
    
    """ from https://gist.github.com/bhawkins/3535131        
    Apply a length-k median filter to a VSP data array x.
    Boundaries are extended by repeating endpoints.

    """
 
    import numpy as np
    
    out = VSP
    
    assert k % 2 == 1, "Median filter length must be odd."

    for b, trace in enumerate(VSP[::1, :]):
        k2 = (k - 1) // 2                
        y = np.zeros ((len (VSP[b,]), k), dtype=VSP.dtype)        
        y[:,k2] = VSP[b,]        
        row=b

        for i in range (k2):
            j = k2 - i
            y[j:,i] = VSP[row,:-j]
            y[:j,i] = VSP[row, 0]
            y[:-j,-(i+1)] = VSP[row,j:]
            y[-j:,-(i+1)] = VSP[row,-1]
            
        out[row,] = np.median(y, axis = 1)

    return out
    
def normalize (Seismic, norm, thead, scal):
    """
    data normalization
    Frobenius norm is the default, L1 norm optional using data.sum
    row_sums = data.sum(axis=1)

    This method takes the whole trace, zeros before TT included so result 
    can be unexpected
    
    I need to create a method for normalizing in a window then apply 
    that norm factor to the whole trace
    stackoverflow.com/questions/8904694/how-to-normalize-a-2-dimensional-numpy-array-in-python-less-verbose
    """
    import numpy as np
    import matplotlib.pyplot as plt
     
#    print("\u0332".join('\nNormalization Stats :'))
    
    rdepth = thead[:,2]
    
    data2 = np.zeros(shape = (Seismic.shape[0], Seismic.shape[1]), dtype=np.float32)    
    data1 = Seismic
    
    if (norm == 'Y') or (norm =='y'):        
        row_sums = np.linalg.norm(data1, axis=1)        
        print (' row_sums shape', row_sums.shape)        
        data2 = (data1 / row_sums[:, np.newaxis])        
        datascaled = data2 * scal

        plt.figure(figsize=(15,7))    
        ax1 = plt.subplot(111)
        ax1.plot(rdepth, row_sums, c = 'red')  # using fftfreq to get x axis    
        ax1.set_title('Norm Factors')   
        ax1.set_xlabel('TVD Depth')    
        ax1.xaxis.grid()    
        ax1.yaxis.grid()
        plt.show()
        
    else:        
        datascaled = data1 * scal
        
                
    return datascaled
    
def normalize2(a, *, datarange = None):
    # not tested
    # from https://community.opengroup.org/osdu/platform/domain-data-mgmt-services
    # /seismic/open-zgy/-/blob/fff61f457d7d3ed90f239cafc5dc74a2d64400f8/python/openzgy/tools/viewzgy.py
    
    a = a.astype(np.float32)
    dead = np.isnan(a)
    amin, amax = datarange or (np.nanmin(a), np.nanmax(a))
    # Zero should be at the center
    if amin * amax < 0:
        x = max(abs(amin), abs(amax))
        amin, amax = (-x, x)
    # NaN and Inf show as smallest number
    a[dead] = amin
    if amin == amax:
        a *= 0
    else:
        # Avoid underflow, because app might have np.seterr(all='raise')
        a = a.astype(np.float64)
        a = (a - amin) / (amax - amin)
    a = (a * 255).astype(np.uint8)
    return a, dead

def cstack(upvsp, thead, fs, cwin, ctrnum, repeat):
    # Generate the a windowed data set and it's stack (corridor) 
    # A  tail mute is first applied at TWT plus cwin to get the corridor of data
    # The bottom ctrnum traces are left untouched
    # Headers must include twt in column 9
    
    import numpy as np
    import warnings

    cwinsamps= int(cwin *1000/fs)
    premtime = thead[:,8]*2*1000/fs
    mtime = (thead[:,8]*2*1000/fs)+cwinsamps
    
    # convert TWT and TWT + window to integer sample numbers 
    premtime = premtime.astype(int)
    mtime = mtime.astype(int)
    
#    corrmute = np.zeros(shape = (upvsp.shape[0], upvsp.shape[1]))#,dtype=np.float32) 
    corrmute = np.empty(shape = (upvsp.shape[0], upvsp.shape[1]))#,dtype=np.float32)
    corrmute = corrmute*np.nan
    
    emuted = corrmute.shape[0]-ctrnum
    snomute = corrmute.shape[0]-ctrnum
    enomute = corrmute.shape[0]
    
    for k in range(0,emuted):
        # copy data array up to break time+window into array of zeros   
        corrmute[k,premtime[k]:mtime[k]] = upvsp[k,premtime[k]:mtime[k]]
    for k in range(snomute,enomute):
        # copy data array from bottom traces   
        corrmute[k,premtime[k]:-1] = upvsp[k,premtime[k]:-1]

    # there will always be nans above the shallowest receiver, calculating the mean or median
    # of all nan samples will generate a warning which we can suppress
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
        corrstk= np.nanmedian(corrmute, axis=0).reshape(-1,1)
#        corrstk= np.nanmean(corrmute, axis=0).reshape(-1,1)
#        corrstk= np.sum(corrmute, axis=0).reshape(-1,1)
    
    print (' np.nanmax(corrmute) :', np.nanmax(corrmute), ' np.nanmax(corrstk) :', np.nanmax(corrstk))

    corrstk = np.repeat(corrstk,repeat,axis=1)

    return corrmute, corrstk.T  