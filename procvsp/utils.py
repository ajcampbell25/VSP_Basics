def depthlimit(VSPdata, headerfile, first, last):
    
    first = first - 1  #python indexing starts at 0, 
    last = last -1

    datanew = VSPdata[first:last,]    
    headernew = headerfile[first:last,]
    
    print("\u0332".join('\nDepthlimit Stats :'))    
    print(' data shape : ', VSPdata.shape, ' data dtype : ', VSPdata.dtype)    
    print (' headers shape :', headerfile.shape)    
    print (' first :', first, ' last : ', last )    
    print (' headers new shape :', headernew.shape)        
    print(' data new shape : ', datanew.shape)
    
    return datanew, headernew

def chosetrace(VSP, thead, num):
    
    theadnew = thead[num:num+1,]
    datanew = VSP[num:num+1,]
    
    return datanew, theadnew
    
def diff(align1, align2):
    import numpy as np
    
    difference  = np.zeros(shape = (align1.shape[0], align1.shape[1]),dtype=np.float32)    
    difference = np.subtract(align1, align2)

    return difference
    
def mute(arr, fs,thead, align):
    '''
    Methods from https://stackoverflow.com/questions/30399534/
    Shift-elements-in-a-numpy-array
    
    fix to get sample index instead of time!
    '''
    import numpy as np
    
    if align == 'owt':        
        mute_time = thead[:,8]*1000/fs 
        
    elif align == 'twt':
        mute_time = thead[:,8]*2*1000/fs
    
    arr_mute = np.zeros(shape = (arr.shape[0], arr.shape[1]),dtype=np.float32) 
    mtime = mute_time.astype(int)
   
    for k in range(0,(arr_mute.shape[0])):        
        arr_mute[k,mtime[k]:-1] = arr[k,mtime[k]:-1]
            
    return arr_mute

def shift(arr, tracehead, align, atime, fs):

    """ Shift VSP to align along direct arrival or to Two Way Time
   
    Adapted from from https://stackoverflow.com/questions/30399534/
                           shift-elements-in-a-numpy-array
    
    arr = VSP data array
    thead = header array to get xshift, the observed time header for indexing
    align = switch for alignment or shits to one-way or two-way time
    atime  = time to be aligned along 
    
    """
    import numpy as np
    
    print("\u0332".join('\nShifting Parameters :'))
    
    newhead = np.copy(tracehead)
    
    xshift = tracehead[:,8] * (fs/1000) # get index number for travel time    
    atime = int(atime *(fs/1000))       # get index number for alignment time
    
    arr2 = np.zeros(shape = (arr.shape[0], arr.shape[1]), dtype=np.float32)            
    print (' fs :', fs,'\n', ' atime :', atime,'\n',' first arr2 shape :', 
           arr2.shape)           
    pad_align = atime      # shallow traces can be cut of if atime > xshift    
    pad_twt = arr.shape[1]-int(np.max(xshift))
                               
    if align == 'up':           # align upgoing (TWT)        
        xshift = xshift.astype(int)        
        arr = np.pad(arr, ((0,0),(0, pad_twt)), 'constant')
        print (' pad twt : ', pad_twt, ' arr shape :', arr.shape)        
#        newhead[:,8] = tracehead[:,8] * 2        
        arr2 = np.zeros(shape = (arr.shape[0], arr.shape[1]),dtype=np.float32)        
        print (' second arr2 shape :', arr2.shape)
               
    elif align == 'down':      # align downgoing(flatten)      
        xshift = xshift - atime        
        xshift = xshift.astype(int) * -1        
        arr = np.pad(arr, ((0,0),(0, pad_align)), 'constant')        
        arr2 = np.zeros(shape = (arr.shape[0], arr.shape[1]),dtype=np.float32)
        
    elif align == 'unalign':  # remove downgoing alignment (back to OWT)       
        xshift = xshift.astype(int) - atime 

    for i, trace in enumerate(arr[::1, :]):        
        arr2[i,] = np.roll(arr[i,],xshift[i]) # careful with input array shape
        
        if xshift[i] > 0:            
            arr2[i,:xshift[i]] = 0            # [1, 4000] need row number 0            
        elif xshift[i] < 0:            
            arr2[i,xshift[i]:] = 0
            
    return arr2, newhead 
