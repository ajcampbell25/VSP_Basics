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
    print("\u0332".join('\nMute Stats :'))    
    print(' data shape : ', arr.shape, ' data dtype : ', arr.dtype)    
    print (' headers shape :', thead.shape)    
    print (' 1WT middle trace :', thead[thead.shape[0]//2,8] )    
    print (' 2WT middle trace :', thead[thead.shape[0]//2,-2] )    

    import numpy as np
    
    if align == 'owt':        
        mute_time = thead[:,8]*fs/1000
        
    elif align == 'twt':
        #mute_time = thead[:,8]*2*fs/1000 # if TWT is not in header array 
        mute_time = thead[:,-2]*fs/1000 # this assumes 2WT is in the header array    
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
    align = switch for alignment or shifts to one-way or two-way time
    atime  = time to be aligned along 
    
    align='up' : data is assumed to be in one way time (OWT). Shifting by  the 
                 measured one-way-time gets 2 way time (2WT). This assumes 
                 vertical ray paths.
    
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
        xshift = (tracehead[:,-2] * (fs/1000))/2# test using geometry TWT
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

def platform_chk():
    r''' check if notebook is running in Jupyter notebook/lab
    or something else like vscode.
    Interactive plots react differently depending on platform.

    See 'https://stackoverflow.com/questions/15411967/how-can-i-check-/
    if-code-is-executed-in-the-ipython-notebook'

    returns : string indicating Jupyter or VS Code

    Useage : 
    # use 'ipympl' interactive graphic cells in jupyter lab??
    # %matplotlib ipympl
    # use 'widget' for interactive graphic cells in vscode

    platform = utilvsp.platform_chk()

    print (platform)
    if (platform =="in VS Code")or(platform =="in VSCODE"):
        %matplotlib widget
    else:
        %matplotlib ipympl
    '''
    import os
    import IPython as ipy

    # add string sources only
    sources = str(os.environ.keys()) + \
            ipy.get_ipython().__class__.__name__

    # make pattern of unique keys
    checks = {'SPYDER': 'Spyder', 'QTIPYTHON': 'qt IPython', 'VSCODE': 
            'VS Code', 'ZMQINTERACTIVEshell': 'Jupyter', }

    results = []
    msg = []

    for k, v in checks.items():
        u = str(k.upper())
        if u in sources.upper():
            results.append(checks[k])

    if not results:
        msg.append("Unknown IDE")
    else:
        msg.append("Program working ")
        while results:
            msg.append(f"in {results.pop()}")
            if results:
                msg.append(' with')

    print(''.join(msg))
    return msg[-1]
