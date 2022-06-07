def write_segyio(data, tracehead,fs, name):
    
    '''Write a segy file using segyio routines
    
    Number of traces can come from binary header or trace headers or
    the shape of the data file

    We use the matrix shape in this example to figure out the number of 
    traces and the number of samples per trace
        
    data - 2d numpy array of VSP data [traces, samples]
    tracehead - 2d array of trace headers [traces, header values]
    fs - sample rate in Hz
    name - output file path
    
    1. Arbitrary byte positions not accessible for headers

    '''
    
    import segyio
#    import numpy as np
    
    ################# print some basic information ###############

    print("\u0332".join('\nWrite Segyio Stats :'))    
    print ('Data shape [0] :', data.shape[0],'Data shape [1] :', data.shape[1])    
    print('Trace header shape', tracehead.shape)    
    print ('time header :', tracehead[0:2,:])
    
    ################### geometry scalars  ########################

    scalcoord = 10
    scaldep = 10
    
    ################# read trace headers ##########################

    FFID = tracehead[:,14].astype(int)
    SRC = tracehead[:,15].astype(int)
    MD = (tracehead[:,1]*scaldep).astype(int)
    TVD = (tracehead[:,2]*scaldep).astype(int)
    RcvX = (tracehead[:,3]*scalcoord).astype(int)
    RcvY = (tracehead[:,4]*scalcoord).astype(int)
    SrcX = (tracehead[:,5]*scalcoord).astype(int)
    SrcY = (tracehead[:,6]*scalcoord).astype(int)
    SrcZ = (tracehead[:,7]*scalcoord).astype(int)
    TVD_SRD = (tracehead[:,9]*scaldep).astype(int)
    TVD_Src = (tracehead[:,10]*scaldep).astype(int)
    SrcZ_SRD = (tracehead[:,11]*scaldep).astype(int)
    Tobs = (tracehead[:,8]*100).astype(int)

    print ('MD shape :', MD.shape, ' MD dtype :' , MD.dtype)
    print (' MD [1:10] :', MD[0:10], ' TVD[0:10] :',TVD[0:10])
    
    ################# create text header ##########################
   
    text_header = {1: 'Synthetic Walkaway VSP', 
    2: 'SEG-Y read-write tests',
    3: 'VSProwess',
    4: ' ',
    5: 'WELL: Allan 1',
    6: ' ',
    7: 'Processed by VSP Consultants',
    8: ' ' ,
    9: 'Reference Elevation:  ft',
    10: ' ',
    11: 'Geophone component: 025-028 Z=1, X=2, Y=3' ,
    12: 'Source Easting:  073-076     Source Elev: 045-048' ,
    13: 'Source Northing: 077-080' ,
    14: 'Receiver Easting: 081-084    Measured Depth: 037-040' ,
    15: 'Receiver Northing: 085-088   Vertical Depth: 041-044' ,
    16: '  ',
    17: 'Uncorrected Pick Time: 107-108' ,
    18: 'TWO-Way Time : 109-110' ,
    19: ' ',
    20: ' ',
    21: 'Units  = Survey Feet' ,
    22: 'Wellhead Easting (ft): 0' ,
    23: 'Wellhead Northing (ft): 0' ,
    24: ' ',
    25: ' ****Processing Steps: ***********' ,
    26: 'Original file from VSProwess' ,
    27: 'Read and write with segyio',
    28: ' ' ,
    29: ' ',
    30: ' ', 
    31: ' ',
    32: ' ',
    33: ' ',
    34: ' ',
    35: ' ',
    36: ' ',
    37: ' ',
    38: 'VSProwess processing system'  }
    
    txt_hed = segyio.tools.create_text_header(text_header)
    
    ################# set some global parameters ###################

    spec = segyio.spec()
    spec.samples = list(range(data.shape[1]))
    spec.sorting = 0 # 0 for unstructured ie. not a cube of data
    spec.format  = 5 #1 = IBM float, 5 = IEEE float
    spec.tracecount = data.shape[0]
    
    print( ' spec.tracecount :', spec.tracecount)
    
    ################# create output file  #########################
    
    with segyio.create(name, spec) as fout:
        fout.trace = data                      # populate traces
        fout.text[0] = txt_hed
        fout.bin.update(hns=len(spec.samples)) # update ample count binary header, su style keys 
        fout.bin.update(format=5)      #1 = IBM float, 5 = IEEE float        
        fout.bin.update(mfeet=2)        
        fout.bin.update(rev=1)                                       
#        fout.bin.update{segyio.BinField.MeasurementSystem: 2} #alternate methods
#        fout.bin.update{3255: 2}      # all binary keys use 3200 + 55 etc for byte positions

    ############## create output trace headers  ####################
        tr=0
        for i, trace in enumerate(data):
#            fout.header[tr] = {segyio.tracefield.TraceField.ReceiverGroupElevation: MD[tr]}
            fout.header[tr] = {1: i+1}    
            fout.header[tr] = {115: data.shape[1]}
            fout.header[tr] = {117:int((1/fs)*1000000)}
            fout.header[tr] = {9: FFID[tr]}
            fout.header[tr] = {13: SRC[tr]}
            fout.header[tr] = {37: MD[tr]}
            fout.header[tr] = {41: TVD[tr]}
            fout.header[tr] = {69: -1*scaldep}
            fout.header[tr] = {71: -1*scalcoord}
            fout.header[tr] = {73: SrcX[tr]}
            fout.header[tr] = {77: SrcY[tr]}
            fout.header[tr] = {49: SrcZ[tr]}
            fout.header[tr] = {81: RcvX[tr]}
            fout.header[tr] = {85: RcvY[tr]}            
            fout.header[tr] = {107: int(Tobs[tr]/100)}
            fout.header[tr] = {197: Tobs[tr]}
                    
            tr +=1
            
def write_segy(data, tracehead,fs, name):
    ''' Use Obspy to write segy
    
    Obspy would need to be installed - but it is big so be careful
        
    This is left in just in case it becomes useful
    '''
    
    print("\u0332".join('\nWrite Segy Stats :'))    
    print ('Data shape [0] :', data.shape[0],'Data shape [1] :', data.shape[1])    
    print('Trace header shape', tracehead.shape)    
    print ('time header :', tracehead[0:2,:])
    
    from obspy.core import read, Trace, AttribDict, Stream, UTCDateTime
    from obspy.io.segy.segy import SEGYBinaryFileHeader

    from obspy.io.segy.segy import SEGYTraceHeader
    from obspy.io.segy.segy import SEGYFile
    
    MD = tracehead[:,1].astype(int)
#    print ('MD shape :', MD.shape, ' MD dtype :' , MD.dtype)
    TVD = tracehead[:,2].astype(int)
    TVD = TVD.reshape(-1)
    RcvX = tracehead[:,3].astype(int)
    RcvY = tracehead[:,4].astype(int)
    SrcX = tracehead[:,5].astype(int)
    SrcY = tracehead[:,6].astype(int)
    SrcZ = tracehead[:,7].astype(int)
    TVD_SRD = tracehead[:,9].astype(int)
    TVD_Src = tracehead[:,10].astype(int)
    SrcZ_SRD = tracehead[:,11].astype(int)
    Tobs = (tracehead[:,8]*100).astype(int)
    IntV = tracehead[:,13].astype(int)
    
    print ('MD shape :', MD.shape, ' MD dtype :' , MD.dtype)
    
    stream = obspy.Stream()
    
    scalcoord = 10
    scaldep = 10
    
    for i, trace in enumerate(data):    
        # Make the trace.
        tr = obspy.Trace(trace)        
        tr.data = np.require(tr.data, dtype=np.float32) 
        # Add required data.
        tr.stats.delta = 1/fs # fs is sample rate in hertz
        # tr.stats.starttime = 0  # Not strictly required.
        # Add yet more to the header (optional).
        tr.stats.segy = {'trace_header': SEGYTraceHeader()}
        
        tr.stats.segy.trace_header.trace_sequence_number_within_line = i + 1
        
        #for reference trhead = np.vstack([mdpth, tvddpth,xrcv, yrcv, xsrc, ysrc, sdpth, ttime])
          
        tr.stats.segy.trace_header.number_of_samples_in_this_trace = data.shape[1]
        tr.stats.segy.trace_header.sample_interval_in_ms_for_this_trace = (1/fs)*1000000
        tr.stats.segy.trace_header.scalar_to_be_applied_to_all_elevations_and_depths = scaldep*-1
        tr.stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group = \
                                int(MD[i]*scaldep)  
        tr.stats.segy.trace_header.receiver_group_elevation = int(TVD[i]*scaldep)
        tr.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = scalcoord *-1
        tr.stats.segy.trace_header.source_coordinate_x = int(SrcX[i]*scalcoord)
        tr.stats.segy.trace_header.source_coordinate_y = int(SrcY[i]*scalcoord)
        tr.stats.segy.trace_header.source_depth_below_surface = int(SrcZ[i]*scaldep)
        tr.stats.segy.trace_header.group_coordinate_x = int(RcvX[i]*scalcoord)
        tr.stats.segy.trace_header.group_coordinate_y = int(RcvY[i]*scalcoord)
        tr.stats.segy.trace_header.lag_time_B = int(Tobs[i]/100)
        tr.stats.segy.trace_header.shotpoint_number = int(Tobs[i])
        tr.stats.segy.trace_header.transduction_units = int(IntV[i])
    
        # Append the trace to the stream.
        stream.append(tr)
        
    stream.stats = AttribDict()
    stream.stats.textual_file_header = b'C01 Synthetic Zero Offset VSP \n\
    C02 Median filtering tests\n\
    C03 VSProwess\n\
    C04\n\
    C05 WELL: Allan 1\n\
    C07 Processed by VSP Consultants\n\
    C08\n\
    C09 Reference Elevation:  ft\n\
    C10\n\
    C11 Geophone component: 025-028 Z=1, X=2, Y=3\n\
    C12 Source Easting:  073-076     Source Elev: 045-048\n\
    C13 Source Northing: 077-080\n\
    C14 Receiver Easting: 081-084    Measured Depth: 037-040\n\
    C15 Receiver Northing: 085-088   Vertical Depth: 041-044\n\
    C16\n\
    C17 Uncorrected Pick Time: 107-108\n\
    C18 TWO-Way Time : 109-110\n\
    C19\n\
    C20\n\
    C21 Units  = Survey Feet\n\
    C22 Wellhead Easting (ft): 0\n\
    C23 Wellhead Northing (ft): 0\n\
    C24\n\
    C25 ****Processing Steps: ***********\n\
    C26 Optimized Stack\n\
    C27 Median noise reduction\n\
    C28 Travel Times picked\n\
    C29\n\
    C30\n\
    C31\n\
    C32\n\
    C33\n\
    C34\n\
    C35\n\
    C36\n\
    C37\n\
    C38 VSProwess processing system\n\ '
    
    print (stream)

    stream.write('%s'%(name), format='SEGY', data_encoding=5)
    
