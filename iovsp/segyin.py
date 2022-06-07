
def readsegyio3(inputfile, file_headers,DF_ASL, SrcElev, SRD_ASL, PTS):
    
    """ Load a segy file using segyio from Equinor

    Number of traces can come from binary header or trace headers or
    the shape of the data file
        
    1. open the segy file.
    2. The numpy data file is created in the first with loop.
    3. The data file shape is used to get number of traces and samples.
    4. The file is closed automatically at end of with...
    5. Then open file again to create the header arrays
        
    inputfile - segy file
    file_headers  - optionally print binary headers to terminal
    DF_ASL - drill floor elevation above Sea Level
    SrcElev - source elevation above sea level
    SRD_ASL - seismic reference datum above sea level
        
    Some useful headers are filled for velocity calculation
        
    For alternatives to my method can use panda dataframes for headers:
        
    https://github.com/equinor/segyio-notebooks/blob/master/notebooks/basic/02_segy_quicklook.ipynb
            
    TRACE_SEQUENCE_LINE= 1
    TRACE_SEQUENCE_FILE= 5
    FieldRecord= 9
    TraceNumber= 13
    EnergySourcePoint= 17
    CDP= 21
    CDP_TRACE= 25
    TraceIdentificationCode= 29
    NSummedTraces= 31
    NStackedTraces= 33
    DataUse= 35
    offset= 37
    ReceiverGroupElevation= 41
    SourceSurfaceElevation= 45
    SourceDepth= 49
    ReceiverDatumElevation= 53
    SourceDatumElevation= 57
    SourceWaterDepth= 61
    GroupWaterDepth= 65
    ElevationScalar= 69
    SourceGroupScalar= 71
    SourceX= 73
    SourceY= 77
    GroupX= 81
    GroupY= 85
    CoordinateUnits= 89
    WeatheringVelocity= 91
    SubWeatheringVelocity= 93
    SourceUpholeTime= 95
    GroupUpholeTime= 97
    SourceStaticCorrection= 99
    GroupStaticCorrection= 101
    TotalStaticApplied= 103
    LagTimeA= 105
    LagTimeB= 107
    DelayRecordingTime= 109
    MuteTimeStart= 111
    MuteTimeEND= 113
    TRACE_SAMPLE_COUNT= 115
    TRACE_SAMPLE_INTERVAL= 117
    GainType= 119
    InstrumentGainConstant= 121
    InstrumentInitialGain= 123
    Correlated= 125
    SweepFrequencyStart= 127
    SweepFrequencyEnd= 129
    SweepLength= 131
    SweepType= 133
    SweepTraceTaperLengthStart= 135
    SweepTraceTaperLengthEnd= 137
    TaperType= 139
    AliasFilterFrequency= 141
    AliasFilterSlope= 143
    NotchFilterFrequency= 145
    NotchFilterSlope= 147
    LowCutFrequency= 149
    HighCutFrequency= 151
    LowCutSlope= 153
    HighCutSlope= 155
    YearDataRecorded= 157
    DayOfYear= 159
    HourOfDay= 161
    MinuteOfHour= 163
    SecondOfMinute= 165
    TimeBaseCode= 167
    TraceWeightingFactor= 169
    GeophoneGroupNumberRoll1= 171
    GeophoneGroupNumberFirstTraceOrigField= 173
    GeophoneGroupNumberLastTraceOrigField= 175
    GapSize= 177
    OverTravel= 179
    CDP_X= 181
    CDP_Y= 185
    INLINE_3D= 189
    CROSSLINE_3D= 193
    ShotPoint= 197
    ShotPointScalar= 201
    TraceValueMeasurementUnit= 203
    TransductionConstantMantissa= 205
    TransductionConstantPower= 209
    TransductionUnit= 211
    TraceIdentifier= 213
    ScalarTraceHeader= 215
    SourceType= 217
    SourceEnergyDirectionMantissa= 219
    SourceEnergyDirectionExponent= 223
    SourceMeasurementMantissa= 225
    SourceMeasurementExponent= 229
    SourceMeasurementUnit= 231
    UnassignedInt1= 233
    UnassignedInt2= 237
    
    * VSProwess puts rcv MD in 37, TVD in 41
    """
    import segyio
    import numpy as np
#    reload(np)
    
    ###### open the segy file and create numpy data array ################
    
    with segyio.open(inputfile, ignore_geometry=True) as f:
        data = segyio.tools.collect(f.trace[:])
        delta_t = segyio.tools.dt(f)

    nrcv, samples = data.shape
    
    ###### open the segy file and read trace headers ######################
    
    with segyio.open(inputfile, ignore_geometry=True) as f:
        trnum = np.array(f.attributes(1))
        ffid = np.array(f.attributes(9))
        src = np.array(f.attributes(17))
        nsamp = np.array(f.attributes(115))
        srate = np.array(f.attributes(117))
        zscale = np.array(f.attributes(69))
        mdpth = np.array(f.attributes(37)/abs(zscale))
        tvddpth = np.array(f.attributes(41)/abs(zscale))
        srd = np.array(f.attributes(53)/abs(zscale))
        scalcoord = np.array(f.attributes(71))
        xsrc = np.array(f.attributes(73)/abs(scalcoord))
        ysrc = np.array(f.attributes(77)/abs(scalcoord))
        sdpth = np.array(f.attributes(49)/abs(zscale))
        xrcv = np.array(f.attributes(81)/abs(scalcoord))
        yrcv = np.array(f.attributes(85)/abs(scalcoord))
        auxtime_ms  = np.array(f.attributes(105))          #trunc. to nearest ms
        ttime_ms  = np.array(f.attributes(107))            # trunc.to nearest ms
        ttime  = np.array(f.attributes(197))/100     # divide if used
        iline = np.array(f.attributes(189))
        
    ###### open the segy file and print binary and text headers ###########
    
    if (file_headers == 'y') or (file_headers == 'Y'):        
        with segyio.open(inputfile, ignore_geometry=True) as f:
            txt_hed = segyio.tools.wrap(f.text[0]) # [1...] are extended headers
            bin_hed = f.bin           # this is a dictionary with keys and values

            for k, v in bin_hed.items():
                keys = str(k)
                value = int(v)
                print ("\t{:<23} {:<10} ".format(keys, value))
           
            print ('\n',txt_hed)
                
    ############## shift elevations to reference datum ######################
    
    SrcZ = (sdpth *0) + SrcElev # careful, comment out if field is correct
    TVD_Src = tvddpth - (DF_ASL - SrcElev)        
    TVD_SRD = tvddpth - (DF_ASL - SRD_ASL)    
    SrcZ_SRD = SrcZ-SRD_ASL
    
    ############## merge arrays and make a nice table ######################
    
    if (PTS == 'y') or (PTS == 'Y'):        
        from tabulate import tabulate        
        thead_tabl = np.vstack((trnum,mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD, iline, ffid,src)).T        
        print (' table header file shape :', thead_tabl.shape)    

        cheaders2 = ["Trc\nNum", "Rcz\nMD", "Rcz\nTVD","Rcv X", 
                    "Rcv Y","Src X", "Src Y",
                    "Src Z", "Obs\nTime","TVD\nSRD", 
                    "TVD \nSrcZ","SrcZ\nSRD",
                    "ILN", "FFID","Src"]
        numfmt = (".0f",".1f", ".1f", ".1f", ".1f",".1f", ".1f",".1f", ".2f",".1f", ".1f",".1f", ".0f",".0f",".0f")                 
        table = tabulate(thead_tabl, headers = cheaders2,  floatfmt = numfmt)#,tablefmt="pretty")

        print(table)
    
    ############## make a file of trace headers #############################
    
    thead = np.vstack((trnum, mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD, iline, ffid,
                      src))
    ############## sort out sample rates ####################################
    
    numsamp = nsamp[0]
    samprate = srate[0]                    # sample interval in microseconds
    samprate_hz = 1000000/samprate
        
    #tindex = np.zeros(shape = (data.shape[0], int(numsamp*(samprate/1000))))
    # 
    # for k in range(0,data.shape[0]):
    #     tindex[k,:] = np.arange(0, numsamp*(samprate/1000) )  # samprate in ms
        
    ########### print some QC data to screen ################################
    
    print("\u0332".join('\nData Loading with segyio Stats :'))
    print ('\n data shape :', data.shape, ' delta t :', delta_t)    
    print (' data type :', data.dtype)    
    print (' trace header file shape :', thead.T.shape)    
    print (' samples :', data.shape[1],' traces :', data.shape[0], \
           ' fs samprate hz : ', samprate_hz, \
           'samprate microseconds : ', samprate, \
           '\n numsamp from headers : ', numsamp)    
    print (' first time header value : ',ttime[0:1], \
           '\n first auxilliary time header value :', auxtime_ms[0:1])    
    print (' source depth from header trace 1 :', sdpth[0:1])

    
    return data, numsamp, samprate, samprate_hz, thead.T
    
def readsegy(inputfile, file_headers,DF_ASL, SrcElev, SRD_ASL, PTS):
    
    """ Load a segy file using obspy read_segy routine

    Obspy would need to be installed - but it is big so be careful
    
    Number of traces can come from binary header or trace headers or
    the shape of the data file
        
    nrcv = section.stats.binary_file_header.number_of_data_traces_per_ens.
    can be useful if filled
    
    We use the matrix shape in this example to figure out the number of 
    traces and the number of samples per trace
    
    inputfile - segy file
    file_headers  - optionally print binary headers to terminal
    DF_ASL - drill floor elevation above Sea Level
    SrcElev - source elevation above sea level
    SRD_ASL - seismic reference datum above sea level
    
    Some useful headers are filled for velocity calculation
    """
    
    from obspy.io.segy.segy import _read_segy
    
    section = _read_segy(inputfile, unpack_headers = 'True')    
    if (file_headers == 'y') or (file_headers == 'Y'):    
        print (section.binary_file_header)    
        x = np.array(list(section.textual_file_header.decode()))        
        print('\n'.join(''.join(row) for row in x.reshape((40, 80))))
            
    data = np.vstack([d.data for d in section.traces])    
    nrcv, samples = data.shape
  
    print("\u0332".join('\nData Loading Stats :'))    
    print (' data shape :', data.shape)    
    print (' samples :', data.shape[1],' traces :', data.shape[0])
    
    ####### get basic info from trace headers and save in 1d arrays #########      
    
    trnum, ffid, src, nsamp, srate, zscale = (np.empty((nrcv,1)),np.empty((nrcv,1)),
                                        np.empty((nrcv,1)),np.empty((nrcv,1)), 
                                        np.empty((nrcv,1)), np.empty((nrcv,1)))

    srd, zrcv, tsr, mdpth, tvddpth = (np.empty((nrcv,1)),np.empty((nrcv,1)),
                                        np.empty((nrcv,1)),np.empty((nrcv,1)),
                                        np.empty((nrcv,1)))

    scalcoord, xsrc, ysrc, sdpth = (np.empty((nrcv,1)),np.empty((nrcv,1)),
                                        np.empty((nrcv,1)),np.empty((nrcv,1))) 

    xrcv, yrcv, auxtime_ms, ttime = (np.empty((nrcv,1)),np.empty((nrcv,1)),
                                        np.empty((nrcv,1)),np.empty((nrcv,1)))

    ttime_ms,iline = (np.empty((nrcv,1)),np.empty((nrcv,1)))
    
    n=0
    for tr in section.traces:
    
        trnum[n] = tr.header.trace_sequence_number_within_line
        ffid[n] = tr.header.original_field_record_number
        src[n] = tr.header.energy_source_point_number
        nsamp[n] = tr.header.number_of_samples_in_this_trace
        srate[n] = tr.header.sample_interval_in_ms_for_this_trace
        zscale[n] = tr.header.scalar_to_be_applied_to_all_elevations_and_depths
        mdpth[n] = tr.header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group/abs(zscale[n])
        tvddpth[n] = tr.header.receiver_group_elevation/abs(zscale[n])
        srd[n] = tr.header.datum_elevation_at_receiver_group/abs(zscale[n])
        scalcoord[n] = tr.header.scalar_to_be_applied_to_all_coordinates
        xsrc[n] = tr.header.source_coordinate_x/abs(scalcoord[n])
        ysrc[n] = tr.header.source_coordinate_y/abs(scalcoord[n])
        sdpth[n] = tr.header.source_depth_below_surface/abs(zscale[n])
        xrcv[n] = tr.header.group_coordinate_x/abs(scalcoord[n])
        yrcv[n] = tr.header.group_coordinate_y/abs(scalcoord[n])
        auxtime_ms[n]  = tr.header.lag_time_A          #trunc. to nearest ms
        ttime_ms[n]  = tr.header.lag_time_B            # trunc.to nearest ms
        ttime[n]  = tr.header.shotpoint_number/100     # divide if used
        iline[n] = tr.header.for_3d_poststack_data_this_field_is_for_in_line_number
        n+=1
    
    ############## shift elevations to reference datum ######################
    
    SrcZ = (sdpth *0) + SrcElev # careful, comment out if field is correct
    TVD_Src = tvddpth - (DF_ASL - SrcElev)        
    TVD_SRD = tvddpth - (DF_ASL - SRD_ASL)    
    SrcZ_SRD = SrcZ-SRD_ASL
    
    ############## make a file of trace headers #############################
    
    thead = np.hstack((trnum, mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD,ttime_ms, iline, ffid,
                      src))
    print (' trace header file shape :', thead.shape)
        
    ############## sort out sample rates ####################################
    
    numsamp = nsamp[0]
    samprate = srate[0]                    # sample interval in microseconds
    samprate_hz = 1000000/samprate
    
    print (' ttime: ',ttime.shape, ' fs samprate hz : ', samprate_hz, \
           'samprate microseconds : ', samprate, \
           '\n numsamp from headers : ', numsamp)
    
    print (' first time header value : ',ttime[0:1], \
           '\n first auxilliary time header value :', auxtime_ms[0:1])
    
    print (' source depth from header trace 1 :', sdpth[0:1])
    
    ############## merge arrays and make a nice table ######################
    
    if (PTS == 'y') or (PTS == 'Y'):
        
        from tabulate import tabulate
        
        thead_tabl = np.hstack((mdpth, tvddpth,xrcv, yrcv,  xsrc, ysrc, sdpth, 
                       ttime,TVD_SRD, TVD_Src, SrcZ_SRD, ttime_ms, ffid))    
        headers = [ "Rcv MD\n\n37-40", "Rcv TVD\nSRD\n41-44","Rcv X\n\n81-84", 
                    "Rcv Y\n\n85-88","Src X\n\n73-76", "Src Y\n\n77-80",
                    "Src Z\n\n49-52", "Obs Time\n\n197-200","TVD Depth\nSRD\nCalc.", 
                    "TVD SrcZ\nSource Z\n Calc.","SrcZ SRD\nSRD\nCalc.","Obs Time\n\n107-108"  "FFID\n\n9-12"]    
        table = tabulate(thead_tabl, headers, tablefmt="pretty")

        print(table)

    return data, numsamp, samprate, samprate_hz, thead    