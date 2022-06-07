def importascii(DF_ASL, SrcElev, SRD_ASL, td_file, dev_file):
    
    """ 
    Used if deviation needs to be loaded. 
    Also useful if travel times are available in an ascii format
    
    """

    f = open(td_file, 'r')

    TD = np.genfromtxt(f, delimiter=',')
    TD = np.delete(TD,slice(0,2),0) # Erases the first 2 rows (i.e. the header)

    f2 = open(dev_file, 'r')

    WellData = np.genfromtxt(f2, delimiter=',')
    
    # take care of well deviation 
    # VSProwess does not apply deviaion if inc is < 5 degrees, so external deviation info may be required
    MD1 = TD[:,12]
    SrcX = TD[:,28] 
    SrcY = TD[:,29] 
    SrcZ = TD[:,30]
    Pick1 = TD[:,19] * 1000
    TVD = TD[:,13] 
    RecvX = TD[:,25]
    RecvY = TD[:,26]
    
    # this example uses an external deviation file from SLB field enginee
#    MD2 = WellData[:,0] 
#    TVD2 = WellData[:,3] 
#    RecvX = WellData[:,1]
#    RecvY = WellData[:,2]

    SrcZ = (SrcZ *0) + SrcElev # field headers were incorrect for this example, comment out if field is correct
    TVD_Src = TVD - (DF_ASL - SrcElev)
    TVD_SRD = TVD - (DF_ASL - SRD_ASL)
    SrcZ_SRD = SrcZ-SRD_ASL
    
    return RecvX, RecvY, TVD_SRD, TVD_Src, SrcX, SrcY, SrcZ_SRD, Pick1

    #DevData = np.delete(WellData,slice(0,2),0) # Erases the first 2 row (i.e. the header)
    
def geovel(trheaders, repvel, PTS):
    '''    
    Calculate source-receiver offset
    Verticalize the times relative to source elevation
    Add in source to SRD time delay
    ''' 
    from tabulate import tabulate
    import numpy as np 

    RcvX = trheaders[:,3]    
    RcvY = trheaders[:,4]    
    SrcX = trheaders[:,5]
    SrcY = trheaders[:,6]    
    TVD_SRD = trheaders[:,9]    
    TVD_Src = trheaders[:,10]    
    SrcZ_SRD = trheaders[:,11]    
    ttime = trheaders[:,8]
    
    a = RcvX - SrcX
    b = RcvY - SrcY
    
    SROffset = np.sqrt(a**2+b**2)
    VertTT_Src = ttime * np.cos(np.arctan(SROffset/TVD_Src))    
    VertTT_SRD = VertTT_Src - ((SrcZ_SRD)/repvel)
    T_Diff = np.ediff1d(VertTT_SRD/1000)    
    D_Diff = np.ediff1d(TVD_SRD)    

    firstTDiff = VertTT_SRD[0]/1000    
    firstDDiff = TVD_SRD[0]    

    TDiff = np.hstack((firstTDiff, T_Diff))    
    DDiff = np.hstack((firstDDiff, D_Diff))
    IntVel = DDiff/TDiff
    FirstVel = TVD_SRD[0]/VertTT_SRD[0]    
    DisplayVelRnd = np.rint(IntVel)             #round and convert to integer    
    DisplayVelInt=DisplayVelRnd.astype(int)    
    IntVel = IntVel.reshape(-1,1)
    
    TWT = VertTT_Src.reshape(-1,1)*2
    
    vheaders = np.hstack((trheaders,TWT, IntVel))
    
    IntVel = IntVel.reshape(-1,)
    
    print("\u0332".join('\nGeovel Stats :'))
    print(' Theader shape', trheaders.shape, 'Vheader shape:', vheaders.shape)    
    print (' TVD_SRD shape :', TVD_SRD.shape, ' IntVel shape :', IntVel.shape)      
    print (' TVDSrc 2 vals  ',TVD_Src[0:2],'\n TVDSRD 2 vals ',TVD_SRD[0:2],
           '\n TTVert 2 vals  ', VertTT_SRD[0:2])
        
    if (PTS == 'y') or (PTS == 'Y'):    
        # first table contains velocity calculation inputs
        pdat = np.vstack((TVD_SRD, VertTT_SRD,TDiff, DDiff, IntVel)).T    
        headers = ["Depth SRD", "Vert Time SRD", "Delta T", "Delta Z","Vp"]    
        table = tabulate(pdat, headers, tablefmt="fancy_grid")
        print(table) 

        # second table is a listing of new theader file
        headers2 = ["Trc\nNum", "Rcv\nMD", "Rcv\nTVD","Rcv X", 
                    "Rcv Y","Src X", "Src Y",
                    "Src Z", "Obs\nTime","TVD\nDpth", 
                    "TVD \nSrcZ","SrcZ\nSRD",
                    "ILN", "FFID","Src","TWT\nSrc", "Int Vel."]
        numfmt = (".0f",".1f", ".1f", ".1f", ".1f",".1f", ".1f",".1f", ".2f",".1f", ".1f",
                  ".1f", ".0f",".0f",".0f",".2f", ".0f")                 
        table2 = tabulate(vheaders, headers = headers2,  floatfmt = numfmt)#,tablefmt="pretty")
        print(table2)
        
    return vheaders

def header_listing(thead, listing):
    
    from tabulate import tabulate

    if listing == 'tz':    
        pdat = thead[:,[0,8,1]]    
        headers = [ " Trace_Num.","Observed Time", "Measured Depth",]    
        table = tabulate(pdat, headers, tablefmt="fancy_grid")
        print(table)     
    