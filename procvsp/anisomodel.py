def aniso_params(Vpz,Vsz, epsilon, delta):
    ''' Estimate eta gamma and sigma
    Convert Thomsen's paramters to Velocity components
    '''
    import numpy as np

    eta = (epsilon-delta)/(1+2*delta)
    gamma = (Vpz/Vsz)
    sigma = gamma**2*(epsilon-delta)
        
    # calculate horizontal and near offset velocities    
    Vpx = Vpz*np.sqrt(1+ 2*epsilon) # horizontal P velocity
    Vpn = Vpz*np.sqrt(1+ 2*delta)   # near offset p velocity
    Vsx = Vsz
    Vsn = Vsz*np.sqrt(1+ 2*sigma)   # near offset S velocity 

    return Vpx,Vpn,Vsx,Vsn
    
def phase_vel(Vpz,Vsz, epsilon, delta):
    ''' Use the exact solution for phase velocity

    Vp^2 = .5*(Vpx^2*sin(theta)^2+Vpz^2*cos(theta)^2 +Vsz^2+
    sqrt{(( Vpx^2-Vsz^2)*sin(theta)^2 +(Vpz^2-Vsz^2)*cos(theta)^2)^2 +
    (Vpz^2-Vsz^2)*(Vpn^2-Vpx^2)*(sin(2*theta))^2)} 

    Inputs
    Vpz: vertical p velocity
    Vsz: vertical s velocity
    epsilon: Thomsen's epsilon
    delta: Thomsen's delta

    Return:
    phase_velp: phase velocity for p waves
    phase_vels: phase velocity for s waves
 
    '''
    import numpy as np
    from math import sqrt, sin, cos, radians
    
    # generate array of wavefront propagation (phase) angles
#     # this version produces same order data as excel spreadsheet, but line plots are not quite correct
#    deg1 = np.arange(0,91,1)
#    deg2 = np.arange(-90,0,1)
#    degs = np.concatenate([deg1,deg2]) # this version produces same order data as
    degs = np.arange(-90,91,1)
    theta = np.radians(degs)    
 
    #print  (' Vpx :',Vpx,' Vpn :',Vpn,' Vsx :',Vsx,' Vsn :',Vsn)    

    Vpx,Vpn,Vsx,Vsn = aniso_params(Vpz,Vsz, epsilon, delta)
    
    phase_velp = np.zeros(shape = theta.shape)
    phase_vels = np.zeros(shape = theta.shape)

    # calculate phase velocities as function of wavefront propagation angle
    for i in range(0,theta.shape[0]):
        phase_velp[i,] = sqrt(.5)*sqrt((Vpx**2*sin(theta[i,])**2+Vpz**2*cos(theta[i,])**2 +Vsz**2)+ 
        sqrt((( Vpx**2-Vsz**2)*sin(theta[i,])**2 +(Vpz**2-Vsz**2)*cos(theta[i,])**2)**2 +
        (Vpz**2-Vsz**2)*(Vpn**2-Vpx**2)*(sin(2*theta[i,]))**2))

        phase_vels[i,] = sqrt(.5)*sqrt((Vpx**2*sin(theta[i,])**2+Vpz**2*cos(theta[i,])**2 +Vsz**2)-
        sqrt((( Vpx**2-Vsz**2)*sin(theta[i,])**2 +(Vpz**2-Vsz**2)*cos(theta[i,])**2)**2 +
        (Vpz**2-Vsz**2)*(Vpn**2-Vpx**2)*(sin(2*theta[i,]))**2))
        
    return phase_velp, phase_vels, theta

def phase_slow(Vpz, Vsz, epsilon, delta):
    '''
    Inputs
    Vpz: vertical p velocity
    Vsz: vertical s velocity
    epsilon: Thomsen's epsilon
    delta: Thomsen's delta

    Return:
    Slowness_vp: phase slowness p waves
    Slowness_vs: phase slowness s waves
    Slowness_xp: horizontal component of Slowness_vp
    Slowness_xs: horizontal component of Slowness_vs
    Slowness_zp: vertical component of Slowness_vp
    Slowness_zs: vertical component of Slowness_vs 
    '''
    import numpy as np

    phase_velp, phase_vels, theta= phase_vel(Vpz, Vsz, epsilon, delta) 
     
     #Convert velocity to slowness and calculate horizontal and vertical components
    Slowness_vp = 1/phase_velp
    Slowness_xp = np.sin(theta)*Slowness_vp
    Slowness_zp = np.cos(theta)*Slowness_vp

    Slowness_vs = 1/phase_vels
    Slowness_xs = np.sin(theta)*Slowness_vs
    Slowness_zs = np.cos(theta)*Slowness_vs
    '''     
    # make some QC tables
    from tabulate import tabulate
    pdat0 = np.vstack((Slowness_xp, Slowness_zp, Slowness_xs, Slowness_zs,Slowness_vp, Slowness_vs)).T    
    headers0 = ["Slowness_xp","Slowness_zp","Slowness_xs", "Slowness_zs", "Slowness_vp", "Slowness_vs"]    
    table0 = tabulate(pdat0, headers0, tablefmt="fancy_grid")
    print(table0)
    '''    
    return Slowness_vp, Slowness_xp,Slowness_zp,Slowness_vs, Slowness_xs,Slowness_zs, theta       
    
def aniso_phase(Vpz, Vsz, epsilon, delta, save):
    ''' Calculate and crossplot the vertical and horizontal Slowness 
    for P and S waves
    Inputs
    Vpz: vertical p velocity
    Vsz: vertical s velocity
    epsilon: Thomsen's epsilon
    delta: Thomsen's delta
        
    '''   
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec

    # get anisotropic response
    _,Slowxp_ani,Slowzp_ani,_,Slowxs_ani,Slowzs_ani,_ = phase_slow(Vpz, Vsz, epsilon, delta)

    # set epsilon and delta to 0 for isotropic       
    epsiloniso=0
    deltaiso=0
    _,Slowxp_iso,Slowzp_iso,_,Slowxs_iso,Slowzs_iso,_ = phase_slow(Vpz, Vsz, epsiloniso, deltaiso)

    # Plot frequency response (in amplitude and dB) and impulse response                        
    fig = plt.figure(figsize=(15,5))    
    gs = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace = .25)
    
    ax1 = plt.subplot(gs[0])    
     
    ax1.plot(Slowxp_ani, Slowzp_ani, label = 'Anisotropic')      
    ax1.plot(Slowxp_iso, Slowzp_iso, label = 'Isotropic',linestyle='--')      
    #ax1.set_xlim(0, f2*2)       # subjective choice
    ax1.set_xlabel('Horizontal Slowness')    
    ax1.set_ylabel('Vertical Slowness')    
    ax1.grid(True)    
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Compressional Phase Slowness',fontsize=12)
    ax1.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)

    ax2 = plt.subplot(gs[1])    
    ax2.plot(Slowxs_ani, Slowzs_ani,label = 'Anisotropic')      
    ax2.plot(Slowxs_iso, Slowzs_iso,label = 'Isotropic',linestyle='--')      
    #ax2.set_xlim(0, f2*4)      # subjective choice
    #ax2.set_ylim(-200,0)       # subjective choice
    ax2.set_xlabel('Horizontal Slowness')    
    ax2.set_ylabel('Vertical Slowness')    
    ax2.grid(True)    
    ax2.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax2.set_title('Shear Phase Slowness')
    ax2.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    ax2.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    plt.show()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\aniso_phase_epsilon%s%%_delta%s%%.png' 
        %(epsilon*100,delta*100), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)

    
def aniso_slowpol(Vpz, Vsz, epsilon, delta, save):
    ''' Calculate and crossplot the vertical and horizontal Slowness 
    for P and S waves. Use the slowness components to calculate the phase angle
    Inputs
    Vpz: vertical p velocity
    Vsz: vertical s velocity
    epsilon: Thomsen's epsilon
    delta: Thomsen's delta
        
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    
    # get anisotropic response slowness components and phase angle vector
    Slowvp_ani,Slowxp_ani,Slowzp_ani,Slowvs_ani,Slowxs_ani,Slowzs_ani, theta = phase_slow(Vpz, Vsz, epsilon, delta)
    
    # set epsilon and delta to 0 for isotropic   
    epsiloniso=0
    deltaiso=0
    Slowvp_iso,Slowxp_iso,Slowzp_iso,Slowvs_iso,Slowxs_iso,Slowzs_iso,_ = phase_slow(Vpz, Vsz, epsiloniso, deltaiso)

    # convert the phase angle vector from radians to degrees(for display)
    degs = np.degrees(theta)
    
    # Plot frequency response (in amplitude and dB) and impulse response                        
    fig = plt.figure(figsize=(15,5))    
    gs = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace = .25)
    
    ax1 = plt.subplot(gs[0])    
    ax1.plot(degs,Slowzp_ani, label = 'Anisotropic')      
    ax1.plot(degs,Slowzp_iso, label = 'Isotropic',linestyle='--')      
    #ax1.set_xlim(0, f2*2)       # subjective choice
    ax1.set_xlabel('Phase Angle')    
    ax1.set_ylabel('Vertical Slowness')    
    ax1.grid(True)    
    ax1.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Compressional Slowness Polarization',fontsize=12)
    ax1.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax2 = plt.subplot(gs[1])    
    ax2.plot(degs, Slowzs_ani,label = 'Anisotropic')      
    ax2.plot(degs, Slowzs_iso, label = 'Isotropic',linestyle='--')      
    #ax2.set_xlim(0, f2*4)      # subjective choice
    #ax2.set_ylim(-200,0)       # subjective choice
    ax2.set_xlabel('Phase Angle')    
    ax2.set_ylabel('Vertical Slowness')    
    ax2.grid(True)    
    ax2.legend(loc='best',borderaxespad=0, fontsize = 8)        
    ax2.set_title('Shear Slowness Polarization')
    ax2.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    ax2.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    plt.show()

    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\aniso_slowpol_epsilon%s%%_delta%s%%.png' 
        %(epsilon*100,delta*100), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)


def aniso_phase_calc(wvsptzoff,tu,Vpz,Vsz,epsilon,delta, save):
    ''' From walkaway VSP time-depth-offset listing calculate Sx and Sz
    and crossplot results.
    Sx is the horizontal slowness
    Sz is vertical slowness
    '''
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import numpy as np
    import warnings
    
    f = open(wvsptzoff, 'r')

    tzoff = np.genfromtxt(f, delimiter=',')
    tzoff = np.delete(tzoff,slice(0,1),0) # Erases the first 1 rows (i.e. the header)
    
    rcvnum = tzoff[:,2] # receiver number
    srcnum = tzoff[:,1] # source number
    MD = tzoff[:,3]
    SrcX = tzoff[:,11] 
    SrcY = tzoff[:,12] 
    SrcZ = tzoff[:,13]
    TTime = tzoff[:,6] 
    TVD = tzoff[:,4] 
    RecvX = tzoff[:,8]
    RecvY = tzoff[:,9]
    
    # ms to s conversion
    if tu == 'ms':
        TTime = TTime/1000
    
    # get number of traces from array shape
    numtraces = rcvnum.shape[0] 

    # Check if the input ascii file is shot or receiver sorted
    # Calculate the number of traces in a shot or receiver gather
    if SrcX[0,] == SrcX[1,]:
        srctest = SrcX[0,]
        print ('Input text sorted by shot ')
        gtype = 's' # gather type is a shot    
        # get the number of traces in a shot gather
        trcgath = np.count_nonzero( SrcX == srctest)
        print (' traces per shot gather :',trcgath)
    else:
        print ('Input text sorted by receiver ')
        gtype = 'r' # gather type is a receiver
        trcgath = np.count_nonzero(rcvnum == 1)
        print (' traces per receiver gather :',trcgath)

    # create output files initialzed to all zeros
    Slowz = np.zeros(shape = (numtraces//trcgath))
    Slowx = np.zeros(shape = (numtraces//trcgath))

    # calculate the vertical and horizontal slownesses        
    if gtype =='s':
        for i in range (0,numtraces,trcgath):
            k = i//trcgath 

            # 3 receiver window for Sz
            Slowz[k,]=1/((MD[i+2,]-MD[i,])/(TTime[i+2,]-TTime[i,])) 
    
            # 3 source window for Sx    
            if (i-trcgath) < 0:
                Slowx[k,] = 1/((SrcX[i,]-SrcX[(i+trcgath),])/(TTime[i,]-TTime[(i+trcgath),]))

            else:
                if (i+trcgath) >= numtraces:
                    Slowx[k,] = 1/((SrcX[i,]-SrcX[(i-trcgath),])/(TTime[i,]-TTime[(i-trcgath),]))

                else:
                    # divide by zero is possible but will generate a warning which we can suppress                   
                    #with warnings.catch_warnings():
                    #    warnings.filterwarnings(action='ignore', category=RuntimeWarning, message='divide by zero encountered in double_scalars')
                    denom=(TTime[i+trcgath,]-TTime[(i-trcgath),]) # if survey is perfectly symmetrical and 1d, this could be zero
                    if denom==0:
                        TTime[i+trcgath,] = TTime[i+trcgath,]+.00001 # add a tiny amount of time to avoid divide by zero                       
                    Slowx[k,] = 1/((SrcX[i+trcgath,]-SrcX[(i-trcgath),])/(TTime[i+trcgath,]-TTime[(i-trcgath),]))
                    
    # temporary fix to getting accurate slowness for edge points
    # just deleted the edges
    Slowx = Slowx[1:-1,]
    Slowz = Slowz[1:-1,]

    # get the modeled slownesses    
    _,Sloxpani,Slozpani,_,_,_,_ = phase_slow(Vpz,Vsz, epsilon, delta)
    epsiso=0
    deltaiso=0
    _,Sloxpiso,Slozpiso,_,_,_,_ = phase_slow(Vpz,Vsz, epsiso, deltaiso)
    

    fig = plt.figure(figsize=(8,5))    
    gs = gridspec.GridSpec(1, 1,  wspace = .25)
    
    ax1 = plt.subplot(gs[0])    
    ax1.scatter(Slowx, Slowz, marker = "+",s=50,color = 'black',label = 'Measured Phase-Slowness')      
    ax1.plot(Sloxpani, Slozpani, label = 'Modeled Anisotropic Phase-Slowness')      
    ax1.plot(Sloxpiso, Slozpiso, label = 'Modeled Isotropic Phase-Slowness',linestyle='--')      
    #ax1.set_xlim(0, f2*2)       # subjective choice
    ax1.set_xlabel('Horizontal Slowness')    
    ax1.set_ylabel('Vertical Slowness')    
    ax1.grid(True)    
    ax1.legend(loc='center',borderaxespad=0, fontsize = 8)        
    ax1.set_title('Measured and Modeled Compressional Phase Slowness',fontsize=12)
    ax1.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\aniso_phcalc_epsilon%s%%_delta%s%%.png' 
        %(epsilon*100,delta*100), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)                  
    plt.show()
    
    
def group_vel(Vpz,Vsz, epsilon, delta,dep):
    ''' This formulation for group velocity comes from Western Geco internal
    anisotropy class.
    Authors: John Markert and Paul Fowler, Western Geco, 2005
    '''
    
    import numpy as np
    from math import sqrt, sin, cos, degrees, tan, atan
    
    # get the horizontal and vertical phase velocity components, using Thomsen's 
    # parameterization
    Vpx,Vpn,Vsx,Vsn = aniso_params(Vpz,Vsz, epsilon, delta)

    # get P and Sv phase velocity    
    Vp, Vs, theta= phase_vel(Vpz, Vsz, epsilon, delta)
     
    Tp=2*(dep/Vpz) # two way time (P) to a fixed depth
    
     # create numpy arrays initialized to zero values
    term0,term1,term2,term3,term4,Vp_grp = [np.zeros(shape = theta.shape) for j in range(6)]
    term3_s,term4_s,Vs_grp = [np.zeros(shape = theta.shape) for i in range(3)]
    Vp_grp_angr, Vp_grp_angd, Vs_grp_angr, Vs_grp_angd = [np.zeros(shape = theta.shape) for i in range(4)]
    Xp_grp, Zp_grp,Xs_grp, Zs_grp = [np.zeros(shape = theta.shape) for i in range(4)]
    Xp_grp_refl,Zp_grp_refl,Xs_grp_refl,Zs_grp_refl = [np.zeros(shape = theta.shape) for i in range(4)]
    
    # group velocity can be calculated from phase velocity and it's components
    for i in range(0,theta.shape[0]):
        term0[i,] = (Vpx**2-Vpz**2)*sin(2*theta[i,])
        term1[i,] = (Vpx**2-Vs[1,]**2)*sin(theta[i,])**2+(Vpz**2-Vsz**2)*cos(theta[i,])**2        
        term2[i,] = (Vpz**2-Vsz**2)*(Vpn**2-Vpx**2)
        term3[i,] = term0[i,]+1/(sqrt(term1[i,]**2+term2[i,]*(sin(2*theta[i,])**2)))*\
        (term0[i,]*(Vpx**2-Vpz**2)*sin(2*theta[i,])+term0[i,]*sin(4*theta[i,]))    
       
        term4[i,] = (1/(4*Vp[i,]))*term3[i,]
        Vp_grp[i,] = sqrt(Vp[i,]**2+term4[i,]**2)

        term3_s[i,] = term0[i,]-1/(sqrt(term1[i,]**2+term2[i,]*(sin(2*theta[i,])**2)))*\
        (term1[i,]*(Vpx**2-Vpz**2)*sin(2*theta[i,])+term2[i,]*sin(4*theta[i,]))    
        
        term4_s[i,] = (1/(4*Vs[i,]))*term3_s[i,]    
        Vs_grp[i,] = sqrt(Vs[i,]**2+term4_s[i,]**2)
        
        # calculate the group angle from phase velocity and phase angle
        Vp_grp_angr[i,] = theta[i,]+atan((1/Vp[i,]) * term4[i,])
        Vs_grp_angr[i,] = theta[i,]+atan((1/Vs[i,]) * term4_s[i,])

        # calculate wavefront x,z compnents
        Xp_grp[i,] = Vp_grp[i,]*Tp*sin(Vp_grp_angr[i,])
        Zp_grp[i,] = Vp_grp[i,]*Tp*cos(Vp_grp_angr[i,])

        Xs_grp[i,] = Vs_grp[i,]*Tp*sin(Vs_grp_angr[i,])
        Zs_grp[i,] = Vs_grp[i,]*Tp*cos(Vs_grp_angr[i,])
        
        # to get the normal moveout curve at fixed reflector depth
        Xp_grp_refl[i,]=dep*tan(Vp_grp_angr[i,])
        Zp_grp_refl[i,]=sqrt((4*dep**2)+Xp_grp_refl[i,]**2)/Vp_grp[i,]
        
        Xs_grp_refl[i,]=dep*tan(Vs_grp_angr[i,])
        Zs_grp_refl[i,]=sqrt((4*dep**2)+Xs_grp_refl[i,]**2)/Vs_grp[i,]

    # convert group angle in radians to group angle in degrees - in case    
    Vp_grp_angd=np.degrees(Vp_grp_angr)
    Vs_grp_angd=np.degrees(Vs_grp_angr)
   
    return   Xp_grp,  Zp_grp, Xs_grp, Zs_grp


def wavefronts(Vpz, Vsz, epsilon, delta,dep,save):
    ''' Plot wave fronts at a constant depth for anisotropic and isotropic models
    
    '''    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec    

    
    # calculate wavefronts in depth and offset for anisotropic model
    Xp_grp, Zp_grp, Xs_grp, Zs_grp = group_vel(Vpz, Vsz, epsilon, delta,dep)
    
    # calculate wavefronts in depth and offset for isotropic model
    epsilon_iso=0
    delta_iso=0
    Xp_grp_iso, Zp_grp_iso, Xs_grp_iso, Zs_grp_iso = group_vel(Vpz, Vsz, epsilon_iso, delta_iso,dep)
    
    '''
    # make some QC tables
    from tabulate import tabulate
    pdat0 = np.vstack((term0, term1, term2, term3,term4, Vp_grp)).T    
    headers0 = ["term0","term1","term2", "term3", "term4", "Vp_grp"]    
    table0 = tabulate(pdat0, headers0, tablefmt="fancy_grid")
    print(table0)
    
    pdat1 = np.vstack((term0, term1, term2, term3_s,term4_s, Vs_grp)).T    
    headers1 = ["term0","term1","term2", "term3_s", "term4_s", "Vs_grp"]    
    table1 = tabulate(pdat1, headers1, tablefmt="fancy_grid")
    print(table1)
    
    pdat2 = np.vstack((Vp_grp_angr, Vs_grp_angr, Xp_grp, Zp_grp,Xs_grp, Zs_grp)).T    
    headers2 = ["Vp_grp_angr","Vs_grp_angr","Xp_grp", "Zp_grp", "Xs_grp", "Zs_grp"]    
    table2 = tabulate(pdat2, headers2, tablefmt="fancy_grid")
    print(table2)
    ''' 
    fig = plt.figure(figsize=(10,5))    
    gs = gridspec.GridSpec(1, 1,  wspace = .25)
    
    ax1 = plt.subplot(gs[0])    
    line1, = ax1.plot(Xp_grp, Zp_grp,  label = 'Anisotropic P ')      
    line2, = ax1.plot(Xs_grp, Zs_grp, label = 'Anisotropic Sv ')
    line3, = ax1.plot(Xp_grp_iso, Zp_grp_iso,  label = 'Isotropic P ',linestyle='--')      
    line4, = ax1.plot(Xs_grp_iso, Zs_grp_iso, label = 'Isotropic Sv ',linestyle='--')

    # Create a legend for the anisotropic curves.
    first_legend = plt.legend(handles=[line1,line2], loc='lower right')
    # Add the legend manually to the current Axes.
    leg1 = plt.gca().add_artist(first_legend)

    # Create another legend for the isotropic curves.
    leg2 = plt.legend(handles=[line3,line4], loc='lower left')

    ax1.set_ylim(max(Zp_grp), min(Zp_grp))  # flip y axis so depth increases from top to bottom
    ax1.set_xlabel('Horizontal Offset')    
    ax1.set_ylabel('Vertical Depth')    
    ax1.grid(True)    
#    ax1.legend(loc='center',borderaxespad=0, fontsize = 8)        
    ax1.set_title('P and S Wavefronts for Fixed Travel Time',fontsize=12)
    ax1.text(.5,.15,"Vpz = %s, Vsz = %s"%(Vpz,Vsz),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(.5,.1,"Epsilon = %s, Delta = %s"%(epsilon,delta),
                  horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    plt.show()
    
    DPI = 200    
    if (save =='Y') or (save =='y'):        
        fig.savefig('data\\wavefronts_epsilon%s%%_delta%s%%.png' 
        %(epsilon*100,delta*100), dpi=DPI, bbox_inches = 'tight', pad_inches = .1)                  

