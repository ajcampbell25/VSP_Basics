def scattering_matrix(vp1, vs1, vp0, vs0, rho1, rho0, theta1):
    '''
    Full Zoeppritz solution, considered the definitive solution.
    Calculates the angle dependent p-wave reflectivity of an interface
    between two mediums.

    Written by: Wes Hamlyn

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.

    :param theta1: A scalar  [degrees].

    :returns: a 4x4 array representing the scattering matrix
                at the incident angle theta1.
                
    The scattering matrix has the form: 
                        ([['PdPu', 'SdPu', 'PuPu', 'SuPu'],
                         ['PdSu', 'SdSu', 'PuSu', 'SuSu'],
                         ['PdPd', 'SdPd', 'PuPd', 'SuPd'],
                         ['PdSd', 'SdSd', 'PuSd', 'SuSd']])

    '''
    import numpy as np
    
    # Make sure theta1 is an array
    # theta1_real = np.radians(np.array(theta1))
    # Complex values allow post-critical calculations
    theta1 = np.radians(theta1).astype(complex)
    
    if theta1.size == 1:
        theta1 = np.expand_dims(theta1, axis=1)

    # Set the ray paramter, p
    p = np.sin(theta1) / vp1  # ray parameter
    #print (' theta1 :',theta1, ' p :',p)

    # Calculate reflection & transmission angles for Zoeppritz
    theta2 = np.arcsin(p * vp0)      # Trans. angle of P-wave
    phi1 = np.arcsin(p * vs1)     # Refl. angle of converted S-wave
    phi2 = np.arcsin(p * vs0)      # Trans. angle of converted S-wave
    #print (' theta2 :',theta2, ' p * vp0 :',p * vp0)

    # Matrix form of Zoeppritz Equations... M & N are matricies
    M = np.array([[-np.sin(theta1), -np.cos(phi1), np.sin(theta2), np.cos(phi2)],
                  [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
                  [2 * rho1 * vs1 * np.sin(phi1) * np.cos(theta1),
                   rho1 * vs1 * (1 - 2 * np.sin(phi1) ** 2),
                   2 * rho0 * vs0 * np.sin(phi2) * np.cos(theta2),
                   rho0 * vs0 * (1 - 2 * np.sin(phi2) ** 2)],
                  [-rho1 * vp1 * (1 - 2 * np.sin(phi1) ** 2),
                   rho1 * vs1 * np.sin(2 * phi1),
                   rho0 * vp0 * (1 - 2 * np.sin(phi2) ** 2),
                   -rho0 * vs0 * np.sin(2 * phi2)]], dtype='complex')#dtype='float')

    N = np.array([[np.sin(theta1), np.cos(phi1), -np.sin(theta2), -np.cos(phi2)],
                  [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
                  [2 * rho1 * vs1 * np.sin(phi1) * np.cos(theta1),
                   rho1 * vs1 * (1 - 2 * np.sin(phi1) ** 2),
                   2 * rho0 * vs0 * np.sin(phi2) * np.cos(theta2),
                   rho0 * vs0 * (1 - 2 * np.sin(phi2) ** 2)],
                  [rho1 * vp1 * (1 - 2 * np.sin(phi1) ** 2),
                   -rho1 * vs1 * np.sin(2 * phi1),
                   - rho0 * vp0 * (1 - 2 * np.sin(phi2) ** 2),
                   rho0 * vs0 * np.sin(2 * phi2)]], dtype='complex')#dtype='float')

    zoep = np.zeros((4, 4, M.shape[-1]))
    for i in range(M.shape[-1]):
        Mi = M[..., i]
        Ni = N[..., i]
        dt = np.dot(np.linalg.inv(Mi), Ni)
        zoep[..., i] = dt.real # get the real part of the complex numbers
    #print (' dt :',dt)
    return zoep

def shuey(vp1, vs1, vp2, vs2, rho1, rho2, theta1):
    """
    Computes the P-wave reflectivity with Shuey (1985) 2 and 3 terms for a 
    two-layerd model.
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    
    Parameters
    ----------
    vp1 : array
        P-wave in the upper layer.
    vs1 : array
        S-wave in the upper layer.
    rho1 : array
        Density in the upper layer.        
    vp2 : array
        P-wave in the lower layer.        
    vs2 : array
        S-wave in the lower layer.      
    rho2 : array
        Density in the lower layer.        
    theta1 : array
        Angles of incidence.

    Returns
    -------
    R0 : array
        Intercept.
    G : array
        Gradient.        
    R2 : array
        Reflection coefficient for the 2-term approximation.
    R3 : array
        Reflection coefficient for the 3-term approximation.   
    """    
    import numpy as np
    
    theta1 = np.radians(theta1)    
    
    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = (vp1+vp2)/2
    vs  = (vs1+vs2)/2
    rho = (rho1+rho2)/2   
    
    R0 = 0.5*(dvp/vp + drho/rho)
    G  = 0.5*(dvp/vp) - 2*(vs**2/vp**2)*(drho/rho+2*(dvs/vs))
    F =  0.5*(dvp/vp)
    
    R2 = R0 + G*np.sin(theta1)**2
    R3 = R0 + G*np.sin(theta1)**2 + F*(np.tan(theta1)**2-np.sin(theta1)**2)
    
    return R0,G,R2, R3

def ruger(vp1, vs1, vp2, vs2, rho1, rho2,d1,  d2, e1,  e2, theta1):

    """
    Coded by Alessandro Amato del Monte and (c) 2016 by him
    https://github.com/aadm/avo_explorer/blob/master/avo_explorer_v2.ipynb
    
    Rüger, A., 1997, P -wave reflection coefficients for transversely
    isotropic models with vertical and horizontal axis of symmetry:
    Geophysics, v. 62, no. 3, p. 713–722.
    
    Provide Vp, Vs, rho, delta, epsilon for the upper and lower intervals,
    and theta, the incidence angle.
    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.
    :param d1: Thomsen's delta of the upper medium.
    :param e1: Thomsen's epsilon of the upper medium.
    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.
    :param d0: Thomsen's delta of the lower medium.
    :param e0: Thomsen's epsilon of the lower medium.
    :param theta: A scalar [degrees].
    :returns: anisotropic reflectivity.
    """
#    print (' d1 :',d1, 'd2 :',d2)
#    print (' e1 :',e1, 'e2 :',e2)
    
    import numpy as np

    #a = np.radians(theta1).astype(complex)
    a = np.radians(theta1)

    vp = np.mean([vp1, vp2])
    vs = np.mean([vs1, vs2])
    z = np.mean([vp1*rho1, vp2*rho2])
    g = np.mean([rho1*vs1**2, rho2*vs2**2])
    dvp = vp2-vp1
    z2, z1 = vp2*rho2, vp1*rho1
    dz = z2-z1
    dg = rho2*vs2**2 - rho1*vs1**2
    dd = d2-d1
    de = e2-e1
    A = 0.5*(dz/z)
    B = 0.5*(dvp/vp - (2*vs/vp)**2 * (dg/g) + dd) * np.sin(a)**2
    C = 0.5*(dvp/vp + de) * np.sin(a)**2 * np.tan(a)**2
    R = A+B+C

    return R



def calc_rc( vp, vs, rho,eps,delta,eqn, maxangle ):#theta1_samp
    """
    Computes the reflection coefficients for given incidence angles in a three layer model.

    Inputs:
    max_angle: maximum incidence angle
    vp: Layer P-wave velocities
    vs: Layer S-wave velocities
    rho: Density of layers
    eqn: Choice of equations - Zoeppritz, Shuey or Rutger

    return:
        theta1_samp : n-samples of incidence angles,
        rc_all : Reflection coefficient at layer1 / layer2,
                and at layer2 / layer3
    """
    import numpy as np
    
    # make an array of incidence angles
    theta1_samp = np.arange(0,maxangle+1,1)
    
#    print (' theta1_samp:',theta1_samp)

    # Calculate refl coeff for top of layer, then bottom of layer
    # rc_1 - refl coeff for top boundary of layer
    # rc_2 - refl coeff for lower boundary of layer
    if eqn=='Zoeppritz':
        zprc_1 = scattering_matrix(vp[0], vs[0], vp[1], vs[1], rho[0], rho[1], theta1_samp)
        zprc_2 = scattering_matrix(vp[1], vs[1], vp[2], vs[2], rho[1], rho[2], theta1_samp)
    
        # isolate the p down, p up refelctivity
        rc_1_pp=zprc_1[0,0,:]
        rc_2_pp=zprc_2[0,0,:]

        # put the 1d reflectivity arrays into a single 2d array
        rc_all=np.vstack((rc_1_pp,rc_2_pp)).T
        
        shuey_attr=np.zeros((4)) # dummy file
        
        ang0 = np.sin(np.radians(theta1_samp))**2

    if eqn=='3-Term-Shuey'or eqn=='2-Term-Shuey':
        
        Int1,Grad1,prc_1_2trm,prc_1_3trm = shuey(vp[0], vs[0], vp[1], vs[1], rho[0], rho[1], theta1_samp)
        Int2,Grad2,prc_2_2trm,prc_2_3trm = shuey(vp[1], vs[1], vp[2], vs[2], rho[1], rho[2], theta1_samp)

        shuey_attr=[Int1,Grad1,Int2,Grad2]
        # put the 1d reflectivity arrays into a single 2d array 
        # 2 layer 1 columns, and 2 layer 2 columns
        rc_all=np.vstack((prc_1_2trm,prc_1_3trm,prc_2_2trm,prc_2_3trm)).T

    if eqn=='Ruger':
        prc_1 = ruger(vp[0], vs[0], vp[1], vs[1], rho[0], rho[1],delta[0], delta[1],eps[0],eps[1], theta1_samp)
        prc_2 = ruger(vp[1], vs[1], vp[2], vs[2], rho[1], rho[2],delta[1], delta[2],eps[1],eps[2], theta1_samp)

        rc_all=np.vstack((prc_1,prc_2)).T
        
        shuey_attr=np.zeros((4)) # dummy file

        ang0 = np.sin(np.radians(theta1_samp))**2
 
    return theta1_samp, rc_all,shuey_attr
    
def layer_times(vp_d,dt,**kwargs):
    ''' Input the model properties Vp, Vs and the sample rate
    Calculate the travel time through each layer, and the 
    number of samples for th times
    Returns:
    tboundaries : time in each layer
    tsampnums : rnning total of samples through each layer
    
    '''    
    
    mod_z=kwargs['mod_z'] # total thickness of 2 bounding layers
    thickness=kwargs['thickness']  # thickness of middle layer
    
    # get p-wave times across layers
    t1=round(2*(mod_z//2)/vp_d[0],3)
    t2=round(2*thickness/vp_d[1],3) # time across layer
    t3=round(2*(mod_z//2)/vp_d[2],3)

    t_totall= t1+t2+t3
    # calculate number of samples needed for each layer (in time)
    l1_numtsamp = int(t1/dt) + 1
    l2_numtsamp = int(t2/dt) + l1_numtsamp
    l3_numtsamp = int(t3/dt) + l2_numtsamp
    numtsamp = l1_numtsamp+l2_numtsamp+l3_numtsamp
    
    interval_tim=[t1,t2,t3]
    tsampnums = [l1_numtsamp,l2_numtsamp,l3_numtsamp]
    
    #print (' \nlayer_times : ')
    #print (' interval_tim :',interval_tim)
    #print (' tsampnums : ',tsampnums)
    
    return interval_tim,tsampnums
    
def model_rc(vp_d,rc_all, angle,dt,**kwargs):
    
    ''' Create a 2d arrray of traces with reflectivity from
    the top and bottom of the middle layer of a 3 layer model
    
    The reflectivity array is then convolved with a Ricker wavelet
    
    Inputs :
    vp_d : 3 element array of p-wave velocites for depth->time
    rc1 : Reflection Coefficient (RC) for top of layer
    rc2 : Reflection Coefficient (RC) for bottom of layer
    maxangle : maximum angle of incidence to model 
    mod_z : total thickness of bounding layers in meters
    thickness : thickness of the middle layer in meters
    dt : time sample rate in seconds
    timlength : length of Ricker wavelet in seconds
    flimits: central frequency of Ricker wavelet
    '''
    import numpy as np
    import scipy.signal as sig
    import procvsp.avawavelets as avawav
    import procvsp.avaprocess as avaproc
    
    #print (' rc_all.shape :',rc_all.shape)
    
    mod_z=kwargs['mod_z'] # total thickness of 2 bounding layers
    thickness=kwargs['thickness']  # thickness of middle layer
    timlength=kwargs['wave_len']   # in seconds, preferably a power of 2    
    flimits=kwargs['wave_freq']    # Central frequency in hertz

    # make an empty array for the combined the reflectivities
    row = int(mod_z+thickness)
    col = angle.shape[0]
    mod_rc = np.zeros((row,col))
    
    anglegather=np.zeros((4,mod_z+thickness,angle.shape[0]))
    
    # Generate a ricker wavelet. An even number of samples means the peak is not 
    # at the middle sample. Generates a one sample shift after convolution
    # with reflectivity, which needs to be corrected for        
    twav,rickwave = avawav.ricker(timlength, dt, flimits)        

    ################ convert model to time and convolve ##################
    # get p-wave times across layers

    timbounds, timsamps=avaproc.layer_times(vp_d,dt,**kwargs)
    numtsamp = int(np.sum(timsamps))
    
    #print ('\nmodel_rc_params : ')
    #print (' timbounds, timsamps : ',timbounds, timsamps)
    
    # get the number of samples to the two boundaries    
    lay1_nsamps = int(timsamps[0])
    lay2_nsamps = int(timsamps[1])
    #print (' lay1_nsamps, lay2_nsamps : ',lay1_nsamps, lay2_nsamps)
    
    row = int(numtsamp)
    col = angle.shape[0]
    tmod_rc = np.zeros((row,col)) # tmod is time model

    # Zoeppritz has 2 columns of RCs, one for each layer
    # Shuey has 4 columns, 2 for each layer (2-term and 3-term reflectivity)
    if rc_all.shape[1]==2: # for array of Zoeppritz RC
        for ang in range(0,np.max(angle)+1,1):
            tmod_rc[lay1_nsamps,ang]=rc_all[ang,0:1]
            tmod_rc[lay2_nsamps,ang]=rc_all[ang,1]
    elif rc_all.shape[1]==4:
        for ang in range(0,np.max(angle)+1,1):
            tmod_rc[lay1_nsamps,ang]=rc_all[ang,1]
            tmod_rc[lay2_nsamps,ang]=rc_all[ang,3]

    #tmod_conv= np.zeros((row,col))
    fft_conv = np.zeros((row,col))
    #print (' fft_conv.shape :',fft_conv.shape,' tmod_rc.shape :',tmod_rc.shape)
 
    for k in range(0,(tmod_rc.shape[1])):        
        #tmod_conv[:,k] = np.convolve(tmod_rc[:,k],rickwave[:-2], mode='same')    
        fft_conv[:,k] = sig.fftconvolve(tmod_rc[:,k], rickwave, mode='same')#[cc_start:cc_end]

    #print (' argmax tmod_rc :',tmod_rc.argmax(axis=0) )
    #print (' argmax tmod_conv before deleting sample :',tmod_conv.argmax(axis=0) )
    #print (' argmax fft_conv :',fft_conv.argmax(axis=0) )
    
    # Convolution always shifted by one sample if number of samples in 
    # wavelet is even. Delete the first sample of output 
    # tmod_conv = np.delete(tmod_conv,0,0) # delete the first sample in every trace
    fft_conv = np.delete(fft_conv,0,0) # delete the first sample in every trace

    #print (' fft_conv after removing first sample :',fft_conv.shape)
                
    return fft_conv
 
def grad_int(rc,angle):
    '''Linear least squares method is used for estimating INTERCEPT and GRADIENT coefficients
    from seismic traces
    '''
    import numpy as np

    # calculate gradient and intercept from the extracted amplitudes
    ang0 = np.sin(np.radians(angle))**2
    #print (' rc.shape, angle.shape, ang0.shape :',rc.shape, angle.shape, ang0.shape )

    Grad1,Int1 = np.polyfit(ang0,rc,1)
    
    return Grad1,Int1

    # the following section is a different way of getting grad, int
    # should try this one day if I figure out the math
    '''    
    #def avo_inv(rc_zoep, ntrc: int, top: list) -> dict:
    #"""
    # from https://github.com/TolaAbiodun/pyavo/blob/master/pyavo/seismodel/tuning_prestack.py
    #AVO inversion for INTERCEPT and GRADIENT from analytic and convolved reflectivity values.
    #Linear least squares method is used for estimating INTERCEPT and GRADIENT coefficients.
    #:param rc_zoep: Zoeppritz reflection coefficient values
    #:param ntrc: Number of traces
    #:param top: Convolved top reflectivity values
    
    #:return: Zoeppritz and Convolved reflectivities
    #"""
    Yzoep = np.array(rc_zoep[:, 0])
    Yzoep = Yzoep.reshape((ntrc, 1))

    Yconv = np.array(top)
    Yconv = Yconv.reshape((ntrc, 1))

    ones = np.ones(ntrc)
    ones = ones.reshape((ntrc, 1))

    sintheta2 = np.sin(np.radians(np.arange(0, ntrc))) ** 2
    sintheta2 = sintheta2.reshape((ntrc, 1))

    X = np.hstack((ones, sintheta2))

    #   ... matrix solution of normal equations
    Azoep = np.dot(np.dot(np.linalg.inv(np.dot(X.T, X)), X.T), Yzoep)
    Aconv = np.dot(np.dot(np.linalg.inv(np.dot(X.T, X)), X.T), Yconv)

    return {'Zoeppritz': (Azoep[0], Azoep[1]),
            'Convolved': (Aconv[0], Aconv[1])
            }
    '''