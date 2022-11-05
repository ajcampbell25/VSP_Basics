def model_logs(vp_d,vs_d,rho_d,eps_d,delta_d,dt,dz, **kwargs):
   
    ''' make time and depth indexed logs of the model properties
    
    Inputs :
    vp_d : 3 element array of p-wave velocites m/s units
    vs_d : 3 element array of s-wave velocites m/s units 
    rho_d : 3 element array of density g/cc units
    
    mod_z : thickness of model, excluding the middle layer
    thickness : thickness of the middle layer
    dt : time sample rate in seconds
    dz = sample rate in depth
    
    Returns:
    tlogs : time indexed model logs
    timsamps : number of time samples for each layer
    tboundaries : cumulative TWT at each boundary
    zlogs : depth indexed model logs
    '''
    import numpy as np
    import procvsp.avaprocess as avaproc
    import plotvsp.avalogs as avalog

    mod_z=kwargs['mod_z'] # total thickness of 2 bounding layers
    thickness=kwargs['thickness']  # thickness of middle layer
    
    ################ convert model to time ##################
    # get p-wave times across layers and cumulative sample numbers

    timbounds, timsamps=avaproc.layer_times(vp_d,dt,**kwargs)
    numtsamp = int(np.sum(timsamps))    
#    print (' \nTime logs :')
#    print (' timbounds, timsamps :',timbounds, timsamps)

    vp_t=np.zeros(shape=int(numtsamp))
    vp_t[0:int(timsamps[0])+1]=vp_d[0]
    vp_t[int(timsamps[0]):int(timsamps[1])]=vp_d[1]          
    vp_t[int(timsamps[1]):]=vp_d[2]

    vs_t=np.zeros(shape=int(numtsamp))
    vs_t[0:int(timsamps[0])]=vs_d[0]
    vs_t[int(timsamps[0]):int(timsamps[1])]=vs_d[1]          
    vs_t[int(timsamps[1]):]=vs_d[2]

    rho_t=np.zeros(shape=int(numtsamp))
    rho_t[0:int(timsamps[0])]=rho_d[0]
    rho_t[int(timsamps[0]):int(timsamps[1])]=rho_d[1]          
    rho_t[int(timsamps[1]):]=rho_d[2]
    
    eps_t=np.zeros(shape=int(numtsamp))
    eps_t[0:int(timsamps[0])]=eps_d[0]
    eps_t[int(timsamps[0]):int(timsamps[1])]=eps_d[1]          
    eps_t[int(timsamps[1]):]=eps_d[2]    
    
    delta_t=np.zeros(shape=int(numtsamp))
    delta_t[0:int(timsamps[0])]=delta_d[0]
    delta_t[int(timsamps[0]):int(timsamps[1])]=delta_d[1]          
    delta_t[int(timsamps[1]):]=delta_d[2]        
    
    tlogs = np.vstack((vp_t,vs_t,rho_t, eps_t,delta_t )).T
    
    t1=timbounds[0]
    t2=t1+timbounds[1]
    t3=t2+timbounds[2]
    tboundaries=[t1,t2,t3]

    zlogs = avalog.depth_logs(vp_d, vs_d, rho_d,eps_d,delta_d,dz,**kwargs)

    return tlogs, timsamps, tboundaries, zlogs

def depth_logs(vp_d,vs_d,rho_d,eps_d,delta_d,dz, **kwargs):
    
    ''' make depth indexed logs of the model properties

    Inputs:
    vp_d : 3 element array of p-wave velocites m/s units
    vs_d : 3 element array of s-wave velocites m/s units 
    rho_d : 3 element array of density g/cc units
    
    mod_z : thickness of model, excluding the middle layer
    thickness : thickness of the middle layer
    dz : depth sample rate in meters
    '''
    import numpy as np
    
    mod_z=kwargs['mod_z'] # total thickness of 2 bounding layers
    thickness=kwargs['thickness']  # thickness of middle layer

    # the properties now will generate a log
    vp1 = np.zeros(mod_z//(2*dz)) + vp_d[0] #m
    vs1 = np.zeros(mod_z//(2*dz)) + vs_d[0]
    rho1 = np.zeros(mod_z//(2*dz)) + rho_d[0] #g/cc
    eps1 = np.zeros(mod_z//(2*dz)) + eps_d[0] #g/cc
    del1 = np.zeros(mod_z//(2*dz)) + delta_d[0] #g/cc
 
    vp2 = np.zeros(thickness//dz) + vp_d[1]
    vs2 = np.zeros(thickness//dz) + vs_d[1] #m
    rho2 = np.zeros(thickness//dz) + rho_d[1] #g/cc
    eps2 = np.zeros(thickness//dz) + eps_d[1] #g/cc
    del2 = np.zeros(thickness//dz) + delta_d[1] #g/cc

    vp3 = np.zeros(mod_z//(2*dz)) + vp_d[2]
    vs3 = np.zeros(mod_z//(2*dz)) + vs_d[2] #m
    rho3 = np.zeros(mod_z//(2*dz)) + rho_d[2] #g/cc
    eps3 = np.zeros(mod_z//(2*dz)) + eps_d[2] #g/cc
    del3 = np.zeros(mod_z//(2*dz)) + delta_d[2] #g/cc
    # Concatenate vp1 then vsp2 then another vp1 to get 3 layer model
    # 
    vp_z=np.concatenate((vp1,vp2,vp3))
    vs_z=np.concatenate((vs1,vs2,vs3))
    rho_z=np.concatenate((rho1,rho2,rho3))
    eps_z=np.concatenate((eps1,eps2,eps3))
    del_z=np.concatenate((del1,del2,del3))
    
    zlogs = np.vstack((vp_z,vs_z,rho_z,eps_z,del_z)).T
    
#    print (' \nDepth logs :')
#    print (' depth logs shape :',zlogs.shape)

    return zlogs