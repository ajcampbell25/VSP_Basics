def model_2d(vp,vs,rho,eps,delta,**kwargs):
    ''' Plot the 2D model for QC of shape 
    '''
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import matplotlib.patches as patches
    import numpy as np

    from matplotlib.patches import ConnectionPatch    

#    print("\u0332".join('\nAVO_offsets - Stats :'))
    
    mod_z=kwargs['mod_z']
    thick=kwargs['thickness']
    save=kwargs['png']

    ################## get model plotting values ###########################
    
    depth1=mod_z/2
    depth2=depth1+thick

    # set up a padding so nothing gets cut off in plot
    pad = (mod_z+thick)*.05
    xmin = 0-pad#sources[0]-pad
    xmax = (mod_z+thick)+pad
    
    # define layers as x,y of bottom left corner, delta x and delta y
    layer1=[xmin-pad, mod_z/2, (xmax-xmin)+2*pad, -1*mod_z/2]
    layer2=[xmin-pad, mod_z/2+thick,(xmax-xmin)+2*pad, -1*thick]
    layer3=[xmin-pad, mod_z,(xmax-xmin)+2*pad, -1*(mod_z/2-thick)]
    
    ################create a well track #####################################
    
    mod_bottom= depth1*2#+thick+depth1
    well_bottom = .75*mod_bottom

    fig = plt.figure(figsize=(13,6) ) # usually needs adjusting    
    gs = gridspec.GridSpec(1, 2,  width_ratios=[1, 1], hspace = .05)
    
    ax1 = plt.subplot(gs[0])    

    # draw arrows for maximum depth and middle layer thickness of model
    # arrow uses x0,y0,dx,dy
    y1=[0,mod_bottom]
    x1=[xmax*.66,xmax*.66]
    y2=[depth1,depth2]
    x2=[xmax*.33,xmax*.33]
    
    # Using FancyArrowPatch keeps arrows the same size when plot is re-scaled
    arrow1=patches.FancyArrowPatch((x1[0], y1[0]), (x1[1],y1[1]), arrowstyle='->', 
              mutation_scale=10, label = 'mod_z',fc='black', ec='black', zorder=20)#shrinkA=1,shrinkB=1,'->'

    arrow2=patches.FancyArrowPatch((x2[0], y2[0]), (x2[1],y2[1]), arrowstyle='->',
              mutation_scale=10, label = 'mod_z',fc='black', ec='black', zorder=20)#shrinkA=1,shrinkB=1,
    
    ax1.add_patch(arrow1)
    ax1.add_patch(arrow2)              
    
    # get slope of arrows to allow parallel annotation  
    dy1 = y1[1] - y1[0]
    dx1 = x1[1] - x1[0]
    angle1 = np.rad2deg(np.arctan2(dy1, dx1))
    
    # annotate arrows with transform_rotates_text to align text and line
    # text uses x, y for left, center or right side of text box
    # rotation angle used to rotate text
    ax1.text(x1[1]-2*pad, (y1[1]-thick-y1[0])/2, 'mod_z', ha='left', va='bottom',
         transform_rotates_text=True, rotation=angle1, rotation_mode='anchor', fontsize=20)
    ax1.text(x2[1]-pad, y2[0], 'thickness', ha='center', va='top',
         transform_rotates_text=True, rotation=angle1, rotation_mode='anchor', fontsize=20)
    
    ax1.set_ylabel('Depth (ft)')
    ax1.set_xlabel('Offset from Wellhead (ft)')

    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(mod_bottom,0)
    ax1.set_title('Simple 3 Layer Model')
    # plot the layers as rectangles - crude
    rect1 = patches.Rectangle((layer1[0], layer1[1]), layer1[2], layer1[3], 
                             linewidth=1, edgecolor='black', facecolor='brown', alpha=.5)
    rect2 = patches.Rectangle((layer2[0], layer2[1]), layer2[2], layer2[3], 
                             linewidth=1, edgecolor='black', facecolor='yellow', alpha=.5)
    rect3 = patches.Rectangle((layer3[0], layer3[1]), layer3[2], layer3[3], 
                             linewidth=1, edgecolor='black', facecolor='orange', alpha=.5)
    ax1.add_patch(rect1)
    ax1.add_patch(rect2)
    ax1.add_patch(rect3)
    
    ax1.text(1.1, 0.75, 'Layer 1: Vp=%s\nVs=%s\nRho=%s\nEpsilon=%s Delta=%s'%(vp[0],vs[0],rho[0], eps[0],delta[0]), 
            fontsize=12, transform=ax1.transAxes ,va='center')
    ax1.text(1.1, 0.5, 'Layer 2: Vp=%s\nVs=%s\nRho=%s\nEpsilon=%s Delta=%s'%(vp[1],vs[1],rho[1], eps[1],delta[1]), 
            fontsize=12, transform=ax1.transAxes, va='center')
    ax1.text(1.1, 0.25, 'Layer 2: Vp=%s\nVs=%s\nRho=%s\nEpsilon=%s Delta=%s'%(vp[2],vs[2],rho[2], eps[2],delta[2]), 
            fontsize=12, transform=ax1.transAxes, va='center' )

    # write out a png file
    DPI = 200    
    if (save =='Y') or (save =='y'):
        
        fig.savefig('C:\\Users\\acampbell45\\Documents\\Python_Scripts'
        '\\AVO\\models\\avo_model.png', dpi=DPI, bbox_inches = 'tight', 
                    pad_inches = .1)

    plt.show()
    
def plot_velmod_t(axvel,vp_t,vs_t,time,tbounds,dt):
    '''Plot velocity logs in time indexing
    Define one track
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    tmin = tbounds[0]
    tmax = tbounds[1]

    velmin = .9*(np.min(vs_t))
    velmax = 1.1*(np.max(vp_t))

    #print (' \nplot velmod_t :')
    #print (' time.shape :',time.shape)
    #print (' vp_t.shape :',vp_t.shape)

    vpvs=vp_t/vs_t

    # twin the y axis at the start so the label positions
    # can be change between ax1 and ax1a
    axvela = axvel.twiny()
    
    #   Plot vp and vs log curves in two-way time
    axvel.plot(vp_t / 1000, time, color='black', lw=2, label = 'Vp',
             drawstyle = 'steps-pre')
    axvel.plot(vs_t / 1000, time, color='red', lw=2, label = 'Vs',
             drawstyle = 'steps-pre')
    axvel.legend(loc='upper right',borderaxespad=0, fontsize = 7)    
    axvel.set_ylim((tmin, tmax))
    axvel.set_xlim(velmin/1000, velmax/1000)
    axvel.invert_yaxis()
    axvel.set_ylabel('TWT (sec)')
    axvel.xaxis.tick_top()
    axvel.xaxis.set_label_position('top')
    axvel.set_xlabel('Velocity\n (km/s)')

    axvel.grid()

    #plot vpvs ratio in two-way time
    axvela.plot(vpvs , time, color='blue', lw=2,linestyle=':', label = 'Vp/Vs Ratio',
              drawstyle = 'steps-pre')
    axvela.legend(loc='lower left',borderaxespad=0, fontsize = 7)    
    axvela.set_xlim(1.5, 3) # VPVS plot limits
    axvela.xaxis.tick_bottom()
    axvela.xaxis.set_label_position('bottom')
    axvela.xaxis.label.set_color('blue')
    axvela.set_xlabel('Vp/Vs ratio', color='blue')
    
def plot_density_t(axrho,rho_t,time,tbounds,dt):
    '''Plot density log in time indexing
    Define one track
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    tmin = tbounds[0]
    tmax = tbounds[1]

    rhomin = .95*(np.min(rho_t))
    rhomax = 1.05*(np.max(rho_t))
    
    axrho.plot(rho_t, time, color='green', lw=2,
             drawstyle = 'steps-pre')
    axrho.set_ylim((tmin, tmax))
    axrho.set_xlim(rhomin, rhomax)
    axrho.invert_yaxis()
    axrho.xaxis.tick_top()
    axrho.xaxis.set_label_position('top')
    axrho.set_xlabel('Density (g/cc)')
    axrho.set_yticklabels([])
    axrho.grid()
    
def plot_aniso_t(axani,eps_t,del_t,time,tbounds,dt):
    '''Plot epsilon and delta logs in time indexing
    Define one track (axis)
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    tmin = tbounds[0]
    tmax = tbounds[1]

    anisomin = (-1.5*np.abs((np.min(del_t))))-.01
    anisomax = 2*(np.max(eps_t))
    
    axania = axani.twiny()
    
    axani.plot(eps_t, time, color='brown', lw=2, label = 'epsilon',
             drawstyle = 'steps-pre')
    axani.legend(loc='lower right',borderaxespad=0, fontsize = 7)    

    axani.set_ylim((tmin, tmax))
    axani.set_xlim(anisomin, anisomax)
    axani.invert_yaxis()
    axani.xaxis.tick_top()
    axani.xaxis.set_label_position('top')
    axani.legend(loc='upper right',borderaxespad=0, fontsize = 7)    
    axani.set_xlabel('Epsilon')
    axani.set_yticklabels([])
    axani.grid()    

    axania.plot(del_t, time, color='orange', lw=2, label = 'delta',
             drawstyle = 'steps-pre')    
    axania.set_xlim(anisomin, anisomax)
    axania.xaxis.tick_bottom()
    axania.xaxis.set_label_position('bottom')
    axania.xaxis.label.set_color('orange')    
    axania.legend(loc='lower right',borderaxespad=0, fontsize = 7)    
    axania.set_xlabel('Delta')
        
    
def plot_model(zlogs, tlogs,tlayers, dz, dt, **kwargs):

    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import numpy as np
    
    vp_t=tlogs[:,0]
    vs_t=tlogs[:,1]
    rho_t=tlogs[:,2]
    eps_t=tlogs[:,3]
    del_t=tlogs[:,4]
    
    vp_d=zlogs[:,0]
    vs_d=zlogs[:,1]
    rho_d=zlogs[:,2]
    eps_d=zlogs[:,3]
    del_d=zlogs[:,4]

    vpvs=vp_t/vs_t
    vpvs_z=vp_d/vs_d

    time=np.arange(0,(tlogs.shape[0])*dt,dt)
    depth=np.arange(0,(zlogs.shape[0])*dz,dz)
    
    # Get a reasonable time window for plot based on time thickness
    # of middle layer
  
    tbounds = t_extents(tlayers)
    
    zmin = kwargs['depth_min']
    zmax = kwargs['depth_max']

    velmin = .5*(np.min(vs_t))
    velmax = 1.1*(np.max(vp_t))
    rhomin = .9*(np.min(rho_t))
    rhomax = 1.1*(np.max(rho_t))
    anisomin = (-1.5*np.abs((np.min(del_t))))-.01
    anisomax = 2*(np.max(eps_t))
    
    fig = plt.figure(figsize=(15,5))    
#    gs = gridspec.GridSpec(1, 3, width_ratios=[.5,.4,2], wspace = .05)
    gs1 = gridspec.GridSpec(1, 7,width_ratios =[.5,.25,.25,.1,.5,.25,.25], wspace = .1)
    gs2 = gridspec.GridSpec(1, 7, width_ratios=[.5,.25,.25,.1,.5,.25,.25], wspace = .1)

    ax1 = plt.subplot(gs1[0])
    ax1a = ax1.twiny()
    #   Plot vp and vs log curves in depth
    ax1.plot(vp_d / 1000, depth, color='black', lw=2, label = 'Vp',
             drawstyle = 'steps-pre')
    ax1.plot(vs_d / 1000, depth, color='red', lw=2, label = 'Vs',
             drawstyle = 'steps-pre')
    ax1.legend(loc='upper right',borderaxespad=0, fontsize = 7)    
    ax1.set_ylim((zmin, zmax))
    ax1.set_xlim(velmin/1000, velmax/1000)
    ax1.invert_yaxis()
    ax1.set_ylabel('Depth (m)')
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax1.set_xlabel('Velocity (km/s)')
    ax1.grid()

    #plot vpvs ratio in depth
    ax1a.plot(vpvs_z , depth, color='blue', lw=2,linestyle=':', label = 'Vp/Vs Ratio',
             drawstyle = 'steps-pre')
    ax1a.legend(loc='lower left',borderaxespad=0, fontsize = 7)    
    ax1a.set_xlim(1.5, 3) # VPVS plot limits
    ax1a.xaxis.tick_bottom()
    ax1a.xaxis.set_label_position('bottom')
    ax1a.xaxis.label.set_color('blue')
    ax1a.set_xlabel('Vp/Vs ratio', color='blue')
    #   Plot log curves in two-way time

    ax2 = plt.subplot(gs1[1])
    
    ax2.plot(rho_d, depth, color='green', lw=2,
             drawstyle = 'steps-pre')
    ax2.set_ylim((zmin, zmax))
    ax2.set_xlim(rhomin, rhomax)
    ax2.invert_yaxis()
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.set_xlabel('Density (g/cc)')
    ax2.set_yticklabels([])
    ax2.grid()    

    ax3 = plt.subplot(gs1[2])
    ax3a = ax3.twiny()
    
    ax3.plot(eps_d, depth, color='brown', lw=2, label = 'epsilon',
             drawstyle = 'steps-pre')
    ax3.legend(loc='lower right',borderaxespad=0, fontsize = 7)    

    ax3.set_ylim((zmin, zmax))
    ax3.set_xlim(anisomin, anisomax)
    ax3.invert_yaxis()
    ax3.xaxis.tick_top()
    ax3.xaxis.set_label_position('top')
    ax3.legend(loc='upper right',borderaxespad=0, fontsize = 7)    

    ax3.set_xlabel('Epsilon')
    ax3.set_yticklabels([])
    ax3.grid()    

    ax3a.plot(del_d, depth, color='orange', lw=2, label = 'delta',
             drawstyle = 'steps-pre')    
    ax3a.set_xlim(anisomin, anisomax)
    ax3a.xaxis.tick_bottom()
    ax3a.xaxis.set_label_position('bottom')
    ax3a.xaxis.label.set_color('orange')    
    ax3a.legend(loc='lower right',borderaxespad=0, fontsize = 7)    
    ax3a.set_xlabel('Delta')

   # Plot the velocity model in two-way time

    ax4 = plt.subplot(gs2[4])    
    plot_velmod_t(ax4,vp_t,vs_t,time,tbounds,dt)

    # Plot the density model in two-way time   
    ax5 = plt.subplot(gs2[5])
    plot_density_t(ax5,rho_t,time,tbounds,dt)  
    
    # Plot the anisotropy model in two-way time   
    ax6 = plt.subplot(gs2[6])
    plot_aniso_t(ax6,eps_t,del_t,time,tbounds,dt)    

def AVO_offsets(vp,vs,rho,**kwargs):
    ''' Plot the 2d model and overlay with straight rays from a
     fixed number of sources (11) to 1 receiver
     -10 dummy receivers are plotted for completeness
 
    Model width = model depth + padding
    '''
    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import matplotlib.patches as patches
    import numpy as np

    from matplotlib.patches import ConnectionPatch    

#    print("\u0332".join('\nAVO_offsets - Stats :'))
    
    max_off=kwargs['max_offset']
    mod_z=kwargs['mod_z']
    thick=kwargs['thickness']
    critang_p=kwargs['critical_angle_p']
    critang_s=kwargs['critical_angle_p']
    rcvz = mod_z/2

    ################### set up the source and receiver arrays #############
    
    numsrcs = 11
    sources=np.linspace(0,max_off,numsrcs)
    numrcvs = 11
    rcvs_z =np.linspace(mod_z/4,mod_z/2,numrcvs)
    rcvs_x=np.zeros(shape=rcvs_z.shape)    

    # calculate the incidence angle at the bottom receiver 
    # The bottom receiver is at the top of the altered layer 
    inc_ang=np.zeros(shape=sources.shape)
    inc_ang= np.degrees(np.arctan(sources/rcvz))

    ################## get model plotting values ###########################
    
    depth1=mod_z/2
    depth2=depth1+thick

    # set up a padding so nothing gets cut off in plot
    pad = (mod_z+thick)*.05
    xmin = sources[0]-pad
    xmax = (mod_z+thick)+pad
    
    # define layers as x,y of bottom left corner, delta x and delta y
    layer1=[xmin-pad, mod_z/2, (xmax-xmin)+2*pad, -1*mod_z/2]
    layer2=[xmin-pad, mod_z/2+thick,(xmax-xmin)+2*pad, -1*thick]
    layer3=[xmin-pad, mod_z,(xmax-xmin)+2*pad, -1*(mod_z/2-thick)]
    
    ################create a well track #####################################
    
    mod_bottom= depth1*2#+thick+depth1
    well_bottom = .75*mod_bottom

    fig = plt.figure(figsize=(14,6) ) # usually needs adjusting    
    gs = gridspec.GridSpec(1, 2,  width_ratios=[1, 1], hspace = .05)
    
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])
    
    ax1.text(0.25, .05, 'Cartoon by Allan Campbell',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax1.transAxes,
        color='black', fontsize=10)
    
    for i, trace in enumerate(sources[:,]):

        ax1.plot(sources[i], 0, marker=11,markersize=15,c='red',label = 'sources', zorder = 15)        
        ax1.arrow(sources[i], 0,  -1*sources[i],rcvz,width = 1,
              length_includes_head = 'False', head_width = 20, head_length = 20,  
              label = 'Source Offset',fc='blue', ec='blue') 

    ax1.set_ylabel('Depth (ft)')
    ax1.set_xlabel('Offset from Wellhead (ft)')

    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(mod_bottom,0)
    ax1.set_title('Walkaway VSP Survey Geometry')
    # plot the layers as rectangles - crude
    rect1 = patches.Rectangle((layer1[0], layer1[1]), layer1[2], layer1[3], 
                             linewidth=1, edgecolor='black', facecolor='brown', alpha=.5)
    rect2 = patches.Rectangle((layer2[0], layer2[1]), layer2[2], layer2[3], 
                             linewidth=1, edgecolor='black', facecolor='yellow', alpha=.5)
    rect3 = patches.Rectangle((layer3[0], layer3[1]), layer3[2], layer3[3], 
                             linewidth=1, edgecolor='black', facecolor='orange', alpha=.5)
    ax1.add_patch(rect1)
    ax1.add_patch(rect2)
    ax1.add_patch(rect3)
    # plot the well and receivers last so they are 'on top'
    ax1.scatter(rcvs_x, rcvs_z, marker='*',s=200,c='red',label = 'receivers', zorder = 20)
    ax1.vlines(x= 0, ymin =0, ymax = well_bottom,
               colors= 'black', linestyle = 'solid', linewidth = 5,zorder = 10 )


    ax2.plot(sources, inc_ang, marker='*',c='red',label = 'Incidence Angle')#, zorder = -1)
    ax2.hlines(y= critang_p, xmin =xmin, xmax = xmax,
               colors= 'black', linestyle = 'solid', linewidth = .3 )  
    ax2.set_ylabel('Incidence Angle (Degrees)')
    ax2.set_xlabel('Offset from Wellhead (ft)')
    ax2.set_title('Angle of Incidence for Walkaway VSP Survey')
    ax2.annotate('Critical Angle', xy=(xmin+pad, critang_p), ha='left', va='center')
    
    plt.show()
    
    
def plot_rc( zoep_rc,shuey_rc,angle,tlogs,tlayers,dt, algtyp,critangp):
    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import numpy as np
    
    vp_t=tlogs[:,0]
    vs_t=tlogs[:,1]
    rho_t=tlogs[:,2]

    vpvs=vp_t/vs_t

    time=np.arange(0,(tlogs.shape[0])*dt,dt)
    
    # Get a reasonable time window for plot based on time thickness
    # of middle layer

    tbounds = t_extents(tlayers)

    fig = plt.figure(figsize=(15,5))    
    gs1 = gridspec.GridSpec(1, 5, width_ratios=[.5,.25,.1,1,.25], wspace = .1)

   # Plot the velocity model in two-way time
    ax1 = plt.subplot(gs1[0])    
    plot_velmod_t(ax1,vp_t,vs_t,time,tbounds,dt)

    # Plot the density model in two-way time   
    ax2 = plt.subplot(gs1[1])
    plot_density_t(ax2,rho_t,time,tbounds,dt)  
    
    ax4 = plt.subplot(gs1[3])
    
    ax4.plot(angle,zoep_rc[:,0],'blue',linewidth=1, label='Zoepp. Top')
    ax4.plot(angle,zoep_rc[:,1],'red',linewidth=1, label='Zoepp. Bottom ')
    ax4.plot(angle,shuey_rc[:,1],'blue',linestyle = '--',linewidth=1, label='%s Top'%(algtyp[2]))
    ax4.plot(angle,shuey_rc[:,3],'red',linestyle = '--',linewidth=1, label='%s Bottom'%(algtyp[2]))
    ax4.plot(angle,shuey_rc[:,0],'blue',linestyle = ':',linewidth=1, label='%s Top'%(algtyp[1]))
    ax4.plot(angle,shuey_rc[:,2],'red',linestyle = ':',linewidth=1, label='%s Bottom'%(algtyp[1]))
    yrange=ax4.get_ylim()
    ax4.vlines(x= critangp, ymin=yrange[0],ymax=yrange[1],colors= 'black', linestyle = 'solid', linewidth = .3 )
    ax4.text(critangp,np.max(zoep_rc[:,0])/2 ,'Critical Angle', rotation =90, ha='left', va='center')
               
    ax4.legend(loc='best',borderaxespad=0, fontsize = 7)    
    ax4.set_xlabel('Angle ($\\theta$)')
    ax4.set_ylabel('Reflection Coefficient')
    ax4.grid()    
    ax4.set_title('AVA for %s and Shuey'%(algtyp[0]))
    plt.show()
    
def plot_rc_aniso( zoep_rc,ruger_rc,angle,tlogs,tlayers,dt, algtyp,critangp):
    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import numpy as np
    
    vp_t=tlogs[:,0]
    vs_t=tlogs[:,1]
    rho_t=tlogs[:,2]
    eps_t=tlogs[:,3]
    del_t=tlogs[:,4]

    vpvs=vp_t/vs_t

    time=np.arange(0,(tlogs.shape[0])*dt,dt)
    
    # Get a reasonable time window for plot based on time thickness
    # of middle layer
   
    tbounds = t_extents(tlayers)
    
    fig = plt.figure(figsize=(15,5))    
    gs1 = gridspec.GridSpec(1, 5, width_ratios=[.5,.25,.25,.1,1], wspace = .1)

    # Plot the velocity model in two-way time
    ax1 = plt.subplot(gs1[0])    
    plot_velmod_t(ax1,vp_t,vs_t,time,tbounds,dt)

    # Plot the density model in two-way time   
    ax2 = plt.subplot(gs1[1])
    plot_density_t(ax2,rho_t,time,tbounds,dt)  
 
    # Plot the anisotropy model in two-way time   
    ax3 = plt.subplot(gs1[2])
    plot_aniso_t(ax3,eps_t,del_t,time,tbounds,dt)
    
    ax4 = plt.subplot(gs1[4])
    
    ax4.plot(angle,zoep_rc[:,0],'blue',linewidth=1, label='Zoepp. Top')
    ax4.plot(angle,zoep_rc[:,1],'red',linewidth=1, label='Zoepp. Bottom ')
    ax4.plot(angle,ruger_rc[:,0],'blue',linestyle = '--',linewidth=1, label='Ruger Top')
    ax4.plot(angle,ruger_rc[:,1],'red',linestyle = '--',linewidth=1, label='Ruger Bottom ')
    yrange=ax4.get_ylim()
    ax4.vlines(x= critangp, ymin=yrange[0],ymax=yrange[1],colors= 'black', linestyle = 'solid', linewidth = .3 )
    ax4.text(critangp,np.max(zoep_rc[:,0])/2 ,'Critical Angle', rotation =90, ha='left', va='center')
    ax4.legend(loc='best',borderaxespad=0, fontsize = 7)    
    ax4.set_xlabel('Angle ($\\theta$)')
    ax4.set_ylabel('Reflection Coefficient')
    ax4.grid()    
    ax4.set_title('AVA for %s and %s'%(algtyp[0],algtyp[3]))
    plt.show()
    

def plot_grad_int( shuey_rc,shuey_attributes,angle,algo):
    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    import numpy as np
    import procvsp.avaprocess as avaproc

    y_sin2 = (np.sin(angle* np.pi / 180))**2

    # calculate gradient and intercept from the extracted amplitudes    
    grad1,int1=avaproc.grad_int(shuey_rc[:,0],angle)
    grad2,int2=avaproc.grad_int(shuey_rc[:,2],angle)
    print (' Top Grad1,Int1  :', grad1,int1 )
    print (' Bottom Grad2,Int2  :', grad2,int2 )
     
    fig = plt.figure(figsize=(15,5))    
    gs = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace = .25)

    ax1 = plt.subplot(gs[0])
    
    ax1.plot(angle,shuey_rc[:,1],'blue',linestyle = '--',linewidth=1, label='%s Top'%(algo[2]))
    ax1.plot(angle,shuey_rc[:,3],'red',linestyle = '--',linewidth=1, label='%s Bottom'%(algo[2]))
    ax1.plot(angle,shuey_rc[:,0],'blue',linestyle = ':',linewidth=1, label='%s Top'%(algo[1]))
    ax1.plot(angle,shuey_rc[:,2],'red',linestyle = ':',linewidth=1, label='%s Bottom'%(algo[1]))
    ax1.legend(loc='best',borderaxespad=0, fontsize = 7)    
    ax1.set_xlabel('Angle ($\\theta$)')
    ax1.set_ylabel('Reflection Coefficient')
    ax1.grid()    
    ax1.set_title('AVA for %s and %s'%(algo[1],algo[2]))
    
    ax2 = plt.subplot(gs[1])
    

    ax2.plot(y_sin2,shuey_rc[:,0],'blue',linestyle = ':',linewidth=1, label='%s Top'%(algo[1]))
    ax2.plot(y_sin2,shuey_rc[:,2],'red',linestyle = ':',linewidth=1, label='%s Bottom'%(algo[1]))
    ax2.legend(loc='best',borderaxespad=0, fontsize = 7)    
    ax2.set_xlabel('sin($\\theta$)^2')
    ax2.set_ylabel('Reflection Coefficient')
    ax2.grid()    
    ax2.set_title('Amplitude from %s Versus sin($\\theta$)^2'%(algo[1]))
    ax2.annotate('Gradient %0.2f Intercept %0.2f'%(grad1,int1), xy=(0, int1,), ha='left', va='center')
    ax2.annotate('Gradient %0.2f Intercept %0.2f'%(grad2,int2), xy=(0, int2,), ha='left', va='center')
    
    
    plt.show()
    
def plot_wig(axw,anglegather,angle,time,cmap3,normalize, scale_factor,tbounds):
    '''Set up a track (axis) for plotting wiggle traces
    in color
    '''

    import matplotlib.pyplot as plt
    from matplotlib  import colors,colorbar

    import numpy as np
    
    anglegatherplt=anglegather*scale_factor #just for the plot
    anglegatherplt[0]=anglegatherplt[0]*scale_factor # Twice just for the plot
    anglegatherplt[1]=anglegatherplt[1]*scale_factor #Twice just for the plot

    for i in range(len(angle)):
        # get the max and min amplitudes in trace for color bar matching
        al3 = np.max(anglegather[:,i])
        al4 = np.min(anglegather[:,i])
        
        axw.plot(i+anglegatherplt[:,i],time,'k',linewidth=1)

        fillx=axw.fill_betweenx(time,anglegatherplt[:,i]+i,i,
                            where=anglegatherplt[:,i]+i>i, linewidth=0,color=cmap3(normalize(al3)))        
        fillx2=axw.fill_betweenx(time,anglegatherplt[:,i]+i,i,
                            where=anglegatherplt[:,i]+i<i, linewidth=0,color=cmap3(normalize(al4)))

    axw.set_ylim(tbounds[0], tbounds[1])
    axw.invert_yaxis()
    axw.yaxis.grid()
    axw.yaxis.tick_right()
    axw.set_yticklabels([])
    axw.yaxis.set_label_position('right')    
    axw.set_xlabel('Angle ($\\theta$)')
    axw.set_ylabel('TWT (sec)')

   
def plot_gath_VI( anglegather,angle,tlogs,tlayers,dt,algo,**kwargs):
    ''' Plot the model and seismic gather
    
    Use:
    plotparams={
          'ang_min':0,
          'scaler':20,
          'time_lines':'y'} # draw time lines for boundaries on wiggles 'y' or 'n'
    avaplt.plot_gath_VI(traces_zoepp,angles,tim_logs,tim_layers,dt,algo[0],**plotparams)
    '''
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    from matplotlib  import colors, colorbar
    #from matplotlib.patches import PathPatch
    import numpy as np
    
    vp_t=tlogs[:-1,0]
    vs_t=tlogs[:-1,1]
    rho_t=tlogs[:-1,2]

    time=np.arange(0,(anglegather.shape[0])*dt,dt)

    # Get a reasonable time window for plot based on time thickness
    # of middle layer
    tbounds = t_extents(tlayers)

    scale_factor=kwargs['scaler']
    tlines=kwargs['time_lines']
 
    fig = plt.figure(figsize=(12,5))    
    # use 2 grids so that the plot spacing can be changed between the model
    # and the traces
    gs1 = gridspec.GridSpec(1, 3, width_ratios=[.5,.3,2], wspace = .1)
 
    # Plot the velocity model in two-way time
    ax1 = plt.subplot(gs1[0])    
    plot_velmod_t(ax1,vp_t,vs_t,time,tbounds,dt)

    # Plot the density model in two-way time   
    ax2 = plt.subplot(gs1[1])
    plot_density_t(ax2,rho_t,time,tbounds,dt)      

    # select  a color map
    cmap3 = plt.get_cmap('bwr')#'RdBu''seismic_r'

    # get the max amplitude in trace plot
    al1 = max(np.max(anglegather), abs(np.min(anglegather)))
    # get normalization to match the trace ampitudes to the colorbar values
    normalize = colors.Normalize(vmin=-1*al1, vmax=al1)

    ax4 = plt.subplot(gs1[2])
    plot_wig(ax4,anglegather,angle,time,cmap3,normalize,scale_factor,tbounds)
    ax4.axhline(tlayers[0]+1*dt, color='blue', lw=2, alpha=0.5)
    ax4.axhline(tlayers[1]+1*dt, color='red', lw=2, alpha=0.5)
    ax4.set_title('%s Synthetic Angle Gather'%(algo))

    # set the color bar location and size
    pad = 0.06    
    width = 0.02   
    pos2 = ax4.get_position()    
    axcol5 = fig.add_axes([pos2.xmax + pad, pos2.ymin, width, 0.9*(pos2.ymax-pos2.ymin) ])
    cb3 = colorbar.ColorbarBase(axcol5, label='Amplitude',cmap=cmap3, norm=normalize, orientation='vertical')

    plt.show()
    
def plot_gath_rc(zoep_rc, anglegather,angle,tlogs,tlayers,dt,algo,**kwargs):
    
    import matplotlib.pyplot as plt
    from matplotlib  import gridspec
    from matplotlib  import colors, colorbar
    #from matplotlib.patches import PathPatch
    import numpy as np

    import procvsp.avaprocess as avaproc
    
    vp_t=tlogs[:-1,0]
    vs_t=tlogs[:-1,1]
    rho_t=tlogs[:-1,2]
    #print (' \nplot gath rc :')
    #print (' anglegather.shape :',anglegather.shape)
    #print (' tlayers :',tlayers)
    #print (' vp_t.shape :',vp_t.shape)
    vpvs=vp_t/vs_t
    
    time=np.arange(0,(anglegather.shape[0])*dt,dt)
    
    # Get a reasonable time window for plot based on time thickness
    # of middle layer
    tbounds = t_extents(tlayers)

    scale_factor=kwargs['scaler']
    tlines=kwargs['time_lines']
    
    velmin = .8*(np.min(vs_t))
    velmax = 1.1*(np.max(vp_t))
    rhomin = .9*(np.min(rho_t))
    rhomax = 1.1*(np.max(rho_t))
    
    ######## calculate gradient and intercept from trace amplitudes ############ 
    
    # Read the number of samples to each boundary for
    # extracting amplitudes from the traces at the boundaries
    line1_samps = int(tlayers[0]/dt)+1
    line2_samps = int(tlayers[1]//dt)+1
    
    # extract the amplitudes at the boundaries
    wig_amp1=np.zeros(shape=angle.shape[0])
    wig_amp2=np.zeros(shape=angle.shape[0])
        
    for j in range(0,anglegather.shape[1]):
        wig_amp1[j]=anglegather[line1_samps,j]
        wig_amp2[j]=anglegather[line2_samps,j]
        
    # calculate gradient and intercept from the extracted amplitudes
    grad1,int1=avaproc.grad_int(wig_amp1,angle)
    grad2,int2=avaproc.grad_int(wig_amp2,angle)
    print (' Top Grad1,Int1  :', grad1,int1 )
    print (' Bottom Grad2,Int2  :', grad2,int2 )

    fig = plt.figure(figsize=(14,5))    

    gs1 = gridspec.GridSpec(1, 5, width_ratios=[.25,.15,1.3,.45,1.1], wspace = .075)

    # Plot the velocity model in two-way time
    ax1 = plt.subplot(gs1[0])    
    plot_velmod_t(ax1,vp_t,vs_t,time,tbounds,dt)

    # Plot the density model in two-way time   
    ax2 = plt.subplot(gs1[1])
    plot_density_t(ax2,rho_t,time,tbounds,dt)  

    # select  a color map
    cmap3 = plt.get_cmap('bwr')#'RdBu''seismic_r'

    # get the max amplitude in trace plot
    al1 = max(np.max(anglegather), abs(np.min(anglegather)))
    # get normalization to match the trace ampitudes to the colorbar values
    normalize = colors.Normalize(vmin=-1*al1, vmax=al1)

    # Plot the wiggle traces
    ax4 = plt.subplot(gs1[2])
    plot_wig(ax4,anglegather,angle,time,cmap3,normalize,scale_factor,tbounds)
    ax4.axhline(tlayers[0]+1*dt, color='blue', lw=2, alpha=0.5)
    ax4.axhline(tlayers[1]+1*dt, color='red', lw=2, alpha=0.5)
    ax4.set_title('%s Synthetic Angle Gather'%(algo))

    # set the color bar location and size
    pad = 0.01    
    width = 0.01    
    pos2 = ax4.get_position()    
    axcol4 = fig.add_axes([pos2.xmax + pad, pos2.ymin, width, 0.9*(pos2.ymax-pos2.ymin) ])
    cb4 = colorbar.ColorbarBase(axcol4, cmap=cmap3, norm=normalize, orientation='vertical')
    
    #Plot the AVA for synthetic traces and Zoeppritz calculation
    ax5=plt.subplot(gs1[4])

    ax5.plot(angle,zoep_rc[:,0],'blue',linewidth=1, label='Zoepp. Top')
    ax5.plot(angle,zoep_rc[:,1],'red',linewidth=1, label='Zoepp. Bottom ')
    ax5.plot(angle,wig_amp1,'blue',linewidth=1, linestyle = '--',label='Seismic Top')
    ax5.plot(angle,wig_amp2,'red',linewidth=1,linestyle = '--', label='Seismic Bottom ')
    ax5.legend(loc='best',borderaxespad=0, fontsize = 7)    
    ax5.set_xlabel('Angle ($\\theta$)')
    ax5.set_ylabel('Amplitude')
    ax5.grid()    
    ax5.set_title('AVA for Zoeppritz and Synthetic')    
    
    plt.show()
    
def t_extents(tlayers):    
    # Get a reasonable time window for plot based on time thickness
    # of middle layer

    import numpy as np

    tlayers_sec=np.array(tlayers)
    t_thick = tlayers_sec[1]-tlayers_sec[0]
    tmin= tlayers_sec[0]-(t_thick*1.5) 
    tmax= tlayers_sec[1]+(t_thick*1.5)
    tbounds=[tmin,tmax]

    return tbounds