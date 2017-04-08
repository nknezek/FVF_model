

##### Display and Save 1D Waves with comparison Spherical Harmonics for found Eigenvalues and vector field map ######
import matplotlib.cm as cm
from numpy import sin
from numpy import cos
import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
import matplotlib as mpl
from matplotlib import gridspec

colors = ['b','g','r','m','y','k','c']

def plot_1D(model,vec,val,ind, dir_name='./',title='1D Wave Plot'):
    E = model.E
    Nk = model.Nk
    Nl = model.Nl
    th = model.th[0,:]
    ## Plot Figures
    fig, axes = plt.subplots(3,1,figsize=(15,8), sharex=True, sharey=True)
    titles = ['Absolute Value','Real','Imaginary']
    for (ax,title,ind) in zip(axes,titles,range(len(axes))):
        ax.set_title(title)
        line = []
        for (var,color) in zip(model.model_variables,colors):
             out=model.get_variable(vec,var)
             if ind == 0:
                  line.append(ax.plot(th*180./np.pi,abs(out.T),color=color))
             elif ind ==1:
                  ax.plot(th*180./np.pi,out.T.real,color=color)
                  ax.grid()
             elif ind==2:
                  ax.plot(th*180./np.pi,out.T.imag,color=color)
                  ax.grid()
        if ind ==0:
             labels = ['ur','uth','uph','br','bth','bph','p']
             ax.legend([x[0] for x in line],labels,loc=0,ncol=4)
             ax.grid()

    plt.suptitle('MAC, Eigenvalue = {0:.5f}, m={1}, l={2}\n Nk={3}, Nl={4}, E={5:.2e})'.format(val,m,l,Nk,Nl,E), size=14)
    plt.savefig('./output/m={1}/MAC_Eig{0:.2f}j_m={1}_l={2}_Nk={3}_Nl={4}_E={5:.2e}.png'.format(val.imag,m,l,Nk,Nl,E))

def plot_1D_noB(model,vec,val,ind, dir_name='./',title='1D Wave Plot'):
    E = model.E
    Nk = model.Nk
    Nl = model.Nl
    th = model.th[0,:]
    ## Plot Figures
    fig, axes = plt.subplots(3,1,figsize=(15,8), sharex=True, sharey=True)
    axtitles = ['Absolute Value','Real','Imaginary']
    for (ax,axtitle,ind) in zip(axes,axtitles,range(len(axes))):
        ax.set_title(axtitle)
        line = []
        for (var,color) in zip(model.model_variables,colors):
             out = model.get_variable(vec,var)
             if ind == 0:
                  line.append(ax.plot(th*180./np.pi,abs(out.T),color=color))
             elif ind ==1:
                  ax.plot(th*180./np.pi,out.T.real,color=color)
                  ax.grid()
             elif ind==2:
                  ax.plot(th*180./np.pi,out.T.imag,color=color)
                  ax.grid()
        if ind ==0:
             labels = ['uth','uph','p']
             ax.legend([x[0] for x in line],labels,loc=0,ncol=4)
             ax.grid()

    fig.suptitle(title, fontsize=14)
    plt.savefig(dir_name+title+'.png')

def plot_mollyweide(model,vec,val,m,l,v_scale=1.0):
   from mpl_toolkits.basemap import Basemap

   E = model.E
   Nk = model.Nk
   Nl = model.Nl
   th = model.th[0,:]
   ## Calculate vector field and contour field for plotting with basemap
   ## Create full vector grid in theta and phi
   u_1D = model.get_variable(vec,'uph')[0][0]
   v_1D = model.get_variable(vec,'uth')[0][0]
   Nph = 2*Nl
   ph = np.linspace(-180.,180.-360./Nph,Nph)
   lon_grid, lat_grid = np.meshgrid(ph,th*180./np.pi-90.,)
   v = ((np.exp(1j*m*lon_grid*np.pi/180.).T*v_1D).T).real
   u = ((np.exp(1j*m*lon_grid*np.pi/180.).T*u_1D).T).real
   absu=u.real**2 + v.real**2
   Nvec = np.floor(Nl/20.)
   ### Plot Mollweide Projection
   plt.figure(figsize=(10,10))
   ## Set up map
   bmap = Basemap(projection='moll',lon_0=0.)
   bmap.drawparallels(np.arange(-90.,90.,15.))
   bmap.drawmeridians(np.arange(0.,360.,15.))
   ## Convert Coordinates to those used by basemap to plot
   lon,lat = bmap(lon_grid,lat_grid)
   bmap.contourf(lon,lat,absu,15,cmap=plt.cm.Reds,alpha=0.5)
   bmap.quiver(lon[::Nvec,::Nvec],lat[::Nvec,::Nvec],u[::Nvec,::Nvec],v[::Nvec,::Nvec], scale=v_scale)
   plt.title('MAC, Mollweide Projection Vector Field for m={0}, l={1}'.format(m,l))
   plt.savefig('./output/m={1}/MAC_MollweideVectorField_m={1}_l={2}.png'.format(val.imag,m,l))

def plot_B_obs(model, vec, m, dir_name='./', title='B-Perturbation at Core Surface'):
   from mpl_toolkits.basemap import Basemap
   Nl = model.Nl
   th = model.th[0,:]

   ## Plot Robinson vector field
   #### Display waves on a Spherical Map Projection
   projtype = 'robin'
   ## Create full vector grid in theta and phi
   Bobsth = model.get_variable(model.BobsMat.tocsr()*vec, 'ur')[-1,:]
   Nph = 2*Nl
   ph = np.linspace(-180.,180.-360./Nph,Nph)
   lon_grid, lat_grid = np.meshgrid(ph,th*180./np.pi-90.,)
   Bobs = (np.exp(1j*m*lon_grid*np.pi/180.).T*Bobsth).T

   ### Plot Robinson Projection
   plt.figure(figsize=(10,10))
   ## Set up map
   bmap = Basemap(projection=projtype,lon_0=0.)
   bmap.drawparallels(np.arange(-90.,90.,15.))
   bmap.drawmeridians(np.arange(0.,360.,15.))
   ## Convert Coordinates to those used by basemap to plot
   lon,lat = bmap(lon_grid,lat_grid)
   bmap.contourf(lon,lat,Bobs.real,15,cmap=plt.cm.RdBu,alpha=0.5)
   plt.title(title)
   plt.savefig(dir_name+title+'.png')

def plot_robinson(model, vec, m, v_scale=1.0, dir_name='./', title='Velocity and Divergence at CMB'):
   from mpl_toolkits.basemap import Basemap
   E = model.E
   Nk = model.Nk
   Nl = model.Nl
   th = model.th[0,:]

   ## Plot Robinson vector field
   #### Display waves on a Spherical Map Projection
   projtype = 'robin'
   ## Create full vector grid in theta and phi
   u_1D = model.get_variable(vec,'uph')[0][-1,:]
   v_1D = model.get_variable(vec,'uth')[0][-1,:]
   Nph = 2*Nl
   ph = np.linspace(-180.,180.-360./Nph,Nph)
   lon_grid, lat_grid = np.meshgrid(ph,th*180./np.pi-90.,)
   v = (np.exp(1j*m*lon_grid*np.pi/180.).T*v_1D).T
   u = (np.exp(1j*m*lon_grid*np.pi/180.).T*u_1D).T
   absu=u**2 + v**2

   div = np.zeros(u.shape)
   for x in range(Nph):
       for y in range(1,model.Nl-1):
           div[y,x] = u[y,x]*m*1j + (v[y+1,x]-v[y-1,x])/(2.*model.dth)
   ### Plot Robinson Projection
   plt.figure(figsize=(10,10))
   ## Set up map
   bmap = Basemap(projection=projtype,lon_0=0.)
   bmap.drawparallels(np.arange(-90.,90.,15.))
   bmap.drawmeridians(np.arange(0.,360.,15.))
   ## Convert Coordinates to those used by basemap to plot
   lon,lat = bmap(lon_grid,lat_grid)
   bmap.contourf(lon,lat,div.real,15,cmap=plt.cm.RdBu,alpha=0.5)
   Nvec = np.floor(Nl/20.)
   bmap.quiver(lon[::Nvec,::Nvec],lat[::Nvec,::Nvec],u[::Nvec,::Nvec].real,v[::Nvec,::Nvec].real, scale=v_scale)
   plt.title(title)
#    plt.show()
   plt.savefig(dir_name+title+'.png')

def plot_A(A,m):
    ### Plot A Matrix (M*l*x = A*x) ###
    plt.figure(figsize=(10,10))
    plt.spy(np.abs(A.todense()))
    plt.grid()
    plt.title('A,m='+str(m)+' matrix all terms (M*l*x = A*x)')
    plt.savefig('./output/m={0}/A_matrix_m={0}.png'.format(m))
    # plt.subplot(3,2,2)
    # plt.spy(np.abs(A.todense().real))
    # plt.grid()
    # plt.title('A'+str(m)+' real')
    #
    # plt.subplot(3,2,3)
    # plt.spy(np.abs(A.todense().imag))
    # plt.grid()
    # plt.title('A'+str(m)+' imaginary')
    #
    # A_pos = np.matrix(A.todense())
    # A_pos[A.todense()<0.] = 0
    # plt.subplot(3,2,4)
    # plt.spy(np.abs(A_pos))
    # plt.grid()
    # plt.title('A'+str(m)+' positive')
    #
    # A_neg = np.matrix(A.todense())
    # A_neg[A.todense()>0.] = 0
    # plt.subplot(3,2,5)
    # plt.spy(np.abs(A_neg))
    # plt.grid()
    # plt.title('A'+str(m)+' negative')
    # plt.tight_layout()
    # plt.savefig('./output/m={0}/A_matrix_m={0}.png'.format(m))

def plot_M(M,m):
    plt.figure(figsize=(10,10))
    plt.title('M Matrix (M*l*x = A*x)')
    plt.spy(np.abs(M.todense()))
    plt.grid()
    plt.savefig('./output/m={0}/M_matrix_m={0}.png'.format(m))

def plot_pcolormesh_rth(model,val,vec,dir_name='./',title='pcolormesh MAC Wave Plot', physical_units = False):
    plt.close('all')
    r_star = model.r_star
    P_star = model.P_star
    B_star = model.B_star
    u_star = model.u_star
    r = np.concatenate([model.rm[:,0], model.rp[-1:,0]],axis=0)*r_star/1e3
    th = np.concatenate([model.thm[0,:], model.thp[0,-1:]],axis=0)*180./np.pi
    rpl, thpl = np.meshgrid(r,th)
    # rpl = np.concatenate([model.rm, model.rp[-1:,:]],axis=0)*r_star/1e3
    # thpl = np.concatenate([model.thm, model.thp[-1:,:]]*180./np.pi
    fig = plt.figure(figsize=(14,14))
    fig.suptitle(title, fontsize=14)
    gs = gridspec.GridSpec(len(model.model_variables), 2, width_ratios=[100, 1])
    axes = []
    gs_data_list = []
    for ind, var in enumerate(model.model_variables):
        gs_data_list.append(gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[ind*2], wspace=0.01))

        var_data = model.get_variable(vec, var)
        axes.append(plt.subplot(gs_data_list[ind][0]))
        axes.append(plt.subplot(gs_data_list[ind][1]))
        axes.append(plt.subplot(gs[ind*2+1]))
        if physical_units:
            if var in ['ur', 'uth', 'uph']:
                var_data = var_data*u_star*31556.926
                axes[ind*3].set_title(var+' real (km/yr)')
                axes[ind*3+1].set_title(var+' imag (km/yr)')
            elif var in ['br', 'bth', 'bph']:
                var_data = var_data*B_star*1e3
                axes[ind*3].set_title(var+' real (mT)')
                axes[ind*3+1].set_title(var+' imag (mT)')
            elif var == 'r_disp':
                var_data = var_data*r_star
                axes[ind*3].set_title(var+' real (m)')
                axes[ind*3+1].set_title(var+' imag (m)')
            elif var == 'p':
                var_data = var_data*P_star
                axes[ind*3].set_title(var+' real (Pa)')
                axes[ind*3+1].set_title(var+' imag (Pa)')
        else:
            axes[ind*3].set_title(var+' real')
            axes[ind*3+1].set_title(var+' imag')
        var_max = np.amax(abs(var_data))
        axes[ind*3].pcolormesh(thpl,rpl,var_data.real.T, cmap='RdBu',vmin=-var_max, vmax=var_max)
        axes[ind*3].set_ylabel('radius (km)')
        p = axes[ind*3+1].pcolormesh(thpl,rpl,var_data.imag.T, cmap='RdBu',vmin=-var_max, vmax=var_max)
        axes[ind*3+1].get_yaxis().set_ticks([])
        plt.colorbar(p, format='%.0e', cax=axes[ind*3+2], ticks=np.linspace(-var_max,var_max,4))
    fig.set_tight_layout(True)
    plt.subplots_adjust(top=0.95)
#    plt.show()
    plt.savefig(dir_name+title+'.png')

def plot_vel_AGU(model,vec,dir_name='./',title='Velocity for AGU', physical_units = False):
    var2plt = ['ur','uth','uph','bth','bph']
    plt.close('all')
    r_star = model.r_star
    P_star = model.P_star
    B_star = model.B_star
    u_star = model.u_star
    rpl = model.r*r_star/1e3
    thpl = model.th*180./np.pi
    ur = model.get_variable(vec, 'ur')*u_star
    ur[:,::2] = ur[:,::2]*-1
    urmax = np.amax(abs(ur))
    uth = model.get_variable(vec, 'uth')*u_star
    uth[:,::2] = uth[:,::2]*-1
    uthmax = np.amax(abs(uth))
    uph = model.get_variable(vec, 'uph')*u_star
    uph[:,::2] = uph[:,::2]*-1
    uphmax = np.amax(abs(uph))
    bth = model.get_variable(vec, 'bth')*B_star
    bth[:,::2] = bth[:,::2]*-1
    bthmax = np.amax(abs(bth))
    bph = model.get_variable(vec, 'bph')*B_star
    bph[:,::2] = bph[:,::2]*-1
    bphmax = np.amax(abs(bph))

    fig = plt.figure(figsize=(10,8))
    fig.suptitle('m=2 Westward Travelling Wave', fontsize=14)

    ax = plt.subplot(511)
    plt.title('Radial Velocity')
    ax.set_xticklabels([])
    p = plt.pcolormesh(thpl,rpl,ur.real, cmap='RdBu',vmin=-urmax, vmax=urmax)
    plt.colorbar(p, format='%.0e', ticks=np.linspace(-urmax,urmax,4))

    ax = plt.subplot(512)
    plt.title('Latiudinal Velocity')
    ax.set_xticklabels([])
    p = plt.pcolormesh(thpl,rpl,uth.real, cmap='RdBu',vmin=-uthmax, vmax=uthmax)
    plt.colorbar(p, format='%.0e', ticks=np.linspace(-uthmax,uthmax,4))

    ax = plt.subplot(513)
    plt.title('Azimuthal Velocity')
    ax.set_xticklabels([])
    p = plt.pcolormesh(thpl,rpl,uph.imag, cmap='RdBu',vmin=-uphmax, vmax=uphmax)
    plt.colorbar(p, format='%.0e', ticks=np.linspace(-uphmax,uphmax,4))

    ax = plt.subplot(514)
    plt.title('Latitudinal Magnetic Field Perturbation')
    ax.set_xticklabels([])
    p = plt.pcolormesh(thpl,rpl,bth.imag, cmap='RdBu',vmin=-bthmax, vmax=bthmax)
    plt.colorbar(p, format='%.0e', ticks=np.linspace(-bthmax,bthmax,4))

    ax = plt.subplot(515)
    plt.title('Azimuthal Magnetic Field Perturbation')
    p = plt.pcolormesh(thpl,rpl,bph.real, cmap='RdBu',vmin=-bphmax, vmax=bphmax)
    plt.colorbar(p, format='%.0e', ticks=np.linspace(-bphmax,bphmax,4))

    fig.set_tight_layout(True)
    plt.subplots_adjust(top=0.9)
    # plt.savefig(dir_name+title+'.png')

def plot_buoy_struct(model, dir_name='./', title='buoyancy_structure'):
    plt.close('all')
    drho_dr = -model.omega_g**2*model.rho/model.g  # density gradient
    fig = plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    drho = np.zeros((model.Nk,1))
    for i in range(0,model.Nk):
        drho[i] = sum(model.dr*model.r_star*drho_dr[:i,model.Nl//2])
    plt.plot(drho,((model.r[:,0]-1)*model.r_star)/1000)
#    plt.plot(drho_dr[:,model.Nl/2]*1e6,(model.r[:,0]-1)*model.r_star/1000)
    plt.title('density perturbation off adiabat')
    plt.ylabel('depth below CMB (km)')
    plt.xlabel('density perturbation off adiabat (kg/m^4)')

    plt.subplot(1,2,2)
    plt.plot(model.omega_g[:,model.Nl//2]*model.t_star,(model.r[:,0]-1)*model.r_star/1000)
    plt.title('buoyancy frequency')
    plt.ylabel('depth below CMB (km)')
    plt.xlabel('buoyancy frequency (omega_g/Omega)')
    fig.set_tight_layout(True)
    plt.savefig(dir_name+title+'.png')

def plot_B(model, dir_name='./', title='B field structure'):
    plt.close('all')
    fig = plt.figure(figsize=(10,10))
    plt.subplot(3,1,1)

    plt.plot(model.th[0,:]*180./np.pi,model.Br[model.Nk//2,:]*model.B_star*1e3)
    plt.title('Br background field')
    plt.ylabel('Br in (10^-3 T)')
    plt.xlabel('colatitude in degrees')
    plt.grid()
    ymin, ymax = plt.ylim()
    if ymin > 0.:
        plt.ylim(ymin=0.)

    plt.subplot(3,1,2)
    plt.plot(model.th[0,:]*180./np.pi,model.Bth[model.Nk//2,:]*model.B_star*1e3)
    plt.title('B_theta background field')
    plt.ylabel('B_theta in (10^-3 T)')
    plt.xlabel('colatitude in degrees')
    plt.grid()
    ymin, ymax = plt.ylim()
    if ymin > 0.:
        plt.ylim(ymin=0.)

    plt.subplot(3,1,3)
    plt.plot(model.th[0,:]*180./np.pi,model.Bph[model.Nk//2,:]*model.B_star*1e3)
    plt.title('B_phi background field')
    plt.ylabel('B_phi in (10^-3 T)')
    plt.xlabel('colatitude in degrees')
    plt.grid()
    ymin, ymax = plt.ylim()
    if ymin > 0.:
        plt.ylim(ymin=0.)

    fig.set_tight_layout(True)
    plt.savefig(dir_name+title+'.png')

def plot_Uphi(model, dir_name='./', title='Uphi structure'):
    plt.close('all')
    fig = plt.figure(figsize=(10,5))
    plt.pcolor(model.th*180./np.pi, (model.r-1)*model.r_star/1000, model.Uphi*model.u_star)
    plt.colorbar()
    plt.title('Uphi background velocity field')
    plt.ylabel('depth below CMB (km)')
    plt.xlabel('colatitude in degrees')
    fig.set_tight_layout(True)
    plt.savefig(dir_name+title+'.png')
























