# #Code used for calculating the hydraulic surface E(psi_leaf, psi_soil) and
# #the SOL used for Carminati & Javaux, Trends in Plant Science, 2020. Code
# #used for Figs.3,4,5. Plant and soil parameters to be adjusted. 
# #Authors: Carminati & Javaux
# #https://www.cell.com/trends/plant-science/fulltext/S1360-1385(20)30118-7
# ##########################################################################
# clear allclose all

# #------------inputs---------------
#   #For varying soil matric potentials hb and transpiration rates E we
# #calculate the leaf water potential hleaf
hb=-seq(10,15010,by=100)#cm bulk soil water potential in cm heads
E=3*(10^-6)*(0:3000) #'#transpiration in cm3 s^-1'

# #root-plant properties
Rroot=0.15*10^7 #resistance#hPa cm^-3 s #for the blue 0.08 # for the green 0.12
r0= 0.05# # cm root radius in cm
r2= 1# #bulk+root+rhizo soil radius in cm
L = 1200# #root length in cm - for the blue 1000 - for the green 10000

## Brooks-Corey parameters for the soil hydraulic properties
h0=-10##this is alfa [cm^-1]
l=1.5#
ths=0.45##not needed in the code, but useful if the decreasing water content needs to be calculated
thr=0.01##not needed in the code, but useful if the decreasing water content needs to be calculated
k0=0.3*10^-2##cm/s - for the blue 0.2
tau=2.5##corresponding to 2+l*(a+2) with l exp for theta(h) and a tortusoity 

# # # radial rhizosphere domain (do not seem to be used here)
r= seq(r0,r2,by=0.01)
# # r=r(:)
# dr=(r2-r0)/(length(r)-1)
# # dr=dr(1)

# #parameters for cavitation if you plot the harmonic mean of 1/Rroot and k_x 
# #you get a conductance that maximum is the one of root and if hx decreases the conductivity drops 
k0_x=1/Rroot# #hPa-1 cm^3 s-1
h0_x=-20000# #in cm (or 10^-4 Mpa) approximately equal p50
tau_x=5
#k_x=k0_x.*(hleaf./h0_x).^-tau_x #this would be the xylem conductance

#------------model---------------
  #--> for given soil and leaf water potential we calculate the transpiration
#rate based on the conductivities. 
#here we prepare the matrix
# hroot=zeros(length(hb),length(E))
hroot = matrix(data = rep(0,length(hb)*length(E)),
               ncol = length(E))
hx=hroot
hx_max=hroot
hleaf=hroot
Emax=hroot
hbmat=hroot
Emat=hroot

# #calculation of hleaf
for (j1 in 1:length(hb))#size of psi_soil
{
  for (j2 in 1:length(E))#size of E
  {
    Emat[j1,j2]=E[j2]
    hbmat[j1,j2]=hb[j1]
    csoil=-2*pi*r0*L* k0/(1-tau)/((h0)^(-tau))/
      (r0/2- r0*r2^2*(log(r2)-log(r0))/(r2^2-r0^2) )
    hroot[j1,j2]=-abs(-E[j2]/csoil+hb[j1]^(1-tau))^(1/(1-tau))            
    if ((hb[j1]-hroot[j1,j2])>200000) break
    #dissipation in the root
    hx[j1,j2]=hroot[j1,j2]-Rroot*E[j2]      # Couvreur model
    #hmax - this is when the system is linear
    hx_max[j1,j2]=hbmat[j1,j2]-Rroot*E[j2]   
    #dissipation in the xylem including cavitation 
    cx=-k0_x/(1-tau_x)*(h0_x^tau_x)
    hleaf[j1,j2]=-abs(E[j2]/cx+hx[j1,j2]^(1-tau_x))^(1/(1-tau_x))
    if ((hx[j1,j2]-hleaf[j1,j2])>30000) break
  }
}

hroot[100,1000]

range(hroot,na.rm=T)
