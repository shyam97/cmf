# -*- coding: utf-8 -*-
"""
Created on Sat May 09 14:05:21 2020

@author: NotCassandra
"""
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as spint
#import BEM_data 





geo = 'wind' ## Choose your rotor. 'wind' for wind turbine, 'prop' for aircraft propeller
N_rotor = 1 ## Choose your number of rotors between 1 and 2
oper = 'land' ## For aircraft propeller. 'cruise' or 'land' condition.
distribution = 'uni' ## either 'uni' or 'cos' for uniform and cosine distribution respectively
N_cp = 20 ## number of control points on a single blade
t_end = 2.5 ## final 'time' for the geometry of the spiralling trailing vortices
N_t = int(t_end*10) ## number of pairs of discretized spiralling trailing vortex filaments. 10 times t_end to get a good spiral
                
## 'L' for the separation distance between the two rotors is specified in the "geometry creation" block


    
''' ================= INDUCED VELOCITY FUNCTION ============================ '''
''' ======================================================================== '''

def velocity_3D_from_vortex_filament(XV1, XV2, XVP,CORE):
   # function to calculate the velocity induced by a straight 3D vortex filament
   # with circulation GAMMA at a point VP1. The geometry of the vortex filament
   # is defined by its edges: the filament starts at XV1 and ends at XV2.
   # the input CORE defines a vortex core radius, inside which the velocity
   # is defined as a solid body rotation.
   # The function is adapted from the algorithm presented in:
   #                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
   #                Vol. 13. Cambridge university press, 2001.

   # read coordinates that define the vortex filament
   X1 = XV1[0] # start point of vortex filament
   Y1 = XV1[1] 
   Z1 = XV1[2] 
   
   X2 = XV2[0] # end point of vortex filament
   Y2 = XV2[1]
   Z2 = XV2[2] 

   # read coordinates of target point where the velocity is calculated
   XP = XVP[0]
   YP = XVP[1]
   ZP = XVP[2]
   
   # calculate geometric relations for integral of the velocity induced by filament
   R1 = np.sqrt((XP-X1)**2 + (YP-Y1)**2 + (ZP-Z1)**2)
   R2 = np.sqrt((XP-X2)**2 + (YP-Y2)**2 + (ZP-Z2)**2)
   
   R1XR2_X = (YP-Y1)*(ZP-Z2)-(ZP-Z1)*(YP-Y2)
   R1XR2_Y = -(XP-X1)*(ZP-Z2)+(ZP-Z1)*(XP-X2)
   R1XR2_Z = (XP-X1)*(YP-Y2)-(YP-Y1)*(XP-X2)
   
   R1XR_SQR = R1XR2_X**2 + R1XR2_Y**2 + R1XR2_Z**2
   R0R1 = (X2-X1)*(XP-X1) + (Y2-Y1)*(YP-Y1) + (Z2-Z1)*(ZP-Z1)
   R0R2 = (X2-X1)*(XP-X2) + (Y2-Y1)*(YP-Y2) + (Z2-Z1)*(ZP-Z2)
   
   # check if target point is in the vortex filament core, and modify to solid body rotation
   if (R1XR_SQR < CORE**2):
       R1XR_SQR = CORE**2
       #print('R1XR_SQR < CORE**2')
       #GAMMA = 0;
   if (R1 < CORE): 
       R1 = CORE
       #print('R1 < CORE**2')
       #GAMMA = 0;
   if (R2 < CORE):
       R2 = CORE
       #print('R2 < CORE**2')
       #GAMMA = 0;
       
   # determine scalar
   K = 1/(4*np.pi*R1XR_SQR) * (R0R1/R1 - R0R2/R2);
   # determine the three velocity components
   U = K*R1XR2_X
   V = K*R1XR2_Y
   W = K*R1XR2_Z
   # output results, vector with the three velocity components
   return U, V, W;


''' =========================== AIRFOIL DATA =============================== '''
''' ======================================================================== '''
alpha_dat = []
cl_dat = []
cd_dat = []
cm_dat = []

if geo == 'wind':
    f = open("DU95W180_polar.txt", 'r')
    lines = f.readlines()
    f.close()
elif geo == 'prop':
    f = open("ARAD8pct_polar.txt", 'r')
    lines = f.readlines()
    f.close()
else:
    print("Check the spelling of the airfoil variable")
#print ("Geometry =  ", geo)
    
for line in lines:
    a = line.split()
    alpha_dat.append(float(a[0]))
    cl_dat.append(float(a[1]))
    cd_dat.append(float(a[2]))
    cm_dat.append(float(a[3]))   

get_Cl = spint.interp1d(alpha_dat, cl_dat) ## INTERPOLATION FUNCTIONS OF Cl FOR THE INTERVAL ALPHA_DAT 
get_Cd = spint.interp1d(alpha_dat, cd_dat) ## INTERPOLATION FUNCTIONS OF Cd FOR THE INTERVAL ALPHA_DAT 
get_Cm = spint.interp1d(alpha_dat, cm_dat) ## INTERPOLATION FUNCTIONS OF Cm FOR THE INTERVAL ALPHA_DAT 


''' =========================================================================== '''
''' =========================================================================== '''
''' ===================== GEOMETRY OF THE LIFTING LINES ======================= '''
''' =========================================================================== '''
''' =========================================================================== '''

## Preliminary geometries
if geo =='wind':
    ## wind turbine
    N_B = 3 ## number of blades
    R = 50
    R0 = 0.2*R     
    D = 2*R
    U0 = 10.0  ## freestream velocity (m/s)
    TSR = 6    
    Omega = TSR*U0/R
    P0 = 101325.0                   ## ISA static pressure
    rho = 1.225                     ## ISA air density at sea-level (kg/m3)
elif geo== 'prop':
    ## aircraft propeller
    N_B = 6                     ## number of blades
    R = 0.7                     ## blade radius (m)
    R0 = 0.25                   ## non-dimensional starting position of the blade (-)
    D = 2*R
    if oper == 'cruise':
        U0 = 60.0                   ## freestream velocity (m/s) (35 or 60 m/s)
    elif oper == 'land':
        U0 = 35.0                   ## freestream velocity (m/s) (35 or 60 m/s)
    else:
        print ("Typo in the variable 'oper'.")
    rpm = 1200.0                ## RPM
    rho = 1.00649               ## ISA air density for altitude of 2000 m (kg/m3) 
    P0 = 79495.22               ## ISA static pressure         
    Omega = rpm*2*np.pi/60      ## angular speed of rotor (rad/s)
else:
    print ('Kindly fix your typo in the variable "geo"')

''' ===================== 1. Control points ======================= '''


X_cp = np.zeros((N_cp, 3)) ## coordinates of the control points. 1st index: index of cp. 2nd index: x, y, or z

## z coordinates of the control points (CHECK OUT THE COORDINATE SYSTEM)
if distribution == 'uni':
    X_cp[:,2] = np.linspace(R0, R, N_cp, endpoint = False) 
    spacing = X_cp[1,2]-X_cp[0,2]
    X_cp[:,2] = X_cp[:,2] + spacing*0.5
elif distribution == 'cos':
    X_cp[:,2] = (0.5*(1-np.cos(np.linspace(0,np.pi, N_cp, endpoint = False))))*(R-R0) + R0  ## (COSINE)
    spacing = np.zeros(N_cp)
    for i in range(N_cp-1):
        spacing[i] = X_cp[i+1,2] - X_cp[i,2] 
    spacing[-1] = R - X_cp[-1,2]
    X_cp[:,2] = X_cp[:,2] + spacing*0.5
else:
    print ("Typo in the distribution variable")

    
#plt.plot( X_cp[:,2], np.ones(N_cp), 'bo')
#plt.grid()
#plt.show()

''' ===================== 2. Bound vortex filaments ======================= '''

X_bv = np.zeros((N_cp, 2, 3)) 
## coordinates of the bound vortices. 1st index: index of cp. 2nd in dex: X1 or X2. 3rd index: x,y,z

X1_bv = X_cp[:,2] - spacing*0.5 ## The coordinate closer to the root
X2_bv = X_cp[:,2] + spacing*0.5 ## The coordinate further from the root
X_bv[:,0,2] =  X1_bv 
X_bv[:,1,2] =  X2_bv

''' =============================================================================== '''
''' ===================== INSERT CHORD AND TWIST BELOW ============================ '''
''' =============================================================================== '''

## For the trailing vortex geometry
## change to the optimized geometry later


if geo == 'wind':
    
    beta = -2.0 ## blade pitch (deg)
    
    R_opt = np.array([0.2, 0.224, 0.24 , 0.256, 0.272, 0.288, 0.304, 0.32 , 0.336,
       0.352, 0.368, 0.384, 0.4  , 0.416, 0.432, 0.448, 0.464, 0.48 ,
       0.496, 0.512, 0.528, 0.544, 0.56 , 0.576, 0.592, 0.608, 0.624,
       0.64 , 0.656, 0.672, 0.688, 0.704, 0.72 , 0.736, 0.752, 0.768,
       0.784, 0.8  , 0.816, 0.832, 0.848, 0.864, 0.88 , 0.896, 0.912,
       0.928, 0.944, 0.96 , 0.976, 1])*R
    
    twist_opt = np.array([12.63387796, 16.2726184 , 16.31120026, 15.61972595, 14.70030495,
       13.71394314, 12.7258278 , 11.76521225, 10.84551051,  9.97244082,
        9.14776302,  8.37114751,  7.64116811,  6.95585153,  6.31298954,
        5.71031809,  5.14561945,  4.61677834,  4.12181004,  3.65887086,
        3.22625736,  2.82239791,  2.44583909,  2.09527158,  1.77004261,
        1.46878615,  1.18992877,  0.93200244,  0.69363484,  0.47353956,
        0.27050611,  0.08338992, -0.08889789, -0.24740105, -0.39312919,
       -0.52706986, -0.6505508 , -0.76539033, -0.87294254, -0.97459137,
       -1.07180702, -1.166224  , -1.25978035, -1.35497588, -1.45538396,
       -1.56676835, -1.69987452, -1.87896498, -2.17849463, -3.07462437])
    
    getTheta_WT = spint.interp1d(R_opt, twist_opt)
    ## CHORD, TWIST FOR COORDINATE "1"
    theta1 = getTheta_WT(X_bv[:,0,2])       ## twist (deg)
    c1 = 4.38*(1- X_bv[:,0,2]/R) + 1.75             ## chord distribution (m)
    
    ## CHORD, TWIST FOR COORDINATE "2" (bound v)
    theta2 = getTheta_WT(X_bv[:,1,2])         ## twist (deg)
    c2 = 4.38*(1- X_bv[:,1,2]/R) + 1.75               ## chord distribution (m)
    
    ## CHORD, TWIST AT CONTROL POINTS
    theta = getTheta_WT(X_cp[:,2])         ## twist (deg)
    c = 4.38*(1- X_cp[:,2]/R) + 1.75               ## chord distribution (m)
elif geo == 'prop':
    ## OPTIMIZED PROPELLER GEOMETRY
    if oper == 'cruise':
        beta_col = 44.7603      ## reference pitch (deg). Cruise: 44.7603. Landing: 31.8266
    else:
        beta_col = 31.8266
    
    theta1 = -37.7777777778*X_bv[:,0,2]/R + 25.0 + beta_col             ## twist/pitch (deg)
    c1 = 0.176315789474 - 0.06*X_bv[:,0,2]/R                         ## chord distribution (-)
    c1 = c1*R                                     ## chord distribution (m)
    
    theta2 = -37.7777777778*X_bv[:,1,2]/R + 25.0 + beta_col             ## twist/pitch (deg)
    c2 = 0.176315789474 - 0.06*X_bv[:,1,2]/R                         ## chord distribution (-)
    c2 = c2*R                                     ## chord distribution (m)
    
    theta = -37.7777777778*X_cp[:,2]/R + 25.0 + beta_col             ## twist/pitch (deg)
    c = 0.176315789474 - 0.06*X_cp[:,2]/R                         ## chord distribution (-)
    c = c*R                                     ## chord distribution (m)
else:
    print 




#plt.plot(X_bv[:,0,2], theta1)
#plt.grid()
#plt.show()

''' ===================== 3. Trailing vortex filaments ======================= '''

t = np.linspace(0,t_end, N_t)

X_tv = np.zeros((N_cp, N_t, 2, 2, 3)) 
## coordinates of the trailing vortices. 
## 1st index: index of cp (control point) 
## 2nd index: trailing vortex ELEMENT
## 3rd index: INDEX OF PART OF PAIR.  INDEX 0: CLOSER TO THE ROOT, INDEX 1: FURTHER FROM THE ROOT
## 4th index: INDEX OF POINT: X1 or X2. INDEX 0: UPSTREAM POINT. INDEX 1: DOWNSTREAM POINT
## 5th index: x,y, or z

## template: X_tv[cp, ELEMENT index, pair, up/downstream POINT, {x,y,z}]

## ===============================================================================
## FIRST: THE TRAILING VORTEX FILAMENT ELEMENTS THAT FOLLOW THE CHORD OF THE BLADE
## ===============================================================================

## z coordinates. upstream points
X_tv[:, 0, 0, 0, 2] = X1_bv ## pair closer to the root, upstream point
X_tv[:, 0, 1, 0, 2] = X2_bv ## pair further from the root, upstream point

## z coordinates. downstream points
X_tv[:, 0, 0, 1, 2] = X1_bv ## pair closer to the root, downstream z 
X_tv[:, 0, 1, 1, 2] = X2_bv ## pair closer to the root, downstream z

## x, y coordinates. downstream point, closer to the root
X_tv[:, 0, 0, 1, 0] = c1*np.sin(theta1*180/np.pi)  ## pair closer to the root, downstream x
X_tv[:, 0, 0, 1, 1] = -c1*np.cos(theta1*180/np.pi) ## pair closer to the root, downstream y 

## x, y coordinates. downstream point, further from the root
X_tv[:, 0, 1, 1, 0] = c2*np.sin(theta2*180/np.pi) ## pair further from the root, downstream x
X_tv[:, 0, 1, 1, 1] = -c2*np.cos(theta2*180/np.pi)## pair further from the root, downstream y 


''' ===============================================================================
                                Initial guess of gamma
    =============================================================================== '''


#gamma = np.matrix(np.linspace(1.0,2.0, N_cp)) ## For wind turbine's current geometry, this works
#gamma = -np.matrix(np.linspace(0.05,0.20, N_cp)) 
#gamma = (0.25*X_cp[:,2]- R*0.5)**2 #+ 0.1
#gamma = np.matrix(0.0*np.ones(N_cp))
#gamma = -(X_cp[:,2] - R0)*(X_cp[:,2] - R)
#gamma = -0.01*np.sqrt(1-(2*X_cp[:,2]/(2*R))**2)



## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Choosing initial gamma to be the mid-span (r = R/2) gamma without induced velocities
## Note that N_cp should be even number below
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

V_p = np.sqrt(U0**2 + (Omega*R/2)**2)
phi = np.arctan2(U0, Omega*R/2)
if geo == 'wind':
    alpha = phi*180.0/np.pi - theta[N_cp/2]
elif geo == 'prop':
    alpha = theta[N_cp/2] - phi*180.0/np.pi
else:
    print ('Typo')
gamma = 0.5*c[N_cp/2]*V_p*get_Cl(alpha)
gamma = gamma*np.ones(N_cp)



## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##      Use BEM data 
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#if geo == 'prop':
#    data = BEM_data.get_prop_data() ## 0 = r/R, 1 = gamma_land, 2 = gamma_cruise
#    if oper == 'cruise':
#        get_gamma = spint.interp1d(data[0], data[2]) 
#    elif oper == 'land':
#        get_gamma = spint.interp1d(data[0], data[1]) 
#    else:
#        print ('Typo in oper')
#        
#gamma = get_gamma(X_cp[:,2]/R)



''' ===============================================================================
                                Solver/start of the loop
    =============================================================================== '''

gamma = np.matrix(gamma)
gamma = np.transpose(gamma)
gamma_hist = [] ## history of gamma
gamma_diff = 3.0#max(np.abs(gamma)) ## some random values, high enough
gamma_tol = 1e-3
N_iter = 50 ## maximum number of iterations
counter = 0

aw = 0.50 ## induced velocity on the blade 
    
u,v,w = np.zeros((3, N_cp, N_cp)) ## Induced velocities

while counter<N_iter and gamma_diff>gamma_tol:
    
    gamma_old = gamma ## store old values for convergence criterion
    
    """ ========================================================= """
    """ ============ below is related to the solver ============= """
    
    Uw = U0*(1-aw) ## convection velocity of the wake
    
    """ ========================================================= """
    """ ========================================================= """
    
    ## ============ 1. TRAILING VORTICES CLOSER TO THE ROOT [-,-,0,-,-] ================
    
    r1 = np.sqrt((X_tv[:,0,0,1,1])**2 + (X_tv[:,0,0,1,2])**2)
    
    for i in range(1,N_t): ## i: looping over all elements, get the upstream point coordinates
        j = i-1 ## for the parameter t, upstream
        X_tv[:,i,0,0,0] = X_tv[:, 0, 0, 1, 0] + Uw*t[j]
        X_tv[:,i,0,0,1] = X_tv[:, 0, 0, 1, 1] + r1*np.sin(Omega*t[j])
        X_tv[:,i,0,0,2] = r1*np.cos(Omega*t[j])
        
    for i in range(1,N_t): ## i: looping over all elements, get the downstream point coordinates
        X_tv[:,i,0,1,0] = X_tv[:, 0, 0, 1, 0] +Uw*t[i]
        X_tv[:,i,0,1,1] = X_tv[:, 0, 0, 1, 1] +r1*np.sin(Omega*t[i])
        X_tv[:,i,0,1,2] = r1*np.cos(Omega*t[i])
    
    ## ============ 2. TRAILING VORTICES FURTHER FROM THE ROOT [-,-,1,-,-] ================
    
    r2 = np.sqrt((X_tv[:,0,1,1,1])**2 + (X_tv[:,0,1,1,2])**2)
    
    for i in range(1,N_t): ## i: looping over all elements, get the upstream point coordinates
        j = i-1
        X_tv[:,i,1,0,0] = X_tv[:, 0, 0, 1, 0] +Uw*t[j]
        X_tv[:,i,1,0,1] = X_tv[:, 0, 0, 1, 1] +r2*np.sin(Omega*t[j])
        X_tv[:,i,1,0,2] = r2*np.cos(Omega*t[j])
        
    for i in range(1,N_t): ## i: looping over all elements, get the downstream point coordinates
        X_tv[:,i,1,1,0] = X_tv[:, 0, 0, 1, 0] +Uw*t[i]
        X_tv[:,i,1,1,1] = X_tv[:, 0, 0, 1, 1] +r2*np.sin(Omega*t[i])
        X_tv[:,i,1,1,2] = r2*np.cos(Omega*t[i])
    
    
    ''' =============================================================================== '''
    ''' =============================================================================== '''
    ''' ================ CREATING THE GEOMETRY OF THE OTHER BLADES ==================== '''
    ''' =============================================================================== '''
    ''' =============================================================================== '''
    
    L = D ## Separation distance between the two rotors. For the assignment: Infty, 5D, 2D and D
    
    angle = 2*np.pi/N_B
    
    X_cp_all = np.zeros((N_B*N_rotor, N_cp, 3))
    X_bv_all = np.zeros((N_B*N_rotor, N_cp, 2, 3))
    X_tv_all = np.zeros((N_B*N_rotor, N_cp, N_t,2,2,3))
    
    for i in range(N_B):
        X_cp_all[i,:,:] = X_cp
        X_bv_all[i,:,:,:] = X_bv
        X_tv_all[i,:,:,:,:,:] = X_tv
        
        # Rotate y: "y' = y*cos(angle) - z*sin(angle)"
        X_cp_all[i,:,1] = X_cp[:,1]*np.cos(i*angle) - X_cp[:,2]*np.sin(i*angle)
        X_bv_all[i,:,:,1] = X_bv[:,:,1]*np.cos(i*angle) - X_bv[:,:,2]*np.sin(i*angle)
        X_tv_all[i,:,:,:,:,1] = X_tv[:,:,:,:,1]*np.cos(i*angle) - X_tv[:,:,:,:,2]*np.sin(i*angle)
    
        # Rotate z: "z' = y*sin(angle) + z*cos(angle)"
        X_cp_all[i,:,2] = X_cp[:,1]*np.sin(i*angle) + X_cp[:,2]*np.cos(i*angle)
        X_bv_all[i,:,:,2] = X_bv[:,:,1]*np.sin(i*angle) + X_bv[:,:,2]*np.cos(i*angle)
        X_tv_all[i,:,:,:,:,2] = X_tv[:,:,:,:,1]*np.sin(i*angle) + X_tv[:,:,:,:,2]*np.cos(i*angle)
        
    if N_rotor!=1:
        for i in range(N_B, N_B*N_rotor):
            X_cp_all[i,:,:] = X_cp_all[i-N_B,:,:]
            X_bv_all[i,:,:,:] = X_bv_all[i-N_B,:,:,:]
            X_tv_all[i,:,:,:,:,:] = X_tv_all[i-N_B,:,:,:,:,:]
            
            ## Shift the y coordinate by L in any direction
            X_cp_all[i,:,1] = X_cp_all[i,:,1] - L
            X_bv_all[i,:,:,1] = X_bv_all[i,:,:,1] - L
            X_tv_all[i,:,:,:,:,1] = X_tv_all[i,:,:,:,:,1] - L
            
    
    
    ''' ========================================================================================== '''
    ''' ================ CREATION OF THE BIOT-SAVART INDUCED VELOCITY MATRICES =================== '''
    ''' ========================================================================================== '''
    
    CORE = 1e-3
    BSM = np.zeros((3,N_cp,N_cp)) #(u, v, w) is the first index
    
#    #print(BSM)
    for i in np.arange(N_cp): #Loops the control point location
        XVP = X_cp_all[0,i,:] ## ONLY THE FIRST BLADE I.E. Z = XCP[0,:,2] = R
        for j in np.arange(N_cp): #Loops the location of influence for the corresponding Gamma 
            for i_b in np.arange(N_B*N_rotor): # Loops over the N_B blades and over N rotors spaced at L distance
                #Bound Vortex
                XV1 = X_bv_all[i_b,j,0,:]
                XV2 = X_bv_all[i_b,j,1,:]
                BSM[:,i,j] += velocity_3D_from_vortex_filament(XV1, XV2, XVP,CORE)
                #Trailing Vortex
                for i_t in np.arange(N_t): #loops over the len of the trailing vortex 
                    #Root Trailing Vortex
                    XV1 = X_tv_all[i_b,j,i_t,0,1,:] #Takes into account sign of circulation as for root trailing vortex circulation sign is opposite of tip TV
                    XV2 = X_tv_all[i_b,j,i_t,0,0,:]
                    BSM[:,i,j] += velocity_3D_from_vortex_filament(XV1, XV2, XVP,CORE)
                    #Tip Trailing Vortex
                    XV1 = X_tv_all[i_b,j,i_t,1,0,:]
                    XV2 = X_tv_all[i_b,j,i_t,1,1,:]
                    BSM[:,i,j] += velocity_3D_from_vortex_filament(XV1, XV2, XVP,CORE)      
    
#    vorBV = velocity_3D_from_vortex_filament(X_bv_all[:,:,0,:], X_bv_all[:,:,1,:], X_cp,CORE)
#    vorBV = np.sum(vorBV,0) ## sum all blades
#    
#    vorTV1 = velocity_3D_from_vortex_filament(X_tv_all[:,:,:,0,1,:], X_tv_all[:,:,:,0,0,:], X_cp,CORE)
#    vorTV1 = np.sum(vorTV1, 1) ## sum all elements
#    vorTV1 = np.sum(vorTV1, 0) ## sum all blades (including all rotor contributions)
#
#    vorTV2 = velocity_3D_from_vortex_filament(X_tv_all[:,:,:,1,0,:], X_tv_all[:,:,:,1,1,:], X_cp,CORE)
#    vorTV2 = np.sum(vorTV2, 1) ## sum all elements
#    vorTV2 = np.sum(vorTV2, 0) ## sum all blades (including all rotor contributions)   
#    
#    BSM = vorBV + vorTV1 + vorTV2
    
    
    
    
    
    ## ONLY THE FIRST BLADE I.E. Z = XCP[0,:,2] = R (line 377)
    
    u = np.array(np.transpose(np.matrix(BSM[0])*gamma))
    v = np.array(np.transpose( np.matrix(BSM[1])*gamma))
    w = np.array(np.transpose(np.matrix(BSM[2])*gamma))
    
    u = u[0]
    v = v[0]
    w = w[0]
    
    # just operate in arrays, much easier
    

    V_ax = U0 + u
    V_tan = Omega*X_cp_all[0,:,2] + (v) ## for the blade aligned with the z axis, this can be checked by inspection
    
    # now the above are column vectors
    
    V_p = np.sqrt(V_ax**2 + V_tan**2)
    phi = np.arctan2(V_ax, V_tan)
    
    if geo == 'wind':    
        alpha = phi*180.0/np.pi - theta 
        alpha[alpha<min(alpha_dat)] = min(alpha_dat)
        alpha[alpha>max(alpha_dat)] = max(alpha_dat)
        
        Cl = get_Cl(alpha)
        Cd = get_Cd(alpha)
        lift = Cl*0.5*c*rho*V_p**2 
        drag = Cd*0.5*c*rho*V_p**2    
        
        F_axial = lift*np.cos(phi) + drag*np.sin(phi) 
        F_tan = lift*np.sin(phi) - drag*np.cos(phi)   
    else:
        alpha = theta - phi*180.0/np.pi  
        alpha[alpha<min(alpha_dat)] = min(alpha_dat)
        alpha[alpha>max(alpha_dat)] = max(alpha_dat)
        
        ## get forces in propeller convention
        Cl = get_Cl(alpha)
        Cd = get_Cd(alpha)    
        lift = Cl*0.5*c*rho*V_p**2 
        drag = Cd*0.5*c*rho*V_p**2    
        
        F_axial = lift*np.cos(phi) - drag*np.sin(phi)  ## note that the rotation is opposite to the wind turbine case
        F_tan = lift*np.sin(phi) + drag*np.cos(phi)    ## note that the rotation is opposite to the wind turbine case
        ## convert the forces into wind turbine convention for plotting the same way as BEM
        F_axial = -F_axial   
        F_tan = -F_tan
    
    gamma = 0.5*c*V_p*Cl
#    aw = 1-np.mean(V_ax)/U0 ## new average convection velocity ratio
    aw = -np.mean(u)/U0
    
    gamma_diff = np.max(np.abs(gamma-np.transpose(gamma_old)))
    gamma_hist.append(gamma_diff)
    
    gamma = np.transpose(np.matrix(gamma)) ## return to matrix form  for next iteration
    
    counter = counter + 1
    
    print (counter, '  ', aw, ' ', np.mean(u))


## THE INITIAL GUESS OF gamma SEEMS TO BE IMPORTANT FOR THE PROPELLER (OTHERWISE ALPHA NOT IN INTERPOLATION RANGE)
## OR WE GET BETTER DATA OF THE AIRFOIL (MORE ALPHA)

## FORCE COEFFICIENTS
plt.figure()
plt.grid()
plt.plot(X_cp[:,2], F_axial/(0.5*rho*U0**2*R), 'r', label = 'C_axial')
plt.plot(X_cp[:,2], F_tan/(0.5*rho*U0**2*R), 'b', label = 'C_tan')
plt.xlabel("r (m)")
plt.ylabel('Force coefficients (-)')
plt.title(geo)
plt.legend(loc = 0)

## CONVERGENCE HISTORY
plt.figure()
plt.grid()
plt.plot(range(len(gamma_hist)), gamma_hist, 'bo')
plt.xlabel('th iteration')
plt.ylabel('Error')

plt.figure()
plt.grid()
plt.plot(X_cp[:,2], gamma)
plt.xlabel('r (m)')
plt.ylabel('$\\Gamma$ (m**2/s)')

plt.show()