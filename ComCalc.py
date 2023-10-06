import numpy as np
AFlow1 = 21.343 # cm2
AFlow3 = 15.241 # cm2
mdot = 2.68E-01 #kg/sec
Mx1 = 3.26E-01 #Station 1 Mach
Mx3 = 2.90E-01 #Station 3 Mach
Vx1 = 110.050 # m/s
Vx3 = 107.497 # m/s
Pt1 = 99804.985 # Pa
Pt3 = 1.70E+05 # Pa
Tt3 = 346.758 # K
Cp1 = 1010.832 # J/(kg*K)
Cp3 = 1014.440 # J/(kg*K)
Gam1 = 1.399
Gam3 = 1.397
Rcons1 = 288.512 # J/(kg*K)
Rcons3 = 288.512 # J/(kg*K)
ht1 = 2.90E+05 # J/kg
ht2 = 3.50E+05 # J/kg

####Design parameters
rh1 = 1.25 # hub radius cm
hB1 = 0.2 # main blade thickness cm
NB1 = 9 # number of main blades of impeller
hB2 = 0.2 # inter-vane thickness cm
NB2 = 9 # number of inter-vanes of impeller 
hB3 = 0.2 #diffuser blade thickness cm
NB3 = 20 #number diffuser blades
RPM = 60000 #Try this and see how blade mach number reacts
rm2 = 4.5 #cm
Vx2 = 108 # m/s assume same axially velocity here
####Design parameters

#Compute inlet tip radius
B = -(hB1*NB1)/np.pi
C = ((hB1*NB1)/np.pi)*rh1 - rh1**2 - AFlow1/np.pi
rt1 = 0.5*(-B+np.sqrt(B**2-4*C))
residual = AFlow1 - (np.pi*(rt1**2-rh1**2) - (rt1-rh1)*hB1*NB1)
print("Residual is: " + str(residual) + " With rt1 of: " + str(rt1) + " cm") # Residual

#Find feasible RPM
omega = RPM*(2*np.pi/60) # rad/sec
SOS1 = Vx1/Mx1 # Speed of sound in m/s
Ut1 = omega*rt1*0.01 #m/s
w1 = np.sqrt(Ut1**2 + Vx1**2)
MachRel1 = w1/SOS1
print("Blade relative Mach: "+ str(MachRel1)+" At RPM of: "+str(RPM)) 

#Compute the hub tip mean impeller inlet blade angles (no inlet swirl)
rm1 = np.sqrt(0.5*(rt1**2+rh1**2)) #cm
Uh1 = omega*rh1*0.01 #m/s
Um1 = omega*rm1*0.01 #m/s
beta1t = np.arctan(-Ut1/Vx1)*(180/np.pi) # deg
beta1h = np.arctan(-Uh1/Vx1)*(180/np.pi) # deg
beta1m = np.arctan(-Um1/Vx1)*(180/np.pi) # deg
print("Beta 1 are (tip,hub,mean): "+ str(beta1t)+" & "+str(beta1h)+" & "+str(beta1m)+" deg") 

#Compute meridional velocity gradient (Vx1 ratio from tip to hub )
print("Meridional velocity gradient is "+str(1))

#Find mean radius at the exist of impeller
dht = ht2 - ht1
Vtheta2 = dht/omega/(rm2*0.01) # m/sec required tangential outlet flow velocity of impeller
Ut2 = omega*rm2*0.01 #m/sec of blade tip speed at station 2 (at station 2 r_tip = r_hub = r_mean for impeller)
beta2 = np.arctan((Vtheta2-Ut2)/Vx2)*(180/np.pi) # deg
alpha2 = np.arctan(Vtheta2/Vx2)*(180/np.pi) #deg
V2 = np.sqrt(Vx2**2 + Vtheta2**2) # m/sec
print("Beta 2 is: "+ str(beta2)+" deg At Sta2 mean radius of: "+str(rm2)+" cm")
print("Alpha 2 is: "+ str(alpha2)+" deg At Sta2 mean radius of: "+str(rm2)+" cm")

#Compute impeller outlet area
Tt2 = Tt3
Gam2 = 0.5*(Gam1+Gam3) #Approximation
Cp2 = 0.5*(Cp1+Cp3) #Approximation J/kg/K
Rcons2 = 0.5*(Rcons1+Rcons3) #Approximation J/kg/K
T2 = Tt2 - (V2**2)/(2*Cp2) # K
P2 = Pt3*((T2/Tt2)**(Gam2/(Gam2-1))) #Approximation using Pt_3 [Pa]
rho2 = P2/(Rcons2*T2) #kg/m^3
AFlow2 = mdot/(Vx2*rho2)*10000 #cm^2
## Estimate the boundary layer thickness
pathHub1_2 = (np.pi/2)*(rm2-rh1) #90degturn
Rex2 = 3e6 #Approxiamte fully turbulent reynolds number (2nd flow)
dDelta2 = ((0.16-0.020)/(Rex2**(1/7)))*pathHub1_2 # cm boundary thickness - displacement thickness  
Hc2 = (AFlow2+2*np.pi*rm2*(2*dDelta2))/(2*np.pi*rm2 - NB1*(hB1+2*dDelta2) - NB2*(hB2+2*dDelta2)) #cm (underpredicted)
print("Station 2 Boundary Layer thickness is: "+str(dDelta2)+" cm")
print("Station 2 channel height is : "+ str(Hc2)+" cm with Vx2 of : "+str(Vx2)+" m/s")

#Compute outlet mach number
SOS2 = np.sqrt(Gam2*Rcons2*T2) #m/sec
W2 = np.sqrt((Vtheta2-Ut2)**2 + Vx2**2)
MachRel2 = W2/SOS2
MachAbs2 = V2/SOS2
print("Station 2 relative Mach is : "+ str(MachRel2)+" & absolute Mach is : "+str(MachAbs2))
print("Station 2 speed of sound is : "+str(SOS2)+" m/sec")

#Compute relative velocity ratio
BVR = W2/Ut2 
print("Blade relative velocity ratio is " + str(BVR))

#Compute station 1 static presssure
P1 = Pt1*(1+0.5*(Gam1-1)*Mx1**2)**(Gam1/(1-Gam1))
print("Station 1 static pressure is:" + str(P1) + " Pa")

#Compute station 2 static presssure
P2 = Pt3*(1+0.5*(Gam2-1)*MachAbs2**2)**(Gam2/(1-Gam2))
print("Station 2 static pressure is:" + str(P2) + " Pa")

#Compute station 3 static pressure (assuming no outlet swirl) (highly inefficient)
P3 = Pt3*(1+0.5*(Gam3-1)*Mx3**2)**(Gam3/(1-Gam3))
print("Station 3 static pressure is:" + str(P3) + " Pa")

#Compute reaction
Reaction = (P2-P1)/(P3-P1)
print("Reaction of the stage is: " + str(Reaction))

#Diffser outlet calculation (2 stage turning)
dDelta3 = 2*dDelta2 #Assume double the streamline path to station 3
rm3 = (AFlow3 + Hc2*NB3*(hB3 + 2*dDelta3))/(2*np.pi*Hc2 - 2*np.pi*(2*dDelta3))
print("Station 3 Boundary Layer thickness is: "+str(dDelta3)+" cm")
print("Station 3 radius is : "+ str(rm3)+" cm with Vx3 of : "+str(Vx3)+" m/s")

print("end")