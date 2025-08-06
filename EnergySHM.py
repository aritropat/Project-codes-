from pylab import *

#Studies energy in simple harmonic motion.

mass=1.0  #the mass has not been mentiond anywhere, I take it to be a unit mass.
g=9.8
l=9.8
gamma=0.25       
alpha_D=0.2
dt=0.01          # differential time step
t_max=50.00     # we can change t_max(in seconds) here and nochange is needed in the rest of the code
total_steps=int(t_max/dt)   #total number of time steps,10000




time_R=[0.0]*total_steps
theta_R=[0.0]*total_steps
omega_R=[0.0]*total_steps
kinetic_energy_R=[0.0]*total_steps
potential_energy_R=[0.0]*total_steps
total_energy_R=[0.0]*total_steps
omega_0=0.95      #driving frequency, taken close to the resonant frequency 

theta_R[0]=0.0      #initial conditions. Of course, we had already initialised with all 0's,but can use any value, keeping \theta within the small angle approximation
potential_energy_R[0]=0.50*mass*g*l*theta_R[0]*theta_R[0]# initial potential energy
omega_R[0]=0.0      #initial \omega
kinetic_energy_R[0]=0.50*mass*l*l*omega_R[0]*omega_R[0] #initial kinetic energy
total_energy_R[0]=kinetic_energy_R[0]+potential_energy_R[0] #initial total energy
time_R[0]= 0.0      #'''''''''''''''''''''''''''' 

for i in range(total_steps-1):
	k_1=dt*omega_R[i]
	l_1=dt*(-theta_R[i]-0.50*omega_R[i]+alpha_D*sin(omega_0*time_R[i]))
	k_2=dt*(omega_R[i] + 0.5*l_1) 
	l_2=dt*(-(theta_R[i]+0.5*k_1)-0.50*(omega_R[i]+0.5*l_1)+alpha_D*sin(omega_0*(time_R[i]     +0.50*dt)))
	k_3=dt*(omega_R[i] + 0.5*k_2)
	l_3=dt*(-(theta_R[i]+0.5*k_2)-0.50*(omega_R[i]+0.5*l_2)+alpha_D*sin(omega_0*(time_R[i]+0.50*dt)))
	k_4=dt*(omega_R[i] + k_3)
	l_4=dt*(-(theta_R[i]+k_3)-0.50*(omega_R[i]+l_3)+alpha_D*sin(omega_0*(time_R[i]+dt)))

	time_R[i+1]= time_R[i] + dt     #updating time
	theta_R[i+1]=theta_R[i] + (1.0/6.0)*(k_1+2*k_2+2*k_3+k_4) #updating \theta
	omega_R[i+1]=omega_R[i] + (1.0/6.0)*(l_1+2*l_2+2*l_3+l_4) #updating \omega
	kinetic_energy_R[i+1]=0.50*mass*l*l*omega_R[i+1]*omega_R[i+1]#updating K.E
	potential_energy_R[i+1]=0.50*mass*g*l*theta_R[i+1]*theta_R[i+1]#updating P.E
	total_energy_R[i+1]=kinetic_energy_R[i+1]+potential_energy_R[i+1]#updating total energy

plot(time_R,potential_energy_R,'r-',label="Potential energy")    #can use this to separately check for the solution with the Runge Kutta
plot(time_R,kinetic_energy_R,'b-',label="Kinetic Energy")
plot(time_R,total_energy_R,'g-',label="Total Energy")
xlabel('Time(Seconds)')
ylabel('Energy(Joule)')
title('Energy variation with time')
legend(bbox_to_anchor=(0.251, 1.015))
grid()
show()	
		
