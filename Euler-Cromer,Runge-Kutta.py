from pylab import *

g=9.8
l=9.8
C=g/l            #I'll use this variable x in the equations, with C=1
# I decided not to write x anywhere and just used unity in place of it, to make things simpler

gamma=0.25       
alpha_D=0.2
dt=0.01          # differential time step
t_max=100.00     # we can change t_max(in seconds) here and nochange is needed in the rest of the code
total_steps=int(t_max/dt)   #total number of time steps,10000


#implementing the Euler Cromer method first:

#time_E=[0.0]*total_steps
#theta_E=[0.0]*total_steps
#omega_E=[0.0]*total_steps
#omega_0=0.90              #driving frequency,taken close to the resonant frequency

#theta_E[0]=0.3      #initial conditions. Of course, we had already initialised with all 0's,but can use any value, keeping \theta within the small angle approximation
#omega_E[0]=0.0      #'''''''''''''''''''''''''''
#time_E[0]= 0.0      #''''''''''''''''''''''''''''

#for i in range(total_steps-1):
#	time_E[i+1]= time_E[i] + dt        #updating time
#	omega_E[i+1]=omega_E[i] +dt*(-theta_E[i]-0.50*omega_E[i]+alpha_D*sin(omega_0*time_E[i]))  #updating omega
#	theta_E[i+1]= theta_E[i] + dt*(omega_E[i+1])  #updating theta, with the Euler Cromer Method

#plot(time_E,theta_E,'r-')    #can use this to separately check for the solution with the Euler Cromer
#ylabel('theta')
#xlabel('time')
#title("Variation of the angle with time in the Euler Cromer Method")
#grid()
#show()


#Now implementing the Runge Kutta 4th order method:

#time_R=[0.0]*total_steps
#theta_R=[0.0]*total_steps
#omega_R=[0.0]*total_steps
#omega_0=0.90       #driving frequency, taken close to the resonant frequency 

#theta_R[0]=0.2      #initial conditions. Of course, we had already initialised with all 0's,but can use any value, keeping \theta within the small angle approximation
#omega_R[0]=0.4      #'''''''''''''''''''''''''''
#time_R[0]= 0.0      #'''''''''''''''''''''''''''' 

#for i in range(total_steps-1):
#	k_1=dt*omega_R[i]
#	l_1=dt*(-theta_R[i]-0.50*omega_R[i]+alpha_D*sin(omega_0*time_R[i]))
#	k_2=dt*(omega_R[i] + 0.5*l_1) 
#	l_2=dt*(-(theta_R[i]+0.5*k_1)-0.50*(omega_R[i]+0.5*l_1)+alpha_D*sin(omega_0*(time_R[i]     +0.50*dt)))
#	k_3=dt*(omega_R[i] + 0.5*k_2)
#	l_3=dt*(-(theta_R[i]+0.5*k_2)-0.50*(omega_R[i]+0.5*l_2)+alpha_D*sin(omega_0*(time_R[i]+0.50*dt)))
#	k_4=dt*(omega_R[i] + k_3)
#	l_4=dt*(-(theta_R[i]+k_3)-0.50*(omega_R[i]+l_3)+alpha_D*sin(omega_0*(time_R[i]+dt)))

#	time_R[i+1]= time_R[i] + dt     #updating time
#	theta_R[i+1]=theta_R[i] + (1.0/6.0)*(k_1+2*k_2+2*k_3+k_4) #updating \theta
#	omega_R[i+1]=omega_R[i] + (1.0/6.0)*(l_1+2*l_2+2*l_3+l_4) #updating \omega

#plot(time_R,omega_R,'r-')    #can use this to separately check for the solution with the Runge Kutta
#ylabel('angular velocity(omega)')
#xlabel('time')
#title("Variation of the angular velocity with time in the Runge Kutta 4th order method")
#grid()
#show()

freq_array=[0.20+0.05*i for i in range(40)]   #array of 40 frequencies close to the resonant frequency. the starting value and the interval can be chosen suitably to show the frequency range that is needed.

def resonant_structure(freq_array):
	u=[0.0]*40                       #array for amplitudes at each of the 40 frequency
	v=[0.0]*40                           #array for phases at each of the 40 frequency
	for m in range(40):
		time_R=[0.0]*total_steps
		theta_R=[0.0]*total_steps
		omega_R=[0.0]*total_steps 

		theta_R[0]=0.2      #initial conditions. Of course, we had already initialised with all 0's,but can use any value, keeping \theta within the small angle approximation
		omega_R[0]=0.4      #'''''''''''''''''''''''''''
		time_R[0]= 0.0      #'''''''''''''''''''''''''''' 

		for i in range(total_steps-1):
			k_1=dt*omega_R[i]
			l_1=dt*(-theta_R[i]-2*gamma*omega_R[i]+alpha_D*sin(freq_array[m]*time_R[i]))                                 #here, freq_array[m] is the m'th driving frequency
			k_2=dt*(omega_R[i] + 0.5*l_1) 
			l_2=dt*(-(theta_R[i]+0.5*k_1)-2*gamma*(omega_R[i]+0.5*l_1)+alpha_D*sin(freq_array[m]*(time_R[i]+0.50*dt)))   #here, freq_array[m] is the m'th driving frequency
			k_3=dt*(omega_R[i] + 0.5*k_2)
			l_3=dt*(-(theta_R[i]+0.5*k_2)-2*gamma*(omega_R[i]+0.5*l_2)+alpha_D*sin(freq_array[m]*(time_R[i]+0.50*dt)))   #here, freq_array[m] is the m'th driving frequency
			k_4=dt*(omega_R[i] + k_3)
			l_4=dt*(-(theta_R[i]+k_3)-2*gamma*(omega_R[i]+l_3)+alpha_D*sin(freq_array[m]*(time_R[i]+dt)))		  #here, freq_array[m] is the m'th driving frequency

			time_R[i+1]= time_R[i] + dt     #updating time
			theta_R[i+1]=theta_R[i] + (1.0/6.0)*(k_1+2*k_2+2*k_3+k_4) #updating \theta
			omega_R[i+1]=omega_R[i] + (1.0/6.0)*(l_1+2*l_2+2*l_3+l_4) #updating \omega	
		
#finding the amplitude for the m'th driving frequency
		w=[]  #array keeping the average  
		j=5000 
		while j<9998:    #taking 5000 as a suitable time, hopefully the homogeneous solution has died down completely by then. the j values run till 9999, and we have (j+2) below, hence used < 9998
			if theta_R[j+1]>theta_R[j] and theta_R[j+2]<theta_R[j+1]:   #if this is true, the maximum definitely lies between j'th time step and j+2'th time step , and we don't know where, so we approximate it by theta_R[j+1]. We find all the maximas and minimas and average them later.
				w.append(theta_R[j+1])
			elif theta_R[j+1]<theta_R[j] and theta_R[j+2]>theta_R[j+1]:  #if this is true, the minimum definitely lies between j'th time step and j+2'th time step,and we approximate the minimum by theta_R[j+1]
				w.append(-(theta_R[j+1]))
			j=j+1
		
		s=0.0
		for x in range(len(w)):
			s=s+w[x]            #summing all the values in w
		u[m]=s/len(w)               #for the i'th driving frequency, the average is taken here, to give the amplitude

		z=1		
		while z<400:  #400 is an arbitrarily large chosen number
			z_3=z*((2*pi)/freq_array[m])   #z is an integer, (2*pi)/freq_array[m]) is the time period for the m'th driving frequency, z_3 is thus an integer multiple of the time period
			if z_3>50:    #taking suitable number of time periods above 5000 time steps(more than 50 seconds), when again, we assume the homogeneous solution has effectively died down to give no contribution
				z_1="%.2f" % z_3            #nearest time step corresponding to the integer number of time steps above 5000
				if (float(z_1)-z_3)>0 :     #finding the other time step, along with the former, within which z*((2*pi)/t[m]) falls
					z_2=float(z_1) -0.01
				elif (float(z_1)-z_3)<0 :   #''''''''''''''''''
					z_2=float(z_1)+0.01
				elif (float(z_1)-z_3)==0 :  #'''''''''''''''''''
					z_2=float(z_1)
				z__1=int(float(z_1)*100)   #convert to an integer time step
				z__2=int(z_2*100)          #convert to an integer time step  
				v[m]=(180.0/(pi))*arcsin((-theta_R[z__1]-theta_R[z__2])/(2.0*(u[m])))   #phase of the steady state solution at the m'th driving frequency. we average the values of the two to find the approximate value of sin(\phi) where \phi is the phase, and hence we require the arcsin
				break
			else:
				z=z+1
	return u,v    #return the array of the amplitude and phase values for the different driving frequencies	

amplitude_array,phase_array=resonant_structure(freq_array)

plot(freq_array,phase_array)
ylabel('phase $\\phi(\\Omega_D)$')
xlabel('frequency $\\Omega_D$(Hz)')
title("Variation of the amplitude of oscillation with driving frequency")
grid()
show()


q=0
while q< len(amplitude_array)-1:   
#	print "a"  #just to check if this while loop is at all running
	if amplitude_array[q+1]< amplitude_array[q]:
		max=amplitude_array[q]
		break
	else:
		q+=1

mid=(max+ amplitude_array[0])/2.0

r=0
while r< len(amplitude_array)-1:
#	print "b"
	if amplitude_array[r]<= mid and amplitude_array[r+1]> mid:
		extreme_1=0.5*(freq_array[r+1]+freq_array[r])   # one end of fwhm in frequency spectrum
		r+=1
	elif amplitude_array[r]>= mid and amplitude_array[r+1]< mid:
		extreme_2=0.5*(freq_array[r+1]+freq_array[r])   #other end of fwhm in frequency spectrum
		r+=1
	else:
		r+=1
FWHM=extreme_2-extreme_1

print "The Full Width for $\\gamma=0.25$ at Half Maximum is: %.2f (Hz)" %FWHM   # it came out to be 0.45, it is taking around 5 minutes to get the result 
        
    



