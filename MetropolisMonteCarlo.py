from pylab import *
import random
#FOR FINDING THE MAGNETISATION
time_0=8*(10**6)  #number of runs initially for the initial value of temperature
time_1=5*(10**4)  # number of runs by when it's assumed that equilibrium has been reached, not veryuseful
time_2=8*(10**6) #  number of runs after equilibriation, till the end
n=50        #length of grid
#n=input("Enter the grid length: ")
N=n*n       #total number of spins
J=1.5       #coupling constant
temp_1=1.5  #lowest temperature value
temp_2=5.5  #highest temperature values
K=300         # number of different temperature values
print ("there will be %d Temperature values" %K)
temp=linspace(temp_1,temp_2,K)  #array of temperature values
E=[] #stores average energy for different temperature values
E_2=[] #stores average of energy squared for different temperature values
M=[] #stores average magnetisation for different temperature values

def periodic(i,j):   #defines the periodic boundary conditions
	if i==0:
		i1=(i+1)   #can move one step right
		i2=(n-1)   #periodic boundary, far right edge taken here
	elif i==(n-1):
		i1=0     #periodic boundary, far left edge taken here 	
		i2=(i-1)   #can move one step left
	else:
		i1=(i+1)   #the two sides
		i2=(i-1)	 #'''''''''''''			

	if j==0:
		j1=(j+1)     #can move one step up
		j2=(n-1)   #periodic boundary, far upper edge taken here
	elif j==(n-1):
		j1=0       #periodic boundary, far bottom edge taken here
		j2=(j-1)     #can move one step down
	else:
		j1=(j+1)    #the two sides
		j2=(j-1)    #'''''''''''''
	p=[i1,i2]         #returns the x coordinates of the adjacent ones
	q=[j1,j2]         #returns the y coordinates of the adjacent ones
	return p,q

def secondmax(a):   #returns the 2nd largest element in an array
	tt=[]    
	for i in range(len(a)):
		tt.append(a[i])
	tt.remove(max(tt))
	return max(tt)

def energy(spin):         # calculates the 
	double_energyy=0.0
	for i in range(n):
		for j in range(n):
			double_energyy=double_energyy-J*spin[i,j]*(spin[periodic(i,j)[0][0],periodic(i,j)[1][0]]+spin[periodic(i,j)[0][0],periodic(i,j)[1][1]]+spin[periodic(i,j)[0][1],periodic(i,j)[1][0]]+spin[periodic(i,j)[0][1],periodic(i,j)[1][1]])
	energyy=0.5*double_energyy    #every spin pair has been counted twice in the earlier step, so needs to be halved here
	return energyy	

def square(x):   #simply returns the square of a number
	y=x*x
	return y

def delta_E(spin,i,j):   #calls the spin grid, and the i,j th coordinate
	a=2*J*spin[i,j]*(spin[periodic(i,j)[0][0],periodic(i,j)[1][0]]+spin[periodic(i,j)[0][0],periodic(i,j)[1][1]]+spin[periodic(i,j)[0][1],periodic(i,j)[1][0]]+spin[periodic(i,j)[0][1],periodic(i,j)[1][1]])  
	return a

def magnetisation(spin):   #returns the initial magnetisation
	summ=0
	for i in range(n):
		for j in range(n):
			summ=summ+spin[i,j]	
	return summ

spin=zeros((n,n)) # defining the array of spins
#a=[-1,1]           #the two distinct spin values
for k in range(n):
	for l in range(n):
		spin[k,l]=1  #assigns spins initially
for m in range(time_0):   #initial relaxation for the first temperature value
	i=random.choice(range(0,n))
	j=random.choice(range(0,n))
	energy_change=delta_E(spin,i,j)   # the change in energy upon flipping the coordinate
	if energy_change<=0:  #flip certainly if you're going to a lower energy state
		spin[i,j]=-spin[i,j]
	elif energy_change>0:  
		r=random.random()
		pp=-(energy_change)/(temp[0])
		lim=exp(pp)
		if r<=lim:     #flip only with this probability. else remain same
			spin[i,j]=-spin[i,j]

for ll in range(len(temp)):  
	print ("the temperature number is:%d",ll)
	energy_array=[]         # stores the energies at the different steps after time_1
	energy_square_array=[]  # stores the energy squared at the different steps after time_1
	magnetisation_array=[]  # stores the magnetisation at the different steps after time_1
	for m in range(time_1):   # begin the Metropolis algorithm,letting it equilibriate till time_1
		i=random.choice(range(0,n))
		j=random.choice(range(0,n))
		energy_change=delta_E(spin,i,j)   # the change in energy upon flipping the coordinate
		if energy_change<=0:  #flip certainly if you're going to a lower energy state
			spin[i,j]=-spin[i,j]
		else:  
			r=random.random()
			pp=-(energy_change)/(temp[ll])
			lim=exp(pp)
			if r<=lim:     #flip only with this probability. else remain same
				spin[i,j]=-spin[i,j]
	mm=magnetisation(spin)
	magnetisation_array.append(mm)
	for m in range(time_2):  #update energy and magnetisation now
		i=random.choice(range(0,n))
		j=random.choice(range(0,n))
		energy_change=delta_E(spin,i,j)   # the change in energy upon flipping the coordinate
		if energy_change<=0:  #flip certainly if you're going to a lower energy state
			spin[i,j]=-spin[i,j]
			mm=mm-2*spin[i,j] #change in magnetisation
			magnetisation_array.append(mm)   #append the new magnetisation value
		else:                    #flip only with a finite probability 
			r=random.random()
			qq=-(energy_change)/(temp[ll])
			lim=exp(qq)
			if r<=lim:     #flip only with this probability. 
				spin[i,j]=-spin[i,j]
				mm=mm-2*spin[i,j] #change in magnetisation
				magnetisation_array.append(mm)   #append the new magnetisation value
			else:
				magnetisation_array.append(mm)   #append the old magnetisation value
	avg_mag=average(magnetisation_array)    #average canonical magnetisation for this temeprature
	M.append(avg_mag) # need the absolute value here.
diff_M=[]
for k in range(K-1):
	diff_M.append(abs(M[k]-M[k+1]))
index_1=diff_M.index(max(diff_M))   #the index in array temp where there's maximum change in temperature
index_2=diff_M.index(secondmax(diff_M)) 
T_1=0.5*(temp[index_1]+temp[index_1+1])
T_2=0.5*(temp[index_2]+temp[index_2+1])
#T_3=0.5*(T_1 + T_2)

#print "this is M:", M
#print "this is diff_M:", diff_M
#print "the first estimate of the critical temperature is:",T_1
#print "the second estimate of the critical temperature is:",T_2
#print "the third estimate of the critical temperature is:",T_3

plot(temp,M)
xlabel('Temperature (K)')
ylabel('Net Magnetisation')
title('Magnetisation versus Temperature in Ising Model, n=50')
grid()
show()

			

