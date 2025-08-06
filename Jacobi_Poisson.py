from pylab import *

# plots the number of iterations by varying the tolerance  for the jacobi algorithm. Part b of the problem

initial_iterations=50   # the number of iterations after which the convergence criterion should be applied
charge=1.0 # (Q/epsilon_{0})*(unit grid area)
R=10	#radius
N=90    #change the number of grid points from here
tolerance= linspace(10**(-5),10**(-6),10) #array of the tolerance values
iterations1=zeros(10)
iteration=array(iterations1,dtype=float)
x=linspace(-R,R,N)
y=linspace(-R,R,N)

a_1=0
a_2=0
b_1=0
b_2=0
for i in range(len(x)):
	if x[i]<-0.3 and x[i+1]>-0.3:
		a_1=i        # we will put the negative charge at -0.3
		a_2=(i+1)    
	if x[i]<0.3 and x[i+1]>0.3:
		b_1=i         
		b_2=(i+1)    # we will put the positive charge at +0.3

def radial(a,b):
	var1=sqrt(a*a+b*b)
	return var1

def distance1(a,b):
	var2=sqrt((a+0.3)**2+b**2)
	return var2

def distance2(a,b):
	var3=sqrt((a-0.3)**2+b**2)
	return var3

def average_array(X):
	summing=0.0
	for i in range(N):
		for j in range(N):
			summing=summing+abs(X[i][j])
	variable= summing/(N*N)
	return variable

for p in range(len(tolerance)):
	V1=zeros(shape=(N,N))
	V=array(V1,dtype=float)    #initialising the potential to 0, initially.... The first coordinate of V gives the y coordinate and the 2nd coordinate gives the x coordinate. This is done so because later while doingthe contour plot, meshgrid does exactly the same thing with its grid.	
	it=0
	initial_average=0.0	
	list_of_averages=[]  #will store all the differences of average values here
	list_of_averages.append(initial_average)
	for m in range(initial_iterations):
		it+=1
		temp1=zeros(shape=(N,N))
		temp=array(temp1,dtype=float)
		for i in range(N):
			for j in range(N):
				if radial(x[i],y[j])<R: #only points within the circle get updated. All other points outside remain 0 at all times.
					if i==a_1 and j==int(0.5*(N-1)):
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i]-charge)	
					elif i==b_2 and j==int(0.5*(N-1)):
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i]+charge)
					else:
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i])
		list_of_averages.append(average_array(temp-V))
		V=temp  # assigning V to the values in temp in this iteration step 

	while (list_of_averages[-1])>tolerance[p]:
		it+=1		
		temp1=zeros(shape=(N,N))
		temp=array(temp1,dtype=float)
		for i in range(N):
			for j in range(N):
				if radial(x[i],y[j])<R: #only points within the circle get updated. All other points outside remain 0 at all times.
					if i==a_1 and j==int(0.5*(N-1)):
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i]-charge)	
					elif i==b_2 and j==int(0.5*(N-1)):
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i]+charge)
					else:
						temp[j][i]=0.25*(V[j][i+1]+ V[j][i-1]+ V[j-1][i]+ V[j+1][i])
						
		list_of_averages.append(average_array(temp-V))
		V=temp  # assigning V to the values in temp in this iteration step 
	iteration[p]=it
	print iteration[p]

plot(tolerance,iteration,'r-^')  #plotting the figure
xlabel('Tolerance values')
ylabel('Number of Iterations required')
title('Plot of Iterations versus tolerance for the Jacobi Algorithm')
grid()
show()













