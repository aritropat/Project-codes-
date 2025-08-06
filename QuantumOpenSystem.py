import numpy as np
from pylab import *
import random

#for n=3 case, the eigenstates are 111,110,101,100,011,010,001,000, similar for larger values of n
n=input("length of array:") #total number of fermionic sites
N=2**n #dimension of the matrices
alpha=1.0
f=1.0
w=complex(0,1)
#w=1
epsilon=1.0
t=0.0 #initial time
delta_t=input("delta t:") #time increment, change from here
T=input("end time:") #end time
rho= (1.0/float(N))*np.identity(N) #initialise rho
print rho
H_0 =[[0 for j in range(N)] for k in range(N)] #initialise the hamiltonian
H=np.matrix(H_0) #in matrix form
D_0 =[[0 for j in range(N)] for k in range(N)] #initialise
D=np.matrix(D_0) #in matrix form 
c_matrices=[]
c_matrices_dagger=[]
trace=[] #trace of rho
trace_0=[]# trace of the lhs of equation
for i in range(n):   #each iteration produces the c_{i} matrices
	temp_matrix1=[[0 for j in range(N)] for k in range(N)] #for the annihilation operators		
	temp_matrix2=[[0 for j in range(N)] for k in range(N)] #for the creation operators
	pi=2**(n-i-1)
#1 means unoccupied, 0 means occupied
	for q in range(N): #bra state
		for r in range(N):#ket state 
			if (r//pi)%2 ==0 and (q//pi)%2 ==0 : #the ket and bra states are both occupied
				temp_matrix1[q][r]=0
			if (r//pi)%2 ==0 and (q//pi)%2 ==1 : #the ket occupied and bra unoccupied
				temp_matrix1[q][r]=1
			if (r//pi)%2 ==1 and (q//pi)%2 ==0 : #the ket unoccupied and bra occupied
				temp_matrix1[q][r]=0
			if (r//pi)%2 ==1 and (q//pi)%2 ==1 : #the ket and bra both unoccupied
				temp_matrix1[q][r]=0
	for y in range(N): #bra state	
		for a in range(N): #ket state 
			if (a//pi)%2 ==0 and (y//pi)%2 ==0 : #the ket and bra states are both occupied
				temp_matrix2[y][a]=0
			if (a//pi)%2 ==0 and (y//pi)%2 ==1 : #the ket occupied and bra unoccupied
				temp_matrix2[y][a]=0
			if (a//pi)%2 ==1 and (y//pi)%2 ==0 : #the ket unoccupied and bra occupied
				temp_matrix2[y][a]=1
			if (a//pi)%2 ==1 and (y//pi)%2 ==1 : #the ket and bra both unoccupied
				temp_matrix2[y][a]=0			
	
	temp1=np.matrix(temp_matrix1) #transform to matrix type	
	temp2=np.matrix(temp_matrix2) #transform to matrix type
	c_matrices.append(temp1)
	c_matrices_dagger.append(temp2)

for i in range(n-1):
	H+=alpha*(c_matrices_dagger[i]*c_matrices[i+1] + c_matrices[i]*c_matrices_dagger[i+1])+ f*(c_matrices_dagger[i]*c_matrices[i]-0.5*np.identity(N))

H+=c_matrices_dagger[n-1]*c_matrices[n-1]
L_1=(sqrt(epsilon))*c_matrices_dagger[0]
L_1dag=(sqrt(epsilon))*c_matrices[0]
L_2=(sqrt(epsilon))*c_matrices[-1]
L_2dag=(sqrt(epsilon))*c_matrices_dagger[-1]

#r=0
while t<T: #and r<10 :
#	if t==20*r:
#		print "progressing",r
#		r+=1
	D=2*L_1*rho*L_1dag-L_1dag*L_1*rho-rho*L_1dag*L_1+ 2*L_2*rho*L_2dag-L_2dag*L_2*rho-rho*L_2dag*L_2
	commutator=H*rho-rho*H
	rho= rho + delta_t*(D-(w*commutator))
	t+=delta_t
	trace.append(np.trace(rho))
	trace_0.append(np.trace(D-(w*commutator)))

final_density=[]
for i in range(n):
	final_density.append(np.trace(c_matrices_dagger[i]*c_matrices[i]*rho))

print c_matrices 
print c_matrices_dagger
print rho
print "the trace of the final rho matrix:", np.trace(rho)

#plot(final_density)
#xlabel('position on the chain')
#ylabel('occupation density')
#grid()
#show()



plot(trace_0)
xlabel("should show 0 on y axis")
show()

#plot(trace)
#xlabel('time')
#ylabel('trace')
#ylim(0.9,1000000000000)
#grid()
#show()
















