#!/usr/bin/env python
# coding: utf-8

# In[1]:


from qiskit import *
from qiskit.circuit.library import QFT
from qiskit.quantum_info import Operator
from qiskit.circuit.library.standard_gates import *
import numpy as np
import cmath


# In[122]:


def R(i,c,qc,qubit):
    qg = QuantumCircuit(1)
    u = np.array([[1,0],[0,np.exp(1j*c*np.pi/(2**i))]])
    u = Operator(u)
    qg.append(u,[0])
    custom = qg.to_gate()
    custom.label = f"R{i,c}"
    qc.append(custom,[qubit])
    return qc
def C_R(i,c,control,target,qc):
    qg = QuantumCircuit(1)
    u = np.array([[1,0],[0,np.exp(1j*c*np.pi/(2**i))]])
    u = Operator(u)
    qg.append(u,[0])
    custom = qg.to_gate()
    custom.label =f'R{i,c}'
    custom_2 = custom.control(1)
    qc.append(custom_2,[control,target])                     
    return qc                     


# In[123]:


def FADD(c,y,z,qc,ny,nz) :
    for i in range(0,ny):
        k=i
        for j in range(0,nz-i):
            C_R(j,c,y[i],z[k],qc)
            k = k+1
    return qc        
def C_FADD(c,y,z,qc,ny,nz,control):
    y1 = QuantumRegister(ny)
    z1 = QuantumRegister(nz)
    qg = QuantumCircuit(y1,z1)
    FADD(c,y1,z1,qg,ny,nz)
    circuit = qg.to_gate()
    circuit.label = f"FADD{c}"
    circuit_2 = circuit.control(1)
    qc.append(circuit_2, [control] + list(y) + list(z))    
    return qc


# In[124]:


def QFMA(x,y,z,nx,ny,nz,a,b,c,d):
    for i in range(0,nx):
      C_FADD((2**(i))*a,y,z,qc,ny,nz,x[i])
    FADD(b,x,z,qc,nx,ny)
    FADD(c,y,z,qc,ny,nz)
    for i in range(0,nz):
        R(i,d,qc,z[i])
    return qc 
def FQNMA(x,y,z,nx,ny,nz,a,b,c,d):
    x1 = QuantumRegister(nx)
    y1 = QuantumRegister(ny)
    z1 = QuantumRegister(nz)
    qg = QuantumCircuit(x1,y1,z1)
    for i in range(0,nx):
      C_FADD((2**(i))*a,y1,z1,qg,ny,nz,x1[i])
    FADD(b,x1,z1,qg,nx,ny)
    FADD(c,y1,z1,qg,ny,nz)
    for i in range(0,nz):
        R(i,d,qg,z1[i])
    iqg = qg.inverse()
    circuit = iqg.to_gate()
    circuit.label="FQNMA"
    qc.append(circuit,list(x)+list(y) + list(z))
    return qc


# In[125]:


def diff_z(qc,nz,z):
    qft =QFT(num_qubits=nz).to_gate()
    qft.label = "QFT"
    iqft = qft.inverse()
    iqft.label = "IQFT"
    t = QuantumCircuit(1)
    t.z(0)
    f = t.to_gate()
    f.label = "Z"
    mcz = f.control(nz-1)
    
    qc.append(iqft,list(z))
    for i in range(0,nz):
        qc.x(z[i])
    qc.append(mcz,[z[i] for i in range(nz-1,-1,-1)])
    for i in range(0,nz):
        qc.x(z[i])
    qc.append(qft,list(z) )   
    return qc    
def diff_xy(nx,ny,x,y,qc):
    t = QuantumCircuit(1)
    t.z(0)
    f = t.to_gate()
    f.label = "Z"
    mcz = f.control(nx+ny-1)
    for i in range(0,nx+ny):
        qc.ry(np.pi/2,i)
    qc.append(mcz,range(0,nx+ny))
    for i in range(0,nx+ny):
        qc.ry(-np.pi/2,i)
    return qc    


# In[126]:


def g_oracle(qc,x,y,z,nx,ny,nz,a,b,c,d):
    qft =QFT(num_qubits=nz).to_gate()
    qft.label = "QFT"
    iqft = qft.inverse()
    iqft.label = "IQFT"
    qc.append(qft,list(z))
    FQNMA(x,y,z,nx,ny,nz,a,b,c,d)
    diff_z(qc,nz,z)
    QFMA(x,y,z,nx,ny,nz,a,b,c,d)
    qc.append(iqft,list(z))
    return qc


# In[127]:


nx = 3
ny = 3
nz = 4
x = QuantumRegister(nx,'x')
y= QuantumRegister(ny,'y')
z = QuantumRegister(nz,'z')
cr = ClassicalRegister(nx+ny)
qc = QuantumCircuit(x,y,z,cr)
for i in range(nx+ny):
    qc.h(i)
a = 6
b=5
c=7
d=5
g_oracle(qc,x,y,z,nx,ny,nz,a,b,c,d)
diff_xy(nx,ny,x,y,qc)
for i in range(0,nx+ny):
    qc.measure([i],cr[i])    


# In[128]:


qc.draw()


# In[112]:


nx = 3
ny = 3
nz = 4
x = QuantumRegister(nx,'x')
y= QuantumRegister(ny,'y')
z = QuantumRegister(nz,'z')
cr1 = ClassicalRegister(nx)
cr2 = ClassicalRegister(ny)
qc = QuantumCircuit(x,y,z,cr1,cr2)
for i in range(nx+ny):
    qc.h(i)
N = 35
S = 3/2 - 1/2*(N%6)
M = (N-S)/6 -1
#Encoding M in z register
#0101
s = 1
qc.x(z[0])
qc.x(z[2])
a = 6
b= 6+s
c= 6+s*S
d= 5+s+s*S
k = (np.pi/4)*(2)**(nx+ny)
for k in range(0,int(k)):
  g_oracle(qc,x,y,z,nx,ny,nz,a,b,c,d)
  diff_xy(nx,ny,x,y,qc)
for i in range(0,nx):
    qc.measure(x[i],cr1[i])  
for i in range(0,nx):
    qc.measure(x[i],cr2[i])      


# In[113]:


simulator = Aer.get_backend('qasm_simulator')
job = execute(qc, backend=simulator, shots=20000)
result = job.result()
counts = result.get_counts()


# In[114]:


print('---')
print(counts)


# In[ ]:





# In[ ]:




