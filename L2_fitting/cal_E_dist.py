import numpy as np
def compute_E(R_DA, R_0):
    E = 1 / (1 + (R_DA / R_0)**6)
    return E

R_0 = 51

data=np.loadtxt("open.input")
dist=data[:,1]
E=compute_E(dist,R_0)
print (f"{E}")

