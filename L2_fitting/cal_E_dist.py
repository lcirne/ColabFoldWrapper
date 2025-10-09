import numpy as np
def compute_E(distances, R_0=51):
    for i in range(len(distances)):
        distances[i] = 1 / (1 + (distances[i] / R_0)**6)
    return distances

R_0 = 51

data=np.loadtxt("open.input")
dist=data[:,1]
E=compute_E(dist,R_0)
print (f"{E}")

