import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ortho_group

n = 10
b = np.random.randint(1000, size=n)
Q = ortho_group.rvs(3)
EigMat = np.diag(np.random.randint(1000, size=n))
A = np.matmul(Q,EigMat,Q.T)
x = np.random.randint(1000, size=n)
r = b-np.dot(A,x)

# print(x_0)
