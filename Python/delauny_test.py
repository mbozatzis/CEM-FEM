import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.spatial import Delaunay

a = 0.76
b = 1.75
V = 1


# Create Geometry
R = np.linspace(a, b, 50)
Theta = np.linspace(0, 2*np.pi, 150)

R_matrix, Theta_matrix = np.meshgrid(R, Theta)
X_n = R_matrix * np.cos(Theta_matrix)
Y_n = R_matrix * np.sin(Theta_matrix)
X = X_n.reshape(X_n.size)
Y = Y_n.reshape(Y_n.size)

points = np.zeros([X.shape[0], 2])
for i in range(1, X.shape[0]):
    points[i, 0] = X[i]
    points[i, 1] = Y[i]
tri = Delaunay(points)
neig = tri.neighbors

intersects = np.where(neig == -1)[0]
tri.simplices = np.delete(tri.simplices, intersects, axis=0)

plt.triplot(X, Y, tri.simplices)
plt.show()

