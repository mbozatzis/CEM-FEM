import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import scipy.sparse as sp

a = 0.76
b = 1.75
V = 1


# Create Geometry
R = np.linspace(a, b, 50)
Theta = np.linspace(0, 2*np.pi, 100)

R_matrix, Theta_matrix = np.meshgrid(R, Theta)
X_n = R_matrix * np.cos(Theta_matrix)
Y_n = R_matrix * np.sin(Theta_matrix)
X = X_n.reshape(X_n.size)
Y = Y_n.reshape(Y_n.size)

triang = mtri.Triangulation(X, Y)
tri_coords = np.zeros((triang.triangles.shape[0], 3, 2))
for i, tri in enumerate(triang.triangles):
    tri_coords[i] = np.column_stack((X[tri], Y[tri]))

mask = np.zeros(triang.triangles.shape[0], dtype=bool)
for i in range(0, triang.triangles.shape[0]):
    t_c = tri_coords[i]
    x1 = t_c[0][0]
    x2 = t_c[1][0]
    x3 = t_c[2][0]
    y1 = t_c[0][1]
    y2 = t_c[1][1]
    y3 = t_c[2][1]
    r1 = np.sqrt(pow(x1, 2) + pow(y1, 2))
    r2 = np.sqrt(pow(x2, 2) + pow(y2, 2))
    r3 = np.sqrt(pow(x3, 2) + pow(y3, 2))
    if r1-a <= pow(10,-6) and r2-a <= pow(10,-6) and r3-a <= pow(10,-6):
        mask[i] = True
triang.set_mask(mask)

nodes = np.column_stack((triang.x, triang.y))
unique_nodes = np.unique(nodes, axis=0)
element_nodes = triang.triangles


# Define known and unknown potentials
node_id = np.ones(nodes.shape[0])
Potentials = np.zeros(nodes.shape[0])
for inode in range(0, nodes.shape[0]):
    radius = np.sqrt(pow(nodes[inode][0],2) + pow(nodes[inode][1], 2))
    if abs(radius-a) <= pow(10,-6):
        node_id[inode] = 0
        Potentials[inode] = V
    elif abs(radius-b) <= pow(10,-6):
        node_id[inode] = 0
        Potentials[inode] = 0

# Index of unkowns
un_index = np.zeros(nodes.size, dtype=int)
counter = 0
for inode in range(0, nodes.shape[0]):
    if node_id[inode] == 1:
        un_index[inode] = counter
        counter = counter + 1
      


# Matrix Calculation
Nf = counter
Sff = sp.lil_matrix((Nf, Nf), dtype = float)
B = np.zeros(Nf)


for ie in range(1, triang.triangles.shape[0]):
    t_c = tri_coords[ie]
    x = np.zeros(3)
    y = np.zeros(3)
    x = [t_c[0][0], t_c[1][0], t_c[2][0]]
    y = [t_c[0][1], t_c[1][1], t_c[2][1]]
    n = [element_nodes[ie][0], element_nodes[ie][1], element_nodes[ie][2]]
        
    A = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
    De = np.linalg.det(A)
    Ae = abs(De/2)
    b = np.zeros(3)
    c = np.zeros(3)
    b = [(y[1]-y[2])/De, (y[2]-y[0])/De, (y[0]-y[1])/De]
    c = [(x[2]-x[1])/De, (x[0]-x[2])/De, (x[1]-x[0])/De]

    Se = np.zeros([3, 3])
    for i in range(0, 3):
        for j in range(0, 3):
            Se[i][j] = (b[i]*b[j] + c[i]*c[j])*Ae

            if node_id[n[i]] == 1:
                if node_id[n[j]] == 1:
                    Sff[un_index[n[i]], un_index[n[j]]] = Sff[un_index[n[i]], un_index[n[j]]] + Se[i][j]
                else:
                    B[un_index[n[i]]] = B[un_index[n[i]]] - Se[i][j]*Potentials[n[j]]


# Unkown Potential Calculation
P = sp.linalg.bicgstab(Sff, B)
pot = P[0]
for inode in range(0, nodes.shape[0]):
    if node_id[inode] == 1:
        Potentials[inode] = pot[un_index[inode]]


# Plot the results
fig, ax = plt.subplots()
ax.triplot(triang)
cax = ax.tripcolor(triang, Potentials, cmap='plasma', shading='flat')
fig.colorbar(cax, ax=ax)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()


