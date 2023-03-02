import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import scipy.sparse as sp

w = 4
t = 0.2
d = 1
V = 100
er = 2.2

A = 2*w
B = w


# Create Geometry
x_n = np.linspace(-A/2, A/2, 100)
y_n = np.linspace(-B/2, B/2, 100)
X_n, Y_n = np.meshgrid(x_n, y_n)

X = X_n.reshape(X_n.size)
Y = Y_n.reshape(Y_n.size)

triang = mtri.Triangulation(X, Y)

tri_coords = np.zeros((triang.triangles.shape[0], 3, 2))
for i, tri in enumerate(triang.triangles):
    tri_coords[i] = np.column_stack((X[tri], Y[tri]))

mask = np.zeros(triang.triangles.shape[0], dtype=bool)
for i in range(0, triang.triangles.shape[0]):
    t_c = tri_coords[i]
    x = [t_c[0][0], t_c[1][0], t_c[2][0]]
    y = [t_c[0][1], t_c[1][1], t_c[2][1]]
    # Capacitor Plates:
    # x: [-w/2, w/2]
    # y: [-d/2-t, -d/2] and [d/2, d/2+t]
    flag = np.zeros(3, dtype = 'bool')
    for j in range(0, 3):
        flag[j] = (x[j] >= -w/2 and x[j] <= w/2) and ((y[j] >= -d/2-t and y[j] <= -d/2) or (y[j] >= d/2 and y[j] <= d/2 + t))
    if flag[0] and flag[1] and flag[2]:
        mask[i] = True
triang.set_mask(mask)
nodes = np.column_stack((triang.x, triang.y))
unique_nodes = np.unique(nodes, axis=0)
element_nodes = triang.triangles


# Define known and unknown potentials
node_id = np.ones(nodes.shape[0])
Potentials = np.zeros(nodes.shape[0])
for inode in range(0, nodes.shape[0]):
    x = nodes[inode][0]
    y = nodes[inode][1]
    if (x >= -w/2 - pow(10,-6) and x <= w/2 + pow(10,-6)):
        if (y >= -d/2- t - pow(10,-6) and y <= -d/2 + pow(10,-6)):
            node_id[inode] = 0
            Potentials[inode] = -V/2
        elif (y >= d/2 - pow(10,-6) and y <= d/2 + t + pow(10,-6)):
            node_id[inode] = 0
            Potentials[inode] = V/2

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

for ie in range(0, triang.triangles.shape[0]):
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
            flag[1] = (x[i] >= -w/2 and x[i] <= w/2) and (y[i] >= -d/2 and y[i] <= d/2 )
            flag[2] = (x[j] >= -w/2 and x[j] <= w/2) and (y[j] >= -d/2 and y[j] <= d/2 )
            if flag[1] and flag[2]:
                Se[i][j] = er*(b[i]*b[j] + c[i]*c[j])*Ae
            else:
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


# Electric Field Calculation
interpolator = mtri.LinearTriInterpolator(triang, Potentials)
[Ex, Ey] = interpolator.gradient(triang.x, triang.y)


 #Plot the results
fig, ax = plt.subplots(1)
ax.triplot(triang)
cax = ax.tripcolor(triang, Potentials, cmap='plasma', shading='flat')
fig.colorbar(cax, ax=ax)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()

