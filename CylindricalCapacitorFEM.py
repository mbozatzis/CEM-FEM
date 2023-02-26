import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

a = 0.76
b = 1.75
V = 1


# Create Geometry
R = np.linspace(a, b, 50)
Theta = np.linspace(0, 2*np.pi, 100)

R_matrix, Theta_matrix = np.meshgrid(R, Theta)
X = R_matrix * np.cos(Theta_matrix)
Y = R_matrix * np.sin(Theta_matrix)
X = X.reshape(X.size)
Y = Y.reshape(Y.size)

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

plt.triplot(triang)
plt.show()

# Define known and unknown potentials
node_id = np.ones(nodes.shape[0])
Potentials = np.zeros(nodes.shape[0])
for inode in range(0, nodes.shape[0]-1):
    radius = np.sqrt(pow(nodes[inode][0],2) + pow(nodes[inode][1], 2))
    if radius-a <= pow(10,-6):
        node_id[inode] = 0
        Potentials[inode] = V
    elif radius-b <= pow(10,-6):
        node_id[inode] = 0
        Potentials[inode] = 0

# Index of unkowns
un_index = np.zeros(nodes.size)
counter = 0
for inode in range(0, nodes.shape[0]-1):
    if node_id[inode] == 1:
        un_index[inode] = counter
        counter = counter + 1
        

