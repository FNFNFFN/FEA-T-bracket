# Libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.sparse
import scipy.linalg
# import sympy as sym

# Shape and Mesh Parameters
L1 = 80*10**-3 
W1 = 40*10**-3
L2 = 20*10**-3
W2 = 100*10**-3
t = 2.5*10**-3
mesh_size = 1*10**-3

# Material Properties
E = 70*10**9
nu = 0.3
P = 5000/2

# Number of Elements on Each Sides
x1 = int(L1/mesh_size)
y1 = int(W1/mesh_size)
x2 = int(L2/mesh_size)
y2 = int(W2/mesh_size)

# Number of Elements, Nodes, DOFs
num_ele_1 = x1*y1
num_ele_2 = x2*y2
num_ele = int(num_ele_1 + num_ele_2)

num_node_1 = (x1 + 1)*(y1 + 1)
num_node_2 = (x2 + 1)*(y2 + 1)
num_node = int(num_node_1 + num_node_2 - (y1 + 1))

num_DOF_1 = 2*num_node_1
num_DOF_2 = 2*num_node_2
num_DOF = int(2*num_node)

# Stiffness, Force, Displacement Matrix Initialization
K = scipy.sparse.lil_matrix((num_DOF, num_DOF))
F = scipy.sparse.lil_matrix((num_DOF, 1))
# K = np.zeros((num_DOF, num_DOF))
# F = np.zeros((num_DOF, 1))
d = np.zeros((num_DOF, 1))

# Strain and Stress Matrix Initialization
epsilon_x = np.zeros((num_ele, 1))
epsilon_y = np.zeros((num_ele, 1))
epsilon_xy = np.zeros((num_ele, 1))

sigma_x = np.zeros((num_ele, 1))
sigma_y = np.zeros((num_ele, 1))
sigma_xy = np.zeros((num_ele, 1))
sigma_vm = np.zeros((num_ele, 1))

# Node Coordinate Matrix Initialization
xcoor = np.zeros((num_node,1))
ycoor = np.zeros((num_node,1))

# Stiffness Matrix for a Single Element
def Stiffness_Matrix_Calculation_Qudra(E,nu,mesh_size,t):
    D = E/(1-nu**2)*np.array([[1,nu,0],
                                [nu,1,0],
                                [0,0,(1-nu)/2]])
    a = mesh_size/2
    xc = np.array([-a, -a, a, a])
    yc = np.array([a, -a, -a, a])
    xyc = np.array([[xc[0], yc[0]],
                    [xc[1], yc[1]],
                    [xc[2], yc[2]],
                    [xc[3], yc[3]]])
  
    xi = np.array([-1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)])
    eta = np.array([1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3)])
    w_xi = np.array(4*[1])
    w_eta = np.array(4*[1])

    K = 0
    H_M = np.array(4*[[]]).tolist()
    for i in range(0,xi.size):
        dN1_dxi = (1/4)*(-1 - eta[i])
        dN1_deta = (1/4)*(1 - xi[i])
        dN2_dxi = (1/4)*(eta[i] - 1)
        dN2_deta = (1/4)*(xi[i] - 1)
        dN3_dxi = (1/4)*(1 - eta[i])
        dN3_deta = (1/4)*(-1 - xi[i])
        dN4_dxi = (1/4)*(1 + eta[i])
        dN4_deta = (1/4)*(1 + xi[i])

        GN_4Q = np.array([[dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi],
                            [dN1_deta, dN2_deta, dN3_deta, dN4_deta]])
        J = np.dot(GN_4Q, xyc)
        detJ = np.linalg.det(J)
        H_star = np.dot(np.linalg.inv(J), GN_4Q)
        H = np.array([[H_star[0][0], 0, H_star[0][1], 0, H_star[0][2], 0, H_star[0][3], 0],
                        [0, H_star[1][0], 0, H_star[1][1], 0, H_star[1][2], 0, H_star[1][3]],
                        [H_star[1][0], H_star[0][0], H_star[1][1], H_star[0][1], H_star[1][2], H_star[0][2], H_star[1][3], H_star[0][3]]])
        H_M[i] = np.copy(H)
        K = K + w_xi[i]*w_eta[i]*H.T@D@H*detJ
    K = t*K

    # H Matrix of Center Location (0,0)
    dN1_dxi = (1/4)*(-1)
    dN1_deta = (1/4)*(1)
    dN2_dxi = (1/4)*(-1)
    dN2_deta = (1/4)*(-1)
    dN3_dxi = (1/4)*(1)
    dN3_deta = (1/4)*(-1)
    dN4_dxi = (1/4)*(1)
    dN4_deta = (1/4)*(1)

    GN_4Q = np.array([[dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi],
                        [dN1_deta, dN2_deta, dN3_deta, dN4_deta]])
    J = np.dot(GN_4Q, xyc)
    detJ = np.linalg.det(J)
    H_star = np.dot(np.linalg.inv(J), GN_4Q)
    Hc = np.array([[H_star[0][0], 0, H_star[0][1], 0, H_star[0][2], 0, H_star[0][3], 0],
                    [0, H_star[1][0], 0, H_star[1][1], 0, H_star[1][2], 0, H_star[1][3]],
                    [H_star[1][0], H_star[0][0], H_star[1][1], H_star[0][1], H_star[1][2], H_star[0][2], H_star[1][3], H_star[0][3]]])

    return K, D, H_M, Hc

K_e, D_e, H_M_e, Hc_e = Stiffness_Matrix_Calculation_Qudra(E,nu,mesh_size,t)

# Coordinate of all nodes, Assembling Stiffness Matrix
node_M = []
DOF_M = []
for i in range(0,x1):
    for j in range(0,y1):
        node = [(y1+1)*i + (j+1), (y1+1)*i + (j+2), (y1+1)*(i+1)+(j+2), (y1+1)*(i+1) + (j+1)]
        node_M.append(node)

        xcoor[(np.array(node) - 1).tolist()] = np.array([i*mesh_size, i*mesh_size, (i+1)*mesh_size, (i+1)*mesh_size]).reshape(4,1)
        ycoor[(np.array(node) - 1).tolist()] = (0.1 - np.array([j*mesh_size, (j+1)*mesh_size, (j+1)*mesh_size, j*mesh_size])).reshape(4,1)

        DOF = [2*node[0] - 1, 2*node[0], 2*node[1] - 1, 2*node[1], 2*node[2] - 1, 2*node[2], 2*node[3] - 1, 2*node[3]]
        DOF_M.append(DOF)
        DOF_row = np.array(8*((np.array(DOF) - 1).tolist())).reshape(8,8).T
        DOF_col = np.array((8*(np.array(DOF) - 1).tolist())).reshape(8,8)
        K[DOF_row,DOF_col] = K[DOF_row,DOF_col] + K_e

for i in range(0,x2):
	for j in range(0,y2):
		node = ((num_node_1 - y1 - 1) + np.array([(y2+1)*i + (j+1), (y2+1)*i + (j+2), (y2+1)*(i+1)+(j+2), (y2+1)*(i+1) + (j+1)])).tolist()
		node_M.append(node)

		xcoor[(np.array(node) - 1).tolist()] = 0.08 + np.array([i*mesh_size, i*mesh_size, (i+1)*mesh_size, (i+1)*mesh_size]).reshape(4,1)
		ycoor[(np.array(node) - 1).tolist()] = (0.1 - np.array([j*mesh_size, (j+1)*mesh_size, (j+1)*mesh_size, j*mesh_size])).reshape(4,1)

		DOF = [2*node[0] - 1, 2*node[0], 2*node[1] - 1, 2*node[1], 2*node[2] - 1, 2*node[2], 2*node[3] - 1, 2*node[3]]
		DOF_M.append(DOF)
		DOF_row = np.array(8*((np.array(DOF) - 1).tolist())).reshape(8,8).T
		DOF_col = np.array((8*(np.array(DOF) - 1).tolist())).reshape(8,8)
		K[DOF_row,DOF_col] = K[DOF_row,DOF_col] + K_e

# Implement Initial Conditions and Restraint
node_F = -P/((x2+1) + (x2+1-2))
F[2*(num_node_1 - y1 + y2) - 1:num_DOF - 1:2*(y2+1)] = 2*node_F
F[2*(num_node_1 - y1 + y2) - 1] = node_F
F[num_DOF - 1] = node_F

dof_rstn_1 = set(np.arange(1,2*(y1+1) + 1,1))
dof_rstn_2 = set(np.arange(num_DOF - 2*(y2+1) + 1, num_DOF + 1, 2))
dof_all = set(np.arange(1,num_DOF + 1,1))
dof_free = np.array(list(dof_all - dof_rstn_1 - dof_rstn_2))

dof_free_row = np.array(len(dof_free)*((np.array(dof_free) - 1).tolist())).reshape(len(dof_free),len(dof_free)).T
dof_free_col = np.array(len(dof_free)*((np.array(dof_free) - 1).tolist())).reshape(len(dof_free),len(dof_free))

# Calculate Displacement of Each Node
# K = sp.sparse.csc_matrix(K)
# d[dof_free - 1] = inv(K[dof_free_row,dof_free_col])@F[dof_free - 1]
d[dof_free - 1] = scipy.linalg.inv(K[dof_free_row,dof_free_col].toarray())@F[dof_free - 1]
d_graph = 15*d

# X and y Coordinates After Deformation
xcoor_de = xcoor + d_graph[0:len(d):2]
ycoor_de = ycoor + d_graph[1:len(d):2]

# Calculate Strain and Stress of Each Element
for i in range(0,len(node_M)):
    node = np.array(node_M[i])

    epsilon_e = np.array(Hc_e)@d[np.array(DOF_M[i]) - 1]
    epsilon_x[i] = epsilon_e[0]
    epsilon_y[i] = epsilon_e[1]
    epsilon_xy[i] = epsilon_e[2]

    sigma_e = D_e@epsilon_e
    sigma_x[i] = sigma_e[0]
    sigma_y[i] = sigma_e[1]
    sigma_xy[i] = sigma_e[2]
    sigma_vm[i] = np.sqrt(0.5*((sigma_e[0] - sigma_e[1])**2 + (sigma_e[1])**2 + (-sigma_e[0])**2) + 3*(sigma_e[2])**2)



# Plot Original Structure
node_top_hor_1 = np.arange(0,num_node_1,y1+1)
node_top_hor_2 = np.arange((num_node_1 - y1 - 1),num_node,y2+1 )
plt.figure(figsize = (12, 12))
plt.plot(xcoor[node_top_hor_1],ycoor[node_top_hor_1],'b')
plt.plot(xcoor[node_top_hor_2],ycoor[node_top_hor_2],'b')
plt.plot(xcoor_de[node_top_hor_1],ycoor_de[node_top_hor_1],'r')
plt.plot(xcoor_de[node_top_hor_2],ycoor_de[node_top_hor_2],'r')

for i in range(0,len(node_M)):
	plt.plot(xcoor[np.array(node_M[i])-1], ycoor[np.array(node_M[i])-1],'b')
	plt.plot(xcoor_de[np.array(node_M[i])-1], ycoor_de[np.array(node_M[i])-1],'r')

plt.axis('equal')
plt.title('Actual and Displaced Shape of the T-Bracket (Exaggerated)')
plt.show()

# Plot Interpolate Colors of Displacement (Patch)
# plt.figure(figsize = (6, 6))

# plt.plot(xcoor[node_top_hor_1],ycoor[node_top_hor_1],'b')
# plt.plot(xcoor[node_top_hor_2],ycoor[node_top_hor_2],'b')
# for i in range(0,len(node_M)):
  # plt.plot(xcoor[np.array(node_M[i])-1], ycoor[np.array(node_M[i])-1],'b')

# plt.pcolormesh(xcoor[np.array(node_M[i])-1], ycoor[np.array(node_M[i])-1], d[0:len(d):2][np.array(node_M[i])-1], edgecolors = None, shading = 'nearest', cmap = 'jet')


# Plot Strain Distribution of Structure
plt.figure(figsize = (72, 72))

plt.subplot(2, 3, 1)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), epsilon_x[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_x), vmax = max(epsilon_x))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), epsilon_x[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_x), vmax = max(epsilon_x))
plt.axis('equal')
plt.colorbar()
plt.title('Strain Distribution of x Direction')

plt.subplot(2, 3, 2)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), epsilon_y[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_y), vmax = max(epsilon_y))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), epsilon_y[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_y), vmax = max(epsilon_y))
plt.axis('equal')
plt.colorbar()
plt.title('Strain Distribution of y Direction')

plt.subplot(2, 3, 3)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), epsilon_xy[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_xy), vmax = max(epsilon_xy))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), epsilon_xy[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(epsilon_xy), vmax = max(epsilon_xy))
plt.axis('equal')
plt.colorbar()
plt.title('Strain Distribution of Shear xy')

# Plot Stress Distribution of Structure
plt.subplot(2, 3, 4)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), sigma_x[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_x), vmax = max(sigma_x))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), sigma_x[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_x), vmax = max(sigma_x))
plt.axis('equal')
plt.colorbar()
plt.title('Stress Distribution of x Direction')

plt.subplot(2, 3, 5)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), sigma_y[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_y), vmax = max(sigma_y))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), sigma_y[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_y), vmax = max(sigma_y))
plt.axis('equal')
plt.colorbar()
plt.title('Stress Distribution of y Direction')

plt.subplot(2, 3, 6)
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), sigma_xy[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_xy), vmax = max(sigma_xy))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), sigma_xy[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_xy), vmax = max(sigma_xy))
plt.axis('equal')
plt.colorbar()
plt.title('Stress Distribution of Shear xy')

plt.figure(figsize = (12, 12))
plt.pcolormesh(np.arange(0, L1+mesh_size, mesh_size), np.arange(W2, W2-W1-mesh_size, -mesh_size), sigma_vm[0:num_ele_1].reshape(x1,y1).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_vm), vmax = max(sigma_vm))
plt.pcolormesh(np.arange(L1, L1+L2, mesh_size), np.arange(W2, 0-mesh_size, -mesh_size), sigma_vm[num_ele_1:num_ele].reshape(x2,y2).T, edgecolors = 'none', shading = 'flat', cmap = 'jet', vmin = min(sigma_vm), vmax = max(sigma_vm))
plt.axis('equal')
plt.colorbar()
plt.title('Von Mises Stress Distribution')
plt.show()