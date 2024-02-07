import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def periodic_boundary(index, num_grid_points):
        if index < num_grid_points/2:
            return 2*index
        else:
            return 2*(num_grid_points-index)-1
        
def index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points):
        return periodic_boundary(index_i, num_i_grid_points) * num_j_grid_points + periodic_boundary(index_j, num_j_grid_points)


def index_periodic_boundary_backward_converter(Periodic_L_index, Periodic_L_indices,num_j_grid_points):
        index_i, index_j = Periodic_L_indices[Periodic_L_index]
        L_index = index_i * num_j_grid_points + index_j
        return L_index

def MatrixCreator(num_i_grid_points, num_j_grid_points, dt, \
                     di, dj, Q_11, Q_22, Q_12, Resvoir, ncols):
        index_plus_i = 0
        index_minus_i = 0
        index_plus_j = 0
        index_minus_j = 0
        index_i = 0
        
        coefficient_matrix = np.zeros((ncols, ncols))
        while index_i < num_i_grid_points:

            if index_i ==0:
                index_minus_i = num_i_grid_points - 1
            else:
                index_minus_i = index_i -1

            if index_i == num_i_grid_points:
                index_plus_i = 0
            else:
                index_plus_i = index_i + 1
        
            index_j = 0
            while index_j < num_j_grid_points:
            
                if index_j ==0:
                    index_minus_j = num_j_grid_points - 1
                else:
                    index_minus_j = index_j -1

                if index_j == num_j_grid_points:
                    index_plus_j = 0
                else:
                    index_plus_j = index_j + 1

                index_current = index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)

                coefficient_matrix[index_current, index_periodic_boundary_forward_converter(index_plus_i, index_j, num_i_grid_points, num_j_grid_points)] = (-Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*di**2)
                coefficient_matrix[index_current, index_periodic_boundary_forward_converter(index_i, index_plus_j, num_i_grid_points, num_j_grid_points)] = (-Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*dj**2)
                coefficient_matrix[index_current, index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] = (Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] *(dt)/(2*di**2) + Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt)/(2*dj**2) + Resvoir[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)]*(dt/2) + 1)
                coefficient_matrix[index_current, index_periodic_boundary_forward_converter(index_minus_i, index_j, num_i_grid_points, num_j_grid_points)] = (-Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*di**2)
                coefficient_matrix[index_current, index_periodic_boundary_forward_converter(index_i, index_minus_j, num_i_grid_points, num_j_grid_points)] = (-Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*dj**2)
                index_j+=1
            index_i+=1
        
        return coefficient_matrix

def phi_generator(num_i_grid_points, num_j_grid_points, dt, di, dj, Q_11, Q_12, Q_22, Resvoir, phi, Temperatures_Solution):
    
        index_plus_i = 0
        index_minus_i = 0
        index_plus_j = 0
        index_minus_j = 0

        index_i = 0
        while index_i <= num_i_grid_points-1:

            if index_i ==0:
                index_minus_i = num_i_grid_points - 1
        
            else:
                index_minus_i = index_i -1

            if index_i == num_i_grid_points:
                index_plus_i = 0
            else:
                index_plus_i = index_i + 1
        
            index_j = 0
            while index_j <= num_j_grid_points -1:
            
                if index_j ==0:
                    index_minus_j = num_j_grid_points - 1
        
                else:
                    index_minus_j = index_j -1

                if index_j == num_j_grid_points:
                    index_plus_j = 0
                else:
                    index_plus_j = index_j + 1

                contribution_0 = ((Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*di**2)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_plus_i, index_j, num_i_grid_points, num_j_grid_points)]
                contribution_1 = ((Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*dj**2)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_plus_j, num_i_grid_points, num_j_grid_points)]
                contribution_2 = -1.0 * (Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] *(dt)/(2*di**2) + Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt)/(2*dj**2) + Resvoir[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)]*(dt/2) - 1) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)]
                contribution_3 = ((Q_11[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*di**2)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_minus_i, index_j, num_i_grid_points, num_j_grid_points)]
                contribution_4 = ((Q_22[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] * (dt))/(4*dj**2)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_minus_j, num_i_grid_points, num_j_grid_points)]

                phi[index_periodic_boundary_forward_converter(index_i, index_j, num_i_grid_points, num_j_grid_points)] = contribution_0 + contribution_1 + contribution_2 + contribution_3 + contribution_4
                index_j +=1
            index_i +=1
        return phi

def Sample_Dataset(Q_11, Q_22,Q_12,num_i_grid_points, num_j_grid_points, Final_Time, Num_Timesteps, Resvoir, Source):

    num_grid_points = num_i_grid_points * num_j_grid_points

    dTheta = (2*np.pi)/(num_i_grid_points)
    dZeta = (2*np.pi)/(num_j_grid_points)
    dTime = (Final_Time)/(Num_Timesteps)

    Temperatures_Solution = np.zeros((num_grid_points,))
    Temperature_Grid = []
    phi = np.zeros((num_grid_points,))
    b = np.zeros((num_grid_points,))

    Periodic_L_indices = np.zeros((num_grid_points, 2), dtype = np.int64)
    index_a = 0
    while index_a <= num_i_grid_points -1:
        index_b = 0
        while index_b <= num_j_grid_points -1:
            Periodic_L_indices[index_periodic_boundary_forward_converter(index_a,index_b, num_i_grid_points, num_j_grid_points)] = (int(index_a),int(index_b))
            index_b +=1
        index_a +=1
    Periodic_L_indices.shape

    time_index = 0
    coefficient_matrix = MatrixCreator(num_i_grid_points, num_j_grid_points, dTime, dTheta, dZeta, Q_11, Q_22, Q_12, Resvoir, num_grid_points)
    coefficient_matrix_sparse = csr_matrix(coefficient_matrix)
    while time_index <= Num_Timesteps-1:

        Temperatures_Solution = np.zeros((num_grid_points,))
        Temperature_Grid_dt = np.zeros((num_grid_points,))
        
        b = phi + 2.0*Source*dTime

        Temperatures_Solution = spsolve(coefficient_matrix_sparse, b)
        
        phi = phi_generator(num_i_grid_points, num_j_grid_points, dTime, dTheta, dZeta, Q_11, Q_12, Q_22, Resvoir, phi, Temperatures_Solution)
        index_i = 0
        while index_i <= num_i_grid_points -1:
            index_j = 0
            while index_j <= num_j_grid_points-1:
                Periodic_L_index = index_i*num_j_grid_points + index_j
                L_index = index_periodic_boundary_backward_converter(Periodic_L_index, Periodic_L_indices, num_j_grid_points)
                Temperature_Grid_dt[L_index] = Temperatures_Solution[Periodic_L_index]
                index_j +=1
            index_i +=1

        Temperature_Grid.append(Temperature_Grid_dt)
        time_index+=1
        
    Temperature_Grid = list(np.array(Temperature_Grid).reshape(-1, num_i_grid_points, num_j_grid_points))
    return Temperature_Grid

    
