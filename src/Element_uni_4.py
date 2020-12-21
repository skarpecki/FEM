import numpy as np
from math import sqrt

from src.Node import Node
from src.Global_Data import GlobalData



class Surface_Uni:
    """pc1, pc2 - tuple, [0] -> ksi, [1] -> eta"""
    def __init__(self, pc, w):
        self.w = w
        self.pc = pc
        self.N_of_pc = np.zeros((len(pc),4))
        for i in range(len(self.pc)):
            self.N_of_pc[i] = self._get_N_in_pc(self.pc[i])

    def _get_N_in_pc(self, pc):
        """returns tuple of 1x4 vector with values of N functions for each pc"""
        arr_N_func = np.zeros((1, 4))
        arr_N_func[0][0] = 1 / 4 * (1 - pc[0]) * (1 - pc[1])
        arr_N_func[0][1] = 1 / 4 * (1 + pc[0]) * (1 - pc[1])
        arr_N_func[0][2] = 1 / 4 * (1 + pc[0]) * (1 + pc[1])
        arr_N_func[0][3] = 1 / 4 * (1 - pc[0]) * (1 + pc[1])
        return arr_N_func

    def _get_det_J(self, node1: Node, node2: Node):
        L = sqrt(pow(node1.x - node2.x, 2) + pow(node1.y - node2.y, 2))
        return L / 2

    def get_H_BC_local(self, node1: Node, node2: Node, alfa):
        detJ = self._get_det_J(node1, node2)
        H_BC_local = np.zeros((4,4))
        for i in range(self.N_of_pc.shape[0]):
            H_BC_local += np.matmul(self.N_of_pc.transpose()[:, i : i+1], self.N_of_pc[i : i+1, :]) * self.w[i]
        H_BC_local = alfa * (H_BC_local * detJ)
        return H_BC_local

    def get_P_local(self, node1: Node, node2: Node, alfa, to):
        detJ = self._get_det_J(node1, node2)
        P = np.zeros((4,1))
        #self.N_of_pc.shape[0] is no. of rows what is equal to no. of npc
        for i in range(self.N_of_pc.shape[0]):
            P += self.N_of_pc.transpose()[:, i:i+1] * self.w[i]
        P = P * -alfa * to * detJ
        return P.transpose()




class Element_Uni_4:
    def __init__(self, npc):
        self.npc = npc
        self.etas = []
        self.ksis = []
        self.wages = [] #list of tuples with wage as (x,y)
        if npc == 2:
          val = 1/sqrt(3)
          surf_wages = [1, 1]
          self.surf_bottom_uni = Surface_Uni([(-val, -1), (val, -1)], surf_wages)
          self.surf_right_uni = Surface_Uni([(1, -val), (1, val)], surf_wages)
          self.surf_top_uni = Surface_Uni([(-val, 1), (val, 1)], surf_wages)
          self.surf_left_uni = Surface_Uni([(-1, -val), (-1, val)], surf_wages)

          for y in range(npc):
              for x in range(npc):
                  self.ksis.append(-val + 2 * x * val)
                  self.etas.append(-val + 2 * y * val)
                  self.wages.append((1,1))

        elif npc == 3:
            val = sqrt(3/5)
            surf_wages = [8/9, 5/9, 8/9]
            self.surf_bottom_uni = Surface_Uni([(-val, -1), (0, -1), (val, -1)], surf_wages)
            self.surf_right_uni = Surface_Uni([(1, -val), (1, 0), (1, val)], surf_wages)
            self.surf_top_uni = Surface_Uni([(-val, 1), (0, 1), (val, 1)], surf_wages)
            self.surf_left_uni = Surface_Uni([(-1, -val), (-1, 0),  (-1, val)], surf_wages)

            for y in range(npc):
                for x in range(npc):
                    self.ksis.append(-val + x * val) #-val + 0val = -val, -val + 1val = 0; -val + 2val = val
                    self.etas.append(-val + y * val)
                    x_wage = 8/9 if x == 1 else 5/9
                    y_wage = 8/9 if y == 1 else 5/9
                    self.wages.append((x_wage, y_wage))

        elif npc == 4:
            val1 = 0.861136
            val2 = 0.339981
            a_wage = 0.347855
            b_wage = 0.652145
            for y in range(npc):
                if y == 0:
                    y_val = -val1
                    y_wage = a_wage
                elif y == 1:
                    y_val = -val2
                    y_wage = b_wage
                elif y == 2:
                    y_val = val2
                    y_wage = b_wage
                elif y == 3:
                    y_val = val1
                    y_wage = a_wage
                for x in range(npc):
                    if x == 0:
                        x_val = -val1
                        x_wage = a_wage
                    elif x == 1:
                        x_val = -val2
                        x_wage = b_wage
                    elif x == 2:
                        x_val = val2
                        x_wage = b_wage
                    elif x == 3:
                        x_val = val1
                        x_wage = a_wage

                    self.ksis.append(x_val)
                    self.etas.append(y_val)
                    self.wages.append((x_wage, y_wage))

        else:
            raise TypeError("Incorrect number of npcs")

        self.N_array = self._local_func_deriv()[0]
        self.N_of_ksi = self._local_func_deriv()[1]
        self.N_of_eta = self._local_func_deriv()[2]




    def _local_func_deriv(self):
        etas = self.etas
        ksis = self.ksis
        dim = pow(self.npc, 2)

        arr_dN_dKsi = np.empty((dim, 4))
        arr_dN_dEta = np.empty((dim,4))
        arr_N_func = np.empty((dim, 4))

        for j in range(dim):
            arr_N_func[j][0] = 1/4 * (1-ksis[j]) * (1-etas[j])
            arr_N_func[j][1] = 1/4 * (1+ksis[j]) * (1-etas[j])
            arr_N_func[j][2] = 1/4 * (1+ksis[j]) * (1+etas[j])
            arr_N_func[j][3] = 1/4 * (1-ksis[j]) * (1+etas[j])

            arr_dN_dKsi[j] = np.array( [ -1/4 * (1-etas[j]),  1/4 * (1-etas[j]), 1/4 * (1+etas[j]), -1/4 * (1+etas[j]) ] )
            arr_dN_dEta[j] = np.array( [ -1/4 * (1-ksis[j]), -1/4 * (1+ksis[j]), 1/4 * (1+ksis[j]),  1/4 * (1-ksis[j]) ] )

        return (arr_N_func, arr_dN_dKsi, arr_dN_dEta)


    def _get_jacobi(self, id, Nodes):
        id = id - 1
        jacobi = np.zeros((2,2))
        for i in range(4):
            jacobi[0][0] += self.N_of_ksi[id][i] * Nodes[i].x
            jacobi[0][1] += self.N_of_ksi[id][i] * Nodes[i].y

            jacobi[1][0] += self.N_of_eta[id][i] * Nodes[i].x
            jacobi[1][1] += self.N_of_eta[id][i] * Nodes[i].y

        return jacobi


    def get_H_C_matrix_for_point(self, Nodes, id, k, c, ro):
        jacobi = self._get_jacobi(id, Nodes)
        jacobi_inverse = np.linalg.inv(jacobi) #inv has 1/detJ in its def, hence no need to additionaly do this later
        jacobi_det = np.linalg.det(jacobi)

        dN_dX = np.zeros((1,4))
        dN_dY = np.zeros((1,4))
        N_array = np.zeros((1,4))

        for i in range(4):
            mat_temp = np.array([[self.N_of_ksi[id-1][i]],
                                 [self.N_of_eta[id-1][i]]])
            mat_temp = np.matmul(jacobi_inverse, mat_temp)
            dN_dX[0][i] = mat_temp[0]
            dN_dY[0][i] = mat_temp[1]

        #H matrix aggregation
        dN_dX_t = dN_dX.transpose()
        dN_dY_t = dN_dY.transpose()
        H = k * jacobi_det * (np.matmul(dN_dX_t, dN_dX ) + np.matmul(dN_dY_t, dN_dY)) *\
            self.wages[id-1][0] * self.wages[id-1][1]

        #C matrix aggregation
        N_array[0] = self.N_array[id - 1]
        N_array_t = N_array.transpose()
        C = c * jacobi_det * ro * (np.matmul(N_array_t, N_array)) * self.wages[id-1][0]  *  self.wages[id-1][1]


        return (H,C)

    def get_H_C_matrix(self, Nodes, globalData: GlobalData):
        #TODO: checking if Nodes are provided in counter clock manner (if possible)
        k = globalData.k
        c = globalData.c
        ro = globalData.ro
        alfa = globalData.alfa
        to = globalData.to

        sumH = self.get_H_C_matrix_for_point(Nodes, 1, k, c, ro)[0]
        sumC = self.get_H_C_matrix_for_point(Nodes, 1, k, c, ro)[1]

        for i in range(2, pow(self.npc,2)+1):
            sumH += self.get_H_C_matrix_for_point(Nodes, i, k, c, ro)[0]
            sumC += self.get_H_C_matrix_for_point(Nodes, i, k, c, ro)[1]

        #Hbc matrix, P vertice
        Hbc = np.zeros((4,4))
        P = np.zeros((1,4))

        if Nodes[0].flag_bc != 0 and Nodes[1].flag_bc != 0:
            Hbc += self.surf_bottom_uni.get_H_BC_local(Nodes[0], Nodes[1], alfa)
            P += self.surf_bottom_uni.get_P_local(Nodes[0], Nodes[1], alfa, to)

        if Nodes[1].flag_bc != 0 and Nodes[2].flag_bc != 0:
            Hbc += self.surf_right_uni.get_H_BC_local(Nodes[1], Nodes[2], alfa)
            P += self.surf_right_uni.get_P_local(Nodes[1], Nodes[2], alfa, to)


        if Nodes[2].flag_bc != 0 and Nodes[3].flag_bc != 0:
            Hbc += self.surf_top_uni.get_H_BC_local(Nodes[2], Nodes[3], alfa)
            P += self.surf_top_uni.get_P_local(Nodes[2], Nodes[3], alfa, to)

        if Nodes[3].flag_bc != 0 and Nodes[0].flag_bc != 0:
            Hbc += self.surf_left_uni.get_H_BC_local(Nodes[3], Nodes[0], alfa)
            P += self.surf_left_uni.get_P_local(Nodes[3], Nodes[0], alfa, to)

        sumH += Hbc

        return (sumH, sumC, P)
