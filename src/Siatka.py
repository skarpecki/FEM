#TODO: sprawdzic dlaczego Siatka generuje niepoprawne wyniki dla nH != nW

from math import sqrt
import numpy as np


#GlobalData def
class GlobalData:
    def __init__(self, path):
        data = self.__get_data__(path)
        if(data == None):
            raise Exception("Incorrect file")
        else:
            self._data_path = path
            self.H = data['H']
            self.W = data['W']
            self.nH = int(data['nH'])
            self.nW = int(data['nW'])
            self.k = int(data['k'])
            self.c = int(data['c'])
            self.ro = int(data['ro'])
            self.to = int(data['to'])
            self.t0 = int(data['t0'])
            self.alfa = int(data['alfa'])

            self.npc = int(data['npc'])
            self.nE = (self.nH - 1)  * (self.nW - 1)
            self.nN = self.nH * self.nW

    def __get_data__(self, path):
        with open(path, "rt") as file:
            data = {}
            lines = file.readlines()
            for line in lines:
                key, value = line.split()
                data[key] = float(value)
            return data
        return None

    def __repr__(self):
        return "{1!s}(path={2!r})".format(self.__class__.__name__, self._data_path)





#Node def
class Node:
    def __init__(self, x, y, t0, flag_bc=0):
        self.x = x
        self.y = y
        self.t0 = t0
        # flag to indicate whether boundary condition is applied to node
        self.flag_bc = flag_bc

    def __str__(self):
        return f"{self.x}, {self.y}"

    def __repr__(self):
        return "{0!s}(x={1!r}, y={2!r})".format(self.__class__.__name__, self.x, self.y)





#Elements def
class Element:
    """IDs must be indexed in counter-clock manner"""
    def __init__(self, *args):
        self.nodes_ids = args

    def __str__(self):
        return f"{self.nodes_ids}"

    def __repr__(self):
        return"{0!s}(*args={1!r}".format(self.__class__.__name__, self.nodes_ids)



class SOE:
    def __init__(self, H, C):
        self.H = H
        self.C = C


class Surface_Uni:
    """pc1, pc2 - tuple, [0] -> ksi, [1] -> eta"""
    def __init__(self, pc1: tuple, pc2: tuple, w1, w2):
        self.w1 = w1
        self.w2 = w2
        self.N_of_pc1 = self._get_N_in_pc(pc1)
        self.N_of_pc2 = self._get_N_in_pc(pc2)

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
        H_BC_local = np.zeros((4, 4))
        H_BC_local = alfa * (((np.matmul(self.N_of_pc1.transpose(), self.N_of_pc1) * self.w1) +
                (np.matmul(self.N_of_pc2.transpose(), self.N_of_pc2) * self.w2)) *
                detJ)
        return H_BC_local

    def get_P_local(self, node1: Node, node2: Node, alfa, to):
        detJ = self._get_det_J(node1, node2)
        P = np.zeros((4,4))
        P = -alfa * to * (self.N_of_pc1.transpose() * self.w1 + self.N_of_pc2.transpose() * self.w2) * detJ
        return P


class Element_Uni_4:
    def __init__(self, npc):
        self.npc = npc
        self.etas = []
        self.ksis = []
        self.wages = [] #list of tuples with wage as (x,y)
        if npc == 2:
          val = 1/sqrt(3)
          for y in range(npc):
              for x in range(npc):
                  self.ksis.append(-val + 2 * x * val)
                  self.etas.append(-val + 2 * y * val)
                  self.wages.append((1,1))

        elif npc == 3:
            val = sqrt(3/5)
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

        val = 1/sqrt(3)
        self.surf_bottom_uni = Surface_Uni((-val, -1), (val, -1), 1, 1)
        self.surf_right_uni = Surface_Uni((1, -val), (1, val), 1, 1)
        self.surf_top_uni = Surface_Uni((-val, 1), (val, 1), 1, 1)
        self.surf_left_uni = Surface_Uni((-1, -val), (-1, val), 1, 1)


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

        #Hbc matrix
        Hbc = np.zeros((4,4))
        if Nodes[0].flag_bc != 0 and Nodes[1].flag_bc != 0:
            Hbc += self.surf_bottom_uni.get_H_BC_local(Nodes[0], Nodes[1], alfa)

        if Nodes[1].flag_bc != 0 and Nodes[2].flag_bc != 0:
            Hbc += self.surf_right_uni.get_H_BC_local(Nodes[1], Nodes[2], alfa)

        if Nodes[2].flag_bc != 0 and Nodes[3].flag_bc != 0:
            Hbc += self.surf_top_uni.get_H_BC_local(Nodes[2], Nodes[3], alfa)

        if Nodes[3].flag_bc != 0 and Nodes[0].flag_bc != 0:
            Hbc += self.surf_left_uni.get_H_BC_local(Nodes[3], Nodes[0], alfa)


        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\local\sum_H.txt", "w") as a_file:
            np.savetxt(a_file, sumH, fmt='%.4f')
            a_file.write("\n")

        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\local\sum_C.txt", "w") as a_file:
            np.savetxt(a_file, sumC, fmt='%.4f')
            a_file.write("\n")

        return (sumH, sumC)


#Siatka
class Siatka:
    def __init__(self, path_to_data):
        self.GlobalData = GlobalData(path_to_data)
        self.Elements = self._fill_elements()
        self.Nodes = self._fill_nodes()

    def _fill_nodes(self):
        nodes = {}
        pos = 1
        dx = self.GlobalData.W / (self.GlobalData.nW - 1)
        dy = self.GlobalData.H / (self.GlobalData.nH - 1)
        for i in range(0,self.GlobalData.nW):
            for j in range(0,self.GlobalData.nH):
                nodes[pos] = Node(i*dx, j*dy, self.GlobalData.t0)
                pos += 1
        return nodes

    def set_bound_cond(self, indexes):
        for key, value in self.Nodes.items():
            if key in outside:
                value.flag_bc = 1



    def list_nodes(self):
        for key, value in self.Nodes.items():
            print(f"{key}: {value}, {value.flag_bc}")

    def _fill_elements(self):
        elements = {}
        j = 1
        for i in range(1, self.GlobalData.nE + 1):
            if j % self.GlobalData.nH == 0:
                j += 1
            elements[i] = Element(j, j + self.GlobalData.nH, j + self.GlobalData.nH + 1, j + 1)
            j += 1
        return elements

    def list_elements(self):
        for key, value in self.Elements.items():
            print(f"{key}: {value}")

    def fill_H_C_global(self):
        el = Element_Uni_4(self.GlobalData.npc)
        H_Global = np.zeros(( len(self.Nodes), len(self.Nodes)))
        C_Global = np.zeros((len(self.Nodes), len(self.Nodes)))
        for i in range(1, self.GlobalData.nE + 1):
            Nodes = []
            for value in self.Elements[i].nodes_ids:
                Nodes.append(self.Nodes[value])
            H = el.get_H_C_matrix(Nodes, self.GlobalData)[0]
            C = el.get_H_C_matrix(Nodes, self.GlobalData)[1]
            for n in range(0, len(self.Elements[i].nodes_ids)):
                for m in range(0, len(self.Elements[i].nodes_ids)):
                    H_Global.itemset((self.Elements[i].nodes_ids[n] - 1,
                                      self.Elements[i].nodes_ids[m] - 1),
                                     H_Global.item(self.Elements[i].nodes_ids[n] - 1, self.Elements[i].nodes_ids[m] - 1) +
                                     H.item(n,m))
                    C_Global.itemset((self.Elements[i].nodes_ids[n] - 1,
                                      self.Elements[i].nodes_ids[m] - 1),
                                     C_Global.item(self.Elements[i].nodes_ids[n] - 1, self.Elements[i].nodes_ids[m] - 1) +
                                     C.item(n, m))
        return (H_Global, C_Global)


####




if __name__ == "__main__":
    if False:
        el = Element_Uni_4(3)
        import numpy as np
        Nodes = [Node(0,0), Node(4,0), Node(4,6), Node(0,6)]
        print(el.N_of_ksi)
        print(el.N_of_eta)
        m = el._get_jacobi(1, Nodes)
        print(m)
        #print(el.get_H_matrix(Nodes, 25))
        #v = el._get_x_y_vertice(1, Nodes)
        #h = el.get_H_matrix(Nodes, k=25)
    if True:
        path = "D:\\DevProjects\\PythonProjects\\MES\\data\\global_data.txt"
        s1 = Siatka(path)
        outside = [1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16]
        s1.set_bound_cond(outside)
        s1.list_nodes()
        val = 1/sqrt(3)
        #powierzchnia = Surface_Uni((-1, -val), (-1, val), 1, 1 )
        #print(powierzchnia.get_H_BC_local(s1.Nodes[1], s1.Nodes[2], 25))
        #print(powierzchnia.get_P_local(s1.Nodes[1], s1.Nodes[2], 25, 1200))
        from pprint import pprint as pp
        #pp(s1.fill_H_global())
        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\h_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
            np.savetxt(a_file, s1.fill_H_C_global()[0], fmt='%.4f')

        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\c_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
            np.savetxt(a_file, s1.fill_H_C_global()[1], fmt='%.4f')
            #s1.fill_H_global()

