#TODO: sprawdzic dlaczego Siatka generuje niepoprawne wyniki dla nH != nW

from math import sqrt
import numpy as np

#Node def
class Node:
    def __init__(self, x, y, t0):
        self.x = x
        self.y = y
        self.t0 = t0

    def __str__(self):
        return f"{self.x}, {self.y}"

    def __repr__(self):
        return "{0!s}(x={1!r}, y={2!r})".format(self.__class__.__name__, self.x, self.y)

#Elements def
class Element:
    def __init__(self, *args):
        self.ids = args

    def __str__(self):
        return f"{self.ids}"

    def __repr__(self):
        return"{0!s}(*args={1!r}".format(self.__class__.__name__, self.ids)






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
            arr_N_func[j][2] = 1/4 * (1+ksis[j]) * (1-etas[j])
            arr_N_func[j][3] = 1/4 * (1-ksis[j]) * (1-etas[j])

            arr_dN_dKsi[j] = np.array( [ -1/4 * (1-etas[j]),  1/4 * (1-etas[j]), 1/4 * (1+etas[j]), -1/4 * (1+etas[j]) ] )
            arr_dN_dEta[j] = np.array( [ -1/4 * (1-ksis[j]), -1/4 * (1+ksis[j]), 1/4 * (1+ksis[j]),  1/4 * (1-ksis[j]) ] )

        return (arr_N_func, arr_dN_dKsi, arr_dN_dEta)


    def _get_jacobi(self, id, Nodes):
        id = id - 1
        row_1 = []
        row_2 = []
        sumx = 0
        sumy = 0
        #i for Ni; id for row in matrix
        for i in range(4):
            sumx += self.N_of_ksi[id][i] * Nodes[i].x
            sumy += self.N_of_ksi[id][i] * Nodes[i].y
        row_1.extend([sumx, sumy])

        sumx = 0
        sumy = 0
        for i in range(4):
            sumx += self.N_of_eta[id][i] * Nodes[i].x
            sumy += self.N_of_eta[id][i] * Nodes[i].y
        row_2.extend([sumx, sumy])
        return np.array([row_1, row_2])


    #calculating {dN/dx}*{dN/dx}T and {dN/dy}*{dN/dy}T and return both as tuples
    def _get_x_y_vertice(self, id, Nodes):
        x_vertice = []
        y_verice = []
        for i in range(4):
            j = self._get_jacobi(id, Nodes)
            jI = np.linalg.inv(j)
            m_x_y = np.array( [ [self.N_of_ksi[id-1][i]] ,
                                [self.N_of_eta[id-1][i]] ] )
            m = np.matmul(jI, m_x_y)
            x_vertice.append(m[0].item(0))
            y_verice.append(m[1].item(0))
        return((x_vertice,y_verice))    


    #TODO: implement built-in transfrom matrixs
    def get_H_matrix_for_point(self, Nodes, id, k, c, ro):

        x_y = self._get_x_y_vertice(id, Nodes)
        x_mat =  [ x_y[0] ]  
        y_mat =  [ x_y[1] ]

        xc_mat = [ self.N_of_ksi[id-1].tolist() ]
        yc_mat = [ self.N_of_eta[id-1].tolist() ]

        x_mat_trans = []
        y_mat_trans = []

        xc_mat_trans = []
        yc_mat_trans = []

        #iterating over to fill trans matrix
        for x_arr in x_mat:
            for x in x_arr:
                x_mat_trans.append([x])

        for y_arr in y_mat:
            for y in y_arr:
                y_mat_trans.append([y])

        for xc_arr in xc_mat:
            for xc in xc_arr:
                xc_mat_trans.append([xc])
        for yc_arr in yc_mat:
            for yc in yc_arr:
                yc_mat_trans.append([yc])

        det = np.linalg.det(self._get_jacobi(id, Nodes))

        x_result = np.array(np.matmul(x_mat_trans, x_mat))
        y_result = np.array(np.matmul(y_mat_trans, y_mat))

        xc_result = np.array(np.matmul(xc_mat_trans, xc_mat))
        yc_result = np.array(np.matmul(yc_mat_trans, yc_mat))
        #doesnt work - why?
        #H = k*(x_result * self.wages[id-1][0] + y_result * self.wages[id-1][1]) * det
        H = k*(x_result + y_result) * det * self.wages[id-1][0] * self.wages[id-1][1]
        C = c*ro*(xc_result + yc_result) * det * self.wages[id-1][0] * self.wages[id-1][1]
        return (H,C)

    def get_H_matrix(self, Nodes, k, c, ro):
        #adding 1 as normallny range = for(i=0;i<x), here I need for(i=1;i<=x)

        sumH = self.get_H_matrix_for_point(Nodes, 1, k, c, ro)[0]
        sumC = self.get_H_matrix_for_point(Nodes, 1, k, c, ro)[1]

        for i in range(2, pow(self.npc,2)+1):
            sumH += self.get_H_matrix_for_point(Nodes, i, k, c, ro)[0]
            sumC += self.get_H_matrix_for_point(Nodes, i, k, c, ro)[1]

        if False:
            H1 = self.get_H_matrix_for_point(Nodes, 1, k)
            H2 = self.get_H_matrix_for_point(Nodes, 2, k)
            H3 = self.get_H_matrix_for_point(Nodes, 3, k)
            H4 = self.get_H_matrix_for_point(Nodes, 4, k)

        with open("sum_H.txt", "a") as a_file:
            np.savetxt(a_file, sumH, fmt='%.4f')
            a_file.write("\n")

        with open("sum_C.txt", "a") as a_file:
            np.savetxt(a_file, sumC, fmt='%.4f')
            a_file.write("\n")

        return (sumH, sumC)
    
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
            self.t0 = int(data['t0'])

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

    def list_nodes(self):
        for key, value in self.Nodes.items():
            print(f"{key}: {value}")

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

    def fill_H_global(self):
        el = Element_Uni_4(self.GlobalData.npc)
        H_Global = np.zeros(( len(self.Nodes), len(self.Nodes)))
        C_Global = np.zeros((len(self.Nodes), len(self.Nodes)))
        for i in range(1, self.GlobalData.nE + 1):
            Nodes = []
            for value in self.Elements[i].ids:
                Nodes.append(self.Nodes[value])
            H = el.get_H_matrix(Nodes, self.GlobalData.k, self.GlobalData.c, self.GlobalData.ro)[0]
            C = el.get_H_matrix(Nodes, self.GlobalData.k, self.GlobalData.c, self.GlobalData.ro)[1]
            for n in range(0, len(self.Elements[i].ids)):
                for m in range(0, len(self.Elements[i].ids)):
                    H_Global.itemset((self.Elements[i].ids[n]-1,
                                    self.Elements[i].ids[m]-1),
                                    H_Global.item(self.Elements[i].ids[n]-1,  self.Elements[i].ids[m]-1) +
                                    H.item(n,m))
                    C_Global.itemset((self.Elements[i].ids[n] - 1,
                                      self.Elements[i].ids[m] - 1),
                                     C_Global.item(self.Elements[i].ids[n] - 1, self.Elements[i].ids[m] - 1) +
                                     C.item(n, m))
        return (H_Global, C_Global)

if False:
    Nodes = [Node(0,0), Node(4,0), Node(4,6), Node(0,6)]
    nodes = [Node(0,0), Node(4,0), Node(4,4), Node(0,5)]
    el = Element_Uni_4()
    m = el._get_jacobi(1, Nodes)
    v = el._get_x_y_vertice(1, Nodes)
    h = el.get_H_matrix(Nodes, k=30)

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
        from pprint import pprint as pp
        #pp(s1.fill_H_global())
        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\h_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
            np.savetxt(a_file, s1.fill_H_global()[0], fmt='%.4f')

        with open(rf"D:\DevProjects\PythonProjects\MES\data\results\c_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
            np.savetxt(a_file, s1.fill_H_global()[1], fmt='%.4f')
            #s1.fill_H_global()

