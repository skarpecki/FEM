#TODO: sprawdzic dlaczego Siatka generuje niepoprawne wyniki dla nH != nW

from math import sqrt
import numpy as np

from src.Element_uni_4 import Element_Uni_4
from src.Node import Node
from src.Global_Data import GlobalData

class Element:
    """IDs must be indexed in counter-clock manner"""
    def __init__(self, *args):
        self.nodes_ids = args

    def __str__(self):
        return f"{self.nodes_ids}"

    def __repr__(self):
        return"{0!s}(*args={1!r}".format(self.__class__.__name__, self.nodes_ids)

class SOE:
    def __init__(self, H, C, P):
        self.H = H
        self.C = C
        self.P = P
        self.H_calc = None

class Siatka:
    iter_to_print = 0
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

    def set_bound_cond(self):
        nH, nW = self.GlobalData.nH, self.GlobalData.nW
        outside = set()
        l_b = 1  # left bottom id
        l_t = 1 + nH - 1  # left top id
        r_b = 1 + nW * (nH - 1)  # left bottom id
        r_t = r_b + nH - 1  # left top id

        for i in range(l_b, l_t + 1):
            outside.add(i)

        for i in range(l_b, r_b + 1, nW):
            outside.add(i)

        for i in range(r_b, r_t + 1):
            outside.add(i)

        for i in range(l_t, r_t + 1, nW):
            outside.add(i)

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
        P_Global = np.zeros((1, len(self.Nodes)))

        for i in range(1, self.GlobalData.nE + 1):
            Nodes = []
            for value in self.Elements[i].nodes_ids:
                Nodes.append(self.Nodes[value])
            H = el.get_H_C_matrix(Nodes, self.GlobalData)[0]
            C = el.get_H_C_matrix(Nodes, self.GlobalData)[1]
            P = el.get_H_C_matrix(Nodes, self.GlobalData)[2]
            for n in range(0, len(self.Elements[i].nodes_ids)):
                P_Global.itemset((0, self.Elements[i].nodes_ids[n] - 1),
                                 P_Global.item(0, self.Elements[i].nodes_ids[n] - 1) +
                                 P.item(0, n))
                for m in range(0, len(self.Elements[i].nodes_ids)):
                    H_Global.itemset( (self.Elements[i].nodes_ids[n] - 1,
                                      self.Elements[i].nodes_ids[m] - 1),
                                     H_Global.item(self.Elements[i].nodes_ids[n] - 1, self.Elements[i].nodes_ids[m] - 1) +
                                     H.item(n,m))
                    C_Global.itemset( (self.Elements[i].nodes_ids[n] - 1,
                                      self.Elements[i].nodes_ids[m] - 1),
                                     C_Global.item(self.Elements[i].nodes_ids[n] - 1, self.Elements[i].nodes_ids[m] - 1) +
                                     C.item(n, m))

        self.soe = SOE(H_Global, C_Global, P_Global)

    def _calculate_H_P(self, t_ver):
        dt = self.GlobalData.dt
        P_Global = self.soe.P
        C_Global = self.soe.C

        if self.soe.H_calc is None:
            self.soe.H_calc = self.soe.H + (C_Global / dt)
        H_Global = self.soe.H_calc

        #t0 pierwsze w matmul poniewaz wektor wierszowy nie kolumnowy
        P_Global = -1 * (-1 * np.matmul(t_ver, (C_Global / dt)) + P_Global)
        return(H_Global, P_Global)

    def calculate_t(self, t_ver):
        H, P = self._calculate_H_P(t_ver)
        t1 = np.linalg.solve(H, P.transpose())
        return t1


    def simulate(self):
        t0_ver = np.zeros((1, len(self.Nodes)))
        for i in range(1, len(self.Nodes) + 1):
            t0_ver.itemset((0, i - 1), self.Nodes[i].t0)
        n = int(self.GlobalData.tk / self.GlobalData.dt)

        t_file = open(rf"D:\DevProjects\PythonProjects\MES\data\results\simulation\t_ver.txt", mode="w")
        open(rf"D:\DevProjects\PythonProjects\MES\data\results\simulation\min_max.txt", mode="w").close()
        with(open(rf"D:\DevProjects\PythonProjects\MES\data\results\simulation\min_max.txt", mode="a")) as file:
            file.write("t, min, max\n")
            for i in range (1, n+1):
                t1_ver = self.calculate_t(t0_ver)
                t0_ver = t1_ver.transpose()
                t_file.write("{}\n\n".format(t0_ver))
                file.write("{}:  {:.4f}, {:.4f}\n".format(self.GlobalData.dt*i, np.min(t0_ver), np.max(t0_ver)))
            t_file.close()


if __name__ == "__main__":
    if True:
        path = "D:\\DevProjects\\PythonProjects\\MES\\data\\global_data.txt"
        s1 = Siatka(path)
        #outside = [1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16]
        s1.set_bound_cond()
        #s1.list_nodes()
        val = 1/sqrt(3)
        s1.fill_H_C_global()
        s1.simulate()
        #powierzchnia = Surface_Uni((-1, -val), (-1, val), 1, 1 )
        #print(powierzchnia.get_H_BC_local(s1.Nodes[1], s1.Nodes[2], 25))
        #print(powierzchnia.get_P_local(s1.Nodes[1], s1.Nodes[2], 25, 1200))
        from pprint import pprint as pp
        #pp(s1.fill_H_global())
        if False:
            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\h_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[0], fmt='%.4f')

            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\c_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[1], fmt='%.4f')
                #s1.fill_H_global()


            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\p_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[2], fmt='%.4f')
                #s1.fill_H_global()




