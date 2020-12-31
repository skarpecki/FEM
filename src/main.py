from src.Mesh import Mesh
from math import sqrt

if __name__ == "__main__":
    if True:
        path = "/src/global_data.txt"
        s1 = Mesh(path)
        s1.set_bound_cond()
        val = 1/sqrt(3)
        s1.fill_H_C_global()
        s1.simulate()
        if False:
            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\h_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[0], fmt='%.4f')

            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\c_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[1], fmt='%.4f')
                #s1.fill_H_global()


            with open(rf"D:\DevProjects\PythonProjects\MES\data\results\p_global\{s1.GlobalData.npc}npc.txt", "w") as a_file:
                np.savetxt(a_file, s1.fill_H_C_global()[2], fmt='%.4f')
                #s1.fill_H_global()
