class GlobalData:
    def __init__(self, path):
        data = self.__get_data__(path)
        if(data == None):
            raise Exception("Incorrect file")
        else:
            self._data_path = path
            self.H = float(data['H'])
            self.W = float(data['W'])
            self.nH = int(data['nH'])
            self.nW = int(data['nW'])
            self.k = float(data['k'])
            self.c = float(data['c'])
            self.ro = float(data['ro'])
            self.to = float(data['to'])
            self.t0 = float(data['t0'])
            self.alfa = int(data['alfa'])
            self.dt = float(data['dt'])
            self.tk = float(data['tk'])

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

