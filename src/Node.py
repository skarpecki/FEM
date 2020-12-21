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