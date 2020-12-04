class SizeFunction:
    def __init__(self, bbox, cell_size, hmin):
        if not isinstance(bbox, tuple):
            raise ValueError("`bbox` must be a tuple")
        self.bbox = bbox
        if not callable(cell_size):
            raise ValueError("`cell_size` must be callable")
        self.cell_size = cell_size
        self.hmin = hmin

    def eval(self, x):
        return self.cell_size(x)
