import numpy as np
import yt

if __name__ == '__main__':
    ds = yt.load("./parthenon.prim.00000.phdf")
    yt.SlicePlot(ds, "x", ("gas", "density"), width=(200.0, "kpc")).save()
    