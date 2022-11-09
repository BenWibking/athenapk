import numpy as np
import yt

if __name__ == '__main__':
    ds = yt.load("./parthenon.prim.final.phdf")
    yt.SlicePlot(ds, "x", ("gas", "density"), width=(200.0, "kpc")).save()
    