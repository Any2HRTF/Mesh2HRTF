import mesh2hrtf as m2h
import pyfar as pf
points = pf.samplings.sph_lebedev(sh_order=10)
m2h.write_evaluation_grid(
    points, "Lebedev_N10", discard=None, show=True)
