import os
from scipy.spatial import Delaunay, ConvexHull
import matplotlib.pyplot as plt
import pyfar as pf


def write_evaluation_grid(
        points, name, start=200000, discard=None, show=False):
    """
    Write evaluation grid for use in Mesh2HRTF.

    Mesh2HRTF evaluation grids consist of the two text files Nodes.txt and
    Elements.txt. Evaluations grids are always triangularized.

    Parameters
    ----------
    points : pyfar Coordinates, numpy array
        pyfar Coordinates object or 2D numpy array containing the cartesian
        points of the evaluation grid in meter. The array must be of shape
        (N, 3) with N > 2.
    name : str
        Name of the folder under which the evaluation grid is saved. If the
        folder does not exist, it is created.
    start : int, optional
        The nodes and elements of the evaluation grid are numbered and the
        first element will have the number `start`. In Mesh2HRTF, each Node
        must have a unique number. The nodes/elements of the mesh for which the
        HRTFs are simulated start at 1. Thus `start` must at least be greater
        than the number of nodes/elements in the evaluation grid.
    discard : "x", "y", "z", optional
        In case all values of the evaluation grid are constant for one
        dimension, this dimension has to be discarded during the
        triangularization. E.g. if all points have a z-value of 0 (or any other
        constant), discarded must be "z". The default ``None`` does not discard
        any dimension.
    show : bool, optional
        If ``True`` the evaluation grid is plotted and the plot is saved to
        the folder given by `name`

    Examples
    --------

    Generate a spherical sampling grid with pyfar and write it to the current
    working directory

    .. plot::

        >>> import mesh2hrtf as m2h
        >>> import pyfar as pf
        >>>
        >>> points = pf.samplings.sph_lebedev(sh_order=10)
        >>> m2h.write_evaluation_grid(
        ...     points, "Lebedev_N10", discard=None, show=True)
    """

    if isinstance(points, pf.Coordinates):
        points = points.get_cart()

    if points.ndim != 2 or points.shape[0] < 3 \
            or points.shape[1] != 3:
        raise ValueError(
            "points must be a 2D array of shape (N, 3) with N > 2")

    # discard dimension
    if discard == "x":
        mask = (1, 2)
    elif discard == "y":
        mask = (0, 2)
    elif discard == "z":
        mask = (0, 1)
    else:
        mask = (0, 1, 2)

    # triangulate
    if discard is None:
        tri = ConvexHull(points[:, mask])
    else:
        tri = Delaunay(points[:, mask])

    # check output directory
    if not os.path.isdir(name):
        os.mkdir(name)

    # write nodes
    N = int(points.shape[0])
    start = int(start)

    nodes = f"{N}\n"
    for nn in range(N):
        nodes += (f"{int(start + nn)} "
                  f"{points[nn, 0]} "
                  f"{points[nn, 1]} "
                  f"{points[nn, 2]}\n")

    with open(os.path.join(name, "Nodes.txt"), "w") as f_id:
        f_id.write(nodes)

    # write elements
    N = int(tri.simplices.shape[0])
    elems = f"{N}\n"
    for nn in range(N):
        elems += (f"{int(start + nn)} "
                  f"{tri.simplices[nn, 0] + start} "
                  f"{tri.simplices[nn, 1] + start} "
                  f"{tri.simplices[nn, 2] + start} "
                  "2 0 1\n")

    with open(os.path.join(name, "Elements.txt"), "w") as f_id:
        f_id.write(elems)

    # plot the evaluation grid
    if show:
        points = pf.Coordinates(
            points[:, 0], points[:, 1], points[:, 2]).show()
        plt.savefig(os.path.join(name, "evaluation_grid.png"), dpi=300)
