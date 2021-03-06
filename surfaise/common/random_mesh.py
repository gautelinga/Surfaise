""" periodic_porous script. """
import numpy as np
import dolfin as df
from meshpy import triangle as tri
import os
#from dolfin import info
from .io import remove_safe, mpi_rank, mpi_is_root, mpi_barrier
from .cmd import info
import matplotlib.pyplot as plt
import h5py


def line_points(x_left, x_right, dx):
    N = int(np.ceil(np.sqrt((x_right[0]-x_left[0])**2 +
                            (x_right[1]-x_left[1])**2)/dx))
    return list(zip(np.linspace(x_left[0], x_right[0], N+1,
                                endpoint=True).flatten(),
                    np.linspace(x_left[1], x_right[1], N+1,
                                endpoint=True).flatten()))


def round_trip_connect(start, end):
    return [(i, i+1) for i in range(start, end)] + [(end, start)]


def rad_points(x_c, rad, dx, theta_start=0., theta_stop=2*np.pi):
    if theta_stop > theta_start:
        arc_length = (theta_stop-theta_start)*rad
    else:
        arc_length = (-theta_stop+theta_start)*rad
    return [(rad * np.cos(theta) + x_c[0], rad * np.sin(theta) + x_c[1])
            for theta in np.linspace(theta_start, theta_stop,
                                     int(np.ceil(arc_length/dx)+1),
                                     endpoint=True)]


def numpy_to_dolfin(nodes, elements):
    """ Convert nodes and elements to a dolfin mesh object. """
    tmpfile = "tmp.h5"

    dim = elements.shape[1]-1
    cell_type = "triangle"
    if dim == 3:
        cell_type = "tetrahedron"

    if mpi_is_root():
        with h5py.File(tmpfile, "w") as h5f:
            cell_indices = h5f.create_dataset(
                "mesh/cell_indices", data=np.arange(len(elements)),
                dtype='int64')
            topology = h5f.create_dataset(
                "mesh/topology", data=elements, dtype='int64')
            coordinates = h5f.create_dataset(
                "mesh/coordinates", data=nodes, dtype='float64')
            topology.attrs["celltype"] = np.string_(cell_type)
            topology.attrs["partition"] = np.array([0], dtype='uint64')

    mpi_barrier()

    mesh = df.Mesh()
    h5f = df.HDF5File(mesh.mpi_comm(), tmpfile, "r")
    h5f.read(mesh, "mesh", False)
    h5f.close()

    mpi_barrier()

    remove_safe(tmpfile)
    return mesh


def place_obstacles(num_obstacles, Lx, Ly, R):
    diam2 = 4*R**2
    pts = np.zeros((num_obstacles, 2))
    for i in range(num_obstacles):
        while True:
            pt = (np.random.rand(2)-0.5) * np.array([Lx, Ly])
            if i == 0:
                break
            dist = pts[:i, :] - np.outer(np.ones(i), pt)
            for j in range(len(dist)):
                if abs(dist[j, 1]) > Ly/2:
                    dist[j, 1] = abs(dist[j, 1])-Ly
                if abs(dist[j, 0]) > Lx/2:
                    dist[j, 0] = abs(dist[j, 0])-Lx
            dist2 = dist[:, 0]**2 + dist[:, 1]**2
            if all(dist2 > diam2):
                break
        pts[i, :] = pt
    pts = pts[pts[:, 0].argsort(), :]
    obstacles = [tuple(row) for row in pts]
    return obstacles


def correct_obstacles(obstacles, rad, x_min, x_max, y_min, y_max):
    shift = None
    for i, x_c in enumerate(obstacles):
        dx_left = x_c[0] - x_min
        dx_right = x_max - x_c[0]
        dy_bottom = x_c[1] - y_min
        dy_top = y_max - x_c[1]
        is_left = dx_left < rad
        is_right = dx_right < rad
        is_bottom = dy_bottom < rad
        is_top = dy_top < rad

        d_topleft = np.sqrt(dx_left**2 + dy_top**2)
        d_topright = np.sqrt(dx_right**2 + dy_top**2)
        d_bottomleft = np.sqrt(dx_left**2 + dy_bottom**2)
        d_bottomright = np.sqrt(dx_right**2 + dy_bottom**2)

        if (is_left and is_top and d_topleft > rad):
            shift = [-dx_left, dy_top]
            break
        elif (is_left and is_bottom and d_bottomleft > rad):
            shift = [-dx_left, -dy_bottom]
            break
        elif (is_right and is_top and d_topright > rad):
            shift = [dx_right, dy_top]
            break
        elif (is_right and is_bottom and d_bottomright > rad):
            shift = [dx_right, -dy_bottom]
            break

    if shift is not None:
        x = np.array(obstacles)
        x[:, 0] += shift[0]
        x[:, 1] += shift[1]
        obstacles = x.tolist()
    return obstacles


def pos(x_c, rad, theta):
    return x_c + rad*np.array([np.cos(theta), np.sin(theta)])


def clamp(pt, x_min, x_max, y_min, y_max):
    tol = 1e-12
    x0 = list(pt)
    if abs(x0[0]-x_min) < tol:
        x0[0] = x_min
    elif abs(x0[0]-x_max) < tol:
        x0[0] = x_max
    if abs(x0[1]-y_min) < tol:
        x0[1] = y_min
    elif abs(x0[1]-y_max) < tol:
        x0[1] = y_max
    return tuple(x0)


def classify_obstacles(obstacles, rad, x_min, x_max, y_min, y_max):
    Lx, Ly = x_max-x_min, y_max-y_min
    interior_obstacles = []
    exterior_obstacles = []
    new_obstacles = []
    for i, x_c in enumerate(obstacles):

        is_left = x_c[0] < x_min + rad
        is_right = x_c[0] > x_max - rad
        is_bottom = x_c[1] < y_min + rad
        is_top = x_c[1] > y_max - rad

        if is_left:
            new_obstacles.append((x_c[0]+Lx, x_c[1]))
        if is_right:
            new_obstacles.append((x_c[0]-Lx, x_c[1]))
        if is_bottom:
            new_obstacles.append((x_c[0], x_c[1]+Ly))
        if is_top:
            new_obstacles.append((x_c[0], x_c[1]-Ly))
        if is_left and is_top:
            new_obstacles.append((x_c[0]+Lx, x_c[1]-Ly))
        if is_left and is_bottom:
            new_obstacles.append((x_c[0]+Lx, x_c[1]+Ly))
        if is_right and is_top:
            new_obstacles.append((x_c[0]-Lx, x_c[1]-Ly))
        if is_right and is_bottom:
            new_obstacles.append((x_c[0]-Lx, x_c[1]+Ly))
        if is_left or is_right or is_top or is_bottom:
            exterior_obstacles.append(x_c)
        else:
            interior_obstacles.append(x_c)

    obst = exterior_obstacles + new_obstacles
    return interior_obstacles, exterior_obstacles, obst


def compute_intersections(obst, rad, x_min, x_max, y_min, y_max):
    theta_low = [[] for _ in range(len(obst))]
    theta_high = [[] for _ in range(len(obst))]

    for i, x_c in enumerate(obst):
        is_left = x_c[0] < x_min + rad
        is_right = x_c[0] > x_max - rad
        is_bottom = x_c[1] < y_min + rad
        is_top = x_c[1] > y_max - rad

        if is_left:
            rx = x_min-x_c[0]
            theta = np.arccos(rx/rad)
            ry = rad*np.sin(theta)
            if rx < 0.:
                theta = np.pi-theta

            x_cross = pos(x_c, rad, theta)
            if x_cross[0] != x_min:
                theta = np.arccos(-(x_cross[0]-x_min)/rad+np.cos(theta))

            y_h_loc = x_c[1] + ry
            y_l_loc = x_c[1] - ry
            if y_l_loc > y_min:
                theta_low[i].append(-theta)
            if y_h_loc < y_max:
                theta_high[i].append(theta)

        if is_right:
            rx = x_max-x_c[0]
            theta = np.arccos(rx/rad)
            ry = rad*np.sin(theta)
            if rx < 0.:
                theta = np.pi-theta

            x_cross = pos(x_c, rad, theta)
            if x_cross[0] != x_max:
                theta = np.arccos(-(x_cross[0]-x_max)/rad+np.cos(theta))

            y_h_loc = x_c[1] + ry
            y_l_loc = x_c[1] - ry
            if y_l_loc > y_min:
                theta_high[i].append(-theta)
            if y_h_loc < y_max:
                theta_low[i].append(theta)

        if is_bottom:
            ry = y_min-x_c[1]
            theta = np.arcsin(ry/rad)
            rx = rad*abs(np.sin(theta))

            x_l_loc = x_c[0] - rx
            x_h_loc = x_c[0] + rx

            if x_l_loc > x_min:
                theta_high[i].append(np.pi-theta)
            if x_h_loc < x_max:
                theta_low[i].append(theta)

        if is_top:
            ry = y_max-x_c[1]
            theta = np.arcsin(ry/rad)
            rx = rad*abs(np.sin(theta))

            x_l_loc = x_c[0] - rx
            x_h_loc = x_c[0] + rx

            if x_l_loc > x_min:
                theta_low[i].append(np.pi-theta)
            if x_l_loc < x_max:
                theta_high[i].append(theta)
    return theta_low, theta_high


def draw_curves(obst, theta_low, theta_high, rad, dx):
    curves = []
    for i, x_c in enumerate(obst):
        for theta_l in theta_low[i]:
            dtheta = (np.array(theta_high[i])-theta_l) % (2*np.pi)
            theta_h = theta_l + dtheta.min()

            rad_pts = rad_points(x_c, rad, dx,
                                 theta_start=theta_l,
                                 theta_stop=theta_h)
            curves.append(rad_pts)
    return curves


def construct_segments(curve_start, curve_stop, x_min, x_max, y_min, y_max):
    cross_stop = np.asarray(list(curve_stop.keys()))
    cross_start = np.asarray(list(curve_start.keys()))

    y_left_h = cross_stop[cross_stop[:, 0] == x_min, 1]
    y_left_l = cross_start[cross_start[:, 0] == x_min, 1]
    y_right_h = cross_start[cross_start[:, 0] == x_max, 1]
    y_right_l = cross_stop[cross_stop[:, 0] == x_max, 1]
    x_top_h = cross_stop[cross_stop[:, 1] == y_max, 0]
    x_top_l = cross_start[cross_start[:, 1] == y_max, 0]
    x_bottom_h = cross_start[cross_start[:, 1] == y_min, 0]
    x_bottom_l = cross_stop[cross_stop[:, 1] == y_min, 0]

    y_left_h_list = sorted(y_left_h.tolist())
    y_left_l_list = sorted(y_left_l.tolist())
    y_right_h_list = sorted(y_right_h.tolist())
    y_right_l_list = sorted(y_right_l.tolist())
    x_top_h_list = sorted(x_top_h.tolist())
    x_top_l_list = sorted(x_top_l.tolist())
    x_bottom_h_list = sorted(x_bottom_h.tolist())
    x_bottom_l_list = sorted(x_bottom_l.tolist())

    empty = not len(y_left_l) or not len(y_left_h)
    if empty or y_left_l.min() < y_left_h.min():
        y_left_h_list = [y_min] + y_left_h_list
        x_bottom_h_list = [x_min] + x_bottom_h_list

    if empty or y_left_h.max() > y_left_l.max():
        y_left_l_list = y_left_l_list + [y_max]
        x_top_h_list = [x_min] + x_top_h_list

    if empty or y_right_l.min() < y_right_h.min():
        y_right_h_list = [y_min] + y_right_h_list
        x_bottom_l_list = x_bottom_l_list + [x_max]

    if empty or y_right_h.max() > y_left_l.max():
        y_right_l_list = y_right_l_list + [y_max]
        x_top_l_list = x_top_l_list + [x_max]

    y_left = list(zip(y_left_h_list, y_left_l_list))
    y_right = list(zip(y_right_h_list, y_right_l_list))
    x_top = list(zip(x_top_h_list, x_top_l_list))
    x_bottom = list(zip(x_bottom_h_list, x_bottom_l_list))

    segments = []
    for y_a, y_b in y_left:
        pt_a = (x_min, y_a)
        pt_b = (x_min, y_b)
        segments.append((pt_a, pt_b))

    for y_a, y_b in y_right:
        pt_a = (x_max, y_a)
        pt_b = (x_max, y_b)
        segments.append((pt_b, pt_a))

    for x_a, x_b in x_top:
        pt_a = (x_a, y_max)
        pt_b = (x_b, y_max)
        segments.append((pt_a, pt_b))

    for x_a, x_b in x_bottom:
        pt_a = (x_a, y_min)
        pt_b = (x_b, y_min)
        segments.append((pt_b, pt_a))

    segments = dict(segments)
    return segments


def get_curve_intersection_points(curves, x_min, x_max, y_min, y_max):
    curve_start = dict()
    curve_stop = dict()
    for i, curve in enumerate(curves):
        curve[0] = clamp(curve[0], x_min, x_max, y_min, y_max)
        curve[-1] = clamp(curve[-1], x_min, x_max, y_min, y_max)
        curve_start[curve[0]] = i
        curve_stop[curve[-1]] = i
    return curve_start, curve_stop


def discretize_loop(pt_start, curve_start, curves, segments, dx):
    pts = []
    pt = pt_start
    while True:
        if pt in curve_start:
            ic = curve_start[pt]
            pt_ = curves[ic][-1]
            pts.extend(curves[ic][1:])
        else:
            pt_ = segments[pt]
            pts.extend(line_points(pt, pt_, dx)[1:])
        pt = pt_
        if pt == pt_start:
            break
    return pts


def method(Lx=4., Ly=4., num_obstacles=25,
           rad=0.25, R=0.3, dx=0.05, seed=123, show=False,
           verbose=False,
           **kwargs):
    x_min, x_max = -Lx/2, Lx/2
    y_min, y_max = -Ly/2, Ly/2

    np.random.seed(seed)
    obstacles = place_obstacles(num_obstacles, Lx, Ly, R)
    obstacles = correct_obstacles(obstacles, rad, x_min, x_max, y_min, y_max)

    interior_obstacles, exterior_obstacles, obst = classify_obstacles(
        obstacles, rad, x_min, x_max, y_min, y_max)

    theta_low, theta_high = compute_intersections(
        obst, rad, x_min, x_max, y_min, y_max)

    curves = draw_curves(obst, theta_low, theta_high, rad, dx)
    if len(curves) > 0:
        curve_start, curve_stop = get_curve_intersection_points(
            curves, x_min, x_max, y_min, y_max)

        segments = construct_segments(curve_start, curve_stop,
                                      x_min, x_max, y_min, y_max)
        pt_start = curves[0][0]
    else:
        curve_start = []
        curve_stop = []
        segments = [((x_min, y_min), (x_max, y_min)),
                    ((x_max, y_min), (x_max, y_max)),
                    ((x_max, y_max), (x_min, y_max)),
                    ((x_min, y_max), (x_min, y_min))]
        segments = dict(segments)
        pt_start = (x_min, y_min)

    if show:
        for x_a, x_b in segments.items():
            x = np.array([x_a[0], x_b[0]])
            y = np.array([x_a[1], x_b[1]])
            plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1],
                       scale_units='xy', angles='xy', scale=1)
        for curve in curves:
            x = np.array(curve)
            plt.quiver(x[:-1, 0], x[:-1, 1],
                       x[1:, 0]-x[:-1, 0],
                       x[1:, 1]-x[:-1, 1],
                       scale_units='xy', angles='xy', scale=1)

    pts = discretize_loop(pt_start, curve_start,
                          curves, segments, dx)[::-1]

    plt.show()

    edges = round_trip_connect(0, len(pts)-1)

    for interior_obstacle in interior_obstacles:
        pts_obstacle = rad_points(interior_obstacle, rad, dx)[1:]
        edges_obstacle = round_trip_connect(len(pts),
                                            len(pts)+len(pts_obstacle)-1)

        pts.extend(pts_obstacle)
        edges.extend(edges_obstacle)

    mi = tri.MeshInfo()
    mi.set_points(pts)
    mi.set_facets(edges)
    mi.set_holes(interior_obstacles)

    max_area = 0.5*dx**2

    mesh = tri.build(mi, max_volume=max_area, min_angle=25,
                     allow_boundary_steiner=False)

    coords = np.array(mesh.points)
    faces = np.array(mesh.elements)

    pp = [tuple(point) for point in mesh.points]
    if verbose:
        info("Number of points:     {}".format(len(pp)))
        info("Number unique points: {}".format(len(set(pp))))

    msh = numpy_to_dolfin(coords, faces)

    # mesh_path = os.path.join(MESHES_DIR,
    #                          "periodic_porous_Lx{}_Ly{}_rad{}_N{}_dx{}".format(
    #                              Lx, Ly, rad, num_obstacles, dx))
    # store_mesh_HDF5(msh, mesh_path)

    # obstacles_path = os.path.join(
    #     MESHES_DIR,
    #     "periodic_porous_Lx{}_Ly{}_rad{}_N{}_dx{}.dat".format(
    #         Lx, Ly, rad, num_obstacles, dx))

    if len(obst) and len(interior_obstacles):
        all_obstacles = np.vstack((np.array(obst),
                                   np.array(interior_obstacles)))
    elif len(interior_obstacles):
        all_obstacles = interior_obstacles
    else:
        all_obstacles = []

    # if len(all_obstacles):
    #     np.savetxt(obstacles_path,
    #                np.hstack((all_obstacles,
    #                           np.ones((len(all_obstacles), 1))*rad)))
    return msh


def unit_square_mesh(dx):
    return rectangle_mesh(1., 1., dx)


def rectangle_mesh(Lx, Ly, dx):
    msh = method(Lx=Lx, Ly=Ly, num_obstacles=0,
                 dx=dx, show=False)
    x = msh.coordinates()
    x[:, 0] += Lx/2
    x[:, 1] += Ly/2
    return msh


if __name__ == "__main__":
    msh = unit_square_mesh(0.025)

    df.plot(msh)
    plt.show()
