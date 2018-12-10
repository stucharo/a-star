from math import pi, sin, cos, atan2, atan

def dist(a, b):
    return ((b[1] - a[1])**2 + (b[0] - a[0])**2)**0.5

def dubins_routes(start, goal, rad, spacing=1):
    paths = get_dubins_paths(start, goal, rad)
    routes = [get_dubins_route_points(path, spacing) for path in paths]
    return routes

def get_dubins_paths(start, goal, rad):
    CSCs = [('l', 'l'), ('r', 'r'), ('l', 'r'), ('r', 'l')]
    paths = []
    for CSC in CSCs:
        p = CSC_path(start, CSC[0], goal, CSC[1], rad)
        if p is not None:
            paths.append(p)
    # CCCs = [('l', 'r', 'l'), ('r', 'l', 'r')]
    return paths

def CSC_path(start, start_dir, goal, goal_dir, rad):
    cs = circle_centre(start, rad, start_dir)
    cg = circle_centre(goal, rad, goal_dir)
    ts = tangents((cs[0], cs[1], rad, start_dir), (cg[0], cg[1], rad, goal_dir))
    if ts is None:
        return ts
    else:
        return [(start, ts[0], cs, start_dir), (ts[1], goal, cg, goal_dir)]

def circle_centre(pt, rad, dir):
    d = {'l': 1,
         'r': -1}
    x = pt[0] - rad * sin(pt[2] + d[dir] * pi/2)
    y = pt[1] - rad * cos(pt[2] + d[dir] * pi/2)
    return x, y

def tangents(p1, p2):
    """ Returns list of outer and inner tangents (in that order).
    Tangents are returned as [x1, y1, x2, y2]
    Modified version of Java algorithm:
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Tangents_between_two_circles

    p1 and p2 are tuples describing the circles and tangent required:
        (cx, cy, rad(>0), d(['l'||'r']))
        cx: float
            x coordinate of circle centre
        cy: float
            y coordinate of circle centre
        rad: float
            circle radius
        d: string ('l' or 'r')
            direction of turning around arc    
    """
    d = dist(p1, p2)
    vx = (p2[0] - p1[0]) / d
    vy = (p2[1] - p1[1]) / d

    signs = {'l': 1,
             'r': -1}

    sign1 = signs[p1[3].lower()]
    sign2 = signs[p2[3].lower()]

    c = (p1[2] - sign1 * p2[2]) / d
    if c**2 > 1:
        return None

    h = (max(0.0, 1.0 - c**2))**0.5

    nx = vx * c - sign2 * h * vy
    ny = vy * c + sign2 * h * vx

    x1 = p1[0] + p1[2] * nx
    y1 = p1[1] + p1[2] * ny
    x2 = p2[0] + sign1 * p2[2] * nx
    y2 = p2[1] + sign1 * p2[2] * ny
    dx = x2-x1
    dy = y2-y1
    if dy == 0:
        if dx > 0:
            t = pi / 2
        else:
            t = 3 * pi / 2
    else:
        t = pi / 2 - atan2(dy, dx)
    return ((x1, y1, t), (x2, y2, t))

def get_dubins_route_points(path, spacing):
    """ Dubins routes are defined by the key locations on each circle. There
    is 2 options to define a Dubins path:
        CSC: two arcs joined by a straight line,
        CCC: two arcs joined by a third arc tangent to both start goal arcs.

    Paths can be defined by arc start, end, centre and direction. 2 or 3
    arcs are required and are passed in as a list of arcs.
    
    """
    if len(path) == 2:
        #we're in a CSC path
        Cs_l, Cs_p = arc_points(path[0], spacing)
        leftover = Cs_l - (len(Cs_p) - 1) * spacing
        S_l, S_p = straight_points(path[0][1], path[1][1], spacing, leftover)
        leftover = (Cs_l + S_l) - (len(Cs_p) + len(S_p) - 1) * spacing
        Cg_l, Cg_p = arc_points(path[1], spacing, leftover)
        p = np.concatenate((Cs_p, S_p, Cg_p), axis=0)
        path_length = Cs_l + S_l + Cg_l
    elif len(path) == 3:
        # We're in a CCC path
        pass
    return path_length, p

def straight_points(p1, p2, spacing, leftover=0):
    path_length = dist(p1, p2)
    start_x = p1[0]
    start_y = p1[1]
    if leftover != 0:
        leftover_l = spacing - leftover
        start_x += leftover_l * sin(p1[2])
        start_y += leftover_l * cos(p1[2])
    num_divs = int(dist((start_x, start_y), p2) / spacing)
    dx = spacing * sin(p1[2])
    dy = spacing * cos(p1[2])
    xs = np.linspace(start_x, start_x + dx*num_divs, num_divs+1)
    ys = np.linspace(start_y, start_y + dy*num_divs, num_divs+1)
    pts = np.stack((xs, ys), axis=-1)
    return path_length, pts

def arc_points(path, spacing, leftover=0):
    p1 = path[0]
    p2 = path[1]
    c = path[2]
    d = path[3]
    if d == 'r':
        start_theta = p1[2] - pi / 2
        end_theta = p2[2] - pi / 2
        theta = end_theta - start_theta
        if theta < 0:
            theta += 2 * pi
    elif d == 'l':
        start_theta = p1[2] + pi / 2
        end_theta = p2[2] + pi / 2
        theta = end_theta - start_theta
        if theta > 0:
            theta -= 2 * pi
    else:
        raise ValueError("d must be 'l' or 'r'.")
    r = dist(p1, c)
    path_length = abs(theta) * r
    theta_inc = spacing / r
    num_thetas = int(theta / theta_inc)
    thetas = np.linspace(start_theta, start_theta+num_thetas*theta_inc, abs(num_thetas)+1)
    xs = c[0] + r * np.sin(thetas)
    ys = c[1] + r * np.cos(thetas)
    pts = np.stack((xs, ys), axis=-1)
    return path_length, pts

def arc_theta(p1, p2, d):
    theta = p2[2] - p1[2]
    if d == "left":
        if theta > 0:
            theta -= 2 * pi
    elif d == "right":
        if theta < 0:
            theta += 2 * pi
    else:
        raise ValueError("d must be 'left' or 'right'.")
    return theta

def CCC_path(start, start_dir, goal, goal_dir, rad):
    # get centre of Cstart
    Cs = circle_centre(start, rad, start_dir)
    Cg = circle_centre(goal, rad, goal_dir)
    return [(0,0)]

import numpy as np
import random
import matplotlib.pyplot as plt

# p1 = (random.randint(0,100), random.randint(0,100), random.uniform(0, 2*pi))
# p2 = (random.randint(0,100), random.randint(0,100), random.uniform(0, 2*pi))
# rad = random.randint(1, 100)

p1 = (50, 50, 0)
p2 = (100, 100, pi/2)
rad = 5

routes = dubins_routes(p1, p2, rad, spacing=0.01)
fig, ax = plt.subplots()
ax.arrow(p1[0],p1[1],5*sin(p1[2]), 5*cos(p1[2]), head_width=0.5, head_length=1)
ax.arrow(p2[0],p2[1],5*sin(p2[2]), 5*cos(p2[2]), head_width=0.5, head_length=1)
for l, p in zip(['LSL', 'RSR', 'LSR', 'RSL'], routes):
    pts = np.array(p[1])
    x = pts[:,0]
    y = pts[:,1]
    ax.plot(x, y, label=l)
plt.axis('equal')
plt.legend()
plt.show()