from math import pi, sin, cos, atan2, atan, acos
import numpy as np

class Point:

    def __init__(self, x, y, theta=None):
        self.x = x
        self.y = y
        self.theta = theta

class Tangent:

    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        if self.p1.theta == self.p2.theta:
            self.theta = self.p1.theta
        else:
            raise ValueError('Both points on tangent must be parallel.')

class Circle:

    def __init__(self, centre, radius, direction):
        self.centre = pc
        self.radius = radius
        if direction.lower() in ['l', 'r']:
            self.direction = direction.lower()
        else:
            raise ValueError("direction must be 'l' or 'r'.")

def dist(a, b):
    return ((b.x - a.x)**2 + (b.y - a.y)**2)**0.5

def dubins_routes(start, goal, min_bend_radius, min_straight_length, spacing):
    paths = get_dubins_paths(start, goal, min_bend_radius, min_straight_length, spacing)
    routes = [get_dubins_route_points(path, spacing) for path in paths]
    return routes


def get_dubins_paths(start, goal, min_bend_radius, min_straight_length, spacing):
    """ These will only be true Dubins routes in the cases where the
    start is straight and the tangent length exceeds the minimum straight
    length.
    
    The shortest path to the goal must satisfy our regular route 
    constraints or it will return artificially cheap paths. These are:
        1. Minimum straight length;
        2. Fixed radius bends;
        3. Minimum bends radius.

    If we are currently on a bend then we must either:
        a. Continue at that bend radius to a tangent point;
        b. Continue at the same heading for the minium straight
           length.

    If we are on a straight section then our shortest path must be to
    turn immediately at the minimum bend radius.

    When we are close to the goal i.e. when a tangent length is less
    than the minimum straight length, we must add in a third circle
    the minimum straight length away of the minimum bend radius. When
    the minimum straight length reduces to zero we will return a
    conventional CCC path.
    """
    ts = CSC_tangents(start, goal, min_bend_radius, spacing, min_straight_length)
    if len(ts) == 0:
        # We have no valid CSC paths because they all violate our minimum straight length
        # constraint. Now we need to compute CSCSC paths. We have far more options here
        # because we have an extra circle. Potential routes are 'RSRSR', 'RSRSL', 'RSLSR',
        # 'RSLSL', 'LSRSR', 'LSRSL', 'LSLSR' and 'LSLSL'. The direction of the mid circle
        # determines it's location between the start and goal.
        ts = CSCSC_tangents(start, goal, min_bend_radius, spacing, min_straight_length)
    
    paths = []
    CSCs = [('r', 'r'), ('l', 'l'), ('r', 'l'), ('l', 'r')]
    for CSC in CSCs:
        p = CSC_path(start, CSC[0], goal, CSC[1], rad)
        if p is not None:
            paths.append(p)
    if dist(start, goal) < (4 * rad):
        CCCs = ['l', 'r']
        for CCC in CCCs:
            p = CCC_path(start, goal, CCC, rad)
            if p is not None:
                paths.append(p)
    return paths

def CSC_tangents(start, goal, min_bend_rad, spacing, min_straight_length):
    """ Compute tangent start and end points for a
    Circle-Straight-Circle path.
    """
    ts = []
    if start.theta > 0:
        # we are already on a bend and we know the start direction is left (and we must
        # continue it to maintain our constant bend radius criteria). This limits us to
        # 2 possible paths: 'LSL' and 'LSR'.
        rs = spacing / abs(start.theta)
        cs = make_circle(start, rs, 'l')
        for d in ['l', 'r']:
            cg = make_circle(goal, min_bend_rad, d)
            t = tangent(cs, cg)
            if t is not None and dist(t.p1, t.p2) > min_straight_length:
                ts.append(t)
    elif start.theta > 0:
        # we are already on a bend and we know the start direction is right (and we must
        # continue it to maintain our constant bend radius criteria). This limits us to
        # 2 possible paths: 'RSL' and 'RSR'.
        rs = spacing / abs(start.theta)
        cs = make_circle(start, rs, 'r')
        for d in ['l', 'r']:
            cg = make_circle(goal, min_bend_rad, d)
            t = tangent(cs, cg)
            if t is not None and dist(t.p1, t.p2) > min_straight_length:
                ts.append(t)
    else:
        # we are on a straight section so we should expect up to 4 tangents to be
        # returned. These are the 'RSR', 'LSL', 'RSL' and 'LSR'. When a tangent is 'None'
        # it can be discarded as another CCC path would be shorter anyway.
        for d in [('r', 'r'), ('l', 'l'), ('r', 'l'), ('l', 'r')]:
            cs = make_circle(start, min_bend_rad, d[0])
            cg = make_circle(goal, min_bend_rad, d[1])
            t = tangent(cs, cg)
            if t is not None and dist(t.p1, t.p2) > min_straight_length:
                ts.append(t)
    return ts
    
def CSCSC_tangents(start, goal, min_bend_rad, spacing, min_straight_length):
    ts = []
    # Find the distances C1 (from the centre of the start circle to the centre of the mid
    # circle) and C2 (from the centre of mid circle to the centre of the goal circle) for
    # the inner and outer tangents the mid and goal circles must always be of minimum bend
    # radius.
    C2_o = min_straight_length
    C2_i = (min_straight_length**2 + (2 * min_bend_rad)**2)**0.5
    # but the start circle radius is based on the radius of the start point
    if start.theta == 0:
        # we are on a straight so the start circle must also be of minimum bend radius
        C1_o = C2_o
        C1_i = C2_i
    else:
        # however, if we're on an arc already then we must maintain the bend radius
        # through the first arc
        rs = spacing / abs(start.theta)
        C1_o = (min_straight_length**2 + (rs - min_bend_rad)**2)**0.5
        C1_i = (min_straight_length**2 + (rs + min_bend_rad)**2)**0.5
    
    CSCSC_paths = [('l', 'r', 'l'), ('l', 'r', 'r'), ('l', 'l', 'l'), ('l', 'l', 'r'),
                   ('r', 'r', 'l'), ('r', 'r', 'r'), ('r', 'l', 'l'), ('r', 'l', 'r')]
    
    # limit the possible paths if we're already on a curve
    if start.theta < 0:
        # we have to keep turning left
        CSCSC_paths = [p for p in CSCSC_paths if p[0] == 'l']
    if start.theta > 0: 
        # we have to keep turning right
        CSCSC_paths = [p for p in CSCSC_paths if p[0] == 'r']

    for CSCSC in CSCSC_paths:
        # find the distance between circle centres for this case
        if CSCSC[0] == CSCSC[1]:
            C1 = C1_o
        else:
            C1 = C1_i
        if CSCSC[1] == CSCSC[2]:
            C2 = C2_o
        else:
            C2 = C2_i
        # get the centre of the start circle
        if start.theta == 0:
            rs = min_bend_rad
        else:
            rs = spacing / abs(start.theta)
        cs = make_circle(start, rs, CSCSC[0])
        # get the centre of the goal circle:
        cg = make_circle(goal, min_bend_rad, CSCSC[2])
        # get the heading from cs to cg
        theta_cs_cg = atan((cg.centre.x-cs.centre.x)/(cg.centre.y-cs.centre.y))
        d = dist(cs.centre, cg.centre)
        # get the angle between vector cs_cg and vector cs_cm
        theta_start_mid = (C1**2 + d**2 - C2**2) / (2 * C1 * d)
        if CSCSC[1] == 'r':
            theta_start_mid *= -1
        # get the heading from cs to cm
        theta_cs_cm = theta_cs_cg + theta_start_mid
        # make a point at the mid circle centre
        mp = Point(cs.centre.x + C1*sin(theta_cs_cm), cs.centre.y + C1*cos(theta_cs_cm))
        cm = Circle(mp, min_bend_rad, CSCSC[1])
        # get tangents
        t1 = tangent(cs, cm)
        t2 = tangent(cm, cg)
        ts.append((t1, t2))
    return ts

def CSC_path(start, start_dir, goal, goal_dir, rad):
    cs = make_circle(start, rad, start_dir)
    cg = make_circle(goal, rad, goal_dir)
    ts = tangent(cs, cg)
    if ts is None:
        return ts
    else:
        return [(start, ts[0], cs, start_dir), (ts[1], goal, cg, goal_dir)]

def make_circle(pt, rad, dir):
    d = {'l': 1,
         'r': -1}
    x = pt[0] - rad * sin(pt[2] + d[dir] * pi/2)
    y = pt[1] - rad * cos(pt[2] + d[dir] * pi/2)
    return Circle(Point(x, y), rad, dir)

def tangent(c1, c2):
    """ Returns the Tangent to c1 and c2, travelling in the appropriate directions.
    Tangents are returned as Tangent objects
    Modified version of Java algorithm:
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Tangents_between_two_circles

    c1 and c2 are tuples describing the circles and tangent required:
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

    signs = {('r', 'r'): (1, 1),
             ('r', 'l'): (-1, 1),
             ('l', 'r'): (-1, -1),
             ('l', 'l'): (1, -1)}

    (sign1, sign2) = signs[(p1[3].lower(), p2[3].lower())]

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
    return Tangent(Point(x1, y1, t), Point(x2, y2, t))

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
        S_l, S_p = straight_points(path[0][1], path[1][0], spacing, leftover)
        leftover = (Cs_l + S_l) - (len(Cs_p) + len(S_p) - 1) * spacing
        Cg_l, Cg_p = arc_points(path[1], spacing, leftover)
        p = np.concatenate((Cs_p, S_p, Cg_p), axis=0)
        path_length = Cs_l + S_l + Cg_l
    elif len(path) == 3:
        # We're in a CCC path
        Cs_l, Cs_p = arc_points(path[0], spacing)
        leftover = Cs_l - (len(Cs_p) - 1) * spacing
        Cm_l, Cm_p = arc_points(path[1], spacing)
        leftover = (Cs_l + Cm_l) - (len(Cs_p) + len(Cm_p) - 1) * spacing
        Cg_l, Cg_p = arc_points(path[2], spacing, leftover)
        p = np.concatenate((Cs_p, Cm_p, Cg_p), axis=0)
        path_length = Cs_l + Cm_l + Cg_l
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
        if len(p2) == 2:
            pass
        else:
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

def CCC_path(start, goal, dir, rad):
    # get centre of Cstart
    Cs = make_circle(start, rad, dir)
    Cg = make_circle(goal, rad, dir)
    dx = Cg[0] - Cs[0]
    dy = Cg[1] - Cs[1]
    d = dist(Cs, Cg)
    if d > 4 * rad:
        return None
    if dir == 'l':
        theta = atan2(dy, dx) - acos(d / 4 / rad)
    elif dir == 'r':
        theta = atan2(dy, dx) + acos(d / 4 / rad)
    Cm = (Cs[0] + 2 * rad * cos(theta), Cs[1] + 2 * rad * sin(theta))
    ts = []
    for c in [Cs, Cg]:
        tx = Cm[0] + (c[0] - Cm[0])/2
        dx = tx - c[0]
        ty = Cm[1] + (c[1] - Cm[1])/2
        dy = ty - c[1]
        theta = (-atan2(dy, dx) + 5 * pi / 2)
        if dir == 'l':
            theta -= pi / 2
        else:
            theta += pi / 2
        ts.append((tx, ty, theta % (2 * pi)))
    if dir == 'l':
        dir_m = 'r'
    else:
        dir_m = 'l'
    
    return [(start, ts[0], Cs, dir), (ts[0], ts[1], Cm, dir_m), (ts[1], goal, Cg, dir)]

if __name__ == '__main__':

    import random
    import matplotlib.pyplot as plt

    p1 = (random.randint(0,100), random.randint(0,100), random.uniform(0, 2*pi))
    p2 = (random.randint(0,100), random.randint(0,100), random.uniform(0, 2*pi))
    rad = random.randint(1, 20)

    # p1 = (50, 50, 0)
    # p2 = (100, 100, pi/2)
    # rad = 20

    routes = dubins_routes(p1, p2, rad, spacing=0.01)
    fig, ax = plt.subplots()
    ax.arrow(p1[0],p1[1],5*sin(p1[2]), 5*cos(p1[2]), head_width=0.5, head_length=1)
    ax.arrow(p2[0],p2[1],5*sin(p2[2]), 5*cos(p2[2]), head_width=0.5, head_length=1)
    for p in routes:
        pts = np.array(p[1])
        x = pts[:,0]
        y = pts[:,1]
        ax.plot(x, y)
    plt.axis('equal')
    plt.show()