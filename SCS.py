import numpy as np
from math import pi, sin, cos, atan2

class Point:

    def __init__(self, x, y, heading):
        self.x = x
        self.y = y
        self.dx = cos(heading)
        self.dy = sin(heading)

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def dist(p1, p2):
    return ((p2.x-p1.x)**2 - (p2.y-p1.y)**2)**0.5

def dist2intersect(p1, p2):
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x-p1.x, p2.y-p1.y])
    return np.linalg.solve(dv, lv)

def angle2intersect(p1, p2):
    return np.arccos(np.clip(np.dot([p1.dx, p1.dy], [p2.dx, p2.dy]), -1.0, 1.0))

def turn_dir(p1, p2):
    return np.cross([p1.dx, p1.dy], [p2.dx, p2.dy])

def dist2tangents(p1, p2, r):
    t_dir = turn_dir(p1, p2)
    rp1 = rotate(p1, t_dir*pi/2)
    rp2 = rotate(p2, t_dir*pi/2)
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x+r*rp2.dx-(p1.x+r*rp1.dx), p2.y+r*rp2.dy-(p1.y+r*rp1.dy)])
    return np.linalg.solve(dv, lv)

def rotate(p, t):
    in_angle = atan2(p.dy, p.dx)
    new_heading = in_angle + t
    if new_heading > pi:
        new_heading -= 2*pi
    elif new_heading < -pi:
        new_heading += 2*pi
    return Point(p.x, p.y, new_heading)


def shortest_route(s, g, min_rad, min_straight, heading_tol, location_tol, spacing):
    """ An SCS path is an SCSCS path where the straight length is 0.
    If this invalidates our minimum straight length then we need to
    increase the bend radius to make the turn in an SCS. 
    """
    if abs(angle2intersect(s, g)) < heading_tol:
        d = dist(s, g)
        if d > min_straight:
            return straight_points(s, g, spacing)
    # find tangent points
    st = Point(s.x+min_straight*s.dx, s.y+min_straight*s.dy, atan2(s.dy, s.dx))
    gt = Point(g.x-min_straight*g.dx, g.y-min_straight*g.dy, atan2(g.dy, g.dx))
    # and the start and goal turn directions
    st_dir = turn_dir(st, gt)
    gt_dir = turn_dir(gt, st)
    # get the centres of both arcs
    sc = rotate(st, st_dir*2/pi)
    sc.x += min_rad*sc.dx
    sc.y += min_rad*sc.dy
    gc = rotate(gt, gt_dir*pi/2)
    gc.x += min_rad*gc.dx
    gc.y += min_rad*gc.dy
    dc = dist(sc, gc)
    # if we can fit a valid straight length in:
    if ((st_dir == gt_dir and dc > min_straight)
        or (st_dir != gt_dir and (dc**2 - (2*min_rad)**2)**0.5 > min_straight)):
            return SCSCS_point(s, st, sc, st_dir, g, gt, gc, gt_dir, min_rad, spacing)
    else:
        # figure out whether we need a SCSCSCS path, or just a SCS path
        pass

def SCSCS_point(s, st, sc, st_dir, g, gt, gc, gt_dir, rad, spacing):
    sc = Point(sc.x, sc.y, atan2(gc.y-sc.y, gc.x-sc.x))
    rv = rotate(sc, -st_dir*acos(2*rad/dist(sc, gc)))
    c1n = angle_c1c2 + st_dir*t
    sc = Point(sc.x, sc.y, atan2(gc.y-sc.y, gc.x-sc.x))
    p = straight_points(s, st, spacing)
    leftover = dist(s, st) - (len(p)-1)*spacing





def straight_points(s, g, spacing, leftover=0):
    if leftover != 0:
        ds = spacing - leftover
        s.x += ds * s.dx
        s.y += ds * s.dy
    d = dist(s, g)
    incs = int(d/spacing)
    xs = np.linspace(s.x, s.x+s.dx*incs, incs+1)
    ys = np.linspace(s.y, s.y+s.dy*incs, incs+1)
    heading = atan2(s.dy, s.dx)
    return [Point(x, y, heading=heading) for x, y in np.stack((xs, ys), axis=-1)]

if __name__ == '__main__':
    
    s = Point(0, 0, pi/2)
    g = Point(10, 10, 0)
    min_rad = 2
    min_straight = 3
    make_SCS(s, g, min_rad, min_straight)


