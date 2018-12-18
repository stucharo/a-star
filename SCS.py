import numpy as np
from math import pi, sin, cos, asin, acos, atan2

class Point:

    def __init__(self, x, y, heading):
        self.x = x
        self.y = y
        self.dx = cos(heading)
        self.dy = sin(heading)

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def dist(p1, p2):
    return ((p2.x-p1.x)**2 + (p2.y-p1.y)**2)**0.5

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
    new_heading = normalize_angle(in_angle + t)
    return Point(p.x, p.y, new_heading)

def normalize_angle(a):
    """ Maintain angle between in range -pi <= a < pi """
    if a > pi:
        return a - 2*pi
    elif a < -pi:
        return a + 2*pi
    else:
        return a


def shortest_route(s, g, min_rad, min_straight, heading_tol, location_tol, spacing):
    """ An SCS path is an SCSCS path where the straight length is 0.
    If this invalidates our minimum straight length then we need to
    increase the bend radius to make the turn in an SCS. 
    """
    if abs(angle2intersect(s, g)) < heading_tol:
        d = dist(s, g)
        if d > min_straight:
            return straight_points(s, g, spacing)[0]
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

def path(s, g, bends=[], spacing=1):
    """ Construct a path from s to g, at increments of spacing, around
    each bend listed in bends. Bends are defined by (bs, bg, bc, d).
    """
    # if we just have a straight path
    if len(bends) == 0:
        return straight_points(s, g, spacing)
    else:
        # build up a path
        p, lo = straight_points(s, bends[0][0], spacing)
        for b in bends:
            pts, lo = bend_points(b, spacing, lo)
            p.extend(pts)
        p.extend(straight_points(s, bends[-1][1], spacing, lo)[0])
        return p

def tangent_points(sc, s_rad, st_dir, gc, g_rad, gt_dir):
    # distance between arc centres
    dc = dist(sc, gc)
    # if we're looking for an outer tangent
    if st_dir*gt_dir > 0:
        da = atan2(gc.y-sc.y, gc.x-sc.y)
        dt = dc
    else:
        da = atan2(gc.y-sc.y, gc.x-sc.y) + st_dir * asin((s_rad+g_rad)/dc)
        dt = (dc**2 - (s_rad+g_rad)**2)**0.5
    da = normalize_angle(da)
    # normal vector points from centre of sc to tp
    vn = normalize_angle(da - st_dir * pi / 2)
    # start point is s_rad along vn from sc, heading is da
    ts = Point(sc.x+s_rad*cos(vn), sc.y+s_rad*sin(vn), da)
    # end point is dt along da from start point, heading is da
    te = Point(ts.x+dt*cos(da), ts.y+dt*sin(da), da)

    return ts, te

def straight_points(s, g, spacing, leftover=0):
    if leftover != 0:
        ds = spacing - leftover
        s.x += ds * s.dx
        s.y += ds * s.dy
    d = dist(s, g)
    incs = int(d/spacing)
    lo = d - incs*spacing 
    xs = np.linspace(s.x, s.x+s.dx*incs, incs+1)
    ys = np.linspace(s.y, s.y+s.dy*incs, incs+1)
    heading = atan2(s.dy, s.dx)
    return [Point(x, y, heading=heading) for x, y in np.stack((xs, ys), axis=-1)], lo

def bend_points(bend, spacing, leftover=0):
    bs, bg, bc, d = bend[0], bend[1], bend[2], bend[3]
    r = dist(bs, bc)
    t = d * spacing / r
    ts = atan2(bs.y-bc.y, bs.x-bc.x)
    tg = atan2(bg.y-bc.y, bg.x-bc.x)
    if leftover != 0:
        ts = normalize_angle(ts + d * (spacing-leftover) / r)
    if ts < 0 and t < 0 and tg > 0:
        # we're going past -pi
        tg = -2*pi + tg
    elif ts > 0 and t > 0 and tg < 0:
        # we're going past pi
        tg = 2*pi + tg
    tc = ts
    return [], 0

if __name__ == '__main__':
    
    g = Point(0, 0, pi/2)
    s = Point(10, 10, 0)
    print(tangent_points(s, 5, -1, g, 5, -1))
    # min_rad = 2
    # min_straight = 3
    # shortest_route(s, g, min_rad, min_straight)


