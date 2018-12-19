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
    rp1 = move(p1, 0, 0, t_dir*pi/2)
    rp2 = move(p2, 0, 0, t_dir*pi/2)
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x+r*rp2.dx-(p1.x+r*rp1.dx), p2.y+r*rp2.dy-(p1.y+r*rp1.dy)])
    return np.linalg.solve(dv, lv)

def move(p, dx, dy, dt):
    in_angle = atan2(p.dy, p.dx)
    new_heading = normalize_angle(in_angle + dt)
    return Point(p.x+dx, p.y+dy, new_heading)

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
    sc = move(st, 0, 0, st_dir*2/pi)
    sc.x += min_rad*sc.dx
    sc.y += min_rad*sc.dy
    gc = move(gt, 0, 0, gt_dir*pi/2)
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
    p = [s]
    # if we just have a straight path
    if len(bends) == 0:
        return p.extend(straight_points(s, g, spacing)[0])
    else:
        # build up a path
        # first straight
        sb = bends[0][0]
        if s.x != sb.x or s.y != sb.y:
            pts, lo = straight_points(s, sb, spacing)
            p.extend(pts)
        # loop through each bend
        for n in range(len(bends) - 1):
            b = bends[n]
            nb = bends[n+1]
            # arc
            pts, lo = bend_points(b, spacing, lo)
            p.extend(pts)
            # then straight
            if b[1].x != nb[0].x or b[1].y != nb[0].y :
                pts, lo = straight_points(b[1], nb[0], spacing, lo)
                p.extend(pts)
        # last bend
        b = bends[-1]
        pts, lo = bend_points(b, spacing, lo)
        p.extend(pts)
        # last straight
        bg = b[1]
        if (g.x != bg.x or g.y != bg.y 
            or g.dx != bg.dx or g.dy != bg.dy):
            p.extend(straight_points(bg, g, spacing, lo)[0])
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
    ds = spacing - leftover
    sp = Point(s.x + ds*s.dx, s.y + ds*s.dy, atan2(s.dy, s.dx))
    d = dist(sp, g)
    incs = int(d/spacing)
    lo = d - incs*spacing 
    xs = np.linspace(sp.x, sp.x+sp.dx*incs, incs+1)
    ys = np.linspace(sp.y, sp.y+sp.dy*incs, incs+1)
    heading = atan2(sp.dy, sp.dx)
    return [Point(x, y, heading=heading) for x, y in np.stack((xs, ys), axis=-1)], lo

def bend_points(bend, spacing, leftover=0):
    bs, bg, bc, d = bend[0], bend[1], bend[2], bend[3]
    r = dist(bs, bc)
    # theta start is normal to bs
    vs = move(bs, bc.x-bs.x, bc.y-bs.y, -d*pi/2)
    vg = move(bg, bc.x-bg.x, bc.y-bg.y, -d*pi/2)
    # move it for leftover
    dt = d * (spacing-leftover) / r
    vs = move(vs, 0, 0, dt)
    # get number of whole increments
    ts = atan2(vs.dy, vs.dx)
    tg = atan2(vg.dy, vg.dx)
    tt = angle_between(vs, vg, d)
    dt = d * spacing / r
    inc = int(tt / dt)
    tl = np.linspace(ts, ts + inc*dt, inc+1)
    p = [Point(bc.x+r*cos(t), bc.y+r*sin(t), t + d*pi/2) for t in tl]
    leftover = abs(abs(tg)-abs(tl[-1])) * r
    return p, leftover

def angle_between(s, g, d):
    # dot product between [x1, y1] and [x2, y2]
    dot = s.dx*g.dx + s.dy*g.dy   
    # determinant   
    det = s.dx*g.dy - s.dy*g.dx   
    # anti-clockwise angle
    angle = atan2(det, dot)
    if d*angle < 0:
        if angle >= 0:
            angle -= 2 * pi
        else:
            angle += 2 * pi
    return angle

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    s = Point(0, -5, pi/2)
    g = Point(15, 15, pi/2)
    bends = [(Point(0,0,pi/2), Point(5,5,0), Point(5,0,0), -1),
             (Point(10, 5, 0), Point(15, 10, pi/2), Point(10, 10, 0), 1)]
    path = path(s, g, bends, 1)
    pts = np.array([[p.x, p.y] for p in path])
    plt.plot(pts[:,0], pts[:,1])
    plt.scatter(pts[:,0], pts[:,1])
    plt.axis('equal')
    plt.show()