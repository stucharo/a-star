import numpy as np
from scipy.optimize import minimize
from math import pi, sin, cos, asin, acos, atan, atan2, tan
import matplotlib.pyplot as plt

class Point():

    def __init__(self, x, y, heading=0, radius=0):
        self.x = x
        self.y = y
        self.dx = normalize_angle(cos(heading))
        self.dy = normalize_angle(sin(heading))
        self.radius = radius

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def normalize_angle(a):
    """ Maintain angle between in range -pi <= a < pi """
    if a > pi:
        return a - 2*pi
    elif a < -pi:
        return a + 2*pi
    else:
        return a

def dist2intersect(p1, p2):
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x-p1.x, p2.y-p1.y])
    try:
        return np.linalg.solve(dv, lv)
    except:
        return np.array([np.inf,np.inf])
 
def dist(p1, p2):
    return ((p2.x-p1.x)**2+(p2.y-p1.y)**2)**0.5

def angle(p1, p2, p3):
    """ Angle between p1 -> p2 -> p3
    Except clause to catch case when all 3 points are on the same line.
    This returns pi e.g. angle of 180 degress
    """
    try:
        return acos((dist(p2,p3)**2 + dist(p1,p2)**2 - dist(p1,p3)**2)
                    / (2 * dist(p2,p3) * dist(p1,p2)))
    except:
        return pi

def v2t(p1, p2, p3, r):
    """ Distance from vertex `p2` to the tangent of circle radius `r`
    through line p1 -> p2 -> p3.
    """
    a = angle(p1, p2, p3)
    if a == pi:
        return 0
    else:
        return r / tan(a / 2)

def arc_length(p1, p2, p3, r):
    a = angle(p1, p2, p3)
    if a == pi:
        return 0
    else:
        return a * r

def turn_dir(p1, p2):
    cp = np.cross([p1.dx, p1.dy], [p2.x-p1.x, p2.y-p1.y])
    if cp < 0:
        return -1
    else:
        return 1

def copy_pt(p, dx, dy, dt):
    in_angle = atan2(p.dy, p.dx)
    new_heading = normalize_angle(in_angle + dt)
    return Point(p.x+dx, p.y+dy, new_heading, radius=p.radius)

def get_shortest_path(s, g, min_bend, min_straight, min_straight_end):

    min_length = np.inf
    bs = []

    d2i = dist2intersect(s, g)
    if 0 < d2i[0] < 10*min_straight and 0 > d2i[1] > -10*min_straight_end:
        # the two lines cross and a triangle is possible
        # lets get a vertex at the intersection point
        ip = Point(s.x+d2i[0]*s.dx, s.y+d2i[0]*s.dy)
        # and the distance between the ip and tp for a minimum radius bend
        d2t = v2t(s, ip, g, min_bend)
        # the minimum required distance depends on whether the start is on a bend
        if s.radius == 0:
            min_d = 0
        else:
            min_d = min_straight
        if d2i[0] > min_d + d2t and -d2i[1] > min_straight_end + d2t: 
            # if we're here we can make a path with at least 1 bend
            # it's length is...
            min_length = d2i[0]-d2t + arc_length(s,ip,g,min_bend) + -1*d2i[1] + d2t
            # and it's tangent points are
            bs = generate_bends([s, ip, g], min_bend)
            
    if s.radius != 0:
        d, bends = get_length_and_bends(s, g, min_bend, min_straight, 'ineq', min_straight, min_straight_end)
        if d < min_length:
            min_length = d
            bs = bends
            
        d, bends = get_length_and_bends(s, g, min_bend, min_straight, 'eq', 0, min_straight_end, s.radius)
        if d < min_length:
            min_length = d
            bs = bends
    else:
        d, bends = get_length_and_bends(s, g, min_bend, min_straight, 'ineq', 0, min_straight_end)
        if d < min_length:
            min_length = d
            bs = bends

    return path(s, g, bs, 1)

def get_length_and_bends(s, g, min_bend, min_straight, con, start_straight, min_straight_end, start_radius=None):

    v = [min_straight, min_straight]

    if start_radius is not None:
        s_bend = abs(start_radius)
    else:
        s_bend = min_bend

    cons = [{'type': con, 'fun': con_len_ab, 'args': (s, g, s_bend, start_straight)},
            {'type': 'ineq', 'fun': con_len_bc, 'args': (s, g, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_len_cd, 'args': (s, g, min_bend, min_straight_end)}]
        
    if start_radius is not None:
        cons.extend([{'type': 'ineq', 'fun': con_turn_dir_b, 'args': (s, g, start_radius)}])
    
    d = minimize(length, v, args=(s,g,min_bend, min_straight), constraints=cons)

    if d.success:
        pts = [s, g]
        pts[1:1] = [copy_pt(s,s.dx*d.x[0],s.dy*d.x[0],0),
                   copy_pt(g,g.dx*-d.x[1],g.dy*-d.x[1],0)]
        return d.fun, generate_bends(pts, min_bend, s_bend)
    else:
        return np.inf, []

def con_len_ab(vars, s, g, r, min_straight):
    a, b, c, _ = unpack(vars, s, g)
    t = v2t(a,b,c,r)
    return vars[0] - t - min_straight

def con_len_bc(vars, s, g, r, min_straight):
    a, b, c, d = unpack(vars, s, g)
    tb = v2t(a,b,c,r)
    tc = v2t(b,c,d,r)
    return dist(b,c) - tb - tc - min_straight

def con_len_cd(vars, s, g, r, min_straight_end):
    _, b, c, d = unpack(vars, s, g)
    t = v2t(b,c,d,r)
    return vars[1] - t - min_straight_end

def con_turn_dir_b(vars, s, g, start_radius):
    a, b, c, _ = unpack(vars, s, g)
    return turn_dir(a, c) * start_radius

def length(vars, s, g, min_bend, min_straight):
    # vars is our minimization variables. these are:
    # vars[0] -> ts
    # vars[1] -> tg
    
    a, b, c, d = unpack(vars, s, g)

    # get distance from vertex to tangent
    tb = v2t(a,b,c,min_bend)
    tc = v2t(b,c,d,min_bend)

    lab = arc_length(a,b,c,min_bend)
    lac = arc_length(b,c,d,min_bend) 

    return dist(a,b) - tb + lab + dist(b,c) - tb - tc + lac + dist(c,d) - tc

def unpack(vars, s, g):
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    return (a, b, c, d)

def graph(path):

    if len(path) > 1:
        ps = np.array([(p.x, p.y) for p in path])
        plt.plot(ps[:,0], ps[:,1])
        plt.axis('equal')
    else:
        cx = (s.x+g.x)/2
        cy = (s.y+g.y)/2
        w = max(abs(g.x-s.x), abs(g.y-s.y))
        plt.axis([cx-w/2-10,cx+w/2+10,cy-w/2-10,cy+w/2+10])
    plt.arrow(s.x, s.y, 5*s.dx, 5*s.dy, head_width=0.5, color='green')
    plt.arrow(g.x, g.y, 5*g.dx, 5*g.dy, head_width=0.5, color='red')
    plt.show()

def generate_bends(pts, min_bend, start_bend=None):
    bends = []
    for i in range(1,len(pts)-1):
        if i == 1:
            if start_bend is None:
                r = min_bend
            else:
                r = start_bend
        else:
            r = min_bend

        t = v2t(pts[i-1], pts[i], pts[i+1],r)
        sp = pts[i-1]
        d1 = dist(pts[i-1], pts[i])
        d2 = dist(pts[i], pts[i+1])

        bs = copy_pt(sp, sp.dx*(d1-t), sp.dy*(d1-t), 0)
        bs = Point(pts[i-1].x, pts[i-1].y)
        bs.dx = (pts[i].x - pts[i-1].x) / d1
        bs.dy = (pts[i].y - pts[i-1].y) / d1
        bs.x += bs.dx * (d1 - t)
        bs.y += bs.dy * (d1 - t)
        d = turn_dir(bs, pts[i+1])
        bs.radius = d*r
        bg = Point(pts[i].x, pts[i].y)
        bg.dx = (pts[i+1].x - pts[i].x) / d2
        bg.dy = (pts[i+1].y - pts[i].y) / d2
        bg.x += bg.dx * t
        bg.y += bg.dy * t
        bg.radius = d*r

        bends.append((bs, bg))
    return bends

def path(s, g, bends=[], spacing=1):
    """ Construct a path from s to g, at increments of spacing, around
    each bend listed in bends. Bends are defined by (bs, bg).
    """
    p = [s]
    # if we just have a straight path
    if len(bends) == 0:
        p.extend(straight_points(s, g, spacing)[0])
        return p
    else:
        # build up a path
        # first straight
        sb = bends[0][0]
        if s.x != sb.x or s.y != sb.y:
            pts, lo = straight_points(s, sb, spacing)
            p.extend(pts)
        else:
            lo = 0
        # loop through each bend
        for n in range(len(bends) - 1):
            b = bends[n]
            nb = bends[n+1]
            # arc
            pts, lo = bend_points(b[0], b[1], spacing, lo)
            p.extend(pts)
            # then straight
            if b[1].x != nb[0].x or b[1].y != nb[0].y :
                pts, lo = straight_points(b[1], nb[0], spacing, lo)
                p.extend(pts)
        # last bend
        b = bends[-1]
        pts, lo = bend_points(b[0], b[1], spacing, lo)
        p.extend(pts)
        # last straight
        bg = b[1]
        if (g.x != bg.x or g.y != bg.y 
            or g.dx != bg.dx or g.dy != bg.dy):
            p.extend(straight_points(bg, g, spacing, lo)[0])
        return p

def straight_points(s, g, spacing, leftover=0):
    """ Generate a list of points along a straight line between the `s`
    and `g` Points. Note that this assumes the previous section 
    included the first point of this list, so this list returns the
    next point along the path. That is, this routine will never return
    a point at `s`.
    """
    ds = spacing - leftover
    sp = Point(s.x + ds*s.dx, s.y + ds*s.dy, atan2(s.dy, s.dx))
    d = dist(sp, g)
    incs = int(d/spacing)
    lo = d - incs*spacing 
    xs = np.linspace(sp.x, sp.x+sp.dx*incs, incs+1)
    ys = np.linspace(sp.y, sp.y+sp.dy*incs, incs+1)
    heading = atan2(sp.dy, sp.dx)
    return [Point(x, y, heading=heading) for x, y in np.stack((xs, ys), axis=-1)], lo

def bend_points(bs, bg, spacing, leftover=0):
    bc = get_centre(bs)
    r = abs(bs.radius)
    d = direction(bs)
    # theta start is normal to bs
    vs = copy_pt(bs, bc.x-bs.x, bc.y-bs.y, -d*pi/2)
    vg = copy_pt(bg, bc.x-bg.x, bc.y-bg.y, -d*pi/2)
    # copy and move it for leftover
    dt = d * (spacing-leftover) / r
    vs = copy_pt(vs, 0, 0, dt)
    # get number of whole increments
    ts = atan2(vs.dy, vs.dx)
    tg = atan2(vg.dy, vg.dx)
    tt = angle_between(vs, vg, d)
    dt = d * spacing / r
    inc = int(tt / dt)
    tl = np.linspace(ts, ts + inc*dt, inc+1)
    p = [Point(bc.x+r*cos(t), bc.y+r*sin(t), t + d*pi/2) for t in tl]
    leftover = abs(abs(tg)-abs(normalize_angle(tl[-1]))) * r
    return p, leftover

def get_centre(p):
    rp = copy_pt(p, 0, 0, -1*direction(p)*pi/2)
    return copy_pt(rp, -1*abs(p.radius)*rp.dx, -1*abs(p.radius)*rp.dy, 0)

def direction(p):
    return p.radius / abs(p.radius)

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
    s = Point(0,0,5*pi/12,radius=-5)
    g = Point(100, 0, -5*pi/12)

    graph(get_shortest_path(s, g, 10, 20, 30))