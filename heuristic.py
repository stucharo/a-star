import numpy as np
from math import pi, sin, cos, asin, acos, atan, atan2


class Point:

    def __init__(self, x, y, heading, cost=0, radius=0, previous=None):
        self.x = x
        self.y = y
        self.dx = cos(heading)
        self.dy = sin(heading)
        self.heading = heading
        self.cost = cost
        self.radius = radius
        self.previous = previous
        
    def __eq__(self, p):
        return dist(self, p) < 1 and abs(self.heading - p.heading) < pi/180

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def dist(p1, p2):
    return ((p2.x-p1.x)**2 + (p2.y-p1.y)**2)**0.5

def dist2intersect(p1, p2):
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x-p1.x, p2.y-p1.y])
    try:
        return np.linalg.solve(dv, lv)
    except:
        return np.array([np.inf,np.inf])
    
def angle2intersect(p1, p2):
    return np.arccos(np.clip(np.dot([p1.dx, p1.dy], [p2.dx, p2.dy]), -1.0, 1.0))

def turn_dir(p1, p2):
    cp = np.cross([p1.dx, p1.dy], [p2.x-p1.x, p2.y-p1.y])
    if cp < 0:
        return -1
    else:
        return 1

def translate_pt(p, t):
    return copy_pt(p, t*cos(p.heading), t*sin(p.heading), 0)

def get_centre(p):
    rp = copy_pt(p, 0, 0, -1*direction(p)*pi/2)
    return copy_pt(rp, -1*abs(p.radius)*rp.dx, -1*abs(p.radius)*rp.dy, 0)

def copy_pt(p, dx, dy, dt):
    in_angle = atan2(p.dy, p.dx)
    new_heading = normalize_angle(in_angle + dt)
    return Point(p.x+dx, p.y+dy, new_heading, radius=p.radius)

def normalize_angle(a):
    """ Maintain angle between in range -pi <= a < pi """
    if a > pi:
        return a - 2*pi
    elif a < -pi:
        return a + 2*pi
    else:
        return a

def S_path(s, g, min_straight, heading_tol, location_tol, spacing=1):
    # are they colinear
    d = dist(s, g)
    if (abs(atan2(g.dy, g.dx) - atan2(s.dy, s.dx)) <= heading_tol
        and dist(copy_pt(s, d*s.dx, d*s.dy, 0), g) <= location_tol):
        # and far enough apart
        if s.radius == 0:
            # they are always at least the minimum straight length apart
            return path(s, g, spacing=spacing)
        else:
            # they have to be the minimum length apart
            if d > min_straight:
                return path(s, g, spacing=spacing)
            else:
                return []
    else:
        return []

def SCS_path(s, g, min_rad, min_straight, spacing):
    d2i = dist2intersect(s, g)
    a2i = angle2intersect(s, g)
    r = min_rad * (a2i/abs(a2i))
    t2i = min_rad * sin(a2i/2)
    bs = translate_pt(s, d2i[0]-t2i)
    bs.radius = r
    be = translate_pt(g, d2i[1]+t2i)
    be.radius = r
    bend = [(bs, be)]
    return path(s, g, bends=bend, spacing=spacing)



def SCSCS_path(s, g, min_rad, min_straight, spacing):
    # first, lets find out where we're going
    gbg = copy_pt(g, -min_straight*g.dx, -min_straight*g.dy, 0)
    # we're really interested in the relative locations of the first bend start and the last bend end
    if s.radius == 0:
        # we're not on a bend so we can start a new bend here
        sbs = copy_pt(s, 0, 0, 0)
    elif s.radius * turn_dir(s, g) > 0:
        # we're on a bend going the right way
        sbs = copy_pt(s, 0, 0, 0)
        sbs.radius = s.radius
    else:
        # we're on a bend going the wrong way.
        # to avoid crossing our old path we must go the minimum straight distance and start a bend there
        sbs = copy_pt(s, min_straight*s.dx, min_straight*s.dy, 0)
        sbs.radius = -1 * direction(s) * min_rad
    if sbs.radius == 0:
        # we need to set it to minimum bend radius but we must calculate the direction to
        # give it the correct sign
        sbs.radius = turn_dir(sbs, gbg) * min_rad
    gbg.radius = turn_dir(gbg, sbs) * min_rad
    t = tangent_points(sbs, gbg)
    if t is None:
        return []
    if turn_dir(gbg, t[0]) != direction(gbg):
        gbg.radius *= -1
        t = tangent_points(sbs, gbg)
    if t is None:
        return []
    if turn_dir(sbs, t[1]) != direction(sbs):
        gbg.radius *= -1
        t = tangent_points(sbs, gbg)
    if t is None:
        return []
    sbg, gbs = t
    if dist(sbg, gbs) > min_straight or (dist(sbg, gbs) == 0 and direction(sbg) == direction(gbs)):
        bends = [(sbs, sbg), (gbs, gbg)]
        return path(s, g, bends=bends, spacing=spacing)
    else:
        return []

def shortest_route(s, g, min_rad, min_straight, heading_tol, location_tol, spacing):
    """ Route options:
        1. We can make it in a straight line.
        2. Stright - curve - straight - curve - striahgt cuts the 
           corner, but only works if we can fit the minimum straight
           lengths.
        3. Straight - curve - straight works if we we can fit the 
           minimum straight lengths
        4. 
    """
    p = S_path(s, g, min_straight, heading_tol, location_tol, spacing)
    if len(p) > 0:
        return p
    scscs_p = SCSCS_path(s, g, min_rad, min_straight, spacing)
    scs_p = SCS_path(s, g, min_rad, min_straight, spacing)
    if len(scscs_p) < 0 and len(scs_p) > 0:
        return scs_p
    elif len(scs_p) < 0 and len(scscs_p) > 0:
        return scscs_p
    elif len(scs_p) > 0 and len(scscs_p) > 0:
        if len(scs_p) < len(scscs_p):
            return scs_p
        else:
            return scscs_p
    else:
        return []

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

def tangent_points(s, g):
    sc = get_centre(s)
    gc = get_centre(g)
    # distance between arc centres
    dc = dist(sc, gc)
    ds = direction(s)
    dg = direction(g)
    if ds * dg > 0:
        # if the turn directions are both the same then we're looking for an outer tangent
        # so the heading of the tangent is parallel to the line between the centre of 1
        # circle and it's tangent on a circle centred on the second circle, with a radius
        # equal to the difference in radii. This defaults to the centre of the second
        # circle when they are equal radius
        if abs(abs(s.radius)-abs(g.radius)) > dc:
            return None
        da = atan2(gc.y-sc.y, gc.x-sc.x) + ds * asin(abs(abs(s.radius)-abs(g.radius))/dc)
        # and the distance
        dt = (dc**2 - abs(abs(s.radius)-abs(g.radius))**2)**0.5
    else:
        # get the angle of the tangent line
        if -1 <= (abs(s.radius)+abs(g.radius))/dc <= 1:
            da = atan2(gc.y-sc.y, gc.x-sc.x) + ds * asin((abs(s.radius)+abs(g.radius))/dc)
            # and the distance
            dt = (dc**2 - (abs(s.radius)+abs(g.radius))**2)**0.5
        else:
            # we can't find an inner tangent for two circles that overlap
            return None
    da = normalize_angle(da)
    # normal vector points from centre of sc to tp
    vn = normalize_angle(da - ds * pi / 2)
    # start point is s_rad along vn from sc, heading is da
    ts = Point(sc.x+abs(s.radius)*cos(vn), sc.y+abs(s.radius)*sin(vn), da)
    ts.radius = s.radius
    # end point is dt along da from start point, heading is da
    te = Point(ts.x+dt*cos(da), ts.y+dt*sin(da), da)
    te.radius = g.radius

    return ts, te

def direction(p):
    return p.radius / abs(p.radius)

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

def bend_point(start, spacing=1, radius=0):
    s = copy_pt(start,0,0,0)
    s.radius = radius
    bc = get_centre(s)
    dt = spacing / radius
    new_x = bc.x+abs(radius)*cos(bc.heading+dt)
    new_y = bc.y+abs(radius)*sin(bc.heading+dt)
    new_t = s.heading+dt
    return Point(new_x, new_y, new_t, radius=radius)


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

def develop_path(s, g, min_straight_length, min_bend_radius):
    # Begin with a straight line between start and goal
    # Check it would arrive in the correct location at the correct heading
    d = dist(s, g)
    if dist(translate_pt(s, d), g) < 1 and abs(s.heading-g.heading) < pi/180:
        # check the path length is valid
        # if s is on a bend then the length must be greater than min straight length
        if s.radius == 0 or d > min_straight_length:
            return [s, g]
    # if we get this far the straight line won't work. we need to bend the straight line
    # to the start and/or goal headings and join them with a tangent. 
    # it might seem obvious which way we have to turn, but because we move the goal by
    # adding the curve, it can subtly change the optimum direction of the start turn and
    # vice versa.
    # lets start by creating a pair of circles at the minimum straight length from the
    # goal sharing a tangent coincident to the goal heading....
    goal_turns = [translate_pt(Point(g.x, g.y, g.heading, radius=min_bend_radius),
                               -min_straight_length),
                  translate_pt(Point(g.x, g.y, g.heading, radius=-min_bend_radius),
                               -min_straight_length)]
    # the start radius is dependent on the start point.
    # if it is straight then we have 2 options from that point.
    # if it is curved then we can continue that curve or have 2 minimum straight bends a
    # minimum straight distance away
    if s.radius == 0:
        start_turns = [Point(s.x, s.y, s.heading, radius=min_bend_radius),
                       Point(s.x, s.y, s.heading, radius=-min_bend_radius)]
    else:
        start_turns = [copy_pt(s, 0, 0, 0)]
        # if the minimum straight requirement means that the two circles pass each other
        # then they cannot work as the path is crossing itself. This is only really a
        # concern with the start length, as there is always an option with no straight
        # length.
        if dist2intersect(s, g)[0] > min_straight_length:
            start_turns.extend([translate_pt(Point(s.x, s.y, s.heading, radius=min_bend_radius),
                                             min_straight_length),
                                translate_pt(Point(s.x, s.y, s.heading, radius=-min_bend_radius),
                                             min_straight_length)])
    # lets evaluate all of the tangents between start and goal arcs
    tangents = []
    for st in start_turns:
        for gt in goal_turns:
            tp = tangent_points(st, gt)
            if tp is not None:
                stt, gtt = tp
                # the tangent is valid if the turn direction of st is equal to the turn
                # direction between st and gtt it st is translated to the location of stt
                # and vice versa. we also need to check the tangent is longer than the
                # minimum straight length.
                if ((turn_dir(copy_pt(st, stt.x-st.x, stt.y-st.y, 0), gtt) * st.radius > 0)
                     and (turn_dir(copy_pt(gt, gtt.x-gt.x, gtt.y-gt.y, 0), stt) * gt.radius > 0)
                     and dist(stt, gtt) > min_straight_length):
                    tangents.append((st, stt, gtt, gt))
    if len(tangents) > 0:
        # we must have at least 1 turn along each valid path so lets format the path
        return [(s,)+t+(g,) for t in tangents]
    return []

def plot_path():
    import matplotlib.pyplot as plt
    s = Point(0, 10, pi/4)
    g = Point(0, 90, -pi/2)
    msl = 10
    mbr = 10

    paths = develop_path(s, g, msl, mbr)

    plt.arrow(s.x, s.y, 5*s.dx, 5*s.dy, head_width=0.5, color='green')
    plt.arrow(g.x, g.y, 5*g.dx, 5*g.dy, head_width=0.5, color='red')
    
    if len(paths) > 0:
        for p in paths:
            pts = np.asarray([(pt.x, pt.y) for pt in p])
            plt.plot(pts[:,0], pts[:,1])
    plt.axis('equal')
    plt.show()



if __name__ == '__main__':
    
    plot_path()

    """
    import matplotlib.pyplot as plt
    import random
    actual = False
    if actual:
        min_rad = 12
        rads = [-15, -14, -13, -12, 0, 12, 13, 14, 15]
        start_rad = 0
        min_straight = 10
        sx = 24.740583304172137
        sy = 50.36969667507735
        sh = 2.840395189256073
        gx = 13.226612622299294
        gy = 58.99167386202944
        gh = -2.771550305600809
    else:
        max_rad = 15
        min_rad = random.randint(2, max_rad)
        rads = []
        rads.extend(range(-max_rad, -min_rad+1))
        rads.extend([0])
        rads.extend(range(min_rad, max_rad+1))
        start_rad = random.choice(rads)
        sx = random.uniform(0, 100)
        sy = random.uniform(0, 100)
        sh = random.uniform(-pi, pi)
        gx = random.uniform(0, 100)
        gy = random.uniform(0, 100)
        gh = random.uniform(-pi, pi)
        min_straight = random.randint(0, 10)
    s = Point(sx, sy, sh, radius=start_rad)
    g = Point(gx, gy, gh)
    heading_tol = pi/180
    location_tol = 1
    spacing = 1
    print(f"min_rad = {min_rad}\nrads = {rads}\nstart_rad = {start_rad}\nmin_straight = {min_straight}")
    print(f"sx = {sx}\nsy = {sy}\nsh = {sh}\ngx = {gx}\ngy = {gy}\ngh = {gh}")
    p = shortest_route(s, g, min_rad, min_straight, heading_tol, location_tol, spacing)
    plt.arrow(s.x, s.y, 5*s.dx, 5*s.dy, head_width=0.5, color='green')
    plt.arrow(g.x, g.y, 5*g.dx, 5*g.dy, head_width=0.5, color='red')
    if len(p) > 0:
        pts = np.asarray([(pt.x, pt.y) for pt in p])
        plt.plot(pts[:,0], pts[:,1])
        plt.scatter(pts[:,0], pts[:,1])
    plt.axis('equal')
    plt.show()
"""