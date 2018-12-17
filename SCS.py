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


def make_SCS(s, g, min_rad, min_straight):
    # criteria for a valid path
    ds = dist2intersect(s, g)
    # intersection point is after start and before goal
    if not (ds[0] >= 0 and ds[1] <= 0):
        return []
    # tangent points are at least a minimum straight length from the vector
    d2tans = dist2tangents(s, g, min_rad)
    if not (d2tans[0]>min_straight and d2tans[1]<-min_straight):
        return []
    # we can fit in a path....what's the max radius bend?
    # find it's centre...
    dt2i = np.min(np.abs(ds)) - min_straight
    print(dt2i)


if __name__ == '__main__':
    
    s = Point(0, 0, pi/2)
    g = Point(10, 10, 0)
    min_rad = 2
    min_straight = 3
    make_SCS(s, g, min_rad, min_straight)


