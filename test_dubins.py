import pytest
from a_star import Point
from dubins import Circle, Tangent, tangent, make_circle, S_route, dist
from math import pi, sin, cos, acos

def test_S_route():
    start = Point(0,0, heading=pi/4)
    goal = Point(10,10, heading=pi/4)
    min_straight_len = 10
    heading_tol = pi/180
    spacing = 1
    actual = S_route(start, goal, min_straight_len, heading_tol, spacing)
    assert len(actual) == 15
    assert actual[0].x == 0
    assert actual[0].y == 0
    assert actual[0].heading == pi/4
    assert actual[-1].x == pytest.approx(9.899, abs=1e-3)
    assert actual[-1].y == pytest.approx(9.899, abs=1e-3)
    assert actual[-1].heading == pi/4

def test_S_route_bad_heading():
    start = Point(0,0, heading=pi/4)
    goal = Point(10,10, heading=0)
    min_straight_len = 10
    heading_tol = pi/180
    spacing = 1
    actual = S_route(start, goal, min_straight_len, heading_tol, spacing)
    assert len(actual) == 0

def test_S_route_too_short():
    start = Point(0,0, heading=pi/4, radius=5)
    goal = Point(10,10, heading=pi/4)
    min_straight_len = 100
    heading_tol = pi/180
    spacing = 1
    actual = S_route(start, goal, min_straight_len, heading_tol, spacing)
    assert len(actual) == 0
    



tangent_cases = [# RR tangents   
                 (Circle(Point(0,0), 2, 'r'), Circle(Point(0,10), 2, 'r'),
                  Tangent(Point(-2,0,heading=0), Point(-2,10,heading=0))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,10), 2, 'r'),
                  Tangent(Point(-2*sin(pi/4),2*cos(pi/4),heading=pi/4), Point(10-2*sin(pi/4),10+2*cos(pi/4),heading=pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,0), 2, 'r'),
                  Tangent(Point(0,2,heading=pi/2), Point(10,2,heading=pi/2))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,-10), 2, 'r'),
                 Tangent(Point(2*sin(pi/4),2*cos(pi/4),heading=3*pi/4), Point(10+2*cos(pi/4),-10+2*sin(pi/4),heading=3*pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(0,-10), 2, 'r'),
                 Tangent(Point(2,0,heading=pi), Point(2,-10,heading=pi))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,-10), 2, 'r'),
                  Tangent(Point(2*sin(pi/4),-2*sin(pi/4),heading=5*pi/4), Point(-10+2*sin(pi/4),-10-2*sin(pi/4),heading=5*pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,0), 2, 'r'),
                  Tangent(Point(0,-2,heading=3*pi/2), Point(-10,-2,heading=3*pi/2))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,10), 2, 'r'),
                  Tangent(Point(-2*sin(pi/4),-2*cos(pi/4),heading=7*pi/4), Point(-10-2*sin(pi/4),10-2*cos(pi/4),heading=7*pi/4))),
                 
                 # LL tangents
                 (Circle(Point(0,0), 2, 'l'), Circle(Point(0,10), 2, 'l'),
                  Tangent(Point(2,0,heading=0), Point(2,10,heading=0))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,10), 2, 'l'),
                  Tangent(Point(2*sin(pi/4),-2*cos(pi/4),heading=pi/4), Point(10+2*sin(pi/4),10-2*cos(pi/4),heading=pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,0), 2, 'l'),
                  Tangent(Point(0,-2,heading=pi/2), Point(10,-2,heading=pi/2))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,-10), 2, 'l'),
                 Tangent(Point(-2*sin(pi/4),-2*cos(pi/4),heading=3*pi/4), Point(10-2*sin(pi/4),-10-2*cos(pi/4),heading=3*pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(0,-10), 2, 'l'),
                 Tangent(Point(-2,0,heading=pi), Point(-2,-10,heading=pi))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,-10), 2, 'l'),
                  Tangent(Point(-2*sin(pi/4),2*sin(pi/4),heading=5*pi/4), Point(-10-2*sin(pi/4),-10+2*sin(pi/4),heading=5*pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,0), 2, 'l'),
                 Tangent(Point(0,2,heading=3*pi/2), Point(-10,2,heading=3*pi/2))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,10), 2, 'l'),
                 Tangent(Point(2*sin(pi/4),2*cos(pi/4),heading=7*pi/4), Point(-10+2*sin(pi/4),10+2*cos(pi/4),heading=7*pi/4)))]

@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p1x(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p1.x == pytest.approx(expected.p1.x, abs=1e-6)

@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p1y(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p1.y == pytest.approx(expected.p1.y, abs=1e-6)

@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p1_heading(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p1.heading == pytest.approx(expected.p1.heading, abs=1e-6)
    
@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p2x(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p2.x == pytest.approx(expected.p2.x, abs=1e-6)
    
@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p2y(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p2.y == pytest.approx(expected.p2.y, abs=1e-6)
    
@pytest.mark.parametrize("c1, c2, expected", tangent_cases)
def test_tangent_p2_heading(c1, c2, expected):
    t = tangent(c1, c2)
    assert t.p2.heading == pytest.approx(expected.p2.heading, abs=1e-6)

test_circles = [(0, 0, 0, 4, 'r', 4, 0),
                (0, 0, 0, 4, 'l', -4, 0),
                (0, 0, pi/4, 4, 'r', 4*sin(pi/4), -4*sin(pi/4)),
                (0, 0, pi/4, 4, 'l', -4*sin(pi/4), 4*sin(pi/4)),
                (0, 0, pi/2, 4, 'r', 0, -4),
                (0, 0, pi/2, 4, 'l', 0, 4),
                (0, 0, 3*pi/4, 4, 'r', -4*sin(pi/4), -4*sin(pi/4)),
                (0, 0, 3*pi/4, 4, 'l', 4*sin(pi/4), 4*sin(pi/4)),
                (0, 0, pi, 4, 'r', -4, 0),
                (0, 0, pi, 4, 'l', 4, 0),
                (0, 0, 5*pi/4, 4, 'r', -4*sin(pi/4), 4*sin(pi/4)),
                (0, 0, 5*pi/4, 4, 'l', 4*sin(pi/4), -4*sin(pi/4)),
                (0, 0, 3*pi/2, 4, 'r', 0, 4),
                (0, 0, 3*pi/2, 4, 'l', 0, -4),
                (0, 0, 7*pi/4, 4, 'r', 4*sin(pi/4), 4*sin(pi/4)),
                (0, 0, 7*pi/4, 4, 'l', -4*sin(pi/4), -4*sin(pi/4))]

@pytest.mark.parametrize("x, y, heading, radius, direction, ex, ey", test_circles)
def test_make_circle_x(x, y, heading, radius, direction, ex, ey):
    a = Point(x, y, heading=heading, radius=radius, direction=direction)
    c = make_circle(a)
    assert c.centre.x == pytest.approx(ex)

@pytest.mark.parametrize("x, y, heading, radius, direction, ex, ey", test_circles)
def test_make_circle_y(x, y, heading, radius, direction, ex, ey):
    a = Point(x, y, heading=heading, radius=radius, direction=direction)
    c = make_circle(a)
    assert c.centre.y == pytest.approx(ey)

