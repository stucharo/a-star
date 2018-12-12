import pytest
from dubins import Circle, Point, Tangent, LinePoint, ArcPoint, tangent, make_circle
from math import pi, sin, cos, acos


tangent_cases = [# RR tangents   
                 (Circle(Point(0,0), 2, 'r'), Circle(Point(0,10), 2, 'r'),
                  Tangent(LinePoint(-2,0,0), LinePoint(-2,10,0))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,10), 2, 'r'),
                  Tangent(LinePoint(-2*sin(pi/4),2*cos(pi/4),pi/4), LinePoint(10-2*sin(pi/4),10+2*cos(pi/4),pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,0), 2, 'r'),
                  Tangent(LinePoint(0,2,pi/2), LinePoint(10,2,pi/2))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(10,-10), 2, 'r'),
                 Tangent(LinePoint(2*sin(pi/4),2*cos(pi/4),3*pi/4), LinePoint(10+2*cos(pi/4),-10+2*sin(pi/4),3*pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(0,-10), 2, 'r'),
                 Tangent(LinePoint(2,0,pi), LinePoint(2,-10,pi))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,-10), 2, 'r'),
                  Tangent(LinePoint(2*sin(pi/4),-2*sin(pi/4),5*pi/4), LinePoint(-10+2*sin(pi/4),-10-2*sin(pi/4),5*pi/4))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,0), 2, 'r'),
                  Tangent(LinePoint(0,-2,3*pi/2), LinePoint(-10,-2,3*pi/2))),

                 (Circle(Point(0,0), 2, 'r'), Circle(Point(-10,10), 2, 'r'),
                  Tangent(LinePoint(-2*sin(pi/4),-2*cos(pi/4),7*pi/4), LinePoint(-10-2*sin(pi/4),10-2*cos(pi/4),7*pi/4))),
                 
                 # LL tangents
                 (Circle(Point(0,0), 2, 'l'), Circle(Point(0,10), 2, 'l'),
                  Tangent(LinePoint(2,0,0), LinePoint(2,10,0))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,10), 2, 'l'),
                  Tangent(LinePoint(2*sin(pi/4),-2*cos(pi/4),pi/4), LinePoint(10+2*sin(pi/4),10-2*cos(pi/4),pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,0), 2, 'l'),
                  Tangent(LinePoint(0,-2,pi/2), LinePoint(10,-2,pi/2))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(10,-10), 2, 'l'),
                 Tangent(LinePoint(-2*sin(pi/4),-2*cos(pi/4),3*pi/4), LinePoint(10-2*sin(pi/4),-10-2*cos(pi/4),3*pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(0,-10), 2, 'l'),
                 Tangent(LinePoint(-2,0,pi), LinePoint(-2,-10,pi))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,-10), 2, 'l'),
                  Tangent(LinePoint(-2*sin(pi/4),2*sin(pi/4),5*pi/4), LinePoint(-10-2*sin(pi/4),-10+2*sin(pi/4),5*pi/4))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,0), 2, 'l'),
                 Tangent(LinePoint(0,2,3*pi/2), LinePoint(-10,2,3*pi/2))),

                 (Circle(Point(0,0), 2, 'l'), Circle(Point(-10,10), 2, 'l'),
                 Tangent(LinePoint(2*sin(pi/4),2*cos(pi/4),7*pi/4), LinePoint(-10+2*sin(pi/4),10+2*cos(pi/4),7*pi/4)))]

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
                (0, 0, pi/4, 4, 'r', 4, 0),
                (0, 0, pi/4, 4, 'l', -4, 0)]

@pytest.mark.parametrize("x, y, heading, radius, direction, ex, ey", test_circles)
def test_make_circle_x(x, y, heading, radius, direction, ex, ey):
    a = ArcPoint(x, y, heading, radius, direction)
    c = make_circle(a)
    assert c.centre.x == pytest.approx(ex)

@pytest.mark.parametrize("x, y, heading, radius, direction, ex, ey", test_circles)
def test_make_circle_y(x, y, heading, radius, direction, ex, ey):
    a = ArcPoint(x, y, heading, radius, direction)
    c = make_circle(a)
    assert c.centre.y == pytest.approx(ey)