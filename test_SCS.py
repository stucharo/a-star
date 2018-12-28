from math import pi
from SCS import Point, tangent_points

def test_tangent_points():
    s = Point(5, 0, 0, radius=-5)
    g = Point(15, 15, 0, radius=5)
    actual = tangent_points(s, g)
    expected = (Point(10, 5, pi/2), Point(10, 10, pi/2))
    assert actual[0].x == expected[0].x
    assert actual[0].y == expected[0].y
    assert actual[0].heading == expected[0].heading
    assert actual[1].x == expected[1].x
    assert actual[1].y == expected[1].y
    assert actual[1].heading == expected[1].heading
