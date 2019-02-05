from math import pi
from heuristic import Point, tangent_points
import pytest

tangent_cases = [(Point(5, 0, 0, radius=5), Point(15, 15, 0, radius=-5),
                  Point(10, 5, pi/2), Point(10, 10, pi/2), 0, 1, 0, 1),
                 (Point(15, 15, pi, radius=5), Point(5, 0, pi, radius=-5),
                  Point(10, 10, -pi/2), Point(10, 5, -pi/2), 0, -1, 0, -1),
                 (Point(5, 0, 0, radius=5), Point(15, 11, 0, radius=-5),
                  Point(10, 5, pi/2), Point(10, 6, pi/2), 0, 1, 0, 1),
                 (Point(5, 0, 0, radius=5), Point(5, 11, pi, radius=5),
                  Point(10, 5, pi/2), Point(10, 6, pi/2), 0, 1, 0, 1)]

@pytest.mark.parametrize("ss, gg, sg, gs, sg_dx, sg_dy, gs_dx, gs_dy", tangent_cases)
def test_tangent_points(ss, gg, sg, gs, sg_dx, sg_dy, gs_dx, gs_dy):
    actual = tangent_points(ss, gg)
    assert actual[0].x == pytest.approx(sg.x)
    assert actual[0].y == pytest.approx(sg.y)
    assert actual[0].dx == pytest.approx(sg_dx)
    assert actual[0].dy == pytest.approx(sg_dy)
    assert actual[1].x == pytest.approx(gs.x)
    assert actual[1].y == pytest.approx(gs.y)
    assert actual[1].dx == pytest.approx(gs_dx)
    assert actual[1].dy == pytest.approx(gs_dy)
