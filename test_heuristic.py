import pytest
import heuristic as h
from math import pi, sin, cos, acos



"""
Testing straight_points function

This function generates a list of Points along a straight line between
the start and the goal with a predefined space. It should return a
tuple containing a list of points and the leftover distance to reach
the goal. 
"""
def test_straight_points_basic():
    s = h.Point(0,0,0)
    g = h.Point(5.5,0,0)
    spacing = 1
    actual = h.straight_points(s, g, spacing)
    expected = ([h.Point(1,0,0), h.Point(2,0,0), h.Point(3,0,0), h.Point(4,0,0), h.Point(5,0,0)], 0.5)

    assert len(actual[0]) == len(expected[0])
    assert actual[1] == expected[1]
    for p in range(len(actual)):
        assert actual[0][p].x == pytest.approx(expected[0][p].x, rel=0.01)
        assert actual[0][p].y == pytest.approx(expected[0][p].y, rel=0.01)
        assert actual[0][p].dx == pytest.approx(expected[0][p].dx, rel=0.01)
        assert actual[0][p].dy == pytest.approx(expected[0][p].dy, rel=0.01)

"""
Testing bend_points function

"""
def test_bend_points_basic():
    s = h.Point(0,0,pi/2,radius=-5)
    g = h.Point(5,5,0,radius=-5)
    spacing = 1
    actual = h.bend_points(s, g, spacing)
    expected = ([h.Point(0.099667111, 0.993346654, 1.370796327),
                 h.Point(0.39469503, 1.947091712, 1.170796327),
                 h.Point(0.873321925, 2.823212367, 0.970796327),
                 h.Point(1.516466453, 3.586780454, 0.770796327),
                 h.Point(2.298488471, 4.207354924, 0.570796327),
                 h.Point(3.188211228, 4.66019543, 0.370796327),
                 h.Point(4.150164285, 4.92724865, 0.170796327)], 0.853981634)

    assert len(actual[0]) == len(expected[0])
    assert actual[1] == pytest.approx(expected[1])
    for p in range(len(actual)):
        assert actual[0][p].x == pytest.approx(expected[0][p].x, rel=0.01)
        assert actual[0][p].y == pytest.approx(expected[0][p].y, rel=0.01)
        assert actual[0][p].dx == pytest.approx(expected[0][p].dx, rel=0.01)
        assert actual[0][p].dy == pytest.approx(expected[0][p].dy, rel=0.01)


"""
Testing path function.

"""

def test_path_basic():
    s = h.Point(0,-5,pi/2)
    g = h.Point(10,5,0)
    bends = [(h.Point(0,0,pi/2,radius=-5), h.Point(5,5,0,radius=-5))]
    spacing = 1

    actual = h.path(s, g, bends, spacing)
    expected = [h.Point(0,-5,pi/2),
                h.Point(0,-4,pi/2),
                h.Point(0,-3,pi/2),
                h.Point(0,-2,pi/2),
                h.Point(0,-1,pi/2),
                h.Point(0,0,pi/2), 
                h.Point(0.099667111, 0.993346654, 1.370796327),
                h.Point(0.39469503, 1.947091712, 1.170796327),
                h.Point(0.873321925, 2.823212367, 0.970796327),
                h.Point(1.516466453, 3.586780454, 0.770796327),
                h.Point(2.298488471, 4.207354924, 0.570796327),
                h.Point(3.188211228, 4.66019543, 0.370796327),
                h.Point(4.150164285, 4.92724865, 0.170796327),
                h.Point(6-0.853981634,-5,pi/2),
                h.Point(7-0.853981634,-5,pi/2),
                h.Point(8-0.853981634,-5,pi/2),
                h.Point(9-0.853981634,-5,pi/2),
                h.Point(10-0.853981634,-5,pi/2)]
    
    assert len(actual) == len(expected)


"""
Testing S_Path function.

This function returns the points along a straight path to the goal.It
should only return a path when the goal is within the heading tolerance
from the start, and when the end of line projected from the start point
by the distance between start and goal ends within the location
tolerance from the goal.
"""

def test_S_Path():
    s = h.Point(0,0,pi/4)
    g = h.Point(5,5,pi/4)
    min_straight = 0
    heading_tol = 1
    location_tol = 1

    actual = h.S_path(s, g, min_straight, heading_tol, location_tol)
    expected = [h.Point(0,0,pi/4), h.Point(0.707,0.707,pi/4), h.Point(1.414,1.414,pi/4),
                h.Point(2.121,2.121,pi/4), h.Point(2.828,2.828,pi/4), h.Point(3.536,3.536,pi/4),
                h.Point(4.243,4.243,pi/4), h.Point(4.95,4.95,pi/4)]
    assert len(actual) == len(expected)
    for p in range(len(actual)):
        assert actual[p].x == pytest.approx(expected[p].x, rel=0.01)
        assert actual[p].y == pytest.approx(expected[p].y, rel=0.01)
        assert actual[p].dx == pytest.approx(expected[p].dx, rel=0.01)
        assert actual[p].dy == pytest.approx(expected[p].dy, rel=0.01)

def test_S_path_bad_heading():
    s = h.Point(0,0,pi/4)
    g = h.Point(5,5,0)
    min_straight = 0
    heading_tol = pi/180
    location_tol = 1

    actual = h.S_path(s, g, min_straight, heading_tol, location_tol)
    expected = []
    assert len(actual) == len(expected)

def test_S_path_bad_goal():
    s = h.Point(0,0,pi/4)
    g = h.Point(5,0,pi/4)
    min_straight = 0
    heading_tol = pi/180
    location_tol = 1

    actual = h.S_path(s, g, min_straight, heading_tol, location_tol)
    expected = []
    assert len(actual) == len(expected)

""" 
Test shortest path functions

"""

def test_shortest_path():
    


"""

Test the tangent calculation to ensure inner and outer tangents are
calculated correctly in every direction and heading


tangent_cases = [(h.Point(5, 0, 0, radius=5), h.Point(15, 15, 0, radius=-5),
                  h.Point(10, 5, pi/2), h.Point(10, 10, pi/2), 0, 1, 0, 1),
                 (h.Point(15, 15, pi, radius=5), h.Point(5, 0, pi, radius=-5),
                  h.Point(10, 10, -pi/2), h.Point(10, 5, -pi/2), 0, -1, 0, -1),
                 (h.Point(5, 0, 0, radius=5), h.Point(15, 11, 0, radius=-5),
                  h.Point(10, 5, pi/2), h.Point(10, 6, pi/2), 0, 1, 0, 1),
                 (h.Point(5, 0, 0, radius=5), h.Point(5, 11, pi, radius=5),
                  h.Point(10, 5, pi/2), h.Point(10, 6, pi/2), 0, 1, 0, 1)]

@pytest.mark.parametrize("ss, gg, sg, gs, sg_dx, sg_dy, gs_dx, gs_dy", tangent_cases)
def test_tangent_points(ss, gg, sg, gs, sg_dx, sg_dy, gs_dx, gs_dy):
    actual = h.tangent_points(ss, gg)
    assert actual[0].x == pytest.approx(sg.x)
    assert actual[0].y == pytest.approx(sg.y)
    assert actual[0].dx == pytest.approx(sg_dx)
    assert actual[0].dy == pytest.approx(sg_dy)
    assert actual[1].x == pytest.approx(gs.x)
    assert actual[1].y == pytest.approx(gs.y)
    assert actual[1].dx == pytest.approx(gs_dx)
    assert actual[1].dy == pytest.approx(gs_dy)
"""