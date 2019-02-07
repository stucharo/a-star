import pytest
from math import pi

import heuristic as h

def test_develop_path_straight():
    s = h.Point(0, 10, pi/2)
    g = h.Point(0, 90, pi/2)
    msl = 100

    actual = h.develop_path(s, g, msl, 0)
    expected = [s, g]

    assert len(actual) == len(expected)
    for i in range(len(actual)):
        assert actual[i].x == pytest.approx(expected[i].x)
        assert actual[i].y == pytest.approx(expected[i].y)
        assert actual[i].dx == pytest.approx(expected[i].dx)
        assert actual[i].dy == pytest.approx(expected[i].dy)

def test_develop_path_straight_too_short():
    s = h.Point(0, 10, pi/2, radius=-10)
    g = h.Point(0, 90, pi/2)
    msl = 100

    actual = h.develop_path(s, g, msl, 0)
    expected = []

    assert len(actual) == len(expected)