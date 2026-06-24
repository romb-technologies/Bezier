"""Smoke + correctness check for the nanobind bindings. Run: pytest bindings/python/tests"""
import numpy as np
import bezier


CP = np.array([[0.0, 0.0], [1.0, 2.0], [3.0, 3.0], [4.0, 0.0]])  # a cubic


def test_construct_and_order():
    c = bezier.Curve(CP)
    assert c.order() == 3
    assert np.allclose(c.controlPoints(), CP)


def test_endpoints_match_value_at():
    c = bezier.Curve(CP)
    assert np.allclose(c.valueAt(0.0), CP[0])
    assert np.allclose(c.valueAt(1.0), CP[-1])


def test_batch_valueat_shape():
    c = bezier.Curve(CP)
    pts = c.valueAt([0.0, 0.5, 1.0])
    assert np.asarray(pts).shape == (3, 2)


def test_length_positive():
    assert bezier.Curve(CP).length() > 0.0


def test_derivative_lowers_order():
    c = bezier.Curve(CP)
    assert c.derivative().order() == 2
    assert c.derivative(2).order() == 1


def test_split_roundtrips_endpoints():
    c = bezier.Curve(CP)
    left, right = c.splitCurve(0.5)
    assert np.allclose(left.valueAt(0.0), CP[0])
    assert np.allclose(right.valueAt(1.0), CP[-1])
    assert np.allclose(left.valueAt(1.0), right.valueAt(0.0))  # C0 at the seam


def test_bounding_box_pair():
    lo, hi = bezier.Curve(CP).boundingBox()
    assert np.all(np.asarray(lo) <= np.asarray(hi))


def test_polycurve_size():
    pc = bezier.PolyCurve([bezier.Curve(CP), bezier.Curve(CP + 4.0)])
    assert pc.size() == 2
    assert len(pc) == 2
    assert pc.curve(0).order() == 3
