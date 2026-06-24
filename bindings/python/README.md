# Bezier — Python bindings

Python bindings for the [Bezier](../../README.md) C++ library, via
[nanobind](https://nanobind.readthedocs.io) and built with
[scikit-build-core](https://scikit-build-core.readthedocs.io).

## Install

Needs a C++17 compiler and Eigen3 on the system.

```bash
pip install .            # from the repo root
```

## Usage

```python
import numpy as np
import bezier

c = bezier.Curve(np.array([[0, 0], [1, 2], [3, 3], [4, 0]], float))
print(c.order())          # 3
print(c.valueAt(0.5))     # point on the curve -> numpy (2,)
print(c.valueAt([0.0, 0.5, 1.0]))   # batched -> numpy (3, 2)
print(c.length())
deriv = c.derivative()    # order-2 curve (borrows from c)

pc = bezier.PolyCurve([c, bezier.Curve(c.controlPoints() + 4)])
print(len(pc))            # 2
```

Method names match the C++ API. Points are numpy arrays of shape `(2,)`;
point/parameter lists are `(N, 2)` / `(N,)`. `boundingBox()` returns a
`(min, max)` tuple of points. The free functions `visvalingamWyatt`,
`solvePolynomial`, and `fitBezier` live under `bezier.Utils`.

## Test

```bash
pip install -e ".[test]" && pytest bindings/python/tests
```
