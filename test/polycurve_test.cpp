#include "test_data.hpp"
#include "test_oracles.hpp"

#include <stdexcept>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
#include "Bezier/polycurve.h"
#include "Bezier/utils.h"

namespace Bezier
{

// Fixture: three C0-chained cubics, each a translated copy of the test cubic
// (the translation is the curve's own start-to-end vector, so endpoints chain).
class PolyCurveTest : public ::testing::Test
{
public:
  PolyCurveTest()
  {
    Eigen::MatrixX2d cp = curvePointsAsMatrix();
    Eigen::RowVector2d delta = cp.row(cp.rows() - 1) - cp.row(0);
    std::deque<Curve> curves;
    for (int k = 0; k < 3; k++)
      curves.emplace_back(Eigen::MatrixX2d(cp.rowwise() + k * delta));
    poly_ = PolyCurve(std::move(curves));
  }

protected:
  PolyCurve poly_;
};

TEST_F(PolyCurveTest, SizeAndAccessors)
{
  EXPECT_EQ(poly_.size(), 3u);
  EXPECT_EQ(poly_.curves().size(), 3u);

  // Subcurves chain C0
  for (unsigned k = 0; k + 1 < poly_.size(); k++)
  {
    Point left_end = poly_.curve(k).endPoints().second;
    Point right_start = poly_.curve(k + 1).endPoints().first;
    EXPECT_NEAR(left_end.x(), right_start.x(), Utils::epsilon);
    EXPECT_NEAR(left_end.y(), right_start.y(), Utils::epsilon);
  }
}

TEST_F(PolyCurveTest, InsertAndRemove)
{
  PolyCurve poly{poly_};
  Curve extra{curvePointsAsMatrix()};

  poly.insertBack(extra);
  EXPECT_EQ(poly.size(), 4u);
  poly.removeBack();
  EXPECT_EQ(poly.size(), 3u);

  poly.insertFront(extra);
  EXPECT_EQ(poly.size(), 4u);
  EXPECT_EQ(poly.curve(0).controlPoints(), extra.controlPoints());
  poly.removeFirst();
  EXPECT_EQ(poly.size(), 3u);

  poly.insertAt(1, extra);
  EXPECT_EQ(poly.size(), 4u);
  EXPECT_EQ(poly.curve(1).controlPoints(), extra.controlPoints());
  poly.removeAt(1);
  EXPECT_EQ(poly.size(), 3u);

  // Back to the original chain
  for (unsigned k = 0; k < poly.size(); k++)
    EXPECT_EQ(poly.curve(k).controlPoints(), poly_.curve(k).controlPoints());
}

TEST_F(PolyCurveTest, CurveIdxBoundaries)
{
  EXPECT_EQ(poly_.curveIdx(0.0), 0u);
  EXPECT_EQ(poly_.curveIdx(0.99), 0u);
  EXPECT_EQ(poly_.curveIdx(1.0), 1u);
  EXPECT_EQ(poly_.curveIdx(1.5), 1u);
  EXPECT_EQ(poly_.curveIdx(2.0), 2u);
  EXPECT_EQ(poly_.curveIdx(2.7), 2u);
  // t == size() maps to the last subcurve
  EXPECT_EQ(poly_.curveIdx(3.0), 2u);
  // out-of-range parameters clamp to the first / last subcurve
  EXPECT_EQ(poly_.curveIdx(-1.0), 0u);
  EXPECT_EQ(poly_.curveIdx(1e9), poly_.size() - 1);
}

TEST(PolyCurveEmptyTest, AccessThrows)
{
  PolyCurve empty;
  EXPECT_THROW(empty.curveIdx(0.0), std::logic_error);
  EXPECT_THROW(empty.valueAt(0.0), std::logic_error);
  EXPECT_THROW(empty.endPoints(), std::logic_error);
  EXPECT_THROW(empty.length(), std::logic_error);
}

TEST_F(PolyCurveTest, ValueAtJointsAndContinuity)
{
  // Integer parameters hit the shared joint points
  for (unsigned k = 1; k < poly_.size(); k++)
  {
    Point joint = poly_.curve(k).endPoints().first;
    Point value = poly_.valueAt(static_cast<double>(k));
    EXPECT_NEAR(value.x(), joint.x(), Utils::epsilon);
    EXPECT_NEAR(value.y(), joint.y(), Utils::epsilon);

    // Approaching the joint from both sides gives the same point
    Point before = poly_.valueAt(k - 1e-9);
    Point after = poly_.valueAt(k + 1e-9);
    EXPECT_NEAR(before.x(), after.x(), 1e-5);
    EXPECT_NEAR(before.y(), after.y(), 1e-5);
  }

  // Endpoints of the global parameter range
  Point start = poly_.valueAt(0.0);
  Point end = poly_.valueAt(3.0);
  EXPECT_NEAR(start.x(), poly_.curve(0).endPoints().first.x(), Utils::epsilon);
  EXPECT_NEAR(end.x(), poly_.curve(2).endPoints().second.x(), Utils::epsilon);
}

TEST_F(PolyCurveTest, LengthAdditivity)
{
  double subcurve_sum{};
  for (unsigned k = 0; k < poly_.size(); k++)
    subcurve_sum += poly_.curve(k).length();
  EXPECT_NEAR(poly_.length(), subcurve_sum, 1e-8);

  // Additivity across one and two joints
  EXPECT_NEAR(poly_.length(0.5, 1.5), poly_.length(0.5, 1.0) + poly_.length(1.0, 1.5), 1e-8);
  EXPECT_NEAR(poly_.length(0.5, 2.5),
              poly_.length(0.5, 1.0) + poly_.length(1.0, 2.0) + poly_.length(2.0, 2.5), 1e-8);

  // Pinned contract: swapped arguments yield the negated length, across both the
  // adjacent-subcurve branch (idx1+1==idx2) and the multi-subcurve accumulate branch
  EXPECT_DOUBLE_EQ(poly_.length(1.5, 0.5), -poly_.length(0.5, 1.5));
  EXPECT_DOUBLE_EQ(poly_.length(2.5, 0.5), -poly_.length(0.5, 2.5));
}

TEST_F(PolyCurveTest, StepWithinAndAcrossJoints)
{
  // Within a single subcurve
  double t_within = poly_.step(0.4, 10.0);
  EXPECT_NEAR(poly_.length(0.4, t_within), 10.0, 1e-6);

  // Across a joint, forward and backward
  double ds_cross = poly_.length(0.9, 1.5);
  EXPECT_NEAR(poly_.step(0.9, ds_cross), 1.5, 1e-6);
  EXPECT_NEAR(poly_.step(1.5, -ds_cross), 0.9, 1e-6);

  // Across two joints
  double ds_cross2 = poly_.length(0.5, 2.5);
  EXPECT_NEAR(poly_.step(0.5, ds_cross2), 2.5, 1e-6);

  // Out-of-range distances clamp to the parameter range [0, size()]
  EXPECT_DOUBLE_EQ(poly_.step(1.5, 1e6), 3.0);
  EXPECT_DOUBLE_EQ(poly_.step(1.5, -1e6), 0.0);
}

TEST_F(PolyCurveTest, PolylineContract)
{
  PointVector polyline = poly_.polyline();
  ParamVector params = poly_.polylineParams();

  ASSERT_EQ(polyline.size(), params.size()) << "polyline and polylineParams must correspond 1:1";
  EXPECT_DOUBLE_EQ(params.front(), 0.0);
  EXPECT_DOUBLE_EQ(params.back(), 3.0);

  // Joint vertices must not be duplicated
  for (size_t i = 1; i < polyline.size(); i++)
    EXPECT_GT((polyline[i] - polyline[i - 1]).norm(), 0.0) << "duplicated vertex at " << i;

  // Every vertex lies on the polycurve at its reported parameter
  for (size_t i = 0; i < params.size(); i++)
  {
    Point on_curve = poly_.valueAt(params[i]);
    EXPECT_NEAR(on_curve.x(), polyline[i].x(), 1e-6) << "vertex " << i << " not on curve";
    EXPECT_NEAR(on_curve.y(), polyline[i].y(), 1e-6) << "vertex " << i << " not on curve";
  }
}

TEST_F(PolyCurveTest, SetControlPointGlobalIndex)
{
  PolyCurve poly{poly_};
  const Point marker{-999, -999};

  // Each cubic contributes 4 control points: global 5 -> subcurve 1 local 1
  poly.setControlPoint(5, marker);
  EXPECT_EQ(poly.curve(1).controlPoint(1), marker);
  EXPECT_EQ(poly.curve(0).controlPoints(), poly_.curve(0).controlPoints());

  // Global 9 -> subcurve 2 local 1
  poly.setControlPoint(9, marker);
  EXPECT_EQ(poly.curve(2).controlPoint(1), marker);
}

TEST_F(PolyCurveTest, ControlPointsEndPointsBoundingBox)
{
  EXPECT_EQ(poly_.controlPoints().size(), 12u); // 3 cubics x 4 points

  auto end_points = poly_.endPoints();
  EXPECT_EQ(end_points.first, poly_.curve(0).endPoints().first);
  EXPECT_EQ(end_points.second, poly_.curve(2).endPoints().second);

  // Bounding box is the union of subcurve bounding boxes
  BoundingBox expected = poly_.curve(0).boundingBox();
  for (unsigned k = 1; k < poly_.size(); k++)
    expected.extend(poly_.curve(k).boundingBox());
  BoundingBox actual = poly_.boundingBox();
  EXPECT_NEAR(actual.min().x(), expected.min().x(), Utils::epsilon);
  EXPECT_NEAR(actual.min().y(), expected.min().y(), Utils::epsilon);
  EXPECT_NEAR(actual.max().x(), expected.max().x(), Utils::epsilon);
  EXPECT_NEAR(actual.max().y(), expected.max().y(), Utils::epsilon);
}

TEST_F(PolyCurveTest, ProjectPointAndDistanceConsistency)
{
  // Probe points near each subcurve land in that subcurve's parameter range
  for (unsigned k = 0; k < poly_.size(); k++)
  {
    Point probe = poly_.valueAt(k + 0.5) + Vector{3, 3};
    double t = poly_.projectPoint(probe);
    EXPECT_GE(t, static_cast<double>(k)) << "projection escaped subcurve " << k;
    EXPECT_LE(t, static_cast<double>(k + 1)) << "projection escaped subcurve " << k;

    // distance must agree with the projection
    EXPECT_NEAR(poly_.distance(probe), (probe - poly_.valueAt(t)).norm(), Utils::epsilon);
  }

  // A point on the polycurve has distance ~0
  EXPECT_NEAR(poly_.distance(poly_.valueAt(1.7)), 0.0, 1e-6);
}

TEST_F(PolyCurveTest, BatchProjectPointAndDistance)
{
  // The vector overloads must agree element-wise with the scalar ones
  PointVector probes;
  for (unsigned k = 0; k < poly_.size(); k++)
    probes.push_back(poly_.valueAt(k + 0.5) + Vector{3, -4});

  ParamVector t_batch = poly_.projectPoint(probes);
  std::vector<double> d_batch = poly_.distance(probes);
  ASSERT_EQ(t_batch.size(), probes.size());
  ASSERT_EQ(d_batch.size(), probes.size());
  for (size_t i = 0; i < probes.size(); i++)
  {
    EXPECT_DOUBLE_EQ(t_batch[i], poly_.projectPoint(probes[i]));
    EXPECT_DOUBLE_EQ(d_batch[i], poly_.distance(probes[i]));
  }
}

TEST_F(PolyCurveTest, IntersectionsWithCurveAndPolyCurve)
{
  // A vertical cubic "line" crossing the middle subcurve once
  Point mid = poly_.valueAt(1.5);
  Eigen::MatrixX2d cross_cp(4, 2);
  cross_cp << mid.x(), mid.y() - 100, mid.x(), mid.y() - 33, mid.x(), mid.y() + 33, mid.x(), mid.y() + 100;
  Curve cross{cross_cp};

  PointVector with_curve = poly_.intersections(cross);
  ASSERT_GE(with_curve.size(), 1u) << "Expected at least one intersection";
  for (const Point& p : with_curve)
  {
    EXPECT_NEAR(poly_.distance(p), 0.0, 1e-4) << "intersection point not on polycurve";
    EXPECT_NEAR(cross.distance(p), 0.0, 1e-4) << "intersection point not on crossing curve";
  }

  PolyCurve cross_poly{std::deque<Curve>{cross}};
  PointVector with_poly = poly_.intersections(cross_poly);
  EXPECT_EQ(with_poly.size(), with_curve.size());
}

} // namespace Bezier
