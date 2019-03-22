#include "bezier.h"

inline double binomial(uint n, uint k) { return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1)); }

namespace Bezier
{

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

void Curve::resetCache()
{
  cached_derivative_.reset();
  cached_ext_points_.reset();
  cached_bounding_box_tight_.reset();
  cached_bounding_box_relaxed_.reset();
  cached_polyline_.reset();
}

Curve::Coeffs Curve::bernsteinCoeffs() const
{
  if (bernstein_coeffs_.find(N_) == bernstein_coeffs_.end())
  {
    bernstein_coeffs_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
    for (uint k = 0; k < N_; k++)
      for (uint i = 0; i <= k; i++)
        bernstein_coeffs_.at(N_)(k, i) = pow(-1, i - k) * binomial(N_ - 1, k) * binomial(k, i);
  }
  return bernstein_coeffs_.at(N_);
}

Curve::Coeffs Curve::splittingCoeffsLeft(double z) const
{
  Curve::Coeffs coeffs(Coeffs::Zero(N_, N_));
  if (z == 0.5)
  {
    if (splitting_coeffs_left_.find(N_) == splitting_coeffs_left_.end())
    {
      splitting_coeffs_left_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
      for (uint k = 0; k < N_; k++)
        splitting_coeffs_left_.at(N_)(k, k) = pow(0.5, k);
      splitting_coeffs_left_.at(N_).noalias() =
          bernsteinCoeffs().inverse() * splitting_coeffs_left_.at(N_) * bernsteinCoeffs();
    }
    coeffs = splitting_coeffs_left_.at(N_);
  }
  else
  {
    for (uint k = 0; k < N_; k++)
      coeffs(k, k) = pow(z, k);
    coeffs.noalias() = bernsteinCoeffs().inverse() * coeffs * bernsteinCoeffs();
  }
  return coeffs;
}

Curve::Coeffs Curve::splittingCoeffsRight(double z) const
{
  Curve::Coeffs coeffs(Coeffs::Zero(N_, N_));
  if (z == 0.5)
  {
    if (splitting_coeffs_right_.find(N_) == splitting_coeffs_right_.end())
    {
      splitting_coeffs_right_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
      for (uint k = 0; k < N_; k++)
        for (uint i = 0; i <= k; i++)
          splitting_coeffs_right_.at(N_)(N_ - 1 - k, N_ - 1 - (k - i)) = splittingCoeffsLeft()(k, i);
    }
    coeffs = splitting_coeffs_right_.at(N_);
  }
  else
  {
    for (uint k = 0; k < N_; k++)
      for (uint i = 0; i <= k; i++)
        coeffs(N_ - 1 - k, N_ - 1 - (k - i)) = splittingCoeffsLeft(z)(k, i);
  }
  return coeffs;
}

Curve::Coeffs Curve::elevateOrderCoeffs(uint n) const
{
  if (elevate_order_coeffs_.find(n) == elevate_order_coeffs_.end())
  {
    elevate_order_coeffs_.insert(std::make_pair(n, Coeffs::Zero(n + 1, n)));
    for (uint k = 0; k < n + 1; k++)
    {
      if (k == 0)
        elevate_order_coeffs_.at(n)(k, 0) = 1;
      else if (k == n)
        elevate_order_coeffs_.at(n)(k, k - 1) = 1;
      else
      {
        elevate_order_coeffs_.at(n)(k, k - 1) = 1. * k / n;
        elevate_order_coeffs_.at(n)(k, k) = 1 - 1. * k / n;
      }
    }
  }
  return elevate_order_coeffs_.at(n);
}

Curve::Coeffs Curve::lowerOrderCoeffs(uint n) const
{
  if (lower_order_coeffs_.find(n) == lower_order_coeffs_.end())
  {
    lower_order_coeffs_.insert(std::make_pair(n, Coeffs::Zero(n - 1, n)));
    lower_order_coeffs_.at(n).noalias() =
        (elevateOrderCoeffs(n - 1).transpose() * elevateOrderCoeffs(n - 1)).inverse() *
        elevateOrderCoeffs(n - 1).transpose();
  }

  return lower_order_coeffs_.at(n);
}

Curve::Curve(const Eigen::MatrixX2d& points)
{
  N_ = static_cast<uint>(points.rows());
  control_points_ = points;
}

Curve::Curve(const PointVector& points)
{
  N_ = static_cast<uint>(points.size());
  control_points_.resize(N_, 2);
  for (uint k = 0; k < N_; k++)
    control_points_.row(k) = points.at(k);
}

uint Curve::getOrder() { return N_ - 1; }

PointVector Curve::getControlPoints() const
{
  PointVector points(N_);
  for (uint k = 0; k < N_; k++)
    points.at(k) = control_points_.row(k);
  return points;
}

std::pair<Point, Point> Curve::getEndPoints() const
{
  return std::make_pair(control_points_.row(0), control_points_.row(N_ - 1));
}

PointVector Curve::getPolyline(double smoothness, double precision) const
{
  if (!cached_polyline_ || smoothness != cached_polyline_smootheness_)
  {
    PointVector* polyline = new PointVector;
    std::vector<std::shared_ptr<Curve>> subcurves;
    subcurves.push_back(std::make_shared<Bezier::Curve>(this->getControlPoints()));
    polyline->push_back(control_points_.row(0));
    while (!subcurves.empty())
    {
      auto c = subcurves.back();
      subcurves.pop_back();
      auto cp = c->getControlPoints();
      double hull = 0;
      for (uint k = 1; k < cp.size(); k++)
        hull += (cp.at(k - 1) - cp.at(k)).norm();
      if (hull < smoothness * (cp.front() - cp.back()).norm() ||
          (cp.front() - cp.back()).norm() / smoothness < precision)
      {
        polyline->push_back(cp.back());
      }
      else
      {
        auto split = c->splitCurve();
        subcurves.push_back(std::make_shared<Bezier::Curve>(split.second));
        subcurves.push_back(std::make_shared<Bezier::Curve>(split.first));
      }
    }
    (const_cast<Curve*>(this))->cached_polyline_smootheness_ = smoothness;
    (const_cast<Curve*>(this))->cached_polyline_.reset(polyline);
  }
  return *cached_polyline_;
}

void Curve::reverse()
{
  control_points_ = control_points_.colwise().reverse().eval();
  resetCache();
}

void Curve::manipulateControlPoint(uint idx, const Point& point)
{
  control_points_.row(idx) = point;
  resetCache();
}

void Curve::manipulateCurvature(double t, const Point& point)
{
  if (N_ < 3 || N_ > 4)
    throw "Only quadratic and cubic curves can be manipulated";

  double r = fabs((pow(t, N_ - 1) + pow(1 - t, N_ - 1) - 1) / (pow(t, N_ - 1) + pow(1 - t, N_ - 1)));
  double u = pow(1 - t, N_ - 1) / (pow(t, N_ - 1) + pow(1 - t, N_ - 1));
  Point C = u * control_points_.row(0) + (1 - u) * control_points_.row(N_ - 1);
  Point B = point;
  Point A = B - (C - B) / r;

  switch (N_)
  {
  case 3:
    control_points_.row(1) = A;
    break;
  case 4:
    Point e1 = control_points_.row(0) * pow(1 - t, 2) + control_points_.row(1) * 2 * t * (1 - t) +
               control_points_.row(2) * pow(t, 2);
    Point e2 = control_points_.row(1) * pow(1 - t, 2) + control_points_.row(2) * 2 * t * (1 - t) +
               control_points_.row(3) * pow(t, 2);
    e1.noalias() = B + e1 - valueAt(t);
    e2.noalias() = B + e2 - valueAt(t);
    Point v1 = A - (A - e1) / (1 - t);
    Point v2 = A + (e2 - A) / t;
    control_points_.row(1).noalias() = control_points_.row(0) + (v1.transpose() - control_points_.row(0)) / t;
    control_points_.row(2).noalias() = control_points_.row(3) - (control_points_.row(3) - v2.transpose()) / (1 - t);
  }
  resetCache();
}

void Curve::elevateOrder()
{
  Eigen::MatrixXd new_points = elevateOrderCoeffs(N_) * control_points_;
  control_points_.resize(++N_, 2);
  control_points_ = new_points;
  resetCache();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw "Cannot further reduce the order of curve.";
  Eigen::MatrixXd new_points = lowerOrderCoeffs(N_) * control_points_;
  control_points_.resize(--N_, 2);
  control_points_ = new_points;
  resetCache();
}

Point Curve::valueAt(double t) const
{
  Eigen::VectorXd power_basis;
  power_basis.resize(N_);
  for (uint k = 0; k < N_; k++)
    power_basis(k) = pow(t, k);

  return (power_basis.transpose() * bernsteinCoeffs() * control_points_).transpose();
}

double Curve::curvatureAt(double t) const
{
  Point d = getDerivative().valueAt(t);
  Point dd = getDerivative().getDerivative().valueAt(t);
  double num = d.x() * dd.y() - d.y() * dd.x();
  double dnm = pow(pow(d.x(), 2) + pow(d.y(), 2), 3 / 2);

  return num / dnm;
}

Vec2 Curve::tangentAt(double t, bool normalize) const
{
  Point p(getDerivative().valueAt(t));
  if (normalize)
    p.normalize();
  return p;
}

Vec2 Curve::normalAt(double t, bool normalize) const
{
  Point tangent = tangentAt(t, normalize);
  return Vec2(-tangent.y(), tangent.x());
}

Curve Curve::getDerivative() const
{
  if (!cached_derivative_)
  {
    Eigen::MatrixX2d new_points;
    new_points.resize(N_ - 1, 2);
    for (uint k = 0; k < N_ - 1; k++)
      new_points.row(k).noalias() = N_ * (control_points_.row(k + 1) - control_points_.row(k));

    (const_cast<Curve*>(this))->cached_derivative_ = std::make_shared<Curve>(new_points);
  }
  return *cached_derivative_;
}

PointVector Curve::getRoots(double step, double epsilon, std::size_t max_iter) const
{
  if (!cached_ext_points_)
  {
    (const_cast<Curve*>(this))->cached_ext_points_ = std::make_shared<PointVector>();
    std::vector<double> added_t;

    // check both axes
    for (long k = 0; k < 2; k++)
    {
      double t = 0;
      while (t <= 1.0)
      {
        double t_current = t;
        std::size_t current_iter = 0;

        // it has to converge in max_iter steps
        while (current_iter < max_iter)
        {
          // Newton-Rhapson: f = f - f' / f''
          double t_new =
              t_current - getDerivative().valueAt(t_current)[k] / getDerivative().getDerivative().valueAt(t_current)[k];

          // if there is no change to t_current
          if (fabs(t_new - t_current) < epsilon)
          {
            // check if between [0, 1]
            if (t_new >= -epsilon && t_new <= 1 + epsilon)
            {
              // check if same value wasn't found before
              if (added_t.end() == std::find_if(added_t.begin(), added_t.end(), [t_new, epsilon](const double& val)
                                                {
                                                  return fabs(val - t_new) < epsilon;
                                                }))
              {
                // add new value and point
                added_t.push_back(t_new);
                (const_cast<Curve*>(this))->cached_ext_points_->push_back(valueAt(t_new));
                break;
              }
            }
          }

          t_current = t_new;
          current_iter++;
        }

        t += step;
      }
    }
  }
  return *cached_ext_points_;
}

BBox Curve::getBBox(bool use_roots) const
{
  if (!(use_roots ? cached_bounding_box_tight_ : cached_bounding_box_relaxed_))
  {
    PointVector extremes;
    if (use_roots)
    {
      extremes = getRoots();
      extremes.push_back(control_points_.row(0));
      extremes.push_back(control_points_.row(N_ - 1));
    }
    else
    {
      for (uint k = 0; k < control_points_.rows(); k++)
        extremes.push_back(control_points_.row(k));
    }

    // find mininum and maximum along each axis
    auto x_extremes = std::minmax_element(extremes.begin(), extremes.end(), [](const Point& lhs, const Point& rhs)
                                          {
                                            return lhs.x() < rhs.x();
                                          });
    auto y_extremes = std::minmax_element(extremes.begin(), extremes.end(), [](const Point& lhs, const Point& rhs)
                                          {
                                            return lhs.y() < rhs.y();
                                          });
    if (use_roots)
      (const_cast<Curve*>(this))->cached_bounding_box_tight_ = std::make_shared<BBox>(
          Point(x_extremes.first->x(), y_extremes.first->y()), Point(x_extremes.second->x(), y_extremes.second->y()));
    else
      (const_cast<Curve*>(this))->cached_bounding_box_relaxed_ = std::make_shared<BBox>(
          Point(x_extremes.first->x(), y_extremes.first->y()), Point(x_extremes.second->x(), y_extremes.second->y()));
  }
  return *(use_roots ? cached_bounding_box_tight_ : cached_bounding_box_relaxed_);
}

std::pair<Curve, Curve> Curve::splitCurve(double z) const
{
  return std::make_pair(Curve(splittingCoeffsLeft(z) * control_points_),
                        Curve(splittingCoeffsRight(z) * control_points_));
}

PointVector Curve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points_of_intersection;
  std::vector<std::pair<ConstCurvePtr, ConstCurvePtr>> subcurve_pairs;

  if (this != &curve)
  {
    // we don't know if shared_ptr of "this" and "curve" exists
    // if not, make temporary copies and use them
    subcurve_pairs.push_back(std::make_pair(std::shared_ptr<Curve>(), std::shared_ptr<Curve>()));
    try // for "this"
    {
      subcurve_pairs.front().first = this->shared_from_this();
    }
    catch (std::bad_weak_ptr const&)
    {
      subcurve_pairs.front().first = std::make_shared<const Curve>(*this);
    }
    try // for "curve"
    {
      subcurve_pairs.front().second = curve.shared_from_this();
    }
    catch (std::bad_weak_ptr const&)
    {
      subcurve_pairs.front().second = std::make_shared<const Curve>(curve);
    }
  }
  else
  {
    // TODO: self-interserction
  }

  while (!subcurve_pairs.empty())
  {
    ConstCurvePtr part_a = std::get<0>(subcurve_pairs.back());
    ConstCurvePtr part_b = std::get<1>(subcurve_pairs.back());
    subcurve_pairs.pop_back();

    BBox bbox1 = part_a->getBBox(false); // very slow with tight BBox (roots)
    BBox bbox2 = part_b->getBBox(false); // very slow with tight BBox (roots)
    if (!bbox1.intersects(bbox2))
    {
      // no intersection
      continue;
    }

    if (bbox1.diagonal().norm() < epsilon && bbox2.diagonal().norm() < epsilon)
    {
      // segments converged, check if not already found and add new
      Point new_point = part_a->valueAt(0.5);
      if (points_of_intersection.end() == std::find_if(points_of_intersection.begin(), points_of_intersection.end(),
                                                       [new_point, epsilon](const Point& point)
                                                       {
                                                         return (point - new_point).norm() < epsilon;
                                                       }))
      {
        points_of_intersection.push_back(new_point);

        // if only first point is needed, stop
        if (stop_at_first)
          return points_of_intersection;
      }
      continue;
    }

    // intersection exists, but segments are still too large
    // divide both segments in half and new pairs
    // LIFO : we want to first discover closest intersection (smallest t on this curve)
    // so it is important which pair of subcurves is inserted first
    std::vector<ConstCurvePtr> subcurves_a;
    std::vector<ConstCurvePtr> subcurves_b;

    if (bbox1.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      subcurves_a.push_back(part_a);
    }
    else
    {
      // divide into two subcurves
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_a.push_back(std::make_shared<Curve>(part_a->splittingCoeffsRight() * part_a->control_points_));
      subcurves_a.push_back(std::make_shared<Curve>(part_a->splittingCoeffsLeft() * part_a->control_points_));
    }

    if (bbox2.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_b.push_back(part_b);
    }
    else
    {
      // divide into two subcurves
      subcurves_b.push_back(std::make_shared<Curve>(part_b->splittingCoeffsRight() * part_b->control_points_));
      subcurves_b.push_back(std::make_shared<Curve>(part_b->splittingCoeffsLeft() * part_b->control_points_));
    }

    // insert all combinations for next iteration
    // last pair is one where both subcurves have smalles t ranges
    for (auto&& subcurve_b : subcurves_b)
      for (auto&& subcurve_a : subcurves_a)
        subcurve_pairs.push_back(std::make_pair(subcurve_a, subcurve_b));
  }
  return points_of_intersection;
}

double Curve::projectPoint(const Point& point, double step) const
{
  double t = 0;
  double t_dist = (valueAt(t) - point).norm();

  // Coarse search
  for (double k = step; k < 1 + step; k += step)
  {
    double new_dist = (valueAt(k) - point).norm();
    if (new_dist < t_dist)
    {
      t_dist = new_dist;
      t = k;
    }
  }

  // fine search
  do
  {
    double new_dist1 = (valueAt(t - step / 2) - point).norm();
    double new_dist2 = (valueAt(t + step / 2) - point).norm();

    if (new_dist1 > t_dist && new_dist2 > t_dist)
    {
      step /= 2;
      continue;
    }
    if ((new_dist1 < new_dist2 ? new_dist1 : new_dist2) < t_dist)
    {
      t_dist = (new_dist1 < new_dist2 ? new_dist1 : new_dist2);
      t += new_dist1 < new_dist2 ? -step / 2 : step / 2;
    }
  } while (step > 0.0001);

  // if closest point is between <1, 1+step>
  if (t > 1)
    t = 1;
  // if closest point is between < -step, 0>
  if (t < 0)
    t = 0;

  return t;
}

}
