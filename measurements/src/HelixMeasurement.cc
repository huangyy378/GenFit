/* HelixMeasurement.cpp */
#include "HelixMeasurement.h"

#include "DetPlane.h"
#include "Exception.h"
#include "HMatrixU.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "RKTrackRep.h"
#include "StateOnPlane.h"

#include "TMath.h"

#include <TClass.h>
#include <cassert>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <stdexcept>
#include <vector>

// ClassImp(genfit::HelixMeasurement)

namespace genfit {

HelixMeasurement::HelixMeasurement(int nDim) : AbsMeasurement(nDim) {
  assert(nDim == 8); // 8维
}

HelixMeasurement::HelixMeasurement(const TVectorD &rawHitCoords,
                                   const TMatrixDSym &rawHitCov, int detId,
                                   int hitId, TrackPoint *trackPoint)
    : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint) {
  // 验证输入维度
  if (rawHitCoords_.GetNrows() != 8) {
    throw Exception("HelixMeasurement requires 8-dimensional rawHitCoords",
                    __LINE__, __FILE__);
  }
}

SharedPlanePtr
HelixMeasurement::constructPlane(const StateOnPlane &state) const {
  // extrapolate to cylinder of radius r
  const AbsTrackRep *rep = state.getRep();
  StateOnPlane st(state);
  rep->extrapolateToCylinder(st, rawHitCoords_(3));
  TVector3 currentPos = rep->getPos(st);

  // find closest point on helix
  auto closest = findClosestPointOnHelix(currentPos);
  const TVector3 &pocaOnHelix = closest.point;
  TVector3 tangent = closest.tangent;
  tangent.SetMag(1.0);

  // get track direction at poca
  TVector3 dirInPoca = rep->getMom(st);
  dirInPoca.SetMag(1.0);

  if (std::fabs(tangent.Angle(dirInPoca)) < 0.01) {
    throw Exception("HelixMeasurement::constructPlane: Direction is parallel "
                    "to helix tangent",
                    __LINE__, __FILE__);
  }

  // construct plane
  TVector3 U = dirInPoca.Cross(tangent);
  if (U.Mag() < 1e-10) {
    U = TVector3(1, 0, 0).Cross(tangent);
    if (U.Mag() < 1e-10) {
      U = TVector3(0, 1, 0).Cross(tangent);
    }
  }
  U.SetMag(1.0);
  auto plane = new DetPlane(pocaOnHelix, U, tangent);

  return SharedPlanePtr(plane);
}

std::vector<MeasurementOnPlane *>
HelixMeasurement::constructMeasurementsOnPlane(
    const StateOnPlane &state) const {
  double d{};
  double V{rawHitCov_(7, 7)};
  return {new MeasurementOnPlane(TVectorD(1, &d), TMatrixDSym(1, &V),
                                 state.getPlane(), state.getRep(),
                                 constructHMatrix(state.getRep()))};
}

const AbsHMatrix *
HelixMeasurement::constructHMatrix(const AbsTrackRep *rep) const {
  if (dynamic_cast<const RKTrackRep *>(rep) == nullptr) {
    throw Exception(
        "HelixMeasurement can only handle state vectors of type RKTrackRep!",
        __LINE__, __FILE__);
  }
  return new HMatrixU();
}

HelixMeasurement::ClosestPointResult
HelixMeasurement::findClosestPointOnHelix(const TVector3 &point) const {
  ClosestPointResult result;

  // get helix parameters
  const TVector3 center(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));
  const double radius = rawHitCoords_(3);
  const double pitchAngle = rawHitCoords_(4);
  const double phi0 = rawHitCoords_(5);
  const double phiTotal = rawHitCoords_(6);

  const auto cosA = std::cos(pitchAngle);
  const auto sinA = std::sin(pitchAngle);
  const auto tanA = sinA / cosA;
  const auto tanAR = radius * tanA;
  const auto zOffset = phiTotal / 2 * tanAR;

  const TVector3 relPos = point - center;

  auto helixPoint = [&](double u) -> TVector3 {
    const double rotatedU = u + phi0;
    return TVector3(radius * std::cos(rotatedU), radius * std::sin(rotatedU),
                    u * tanAR - zOffset);
  };

  // 距离平方函数
  auto distanceSq = [&](double u) -> double {
    const TVector3 helixPt = helixPoint(u);
    return (helixPt.X() - relPos.X()) * (helixPt.X() - relPos.X()) +
           (helixPt.Y() - relPos.Y()) * (helixPt.Y() - relPos.Y()) +
           (helixPt.Z() - relPos.Z()) * (helixPt.Z() - relPos.Z());
  };

  // 一阶导数函数
  auto derivative = [&](double u) -> double {
    const double u1 = u + phi0;
    const TVector3 helixPt = helixPoint(u);

    // 螺旋线导数分量
    const double dx_du = -radius * std::sin(u1);
    const double dy_du = radius * std::cos(u1);
    const double dz_du = tanAR;

    return 2 * (helixPt.X() - relPos.X()) * dx_du +
           2 * (helixPt.Y() - relPos.Y()) * dy_du +
           2 * (helixPt.Z() - relPos.Z()) * dz_du;
  };

  // 二阶导数函数
  auto secondDerivative = [&](double u) -> double {
    const double u1 = u + phi0;
    const TVector3 helixPt = helixPoint(u);

    // 螺旋线一阶导数分量
    const double dx_du = -radius * std::sin(u1);
    const double dy_du = radius * std::cos(u1);
    const double dz_du = tanAR;

    // 螺旋线二阶导数分量
    const double d2x_du2 = -radius * std::cos(u1);
    const double d2y_du2 = -radius * std::sin(u1);
    const double d2z_du2 = 0;

    return 2 * (dx_du * dx_du + (helixPt.X() - relPos.X()) * d2x_du2) +
           2 * (dy_du * dy_du + (helixPt.Y() - relPos.Y()) * d2y_du2) +
           2 * (dz_du * dz_du + (helixPt.Z() - relPos.Z()) * d2z_du2);
  };

  double bestU = 0.0;
  double minDistSq = std::numeric_limits<double>::max();
  const int numSamples = 100;
  const double startU = 0;
  const double endU = phiTotal;

  for (int i = 0; i <= numSamples; ++i) {
    const double u = startU + i * (endU - startU) / numSamples;
    const double distSq = distanceSq(u);
    if (distSq < minDistSq) {
      minDistSq = distSq;
      bestU = u;
    }
  }

  constexpr int maxIterations = 50;
  constexpr double tolerance = 1e-10;
  double u = bestU;

  for (int i = 0; i < maxIterations; ++i) {
    const double d1 = derivative(u);
    const double d2 = secondDerivative(u);

    if (std::abs(d2) < 1e-15)
      break;

    const double delta = -d1 / d2;
    u += delta;

    if (std::abs(delta) < tolerance)
      break;
  }

  const TVector3 closestPoint = helixPoint(u);
  result.point = center + closestPoint;

  const double u1 = u + phi0;
  result.tangent.SetXYZ(-radius * std::sin(u1), radius * std::cos(u1), tanAR);

  result.tangent = result.tangent.Unit();

  return result;
}

} // namespace genfit