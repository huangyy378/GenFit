/* HelixMeasurement.cpp */
#include "HelixMeasurement.h"

#include "DetPlane.h"
#include "Exception.h"
#include "HMatrixU.h"
#include "RKTrackRep.h"
#include "StateOnPlane.h"
#include "TMath.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <TClass.h>
#include <vector>

namespace genfit {

HelixMeasurement::HelixMeasurement(int nDim) :
    AbsMeasurement(nDim),
    maxDistance_(2.0),
    leftRight_(0) {
    assert(nDim == 8);
}

HelixMeasurement::HelixMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov,
                                   int detId, int hitId, TrackPoint* trackPoint) :
    AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint),
    maxDistance_(2.0),
    leftRight_(0) {

    if (rawHitCoords_.GetNrows() != 8) {
        throw Exception("HelixMeasurement requires 7-dimensional rawHitCoords", __LINE__, __FILE__);
    }
}

SharedPlanePtr HelixMeasurement::constructPlane(const StateOnPlane& state) const {
    // get state
    const AbsTrackRep* rep = state.getRep();
    TVector3 currentPos = rep->getPos(state);

    // find closest point on helix
    auto closest = findClosestPointOnHelix(currentPos);
    const TVector3& pocaOnHelix = closest.point;
    TVector3 tangent = closest.tangent;
    tangent.SetMag(1.0);

    // get track direction
    TVector3 dirInPoca = rep->getMom(state);
    dirInPoca.SetMag(1.0);

    // check for parallelism
    if (fabs(tangent.Angle(dirInPoca)) < 0.01) {
        throw Exception(
            "HelixMeasurement::constructPlane: Direction is parallel to helix tangent",
            __LINE__, __FILE__);
    }

    // construct orthogonal vector
    TVector3 U = dirInPoca.Cross(tangent);
    U.SetMag(1.0);

    // construct plane
    return SharedPlanePtr(new DetPlane(pocaOnHelix, U, tangent));
}

std::vector<MeasurementOnPlane*> HelixMeasurement::constructMeasurementsOnPlane(
    const StateOnPlane& state) const {
    double mR = rawHitCoords_(7);
    double mL = -mR;
    double V = rawHitCov_(7, 7);

    MeasurementOnPlane* mopL = new MeasurementOnPlane(
        TVectorD(1, &mL),
        TMatrixDSym(1, &V),
        state.getPlane(),
        state.getRep(),
        constructHMatrix(state.getRep()));

    MeasurementOnPlane* mopR = new MeasurementOnPlane(
        TVectorD(1, &mR),
        TMatrixDSym(1, &V),
        state.getPlane(),
        state.getRep(),
        constructHMatrix(state.getRep()));

    // set weights
    if (leftRight_ < 0) {
        mopL->setWeight(1);
        mopR->setWeight(0);
    } else if (leftRight_ > 0) {
        mopL->setWeight(0);
        mopR->setWeight(1);
    } else {
        double val = 0.5 * pow(std::max(0., 1 - mR / maxDistance_), 2.);
        mopL->setWeight(val);
        mopR->setWeight(val);
    }

    std::vector<MeasurementOnPlane*> retVal;
    retVal.push_back(mopL);
    retVal.push_back(mopR);
    return retVal;
}

const AbsHMatrix* HelixMeasurement::constructHMatrix(const AbsTrackRep* rep) const {
    if (dynamic_cast<const RKTrackRep*>(rep) == nullptr) {
        throw Exception("HelixMeasurement can only handle state vectors of type RKTrackRep!",
                        __LINE__, __FILE__);
    }
    return new HMatrixU();
}

void HelixMeasurement::setLeftRightResolution(int lr) {
    if (lr == 0)
        leftRight_ = 0;
    else if (lr < 0)
        leftRight_ = -1;
    else
        leftRight_ = 1;
}

HelixMeasurement::ClosestPointResult
HelixMeasurement::findClosestPointOnHelix(const TVector3& point) const {
    ClosestPointResult result;

    // get helix parameters
    const TVector3 center(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));
    const double radius = rawHitCoords_(3);
    const double pitch = rawHitCoords_(4);
    const double phi0 = rawHitCoords_(5);
    const double phiTotal = rawHitCoords_(6);
    double phi1 = (phi0 + phiTotal) / 2;
    const TVector3 relPos = point - center;
    const double k = pitch / (2 * TMath::Pi());

    // construct functions
    auto distanceSq = [&](double phi) -> double {
        const double deltaPhi = phi - phi1;
        const double x = radius * std::cos(deltaPhi);
        const double y = radius * std::sin(deltaPhi);
        const double z = k * phi;

        return (x - relPos.X()) * (x - relPos.X()) +
               (y - relPos.Y()) * (y - relPos.Y()) +
               (z - relPos.Z()) * (z - relPos.Z());
    };

    auto derivative = [&](double phi) -> double {
        const double deltaPhi = phi - phi1;
        const double x = radius * std::cos(deltaPhi);
        const double y = radius * std::sin(deltaPhi);
        const double z = k * phi;

        const double dx_dphi = -radius * std::sin(deltaPhi);
        const double dy_dphi = radius * std::cos(deltaPhi);
        const double dz_dphi = k;

        return 2 * (x - relPos.X()) * dx_dphi +
               2 * (y - relPos.Y()) * dy_dphi +
               2 * (z - relPos.Z()) * dz_dphi;
    };

    // second derivative
    auto secondDerivative = [&](double phi) -> double {
        const double deltaPhi = phi - phi1;
        const double x = radius * std::cos(deltaPhi);
        const double y = radius * std::sin(deltaPhi);

        const double dx_dphi = -radius * std::sin(deltaPhi);
        const double dy_dphi = radius * std::cos(deltaPhi);

        const double d2x_dphi2 = -radius * std::cos(deltaPhi);
        const double d2y_dphi2 = -radius * std::sin(deltaPhi);

        return 2 * (dx_dphi * dx_dphi + (x - relPos.X()) * d2x_dphi2) +
               2 * (dy_dphi * dy_dphi + (y - relPos.Y()) * d2y_dphi2) +
               2 * k * k;
    };

    double bestPhi = 0.0;
    double minDistSq = std::numeric_limits<double>::max();
    const int numSamples = 100;
    const double zGuess = (std::abs(k) > 1e-9) ? relPos.Z() / k : 0.0;

    // find good starting point
    const double startPhi = phi0;
    const double endPhi = phi0 + phiTotal;

    for (int i = 0; i <= numSamples; ++i) {
        const double phi = startPhi + i * (endPhi - startPhi) / numSamples;
        const double distSq = distanceSq(phi);
        if (distSq < minDistSq) {
            minDistSq = distSq;
            bestPhi = phi;
        }
    }

    // refine using Newton's method
    constexpr int maxIterations = 50;
    constexpr double tolerance = 1e-10;
    double phi = bestPhi;

    for (int i = 0; i < maxIterations; ++i) {
        const double d1 = derivative(phi);
        const double d2 = secondDerivative(phi);

        if (std::abs(d2) < 1e-15) break;

        const double delta = -d1 / d2;
        phi += delta;

        if (std::abs(delta) < tolerance) break;
    }

    const double deltaPhi = phi - phi1;
    const double normalizedPhi = deltaPhi - 2 * TMath::Pi() * std::floor(deltaPhi / (2 * TMath::Pi()));

    result.point.SetXYZ(
        center.X() + radius * std::cos(normalizedPhi),
        center.Y() + radius * std::sin(normalizedPhi),
        center.Z() + k * phi);

    result.tangent.SetXYZ(
        -radius * std::sin(normalizedPhi),
        radius * std::cos(normalizedPhi),
        k);

    // unitize tangent
    result.tangent = result.tangent.Unit();

    return result;
}

} // namespace genfit
