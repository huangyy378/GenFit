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

// ClassImp(genfit::HelixMeasurement)

namespace genfit {

HelixMeasurement::HelixMeasurement(int nDim) :
    AbsMeasurement(nDim),
    maxDistance_(2.0),
    leftRight_(0) {
    assert(nDim == 7); // 必须为7维
}

HelixMeasurement::HelixMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov,
                                   int detId, int hitId, TrackPoint* trackPoint) :
    AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint),
    maxDistance_(2.0),
    leftRight_(0) {
    // 验证输入维度
    if (rawHitCoords_.GetNrows() != 7) {
        throw Exception("HelixMeasurement requires 7-dimensional rawHitCoords", __LINE__, __FILE__);
    }
}

SharedPlanePtr HelixMeasurement::constructPlane(const StateOnPlane& state) const {
    // 获取当前状态的位置
    const AbsTrackRep* rep = state.getRep();
    TVector3 currentPos = rep->getPos(state);

    // 找到螺旋线上距离当前点最近的点及切线方向
    auto closest = findClosestPointOnHelix(currentPos);
    const TVector3& pocaOnHelix = closest.point;
    TVector3 tangent = closest.tangent;
    tangent.SetMag(1.0); // 确保单位向量

    // 获取轨迹在当前位置的方向
    TVector3 dirInPoca = rep->getMom(state);
    dirInPoca.SetMag(1.0);

    // 检查方向是否与切线平行
    if (fabs(tangent.Angle(dirInPoca)) < 0.01) {
        throw Exception(
            "HelixMeasurement::constructPlane: Direction is parallel to helix tangent",
            __LINE__, __FILE__);
    }

    // 构造正交向量 U = 轨迹方向 × 切线方向
    TVector3 U = dirInPoca.Cross(tangent);
    if (U.Mag() < 1e-10) {
        // 防止零向量：使用备选方法
        U = TVector3(1, 0, 0).Cross(tangent);
        if (U.Mag() < 1e-10) {
            U = TVector3(0, 1, 0).Cross(tangent);
        }
    }
    U.SetMag(1.0);

    // 创建平面：原点为螺旋线上最近点，U轴为法向，V轴为切线方向
    return SharedPlanePtr(new DetPlane(pocaOnHelix, U, tangent));
}

std::vector<MeasurementOnPlane*> HelixMeasurement::constructMeasurementsOnPlane(
    const StateOnPlane& state) const {
    double mR = rawHitCoords_(6); // 漂移距离
    double mL = -mR;
    double V = rawHitCov_(6, 6); // 漂移距离方差

    // 创建左右测量
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

    // 设置左右权重
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

    // 1. 获取螺旋线参数
    const TVector3 center(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));
    const double radius = rawHitCoords_(3);
    const double pitch = rawHitCoords_(4);
    const double phi0 = rawHitCoords_(5);
    const TVector3 relPos = point - center;
    const double k = pitch / (2 * TMath::Pi());

    // 2. 定义连续的目标函数（无周期性折叠）
    auto distanceSq = [&](double phi) -> double {
        const double deltaPhi = phi - phi0;
        const double x = radius * std::cos(deltaPhi);
        const double y = radius * std::sin(deltaPhi);
        const double z = k * phi; // 连续z坐标

        return (x - relPos.X()) * (x - relPos.X()) +
               (y - relPos.Y()) * (y - relPos.Y()) +
               (z - relPos.Z()) * (z - relPos.Z());
    };

    // 3. 定义连续导数
    auto derivative = [&](double phi) -> double {
        const double deltaPhi = phi - phi0;
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

    // 4. 定义连续二阶导数
    auto secondDerivative = [&](double phi) -> double {
        const double deltaPhi = phi - phi0;
        const double x = radius * std::cos(deltaPhi);
        const double y = radius * std::sin(deltaPhi);

        const double dx_dphi = -radius * std::sin(deltaPhi);
        const double dy_dphi = radius * std::cos(deltaPhi);

        const double d2x_dphi2 = -radius * std::cos(deltaPhi);
        const double d2y_dphi2 = -radius * std::sin(deltaPhi);

        return 2 * (dx_dphi * dx_dphi + (x - relPos.X()) * d2x_dphi2) +
               2 * (dy_dphi * dy_dphi + (y - relPos.Y()) * d2y_dphi2) +
               2 * k * k; // z分量二阶导为0，但一阶导平方项保留
    };

    // 5. 粗搜索初始化（关键改进）
    double bestPhi = 0.0;
    double minDistSq = std::numeric_limits<double>::max();
    const int numSamples = 100;
    const double zGuess = (std::abs(k) > 1e-9) ? relPos.Z() / k : 0.0;

    // 搜索范围：zGuess附近±5个螺距
    const double startPhi = zGuess - 5 * 2 * TMath::Pi();
    const double endPhi = zGuess + 5 * 2 * TMath::Pi();

    for (int i = 0; i <= numSamples; ++i) {
        const double phi = startPhi + i * (endPhi - startPhi) / numSamples;
        const double distSq = distanceSq(phi);
        if (distSq < minDistSq) {
            minDistSq = distSq;
            bestPhi = phi;
        }
    }

    // 6. 牛顿迭代法
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

    // 7. 计算结果（保持角度在[0, 2π)）
    const double deltaPhi = phi - phi0;
    const double normalizedPhi = deltaPhi - 2 * TMath::Pi() * std::floor(deltaPhi / (2 * TMath::Pi()));

    result.point.SetXYZ(
        center.X() + radius * std::cos(normalizedPhi),
        center.Y() + radius * std::sin(normalizedPhi),
        center.Z() + k * phi);

    result.tangent.SetXYZ(
        -radius * std::sin(normalizedPhi),
        radius * std::cos(normalizedPhi),
        k);

    // 归一化切线
    result.tangent = result.tangent.Unit();

    return result;
}

} // namespace genfit