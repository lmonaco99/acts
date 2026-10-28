// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <format>

namespace Acts::Experimental::detail {

/// Convert an angle to degrees (for printouts)
constexpr double inDeg(double angle) {
  return angle / UnitConstants::degree;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
PatternStateAux<Bucket_t, Sector_t, nGroups>::PatternStateAux(
    const OrderedHit& seed, const typename Sector_t::Index_t sectorIdx,
    const Config& config,
    const Logger& logger)
    : m_cfg{std::addressof(config)},
      m_logger{std::addressof(logger)},
      lastInsertedHit{seed},
      prevLayerHit{seed},
      lineAnchorHit{seed},
      patTheta{VectorHelpers::theta(seed->globalPosition())},
      sector{sectorIdx} {
  /// Add the new hit
  hitsPerGroup[seed->groupIdx()].push_back(seed);
  /// Initialize the bending plane
  setBendPlaneFromSector();
  /// Update the hit counts
  if (seed->isPrecision()) {
    nPrecisionLayers++;
  } else {
    nTriggerLayers++;
  }
  if (seed->spacePoint().measuresLoc0()) {
    setStateVars({PatternState::eBendPlane});
  }
  setStateVars({PatternState::eLine});
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::addHit(
    const OrderedHit& hit, const LineTestRes& testRes) {
  /// Add the new hit
  hitsPerGroup[hit->groupIdx()].push_back(hit);
  /// Update the hit counts
  if (hit->isPrecision()) {
    nPrecisionLayers++;
  } else {
    nTriggerLayers++;
  }
  if (hit->spacePoint().measuresLoc0()) {
    setStateVars({PatternState::eBendPlane});
  }
  /// Update the pointers to previous layer hit
  const bool isNewSGroup{hit->groupIdx() !=
                         lastInsertedHit->groupIdx()};
  prevLayerHit = lastInsertedHit;
  lastInsertedHit = hit;
  /// Update the residual
  sumNormResidual2 += Acts::square(testRes.residual / testRes.sigma);
  nDoF += testRes.nDof;
  lastHitTestRes = testRes;
  /// If the new compatible hit is in a different group, update the line anchor
  if (isNewSGroup) {
    setStateVars({PatternState::eLineAnchor});
  }
  setStateVars({PatternState::eLine});
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::overWriteHit(
    const OrderedHit& newHit, const LineTestRes& newTestRes) {
  const auto group{newHit->groupIdx()};
  if (group != lastInsertedHit->groupIdx() ||
      lastInsertedHit.globLayer != newHit.globLayer) {
    throw std::runtime_error(std::format(
        "Trying to overwrite a hit in group/layer {}/{} with another one from "
        "group/layer {}/{}",
        static_cast<int>(lastInsertedHit->groupIdx()),
        lastInsertedHit.globLayer, static_cast<int>(group), newHit.globLayer));
  }
  /// We expect to overwrite hits of the same type (precision/trigger), since we
  /// only branch when we have compatible hits in the same layer, except for
  /// sTGC hits, where we have pad and strips in the same layer
  if (lastInsertedHit->isPrecision() != newHit->isPrecision()) {
    if (newHit->isPrecision()) {
      nPrecisionLayers++;
      nTriggerLayers--;
    } else {
      nPrecisionLayers--;
      nTriggerLayers++;
    }
  }
  /// Update the phi counts
  if (lastInsertedHit->spacePoint().measuresLoc0() ||
      newHit->spacePoint().measuresLoc0()) {
    setStateVars({PatternState::eBendPlane});
  }
  /// Update the residual
  sumNormResidual2 +=
      Acts::square(newTestRes.residual / newTestRes.sigma) -
      Acts::square(lastHitTestRes.residual / lastHitTestRes.sigma);
  nDoF += newTestRes.nDof - lastHitTestRes.nDof;
  lastHitTestRes = newTestRes;
  /// Remove the last inserted hit
  auto& stHits{hitsPerGroup[group]};
  if (stHits.back() != lastInsertedHit) {
    std::stringstream ss{};
    ss << "Trying to overwrite a hit that is not the last inserted hit in "
          "group/layer "
       << static_cast<int>(group) << "/" << lastInsertedHit.globLayer << "\n";
    ss << "Last inserted hit: " << Acts::toString(lastInsertedHit->spacePoint())
       << "\n";
    ss << "Last hit in group: " << Acts::toString(stHits.back()->spacePoint());
    throw std::runtime_error(ss.str());
  }
  stHits.pop_back();
  /// Add the new hit
  stHits.push_back(newHit);
  lastInsertedHit = newHit;
  setStateVars({PatternState::eLine});
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
double PatternStateAux<Bucket_t, Sector_t, nGroups>::covarianceContraction(
    const Hit_t& hit,
    const Vector3& contractionVector, bool isProjected) {
  const typename Hit_t::SpacePoint_t& sp {hit.spacePoint()};
  using CovIdx = CompSpacePointAuxiliaries::ResidualIdx;
  if (sp.isStraw()) {
    /// For straw hits, the math simplifies according to whether the hit is
    /// projected or not.
    const double discVar{Acts::square(sp.driftRadius()) +
                         sp.covariance()[Acts::toUnderlying(CovIdx::bending)]};
    if (isProjected) {
      /// If the hit is projected, the contraction vector is J^T *
      /// residualDirection, where J is the Jacobian of the projection. The phi
      /// measurement if available is cancelled by the projection.
      return discVar * contractionVector.squaredNorm();
    }
    const double vDotRsq{
        Acts::square(hit.globalSensorDirection().dot(contractionVector))};
    return discVar * (contractionVector.squaredNorm() - vDotRsq) +
           vDotRsq * sp.covariance()[Acts::toUnderlying(CovIdx::nonBending)];
  }
  /// Strip measurements: by definition, the globalSensorDirection is along the
  /// primary (eta) strip, while the globalToNextSensor is either along the secondary 
  /// (phi) strip for combined measurements or orthogonal to the primary strip for 
  /// single-strip measurements.
  const bool phiOnly{!sp.measuresLoc1()};
  const Vector3 etaStripDir{phiOnly ? hit.globalToNextSensor() : hit.globalSensorDirection()};
  const Vector3 phiStripDir{phiOnly ? hit.globalSensorDirection() : hit.globalToNextSensor()};
  /// An error of one measurement moves the point along the other strip, by the error
  /// divided by the sine of the angle between the strips (1 for orthogonal strips).
  /// Only squares enter, so 1/sin^2 = 1/(1 - cos^2) of the angle between the strips is enough
  const double invSin2{1. / (1. - Acts::square(etaStripDir.dot(phiStripDir)))};
  /// Helper method to compute the contribution of a 1D measurement to the
  /// contracted covariance, given the strip along which it moves the point
  auto oneDimContribution = [&](CovIdx idx, const Vector3& stripDir) -> double {
    return invSin2 * sp.covariance()[Acts::toUnderlying(idx)] *
           Acts::square(stripDir.dot(contractionVector));
  };
  /// The projection runs along the sensor direction, i.e. the eta strip (phi-only
  /// hits are never projected): the contribution along it is cancelled
  if (isProjected) {
    assert(!phiOnly);
    return oneDimContribution(CovIdx::bending, phiStripDir);
  }
  return oneDimContribution(CovIdx::bending, phiStripDir) + 
         oneDimContribution(CovIdx::nonBending, etaStripDir);
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
typename PatternStateAux<Bucket_t, Sector_t, nGroups>::LineTestRes
PatternStateAux<Bucket_t, Sector_t, nGroups>::computeLineResidual(
    const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const {
  LineTestRes res{};
  /// We project the test hit onto the phi plane only when the test hit does
  ///  not measure phi or when we have no phi layers, otherwise we do not
  ///  project so the residual includes the error in the phi direction.
  const bool projectTestHit{!testHit->spacePoint().measuresLoc0() ||
                            nPhiLayers == 0u};
  const Vector3 testPos{projectTestHit ? projToPhiPlane(*testHit)
                                       : testHit->globalPosition()};
  const Vector3 K{testPos - linePos};
  const double KdotD{K.dot(lineDir)};

  const Vector3 R{K - KdotD * lineDir};  // Residual vector
  res.residual = R.norm();
  if (res.residual < Acts::s_epsilon) {
    /// If the residual is too small, it's likely due to bad topology, reject.
    res.residual = std::numeric_limits<double>::max();
    return res;
  }
  const Vector3 resDir{R / res.residual};  // Residual direction
  /// Alpha represents the extrapolation distance along the pattern line
  const double alpha{KdotD / leverArm};
  /// Unit vector along the bending plane
  const Vector3 eR{bendPlaneDir()};

  /// Accumulate the derivatives of the scalar residual with respect to the
  /// common phi-plane angle.
  double phiPlaneDerivativeAcc{0.};
  /// Accumulate the derivatives of the scalar residual with respect to the
  /// common phi-plane displacement.
  double phiDisplDerivativeAcc{0.};
  /// Accumulate the covariance contributions of the hits to the residual
  double residualCovAcc{0.};

  /// @brief Accumulate the covariance contributions of one hit to the residual.
  ///         Each hit contributes with its intrinsic covariance and, if
  /// projected, with the uncertainty in the common phi-plane angle of the
  /// pattern.
  ///
  ///  Intrinsic covariance: for a projected hit, the residual direction is
  ///  transformed with the projection jacobian J^T * dir.
  ///  dir^T * J * cov * J^T * dir = (J^T * dir)^T * cov * (J^T * dir)
  ///
  ///  Phi plane uncertainty: for a projected hit, the derivative of the scalar
  ///  residual with respect to the common phi-plane angle for the projected
  /// hits.
  ///
  ///  @param hit: the hit for which to compute the covariance
  ///  @param pos: projected position (needed for the phi-plane derivative of projected hits)
  ///  @param preFactor: factor, function of alpha, to scale the covariance contributions
  ///  @param isProjected: flag indicating if the hit is projected
  auto covarianceTerm = [&](const Hit_t& hit, const Vector3& pos,
                            const double preFactor, bool isProjected) -> void {
    if (!isProjected) {
      residualCovAcc +=
          Acts::square(preFactor) *
          covarianceContraction(hit, resDir, /*isProjected=*/false);
      return;
    }
    const Vector3 sensorDir{hit.globalSensorDirection()};
    const double projFactor{sensorDir.dot(resDir) /
                            sensorDir.dot(bendPlaneNorm)};
    const Vector3 trfDir{resDir - projFactor * bendPlaneNorm};

    residualCovAcc +=
        Acts::square(preFactor) *
        covarianceContraction(hit, trfDir, /*isProjected=*/true);
    phiPlaneDerivativeAcc += preFactor * projFactor * pos.dot(eR);
    phiDisplDerivativeAcc += preFactor * projFactor;
  };

  /// Compute the covariance contributions of the first line point
  if (useBeamspot) {
    const double covS1 =
        Acts::square(beamSpot.length) * Acts::square(resDir.z()) +
        Acts::square(beamSpot.radius) * (1 - Acts::square(resDir.z()));
    residualCovAcc += Acts::square(alpha - 1) * covS1;
  } else {
    covarianceTerm(*lineAnchorHit, linePos, alpha - 1, /*isProjected=*/true);
  }

  /// Compute the covariance contributions of the second line point
  const Vector3 pos2{linePos + leverArm * lineDir};
  covarianceTerm(*lineEndHit, pos2, -alpha, /*isProjected=*/true);

  /// Compute the covariance contributions of the test hit
  covarianceTerm(*testHit, testPos, 1., projectTestHit);

  /// Contribution to the residual given by the projection
  const double projSigmaTerm{
      Acts::square(phiPlaneDerivativeAcc) *
          bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] +
      Acts::square(phiDisplDerivativeAcc) *
          bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eOffsetOffset)] +
      2 * phiPlaneDerivativeAcc * phiDisplDerivativeAcc *
          bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiOffset)]};

  res.sigma = std::sqrt(residualCovAcc + projSigmaTerm);
  res.nDof = projectTestHit ? 1u : 2u;

  ACTS_VERBOSE(__func__ << "() Use beamspot: " << useBeamspot << ", alpha: "
                        << alpha << ", Residual: " << res.residual << " +- "
                        << res.sigma << ", nDof: " << static_cast<int>(res.nDof)
                        << ", pull: " << res.residual / res.sigma << std::endl
                        << "linePos: " << toString(linePos)
                        << ", lineDir: " << toString(lineDir) << ", testPos: "
                        << toString(testPos) << ", resDir: " << toString(resDir)
                        << ", hit pos sigma: " << std::sqrt(residualCovAcc)
                        << ", phi plane sigma: " << std::sqrt(projSigmaTerm));
  return res;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
Vector3 PatternStateAux<Bucket_t, Sector_t, nGroups>::projToPhiPlane(
    const Hit_t& hit) const {
  return Acts::PlanarHelper::intersectPlane(
             hit.globalPosition(), hit.globalSensorDirection(),
             bendPlaneNorm, patPhiOffset * bendPlaneNorm)
      .position();
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::updateLineParameters(
    const BeamspotInfo& beamSpot) {
  needUpdate.reset(toUnderlying(PatternState::eLine));
  lineEndHit = lastInsertedHit;
  Vector3 pos1{projToPhiPlane(*lineAnchorHit)};
  Vector3 pos2{projToPhiPlane(*lineEndHit)};
  Vector3 d{pos2 - pos1};
  leverArm = d.norm();

  /// Check whether we have to use the beamspot instead of the anchor hit to
  /// draw the line.
  useBeamspot = leverArm < m_cfg->minHitDistance4Line;
  if (useBeamspot) {
    linePos = beamSpot.position;
    d = pos2 - beamSpot.position;
    leverArm = d.norm();
  } else {
    linePos = pos1;
  }
  lineDir = d / leverArm;

  ACTS_VERBOSE(__func__ << "() Updated --> linePos: " << toString(linePos)
                        << ", lineDir: " << toString(lineDir) << ", LeverArm: "
                        << leverArm << ", Use beamspot: " << useBeamspot);
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::updateBendPlane(
    const BeamspotInfo& beamSpot) {
  using namespace Acts::UnitLiterals;
  needUpdate.reset(toUnderlying(PatternState::eBendPlane));
  needUpdate.set(toUnderlying(PatternState::eLine));
  /// Helper method to align the final solution with the previous one.
  auto alignPhiAndStore = [&](double newPatphi) -> void {
    const double refPhi{sector.phi()};
    while (newPatphi - refPhi > 90._degree) {
      newPatphi -= 180._degree;
    }
    while (newPatphi - refPhi < -90._degree) {
      newPatphi += 180._degree;
    }
    patPhi = newPatphi;
    bendPlaneNorm = Vector3{-std::sin(patPhi), std::cos(patPhi), 0.};
  };
  /// Count the phi hits and accumulate the regression moments. The only approx
  /// is that weights are frozen at the previous plane estimate. Positions are
  /// taken relative to the beam spot: this avoids the cancellation in the
  /// centred moments and puts the beam spot at the origin, so it only adds its
  /// weight. We exploit that the sectorization is in phi, so everything is
  /// computed in the transverse plane.
  const Vector2 bsPos{beamSpot.position.template head<2>()};
  nPhiLayers = 0u;
  double weightSum{0.}, Sxx{0.}, Sxy{0.}, Syy{0.};
  Vector2 centroid{Vector2::Zero()};  // Sx, Sy

  auto accumulate = [&](const Hit_t& hit, const Vector3& pos) -> void {
    if (!hit.spacePoint().measuresLoc0()) {
      return;
    }
    nPhiLayers++;
    const double variance{
        covarianceContraction(hit, bendPlaneNorm, /*isProjected=*/false)};
    if (variance < Acts::s_epsilon) {
      ACTS_ERROR("updateBendPlane() Phi hit with zero variance, hit:"
                 << std::endl
                 << Acts::toString(hit.spacePoint()));
      return;
    }
    double weight{1. / variance};
    const Vector2 d{pos.template head<2>() - bsPos};
    const Vector2 wd{weight * d};

    weightSum += weight;
    centroid += wd;
    Sxx += wd.x() * d.x();
    Sxy += wd.x() * d.y();
    Syy += wd.y() * d.y();
  };
  std::ranges::for_each(hitsPerGroup, [&](const std::vector<OrderedHit>& hits) {
    std::ranges::for_each(hits, [&](const OrderedHit& hit) {
      /// Hits measuring eta and phi are already constrained, so we use directly
      /// the global position
      accumulate(*hit, hit->globalPosition());
    });
  });
  std::ranges::for_each(phiOnlyHits, [&](const PhiOnlyHit& phiHit) {
    /// Phi-only hits are not constrained along the bending coordinate,
    /// so we use the position on the strip crossed by the pattern line.
    accumulate(phiHit.hit, phiHit.crossPos);
  });

  if (nPhiLayers == 0) {
    setBendPlaneFromSector();
    return;
  }

  /// Beam spot as an additional point measurement, sitting at the origin of the
  /// shifted frame
  weightSum += 1. / Acts::square(beamSpot.radius);

  if (weightSum < Acts::s_epsilon) {
    throw std::runtime_error(
        std::format("Unexpected to have zero total weight when updating "
                    "pattern phi with {} phi hits and beamspot",
                    nPhiLayers));
  }
  const double weightSumInv{1. / weightSum};
  centroid *= weightSumInv;

  /// Weighted centered scatter matrix (translation invariant).
  const double Mxx{Sxx - weightSum * Acts::square(centroid.x())};
  const double Mxy{Sxy - weightSum * centroid.x() * centroid.y()};
  const double Myy{Syy - weightSum * Acts::square(centroid.y())};

  const double eigenvalSeparation{Acts::fastHypot(Mxx - Myy, 2. * Mxy)};
  if (eigenvalSeparation < Acts::s_epsilon) {
    throw std::runtime_error(
        std::format("Unexpected to have zero eigenvalue separation when "
                    "updating pattern phi with {} phi hits and beamspot",
                    nPhiLayers));
  }
  const double eigenValSeparationInv{1. / eigenvalSeparation};

  alignPhiAndStore(0.5 * std::atan2(2. * Mxy, Mxx - Myy));

  /// Back to the original frame
  centroid += bsPos;
  patPhiOffset = bendPlaneNorm.template head<2>().dot(centroid);
  /// Centroid coordinate along the line
  const Vector2 eR{bendPlaneDir().template head<2>()};
  const double rhoCentroid{centroid.dot(eR)};

  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] =
      eigenValSeparationInv;
  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiOffset)] =
      -rhoCentroid * eigenValSeparationInv;
  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eOffsetOffset)] =
      Acts::square(rhoCentroid) * eigenValSeparationInv + weightSumInv;

  ACTS_VERBOSE(__func__ << "() Updated pattern phi " << inDeg(patPhi)
                        << ", offset " << patPhiOffset << ", BendPlaneCov "
                        << print(bendPlaneCov) << ", with "
                        << static_cast<int>(nPhiLayers)
                        << " phi hits and beamspot");
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::setBendPlaneFromSector() {
  /// If there are no phi hits, we just use the central phi of the sector,
  ///  with a standard deviation based on the expanded sector size. The plane
  /// passes through the origin by construction, so s = 0 exactly.
  patPhi = sector.phi();
  bendPlaneNorm = Vector3{-std::sin(patPhi), std::cos(patPhi), 0.};
  patPhiOffset = 0.;
  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] =
      Acts::square(sector.sectorSize()) / 3.;
  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eOffsetOffset)] = 0.;
  bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiOffset)] = 0.;
  ACTS_VERBOSE(__func__ << "() No phi hits in the pattern, set pattern phi "
                        << inDeg(patPhi) << ", offset " << patPhiOffset
                        << ", BendPlaneCov " << print(bendPlaneCov));
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
double PatternStateAux<Bucket_t, Sector_t, nGroups>::angleToBendPlane(
    const Vector3& pos) const {
  /// A point X has distance = n.X - s and signed coordinate along the plane
  /// rho  = e_r.X, with e_r = (n_y, -n_x). The atan of these gives the angle
  /// between the point and the plane, as seen from the origin. The angle is in
  /// (-pi, pi]
  const double dist{bendPlaneNorm.dot(pos) - patPhiOffset};
  const double rho{pos.dot(bendPlaneDir())};
  return std::atan2(dist, rho);
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
typename PatternStateAux<Bucket_t, Sector_t, nGroups>::LineTestRes
PatternStateAux<Bucket_t, Sector_t, nGroups>::phiHitResidual(
    const Hit_t& phiHit, const Vector3& pos) const {
  /// The distance from the plane is d = n.X - s. A change of the plane
  /// parameters moves it by -rho dphi - ds, with rho = e_r.X the coordinate of
  /// the point along the plane
  const double rho{pos.dot(bendPlaneDir())};
  const double planeVar{
      Acts::square(rho) * bendPlaneCov[toUnderlying(BendPlaneCov::ePhiPhi)] +
      2. * rho * bendPlaneCov[toUnderlying(BendPlaneCov::ePhiOffset)] +
      bendPlaneCov[toUnderlying(BendPlaneCov::eOffsetOffset)]};
  LineTestRes res{};
  res.residual = bendPlaneNorm.dot(pos) - patPhiOffset;
  res.sigma = std::sqrt(
      covarianceContraction(phiHit, bendPlaneNorm, /*isProjected=*/false) +
      planeVar);
  res.nDof = 1u;
  return res;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
double PatternStateAux<Bucket_t, Sector_t, nGroups>::phiAtRadius(
    const double radius) const {
  /// Approximate azimuth of the pattern's line at transverse distance R from
  /// the beam line (phi + atan(s/R) ~ phi + s/R
  return patPhi + patPhiOffset / radius;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
Vector3 PatternStateAux<Bucket_t, Sector_t, nGroups>::bendPlaneDir() const {
  /// Here we heavily exploit the fact that the sectorization is in phi, so
  /// the normal and the direction of the bending plane are always in the
  /// transverse plane.
  return Vector3{bendPlaneNorm.y(), -bendPlaneNorm.x(), 0.};
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
bool PatternStateAux<Bucket_t, Sector_t, nGroups>::defineLineAroundHit(
    const OrderedHit& refHit) {
  const auto& groupHits{hitsPerGroup[refHit->groupIdx()]};

  /// We use useBeamspot as a flag to indicate whether the pattern line
  ///  has been determined with enough lever arm
  useBeamspot = true;

  if (groupHits.size() > 1) {
    /// if we have >= 2 eta hits in the group, we use the furthest ones to
    /// define the pattern line
    const auto [minIt, maxIt]{std::ranges::minmax_element(
        groupHits, {}, [](const OrderedHit& c) { return c.globLayer; })};
    lineAnchorHit = *minIt;
    lastInsertedHit = *maxIt;
    updateLineParameters(BeamspotInfo{});
  }
  if (useBeamspot && !groupHits.empty()) {
    /// if we have only one eta hit or the layer separation is too small, to
    /// find the second hit we use the functionality of anchor hit
    resetLineAnchorHit(groupHits.front());
    lastInsertedHit =
        *std::ranges::max_element(groupHits, {}, [&](const OrderedHit& c) {
          return layerDistance(c.globLayer, lineAnchorHit.globLayer);
        });
    updateLineParameters(BeamspotInfo{});
  }
  if (useBeamspot) {
    /// Define the closest hits in the upper and lower groups
    const OrderedHit *lowerHit{nullptr}, *upperHit{nullptr};
    for (const auto& [group, hits] : Acts::enumerate(hitsPerGroup)) {
      if (hits.empty() || group == refHit->groupIdx()) {
        continue;
      }
      assert(hits.front().globLayer != refHit.globLayer);
      if (hits.front().globLayer < refHit.globLayer) {
        const OrderedHit& lowerInGroup{*std::ranges::max_element(
            hits, {}, [](const OrderedHit& c) { return c.globLayer; })};
        if (lowerHit == nullptr ||
            lowerInGroup.globLayer > lowerHit->globLayer) {
          lowerHit = std::addressof(lowerInGroup);
        }
      } else {
        const OrderedHit& upperInGroup{*std::ranges::min_element(
            hits, {}, [](const OrderedHit& c) { return c.globLayer; })};
        if (upperHit == nullptr ||
            upperInGroup.globLayer < upperHit->globLayer) {
          upperHit = std::addressof(upperInGroup);
        }
      }
    }
    if (lowerHit == nullptr || upperHit == nullptr) {
      return false;
    }
    lineAnchorHit = *lowerHit;
    lastInsertedHit = *upperHit;
    updateLineParameters(BeamspotInfo{});
  }
  return !useBeamspot;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
std::pair<Vector3, double>
PatternStateAux<Bucket_t, Sector_t, nGroups>::phiStripCrossing(
    const Hit_t& phiHit) const {
  using CovIdx = CompSpacePointAuxiliaries::ResidualIdx;
  const Vector3 pos{phiHit.globalPosition()};
  /// The half length of a 1D phi strip is stored as the covariance of the
  /// unmeasured (bending) coordinate
  const double halfLength{std::sqrt(
      phiHit.spacePoint().covariance()[toUnderlying(CovIdx::bending)])};
  const Vector3 sensorDir{phiHit.globalSensorDirection()};
  auto projOntoPhiPlane = [this](const Vector3& p) -> Vector3 {
    return p - (p.dot(bendPlaneNorm) - patPhiOffset) * bendPlaneNorm;
  };
  const Vector3 stripLow{projOntoPhiPlane(pos - halfLength * sensorDir)};
  const Vector3 stripHigh{projOntoPhiPlane(pos + halfLength * sensorDir)};
  const double projLength{(stripHigh - stripLow).norm()};
  if (projLength < Acts::s_epsilon) {
    return {pos, 0.};
  }
  const auto intersection{Acts::detail::LineHelper::lineIntersect<3>(
      linePos, lineDir, stripLow,
      Vector3{(stripHigh - stripLow) / projLength})};
  /// The pattern line is parallel to the projected strip: not compatible
  if (!intersection.isValid()) {
    return {pos, std::numeric_limits<double>::infinity()};
  }
  const double intersect{intersection.pathLength()};
  /// The projection onto the bending plane is affine along the strip, so the
  /// fraction of the projected strip length maps directly onto the strip
  return {pos + (intersect / projLength - 0.5) * 2. * halfLength * sensorDir,
          std::max(-intersect, intersect - projLength)};
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::resetLineAnchorHit(
    const OrderedHit& refHit) {
  needUpdate.reset(toUnderlying(PatternState::eLineAnchor));
  needUpdate.set(toUnderlying(PatternState::eLine));
  /// Treat first the special case where we have only one group
  if (countGroups() < 2) {
    const auto& hits{*std::ranges::find_if(
        hitsPerGroup, [](const auto& h) { return !h.empty(); })};
    lineAnchorHit = *std::ranges::max_element(
        hits, std::ranges::less{}, [&refHit](const OrderedHit& hit) {
          return layerDistance(hit.globLayer, refHit.globLayer);
        });
    return;
  }
  const auto refGroup{refHit->groupIdx()};
  unsigned bestDist{std::numeric_limits<unsigned>::max()};
  for (const auto& hits : hitsPerGroup) {
    if (hits.empty() || hits.front()->groupIdx() == refGroup) {
      continue;
    }
    for (const auto& hit : hits) {
      const unsigned dist = layerDistance(hit.globLayer, refHit.globLayer);
      if (dist < bestDist) {
        bestDist = dist;
        lineAnchorHit = hit;
      }
      if (dist == 1u) {
        return;
      }
    }
  }
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::syncState(
    const BeamspotInfo& beamSpot) {
  if (needUpdate[toUnderlying(PatternState::eBendPlane)]) {
    updateBendPlane(beamSpot);
  }
  if (needUpdate[toUnderlying(PatternState::eLineAnchor)]) {
    /// eLineAnchor is only set by addHit when a new group is opened: the
    /// closest hit in another group is the previously added one
    lineAnchorHit = prevLayerHit;
    needUpdate.reset(toUnderlying(PatternState::eLineAnchor));
    needUpdate.set(toUnderlying(PatternState::eLine));
  }
  /// Searching in a group not yet reached: extrapolate with the line spanned by
  /// the group of the last inserted hit (anchor = furthest hit from it). Once a
  /// hit in the new group is accepted, addHit flags eLineAnchor and we return
  /// to the default anchor.
  const auto lastHitGroup{lastInsertedHit->groupIdx()};
  if (lastHitGroup != lineAnchorHit->groupIdx()) {
    OrderedHit furthest{*std::ranges::max_element(
        hitsPerGroup[lastHitGroup], {}, [this](const OrderedHit& h) {
          return layerDistance(h.globLayer, lastInsertedHit.globLayer);
        })};
    if (lastInsertedHit != furthest &&
        (projToPhiPlane(*lastInsertedHit) -
         projToPhiPlane(*furthest))
                .squaredNorm() >= square(m_cfg->minHitDistance4Line)) {
      lineAnchorHit = furthest;
      needUpdate.set(toUnderlying(PatternState::eLine));
    }
  }
  if (needUpdate[toUnderlying(PatternState::eLine)]) {
    updateLineParameters(beamSpot);
  }
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::setStateVars(
    const std::initializer_list<PatternState>& vars) {
  for (const PatternState var : vars) {
    needUpdate.set(toUnderlying(var));
  }
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
bool PatternStateAux<Bucket_t, Sector_t, nGroups>::isInPattern(
    const Hit_t& hit) const {
  const auto& hits{hitsPerGroup[hit.groupIdx()]};
  return std::ranges::find(hits, std::addressof(hit), &OrderedHit::hitPtr) != hits.end();
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
typename PatternStateAux<Bucket_t, Sector_t, nGroups>::BucketRefs_t
PatternStateAux<Bucket_t, Sector_t, nGroups>::bucketRefs() const {
  BucketRefs_t bucketRefs{};
  for (const auto& [group, hits] : Acts::enumerate(hitsPerGroup)) {
    for (const OrderedHit& hit : hits) {

      const bool alreadyIn {std::ranges::any_of(bucketRefs[group], [&hit](const OrderedHit* ref) {
        return std::addressof((*ref)->parentBucket()) == std::addressof(hit->parentBucket());
      })};
      if (!alreadyIn) {
        bucketRefs[group].push_back(std::addressof(hit));
      }
    }
  }
  return bucketRefs;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
std::uint8_t PatternStateAux<Bucket_t, Sector_t, nGroups>::countGroups() const {
  return std::ranges::count_if(hitsPerGroup,
                               [](const auto& hits) { return !hits.empty(); });
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
double PatternStateAux<Bucket_t, Sector_t, nGroups>::getMeanResidual2() const {
  return sumNormResidual2 / std::max(1., static_cast<double>(nDoF));
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
std::uint8_t PatternStateAux<Bucket_t, Sector_t, nGroups>::nBendingLayers()
    const {
  return nPrecisionLayers + nTriggerLayers;
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
typename PatternStateAux<Bucket_t, Sector_t, nGroups>::PatternPrintView
PatternStateAux<Bucket_t, Sector_t, nGroups>::brief(const PatternStateAux& p) {
  return {p, /*detailed=*/false};
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
typename PatternStateAux<Bucket_t, Sector_t, nGroups>::PatternPrintView
PatternStateAux<Bucket_t, Sector_t, nGroups>::detailed(
    const PatternStateAux& p) {
  return {p, /*detailed=*/true};
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::print(std::ostream& ostr,
                                                         bool detailed) const {
  ostr << "PatternStateAux Exp Sector: " << static_cast<int>(sector.sector())
       << ", Theta: " << inDeg(patTheta) << ", Phi " << inDeg(patPhi)
       << ", phiOffset: " << patPhiOffset;
  ostr << ", nPrec: " << static_cast<int>(nPrecisionLayers)
       << ", nEtaTrig: " << static_cast<int>(nTriggerLayers)
       << ", nPhi: " << static_cast<int>(nPhiLayers);
  ostr << ", mean norma res sq: " << getMeanResidual2()
       << ", nDof: " << static_cast<int>(nDoF)
       << ", BendPlane Cov: " << print(bendPlaneCov);
  ostr << ", Hit per group: \n";
  for (std::uint8_t g = 0u; g < nGroups; ++g) {
    const auto& hits{hitsPerGroup[g]};
    if (hits.empty()) {
      continue;
    }
    ostr << "  Group " << static_cast<int>(g) << " has " << hits.size()
         << " hits ";
    if (detailed) {
      ostr << "\n";
      for (const auto& hit : hits) {
        ostr << "    " << Acts::toString(hit->spacePoint())
             << ", globLay: " << static_cast<int>(hit.globLayer) << "\n";
      }
    }
  }
  if (!detailed) {
    ostr << "\n";
    ostr << "  Last inserted hit: " << lastInsertedHit << "\n";
    ostr << "  Previous layer hit: " << prevLayerHit << "\n";
    ostr << "  Line anchor hit: " << lineAnchorHit << "\n";
  }
}

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
void PatternStateAux<Bucket_t, Sector_t, nGroups>::OrderedHit::print(
    std::ostream& ostr) const {
  ostr << Acts::toString(hitPtr->spacePoint())
       << ", group: " << static_cast<int>(hitPtr->groupIdx())
       << ", globLay: " << static_cast<int>(globLayer);
}

}  // namespace Acts::Experimental::detail
