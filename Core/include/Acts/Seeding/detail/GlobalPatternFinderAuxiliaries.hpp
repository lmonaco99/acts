// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <bitset>
#include <cstdint>
#include <format>
#include <sstream>

/// Auxiliary data structures and concepts used by the Global Pattern Finder.
///
/// This file defines the hit object used by the Global Pattern Finder, the
/// concept used to abstract the detector-specific sectorization, together
/// with the main PatternStateAux object used internally by the Global Pattern
/// Finder.

namespace Acts::Experimental::detail {

/// @brief Sectorization of the detector in the non-bending coordinate.
///  In this first implementation of Global Pattern Finder with toroidal
///  azimuthal field, the sectors are defined along the phi direction.
template <typename Sector_t>
concept SectorType =
    requires(const Sector_t& sector, const Sector_t& other, const double phi) {
      /// Type definition for the sector index type.
      typename Sector_t::Index_t;
      /// The sector type must be constructible from the sector index.
      requires std::constructible_from<Sector_t, typename Sector_t::Index_t>;
      /// Returns the sector index of the sector
      { sector.sector() } -> std::same_as<typename Sector_t::Index_t>;
      /// Returns the phi coordinate of the sector center in [rad]
      { sector.phi() } -> std::same_as<double>;
      /// Returns the angular size of the sector in [rad]
      { sector.sectorSize() } -> std::same_as<double>;
      /// Returns whether the sector is a neighbour of another sector.
      /// This helper incorporates the wrap-around of the sectorization.
      { sector.isNeighbour(other) } -> std::same_as<bool>;
      /// Returns whether a given phi coordinate is inside the sector.
      { sector.insideSector(phi) } -> std::same_as<bool>;
    };
/// @brief Pattern state object storing pattern information during construction
template <CompositeSpacePoint Hit_t, SectorType Sector_t,
          std::uint8_t nGroups>
struct PatternStateAux {
  /// @brief Configuration object of the residual calculator
  struct Config {
    /// @brief Minimum distance between two hits for being used to
    ///         compute a reliable pattern line. Use the beamspot otherwise.
    double minHitDistance4Line{200. * UnitConstants::mm};
  };
  /// @brief Beamspot information, used as an additional constraint on the
  ///  pattern line when the lever arm of the hits is too small. The beamspot
  ///  is also used to compute the pattern bending plane through linear
  ///  regression in the transverse plane. The beamspot is modelled as a
  ///  cylinder along the beam axis (z), whose radius and half-length act as
  ///  the transverse and longitudinal uncertainties of the constraint. They
  ///  should therefore be large enough to also cover the deviation of the
  ///  straight pattern line from the true trajectory in the transverse plane
  ///  close to the beam, e.g. due to the bending in the solenoidal field of the
  ///  inner detector.
  struct BeamspotInfo {
    /// @brief Beamspot position
    Vector3 position{Vector3::Zero()};
    /// @brief Beamspot radius
    double radius{30. * UnitConstants::cm};
    /// @brief Beamspot (half-)length
    double length{2. * UnitConstants::m};
  };
  /// @brief Small wrapper for hits used to build patterns. In general, due to detector
  ///         complexity, the layer number cannot be defined globally, but it
  ///         can
  /// be computed given a set of hits.
  struct OrderedHit {
    /// @brief Pointer to the underlying hit
    const Hit_t* hitPtr{nullptr};
    /// @brief Global measurement layer number
    std::uint8_t globLayer{0u};
    /// @brief Group index of the hit in the pattern topology (Config::groupIndex)
    std::uint8_t groupIdx{0u};
    /// @brief Whether the hit is a precision measurement (Config::isPrecision)
    bool isPrecision{false};
    /// @brief Access to the underlying hit methods
    const Hit_t* operator->() const { return hitPtr; }
    /// @brief Access to the underlying hit
    const Hit_t& operator*() const { return *hitPtr; }
    /// @brief Check if two OrderedHit objects are equal: they are equal
    ///  if they point to the same underlying hit
    bool operator==(const OrderedHit& other) const {
      return hitPtr == other.hitPtr;
    }
    /// @brief Stream operator
    friend std::ostream& operator<<(std::ostream& ostr, const OrderedHit& c) {
      c.print(ostr);
      return ostr;
    }
    /// @brief Print the OrderedHit information to a stream
    void print(std::ostream& ostr) const;
  };
  /// @brief Phi-only hit attached to the pattern, together with the point on
  ///        the strip crossed by the pattern line. The strip position is in
  ///        general not the crossing point (e.g. the strip centre), so the
  ///        crossing point is used for the phi of the hit.
  struct PhiOnlyHit {
    /// @brief The phi-only hit, owned by the user
    const Hit_t* hit{nullptr};
    /// @brief Point on the strip crossed by the pattern line
    Vector3 crossPos{Vector3::Zero()};
  };
  /// @brief Enum for the pattern state variables, needed when updating the pattern state
  enum class PatternState : std::int8_t {
    eBendPlane,   ///< Bending plane
    eLineAnchor,  ///< Line anchor
    eLine,        ///< Line parameters
    NVars         ///< Number of pattern state variables
  };
  /// @brief Small struct to encapsulate the result of the line compatibility test
  struct LineTestRes {
    /// @brief Distance of the test hit from the pattern line
    double residual{0.};
    /// @brief Uncertainty on the residual
    double sigma{0.};
    /// @brief Number of degrees of freedom of the residual (1 or 2)
    std::uint8_t nDof{0u};
  };
  /// @brief Enum for the bending plane parameter covariance
  enum class BendPlaneCov : std::int8_t {
    ePhiPhi,        ///< Variance of the plane azimuth phi
    eOffsetOffset,  ///< Variance of the plane offset from the beam axis
    ePhiOffset      ///< Covariance of the plane azimuth and offset
  };
  using BendPlaneCov_t = std::array<double, 3>;
  /// @brief Constructor taking the seed information
  ///  @param seed: seed hit
  ///  @param sectorIdx: sector index of the seed hit
  ///  @param config: configuration object
  ///  @param logger: messaging object
  explicit PatternStateAux(const OrderedHit& seed,
                           const typename Sector_t::Index_t sectorIdx,
                           const Config& config, const Logger& logger);
  /// @brief Move constructor
  ///  @param other: other pattern state to move from
  PatternStateAux(PatternStateAux&& other) noexcept = default;
  /// @brief Move assignment operator
  ///  @param other: other pattern state to move from
  PatternStateAux& operator=(PatternStateAux&& other) noexcept = default;
  /// @brief Copy constructor
  ///  @param other: other pattern state to copy from
  PatternStateAux(const PatternStateAux& other) = default;
  /// @brief Copy assignment operator
  ///  @param other: other pattern state to copy from
  PatternStateAux& operator=(const PatternStateAux& other) = default;
  /// @brief Default destructor
  ~PatternStateAux() = default;
  
  /// ********* Add and remove hits from the pattern *********
  /// @brief Add a hit to the pattern and update the internal state
  ///  @param hit: hit to be added
  ///  @param testRes: line compatibility test result of the hit
  void addHit(const OrderedHit& hit, const LineTestRes& testRes);
  /// @brief Overwrite the last inserted hit with the new one on the same layer
  ///  @param newHit: new hit to replace with
  ///  @param newTestRes: line compatibility test result of the new hit
  void overWriteHit(const OrderedHit& newHit, const LineTestRes& newTestRes);
  /// ********* Residual computation *********
  /// @brief Contract the intrinsic covariance of the hit with a vector, i.e.
  ///        compute the variance of the hit position along that vector. It is
  ///        used both for the residual uncertainty and for the weight of the
  ///        hit in the bending plane regression.
  ///  @param hit: hit for which to compute the contracted covariance
  ///  @param contractionVector: vector used to contract the covariance
  ///  @param isProjected: flag indicating if the hit is projected onto the
  ///         bending plane. Only used in the residual computation, false
  ///         otherwise.
  ///  @return: the contracted covariance of the hit
  static double covarianceContraction(const Hit_t& hit,
                                      const Vector3& contractionVector,
                                      bool isProjected = false);
  /// @brief Method to compute the residual of a test hit against the pattern line
  ///  @param testHit: test hit information
  ///  @param beamSpot: position of the beam spot
  ///  @return: Test result holding the residual, its uncertainty, and number of degrees of freedom.
  LineTestRes computeLineResidual(const OrderedHit& testHit,
                                  const BeamspotInfo& beamSpot) const;
  /// ********* Bending plane computation and hit projection *********
  /// @brief Project a certain hit position onto the bending plane where the pattern is defined.
  ///         The hit is moved along the sensor direction.
  ///  @param hit: hit whose position is to be projected
  ///  @return: projected position
  Vector3 projToPhiPlane(const Hit_t& hit) const;
  /// @brief Update the line parameters based on the current hits
  ///  @param beamSpot: beam spot information, needed when there are not enough hits
  void updateLineParameters(const BeamspotInfo& beamSpot);
  /// @brief Helper method to update the pattern phi and bending plane normal
  ///  @param beamSpot: position of the beam spot
  void updateBendPlane(const BeamspotInfo& beamSpot);
  /// @brief Helper method to set the bending plane from the sector
  void setBendPlaneFromSector();
  /// @brief Signed angle under which the pattern sees a point, relative to its bending plane.
  ///  @param pos: position to check
  ///  @return: the angle between the position and the bending plane
  double angleToBendPlane(const Vector3& pos) const;
  /// @brief Return the pattern's phi at a given radius
  double phiAtRadius(const double R) const;
  /// @brief Return the unit vector along the bending plane
  Vector3 bendPlaneDir() const;
  /// ********** Phi-only strip utilities **********
  /// @brief Define the pattern line locally around a reference hit, using eta hits in its group
  ///         or, failing that, the closest hits in the groups below and above.
  ///  @param refHit: reference hit
  ///  @return true if a line with sufficient lever arm was found
  bool defineLineAroundHit(const OrderedHit& refHit);
  /// @brief Compute where the pattern line crosses a phi-only strip. The
  ///        strip ends are projected onto the bending plane of the pattern
  ///        and the projected strip is intersected with the pattern line.
  ///  @param phiHit: phi-only hit
  ///  @return: pair of the point on the strip crossed by the pattern line and
  ///           the distance of this point from the closest strip edge, along
  ///           the strip projected onto the bending plane (negative if the
  ///           pattern line crosses the strip within its edges)
  std::pair<Vector3, double> phiStripCrossing(const Hit_t& phiHit) const;
  /// @brief Residual of a phi hit from the bending plane: signed distance of the
  ///        hit position from the plane, with the uncertainty of the hit across
  ///        the plane and of the plane parameters at the hit position.
  ///  @param phiHit: hit measuring phi
  ///  @param pos: position of the hit, e.g. the point on a phi-only strip
  ///         crossed by the pattern line
  ///  @return: the residual, its uncertainty and one degree of freedom
  LineTestRes phiHitResidual(const Hit_t& phiHit, const Vector3& pos) const;
  /// *********** Pattern state variables *********
  /// @brief Move the line anchor hit to the closest hit in the closest (but different) group to the reference hit.
  ///  @param refHit: reference hit
  void resetLineAnchorHit(const OrderedHit& refHit);
  /// @brief Synchronize the pattern with the current state
  ///  @param beamSpot: position of the beam spot
  void syncState(const BeamspotInfo& beamSpot);
  /// @brief Book pattern state variables that need to be updated when syncState is called
  ///  @param vars: list of variables to update
  void setStateVars(const std::initializer_list<PatternState>& vars);
  /// ************ Pattern information getters ************
  /// @brief Check whether a hit is present in the pattern
  ///  @param hit: hit to be checked
  ///  @param group: group index of the hit
  ///  @return: true if the hit is in the pattern
  bool isInPattern(const Hit_t& hit, std::uint8_t group) const;
  /// @brief Method returning the number of groups
  ///  @return: number of groups in the pattern
  std::uint8_t countGroups() const;
  /// @brief Return the mean normalized residual squared
  double getMeanResidual2() const;
  /// @brief Return the number of layers in bending coordinate
  ///  @return: number of layers in the bending coordinate
  std::uint8_t nBendingLayers() const;
  /// ************ Retrieve the logger **********
  /// @brief Return the logger
  const Logger& logger() const { return *m_logger; }

  /// @brief Configuration options
  const Config* m_cfg{nullptr};
  /// @brief Logger
  const Logger* m_logger{nullptr};
  /// @brief Last inserted hit. Needed to speed-up lookup
  OrderedHit lastInsertedHit{};
  /// @brief Last hit in the second-to-last layer
  OrderedHit prevLayerHit{};
  /// @brief Line anchor hit, used as first point for the pattern line
  OrderedHit lineAnchorHit{};
  /// @brief Hit used as second point of the pattern line
  OrderedHit lineEndHit{};
  /// @brief Normal vector to the bending plane where the pattern lies
  Vector3 bendPlaneNorm{Vector3::Zero()};
  /// @brief Position and direction of the pattern line. Both are constructed to be
  ///         within the bending plane of the pattern. We store them as vectors
  /// to facilitate vector operations and avoid constructing them repeatedly.
  Vector3 linePos{Vector3::Zero()};
  Vector3 lineDir{Vector3::Zero()};
  /// @brief Distance between the two points defining the pattern line
  double leverArm{0.};
  /// @brief Sum over eta hits of the square of their residual divided by residual uncertainty
  double sumNormResidual2{0.};
  /// @brief Pattern phi, which is the phi of the bending plane where the pattern lies
  double patPhi{0.};
  /// @brief Signed offset of the bending plane from the origin
  double patPhiOffset{0.};
  /// @brief Covariance of the bending plane
  BendPlaneCov_t bendPlaneCov{};
  /// @brief Pattern theta, which is the value of the seed hit
  double patTheta{0.};
  /// @brief Total number of degrees of freedom of the pattern
  std::size_t nDoF{0u};
  /// @brief Bookkeeping of the pattern state variables that need to be
  ///        updated, one bit per PatternState. Placed among the 8-byte
  ///        aligned members, as std::bitset occupies at least one word.
  std::bitset<toUnderlying(PatternState::NVars)> needUpdate{};
  /// @brief Residual and its uncertainty of the last hit tested (needed when replacing a hit)
  LineTestRes lastHitTestRes{};
  /// @brief Sector
  Sector_t sector{static_cast<typename Sector_t::Index_t>(0u)};
  /// @brief Counts of precision / non-precision / phi layers
  std::uint8_t nPrecisionLayers{0u};
  std::uint8_t nTriggerLayers{0u};
  std::uint8_t nPhiLayers{0u};
  /// @brief Flag to indicate if the pattern is overlapping with another one, used during overlap removal
  bool isOverlap{false};
  /// @brief Whether we used the beamspot to compute the line parameters
  bool useBeamspot{false};

  /// @brief Collection of hits per group. A pattern is determined by the hits belonging to it.
  std::array<std::vector<OrderedHit>, nGroups> hitsPerGroup{};
  /// @brief Phi-only hits attached to the pattern
  std::vector<PhiOnlyHit> phiOnlyHits{};

  /// @brief Print the covariance matrix
  static std::string print(const BendPlaneCov_t& cov) {
    std::ostringstream os;
    os << std::format("[varPhi: {}, varOffset: {}, covPhiOffset: {}]",
                      cov[toUnderlying(BendPlaneCov::ePhiPhi)],
                      cov[toUnderlying(BendPlaneCov::eOffsetOffset)],
                      cov[toUnderlying(BendPlaneCov::ePhiOffset)]);
    return os.str();
  }

  /// @brief Patterns are considered identical if they have the same hit content. However, map comparison is very expensive
  bool operator==(const PatternStateAux& other) const = delete;
  /// @brief Print the pattern candidate
  void print(std::ostream& ostr, bool detailed) const;
  /// @brief A view of the pattern state for printing purposes
  struct PatternPrintView {
    /// @brief The pattern state to be printed
    const PatternStateAux& pat;
    /// @brief Whether to print detailed information
    bool detailed = false;

    /// @brief Print the pattern print view
    friend std::ostream& operator<<(std::ostream& os,
                                    const PatternPrintView& v) {
      v.pat.print(os, v.detailed);
      return os;
    }
  };
  /// @brief Print the pattern state with brief information
  static PatternPrintView brief(const PatternStateAux& p);
  /// @brief Print the pattern state with detailed information
  static PatternPrintView detailed(const PatternStateAux& p);
};

/// Absolute difference of two layer indices
constexpr unsigned layerDistance(const std::uint8_t a, const std::uint8_t b) {
  return static_cast<unsigned>(
      Acts::abs(static_cast<int>(a) - static_cast<int>(b)));
}

}  // namespace Acts::Experimental::detail
#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.ipp"
