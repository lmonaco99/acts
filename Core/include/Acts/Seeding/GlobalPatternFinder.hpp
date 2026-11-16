// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <compare>
#include <span>

/// The Global Pattern Finder (GPF) is a fast reconstruction algorithm whose
/// primary purpose is to identify and collect detector hits belonging to the
/// same particle trajectory. In contrast to conventional tracking algorithms,
/// the GPF does not attempt to estimate track parameters; instead, it focuses
/// on efficiently identifying compatible hit patterns that can subsequently be
/// used by downstream reconstruction algorithms.
///
///  The basic assumption of the GPF is that, to first order, the trajectory can
/// be described by a bending and a non-bending coordinate, where the bending
/// coordinate describes the direction in which the trajectory changes due to
/// the magnetic field, while the non-bending coordinate remains approximately constant
/// over the relevant detector region. As a consequence, the detector can be
/// sectorized in the non-bending coordinate. This sectorization reduces the
/// search space and allows the GPF to efficiently identify compatible hits
/// across the detector layers.
///
///  The detector is assumed to be organized in layers, which are grouped
/// into a number of groups. The layers are ordered along the trajectory by a
/// user-provided ordering of the hits (Config::layerOrder). This allows the
/// pattern-finding algorithm to operate independently of the actual detector
/// geometry.
///
///  The precise definition of the bending and non-bending coordinates, as well
/// as the corresponding sectorization, is detector-dependent. The pattern-finding
/// strategy itself is independent of these specific coordinate choices.
///  **** This first implementation is developed for a toroidal magnetic field
/// where the bending coordinate is the polar angle theta, while the non-bending
/// one is the azimuthal angle phi. ****
///
///  The hits are space points satisfying the @ref CompositeSpacePoint concept,
/// expressed in the global frame: their localPosition(), sensorDirection(),
/// toNextSensor() and planeNormal() are global. The hit type declares it with
/// static constexpr bool isGlobalFrame = true and provides its group index and
/// precision flag. The detector geometry stays with the experiment: the layer 
/// ordering, the hit overlap and the phi-only hits that can be attached to a 
/// pattern hit are provided by the user. All the hits are owned by the user and 
/// must outlive the output patterns.
///
///  The definition of the sectorization concept, together
///  with the main @ref PatternStateAux object used by the GPF, are provided in the GlobalPatternFinderAuxiliaries.
namespace Acts::Experimental::detail {
template <CompositeSpacePoint Hit_t, SectorType Sector_t,
          std::uint8_t nGroups>
class GlobalPatternFinder {
 public:
  /// @brief Pattern state auxiliaries used during the pattern building
  using PatternStateAux = detail::PatternStateAux<Hit_t, Sector_t, nGroups>;
  /// @brief Beamspot information
  using BeamspotInfo = PatternStateAux::BeamspotInfo;
  /// The hits must be expressed in the global frame: the hit type declares it
  /// with static constexpr bool isGlobalFrame = true
  static_assert(requires { requires Hit_t::isGlobalFrame; },
                "GlobalPatternFinder: the hits must be expressed in the global frame "
                "(static constexpr bool isGlobalFrame = true in the hit class)");
  /// Whether two hits represent overlapping detector information. This
  /// relation may be less restrictive than requiring equality and is used
  /// when determining whether reconstructed patterns share detector
  /// information.
  using HitOverlap_t = Delegate<bool(const Hit_t&, const Hit_t&)>;
  /// @brief Delegate to order two hits by their measurement layer along the trajectory
  using LayerOrder_t = Delegate<std::weak_ordering(const Hit_t&, const Hit_t&)>;
  /// @brief Delegate returning the phi-only hits that can be attached to a pattern hit
  using PhiOnlyHits_t = Delegate<std::span<const Hit_t>(const Hit_t&)>;
  /// @brief Delegate returning the group index of a hit in the pattern topology
  using GroupIndex_t = Delegate<std::uint8_t(const Hit_t&)>;
  /// @brief Delegate returning whether a hit is a precision measurement
  using IsPrecision_t = Delegate<bool(const Hit_t&)>;
  /// @brief Delegate deciding whether a hit seeds a pattern
  using SeedSelection_t = Delegate<bool(const Hit_t&)>;
  /// @brief Delegate returning the theta search window, per group, around a seed hit
  using ThetaSearchWindow_t = Delegate<const std::array<double, nGroups>&(const Hit_t&)>;
  /// @brief Configuration object for the pattern finder
  struct Config : PatternStateAux::Config {
    /// Pattern building acceptance
    /// @brief Number of standard deviations to consider for residual acceptance
    double nResidualSigma{3.0};
    /// @brief Maximum phi difference [rad] allowed between two hits belonging to the same pattern
    double phiTolerance{3.5 * UnitConstants::degree};
    /// @brief Maximum theta difference [rad] allowed between two patterns to be considered overlapping
    double maxThetaSizeOverlap{6. * UnitConstants::degree};
    /// @brief Residual uncertainty to consider the hit as low confidence
    double lowConfidenceResSigma{50.};
    /// @brief Factor by which the acceptance window is loosened when using the beamspot or searching in a new group
    double loosenAcceptanceFactor{2.5};
    /// Pile-up & Fake rate suppression
    /// @brief Minimum number of strip layers in the bending direction required to accept a pattern
    unsigned int minStripEtaLayers{2};
    /// @brief Minimum number of precision layers in the bending direction required to accept a pattern
    unsigned int minPrecisionLayers{8};
    /// @brief Minimum number of phi layers required to accept a pattern
    unsigned int minPhiLayers{1};
    /// @brief Minimum number of layer groups required to accept a pattern
    unsigned int minGroups{2};
    /// @brief Minimum number of layers in a group to be considered a good group
    unsigned int minGroupLayers{4};
    /// @brief Quality cut on the pattern's mean squared normalized residual
    double meanNormRes2Cut{3.5};
    /// Pattern recovery power
    /// @brief Maximum number of attempts to build a pattern from hits already used in existing patterns
    unsigned int maxSeedAttempts{3};
    /// @brief Maximum number of missed candidate hits in different measurement layers in a group
    unsigned int maxMissLayersInGroup{3};
    /// @brief Maximum number of active pattern branches kept after each extension step
    ///        (0: no limit). The best branches according to isBetterBranch are
    ///        kept.
    unsigned int maxActivePatterns{0};
    /// @brief (optional) Delegate to determine whether two hits represent overlapping detector information.
    ///        If not provided, only a hit with itself is considered overlapping.
    HitOverlap_t hitOverlap{};
    /// @brief Delegate ordering two hits by their measurement layer along the trajectory.
    ///        The hits on the same layer keep the order in which the user provided them.
    LayerOrder_t layerOrder{};
    /// @brief (optional) Delegate returning the phi-only hits that can be attached to a pattern hit,
    ///        i.e. those of the same detector region. Each region is searched only once.
    PhiOnlyHits_t phiOnlyHits{};
    /// @brief Delegate returning the group index of a hit in the pattern topology.
    GroupIndex_t groupIndex{};
    /// @brief Delegate returning whether a hit is a precision measurement.
    IsPrecision_t isPrecision{};
    /// @brief Delegate deciding whether a hit seeds a pattern.
    SeedSelection_t goodForSeeding{};
    /// @brief Delegate returning the (half) search window in the theta coordinate, per group of the
    ///        pattern topology, of the hits compatible with a pattern seeded by a hit.
    ThetaSearchWindow_t thetaSearchWindow{};
  };
  /// @brief Entry of the hit search: a hit with the sector and the theta coordinate in which it is
  ///        searched. A hit can have several entries, e.g. in the overlap of two sectors.
  struct SearchEntry {
    /// @brief Sector index
    typename Sector_t::Index_t sector{};
    /// @brief Theta coordinate [rad]
    double theta{0.};
    /// @brief The hit, owned by the user
    const Hit_t* hit{nullptr};
  };
  /// @brief Entries of the hit search
  using SearchEntries_t = std::vector<SearchEntry>;

  /// @brief Structure to hold the pattern result
  struct OutputPattern {
    explicit OutputPattern(typename Sector_t::Index_t sectorIdx)
        : sector{sectorIdx} {}
    /// @brief Hits of the pattern per group, owned by the user, including phi-only hits
    std::array<std::vector<const Hit_t*>, nGroups> hitsPerGroup{};
    /// @brief Mean over eta hits of the square of their residual divided by residual uncertainty
    double meanNormResidual2{0.};
    /// @brief Pattern phi, which is the phi of the bending plane where the pattern lies
    double patPhi{0.};
    /// @brief Signed offset of the bending plane from the origin
    double patPhiOffset{0.};
    /// @brief Covariance of the bending plane
    PatternStateAux::BendPlaneCov_t bendPlaneCov{};
    /// @brief Pattern theta, which is the value of the seed hit
    double patTheta{0.};
    /// @brief Sector
    Sector_t sector{static_cast<typename Sector_t::Index_t>(0u)};
    /// @brief Counts of precision / non-precision / phi layers
    std::uint8_t nPrecisionLayers{0u};
    std::uint8_t nTriggerLayers{0u};
    std::uint8_t nPhiLayers{0u};
  };

  /// @brief Standard constructor
  ///  @param config: Configuration parameters
  ///  @param outputLogger: Logger instance used for the messaging
  explicit GlobalPatternFinder(Config&& config,
                               std::unique_ptr<const Logger> outputLogger =
                                   getDefaultLogger("GlobalPatternFinder",
                                                    Logging::Level::INFO));

  /// @brief Main method steering the pattern finding. Given the search entries
  ///         of the hits, it builds patterns in eta, attaches the compatible
  ///         phi-only hits if Config::phiOnlyHits is connected, and converts the
  ///         patterns into OutputPatterns. The entries are sorted by sector and
  ///         theta (stable: equal entries keep their order); the seeds are tried
  ///         in this order and the candidate hits of a seed come in theta order.
  ///  @param searchEntries: Search entries of the hits (any order)
  ///  @param beamspotInfo: Beamspot information
  ///  @return: Vector of found patterns
  std::vector<OutputPattern> findPatterns(
      SearchEntries_t searchEntries,
      const BeamspotInfo& beamspotInfo) const;

 private:
  using PatternStateAuxVec = std::vector<PatternStateAux>;
  using OrderedHit = typename PatternStateAux::OrderedHit;
  using LineTestRes = typename PatternStateAux::LineTestRes;

  /// @brief Enum for possible outcomes of pattern line compatibility test
  enum class LineTestDecision : std::uint8_t {
    /// @brief Test successful, add hit to pattern
    eAddHit,
    /// @brief Test successful with multiple pattern hits on same layer, branch the pattern
    eBranchPattern,
    /// @brief Test successful and the new hit replaces the last inserted hit
    eReplaceHit,
    /// @brief Test failed, discard the hit
    eRejectHit
  };
  using LineTestResult_t = std::pair<LineTestDecision, LineTestRes>;

  /// @brief Method steering the global pattern building in the bending plane.
  ///  @param searchEntries: Search entries of the hits, sorted by sector and theta
  ///  @param beamspotInfo: Information about the beam spot
  ///  @return: resulting vector of PatternStateAuxs successfully built
  PatternStateAuxVec findPatternsInEta(const SearchEntries_t& searchEntries,
                                       const BeamspotInfo& beamspotInfo) const;
  /// @brief Main function controlling the development of patterns, including pattern branching when necessary.
  ///         It tests pattern compatibility of a set of active patterns
  /// (patterns produced from the same seed hit) against one test hit. At the
  /// end, endPatterns contains the surviving patterns.
  ///  @param startPatterns: Vector of active patterns to be extended
  ///  @param endPatterns: Vector to store the surviving patterns after testing against the test hit.
  ///  @param testHit: Hit to be tested against the patterns
  ///  @param beamSpotInfo: Information about the beam spot, needed when the pattern line cannot be reliably defined from the pattern hits
  void extendPatterns(PatternStateAuxVec& startPatterns,
                      PatternStateAuxVec& endPatterns,
                      const OrderedHit& testHit,
                      const BeamspotInfo& beamSpotInfo) const;
  /// @brief Method checking line compatibility of a test hit against the pattern
  ///  @param pat: pattern to be tested, synchronized when the test hit is on a new layer
  ///  @param testHit: test hit information
  ///  @param beamSpot: Beam spot information, needed to update the pattern line
  ///  @return: result of the test, including the computed line residual and acceptance window
  LineTestResult_t checkLineCompatibility(PatternStateAux& pat,
                                          const OrderedHit& testHit,
                                          const BeamspotInfo& beamSpot) const;
  /// @brief Method to check the phi compatibility of a test hit with the pattern
  ///  @param pat: pattern to be checked
  ///  @param pos: global position at which the phi of the test hit is
  ///         evaluated
  ///  @return: true if the test hit is phi compatible with the pattern
  bool isPhiCompatible(const PatternStateAux& pat, const Vector3& pos) const;
  /// @brief Check if a pattern passes the quality cuts
  ///  @param pat: Pattern to be checked
  ///  @return: true if the pattern passes the cuts, false otherwise
  bool passPatternCuts(const PatternStateAux& pat) const;
  /// @brief Check if a pattern has the minimum number of phi layers
  ///  @param pat: Pattern to be checked
  ///  @return: true if the pattern passes the cut, false otherwise
  bool passPhiLayerCut(const PatternStateAux& pat) const;
  /// @brief Method to insert a pattern into a vector if it is better than the existing patterns
  ///  @param patterns: vector of patterns to insert into
  ///  @param candidate: candidate pattern to be inserted
  ///  @return: true if the candidate was inserted, false otherwise
  bool insertIfBest(PatternStateAuxVec& patterns,
                    PatternStateAux&& candidate) const;
  /// @brief Check if two patterns overlap in space or in hit content
  ///  @param a: first pattern
  ///  @param b: second pattern
  ///  @return: true if the patterns overlap, false otherwise
  bool areOverlapping(const PatternStateAux& a, const PatternStateAux& b) const;
  /// @brief Overlap resolution: Compare two patterns and define which one is better.
  ///  @param a: first pattern
  ///  @param b: second pattern
  ///  @return: true if pattern a is better than pattern b, false otherwise
  bool isBetterOverlap(const PatternStateAux& a,
                       const PatternStateAux& b) const;
  /// @brief Ranking of the pattern branches of one seed, used by the pruning
  ///        of branches with the same last hit and by the cap on active
  ///        branches.
  ///  @param a: first pattern
  ///  @param b: second pattern
  ///  @return: true if pattern a is ranked before pattern b
  static bool isBetterBranch(const PatternStateAux& a,
                             const PatternStateAux& b);
  /// @brief Default hit overlap: only a hit with itself
  ///  @param a: first hit
  ///  @param b: second hit
  ///  @return: true if the hits are overlapping, false otherwise
  static bool defaultHitOverlap(const Hit_t& a, 
                                const Hit_t& b);
  /// @brief Count the number of good groups in a pattern
  ///  @param pat: pattern to be checked
  ///  @return: number of good groups in the pattern
  std::uint8_t countGoodGroups(const PatternStateAux& pat) const;
  /// @brief Method to add phi-only measurements to existing PatternStateAuxs. The
  ///         phi-only hits are searched in the parent buckets of the
  ///         pattern hits. Patterns with too few phi layers afterwards are
  ///         removed.
  ///  @param beamSpot: Beam spot information, needed to update the bending plane
  ///  @param patterns: Vector of pattern states to which to add phi-only hits
  void addPhiOnlyHits(const BeamspotInfo& beamSpot,
                      PatternStateAuxVec& patterns) const;
  /// @brief Count the number of patterns that contain a specific hit
  ///  @param patterns: Vector of pattern states to check
  ///  @param seed: The seed hit to check
  ///  @param seedEntry: The search entry of the seed hit
  ///  @return: The number of patterns that contain the hit
  unsigned countPatterns(const PatternStateAuxVec& patterns, const Hit_t& seed,
                         const SearchEntry& seedEntry) const;
  /// @brief Convert a PatternStateAux into an OutputPattern object
  ///  @param candidate: PatternStateAux to be converted
  ///  @return: Converted OutputPattern
  OutputPattern convertToPattern(PatternStateAux&& candidate) const;
  /// @brief Convert a vector of PatternStateAuxs into OutputPattern objects
  ///  @param candidates: PatternStateAuxs to be converted
  ///  @return: Vector of converted OutputPatterns
  std::vector<OutputPattern> convertToPattern(
      PatternStateAuxVec&& candidates) const;

  /// @brief Global Pattern Recognition configuration
  Config m_cfg;
  /// @brief Cut to the residual when we have 2 dof
  double m_cut2Dof{std::sqrt(
      -2. * std::log(std::erfc(m_cfg.nResidualSigma / std::numbers::sqrt2)))};
  /// @brief Logger for the Global Pattern Finder
  std::unique_ptr<const Acts::Logger> m_logger{};
  /// Reference to the logger object
  const Logger& logger() const { return *m_logger; }
};
}  // namespace Acts::Experimental::detail

#include "Acts/Seeding/GlobalPatternFinder.ipp"
