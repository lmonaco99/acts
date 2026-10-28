// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/KDTree.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <compare>
#include <optional>

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
/// bending, while the non-bending coordinate remains approximately constant
/// over the relevant detector region. As a consequence, the detector can be
/// sectorized in the non-bending coordinate. This sectorization reduces the
/// search space and allows the GPF to efficiently identify compatible hits
/// across the detector layers.
///
///  The detector is assumed to be organized in layers, which are grouped
/// into a number of groups. Within a local frame, the layers are ordered by
/// the layer index of the hits. Hits in different local frames are ordered by
/// an optional, user-provided fixed ordering of the groups and, where the
/// ordering is not fixed, by their depth in the detector. This allows the
/// pattern-finding algorithm to operate independently of the actual detector
/// geometry.
///
///  The precise definition of the bending and non-bending coordinates, as well
/// as the corresponding sectorization, is detector-dependent. The
/// pattern-finding strategy itself is independent of these specific coordinate
/// choices.
///  **** This first implementation is developed for a toroidal magnetic field
/// where the bending coordinate is the polar angle theta, while the non-bending
/// one is the azimuthal angle phi. ****
///
///  The hits are space points of containers, referred to as buckets, satisfying the
/// @ref CompositeSpacePointContainer concept. All the space points of a
/// bucket are expressed in the same local frame. The phi-only hits are
/// searched in the parent buckets of the pattern hits.
///
///  The definition of the hit type, of the sectorization concept, together
///  with the main @ref PatternStateAux object used by the GPF, are provided in the GlobalPatternFinderAuxiliaries.
namespace Acts::Experimental::detail {
/// @brief Concept for a pattern seed selector
template <typename Selector_t, typename Bucket_t, std::uint8_t nGroups>
concept PatternSeedSelector =
    requires(const Selector_t& selector, const GlobPatFinderHit<Bucket_t>& hit) {
      /// Whether the hit is good for seeding a pattern
      { selector.goodForSeeding(hit) } -> std::same_as<bool>;
      /// The (half) search window in the theta coordinate for compatible hits
      /// with a pattern for a given seed hit. The search window is defined per
      /// group of the pattern topology.
      {
        selector.thetaSearchWindow(hit)
      } -> std::convertible_to<std::array<double, nGroups>>;
    };

template <CompositeSpacePointContainer Bucket_t, SectorType Sector_t,
          std::uint8_t nGroups>
class GlobalPatternFinder {
 public:
  /// @brief Pattern state auxiliaries used during the pattern building
  using PatternStateAux = detail::PatternStateAux<Bucket_t, Sector_t, nGroups>;
  /// @brief Beamspot information
  using BeamspotInfo = PatternStateAux::BeamspotInfo;
  /// @brief Hit type used by the pattern finder
  using Hit_t = GlobPatFinderHit<Bucket_t>;
  /// @brief Space point type of the hits
  using SpacePoint_t = typename Hit_t::SpacePoint_t;
  /// Whether two hits represent overlapping detector information. This
  /// relation may be less restrictive than requiring equality and is used
  /// when determining whether reconstructed patterns share detector
  /// information.
  using HitOverlap_t = Delegate<bool(const Hit_t&, const Hit_t&)>;
  /// @brief Delegate to determine the layer index in the local frame of phi-only hits
  using PhiOnlyLayerIdx_t = Delegate<std::uint8_t(const SpacePoint_t&)>;
  /// @brief Delegate to fix the ordering of two different groups along the trajectory
  using GroupOrder_t = Delegate<std::partial_ordering(std::uint8_t, std::uint8_t)>;
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
    /// @brief Whether to search for phi-only hits in the parent buckets of the pattern hits
    bool attachPhiOnlyHits{true};
    /// @brief (optional) Delegate to determine whether two hits represent overlapping detector information.
    ///        If not provided, hits are considered overlapping if they are on the same layer of the
    ///        same local frame and their measured bending coordinate cannot be distinguished within
    ///        the measurement uncertainty. The same measurement stored in different local frames is
    ///        therefore not recognised as overlapping by the default.
    HitOverlap_t hitOverlap{};
    /// @brief Delegate to determine the layer index of phi-only hits in the local frame. It must
    ///        follow the same numbering as the layer index of the pattern hits (increasing along
    ///        the local z-axis). Required only if attachPhiOnlyHits is set.
    PhiOnlyLayerIdx_t phiOnlyLayerIdx{};
    /// @brief (optional) Fixed ordering of two different groups a and b along the trajectory:
    ///        less if group a comes before group b, greater if after, and unordered if the
    ///        order is not fixed by the topology. Must be antisymmetric and transitive. Pairs
    ///        without a fixed order are ordered by geometry.
    ///        If not connected, hits in different local frames are ordered by geometry only.
    GroupOrder_t groupOrder{};
  };
  /// @brief Abbreviation of the seed coordinates
  enum class HitCoords : std::uint8_t {
    /// Sector coordinate of the associated spectrometer sector
    eSector = 0,
    /// Global Theta
    eTheta = 1,
    /// Number of coordinates
    eNCoords = 2
  };

  /// @brief Definition of the search tree class
  using SearchTree_t = KDTree<toUnderlying(HitCoords::eNCoords), const Hit_t*,
                              double, std::array, 50>;

  /// @brief Structure to hold the pattern result
  struct OutputPattern {
    explicit OutputPattern(typename Sector_t::Index_t sectorIdx)
        : sector{sectorIdx} {}
    /// @brief Vector of hits in the pattern
    std::array<std::vector<const Hit_t*>, nGroups> hitsPerGroup{};
    /// @brief Vector of phi-only hits
    std::vector<Hit_t> phiOnlyHits{};
    /// @brief Distinct parent buckets of the pattern hits
    std::vector<const Bucket_t*> parentBuckets{};
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

  /// @brief Main method steering the pattern finding. Given the search tree of
  ///         hits, it builds patterns in eta, attaches the compatible phi-only
  ///         hits if Config::attachPhiOnlyHits is set, and converts the
  ///         patterns into OutputPatterns.
  ///  @param treeData: Data for the search tree
  ///  @param seedSelector: Selector for the seed
  ///  @param beamspotInfo: Beamspot information
  ///  @return: Vector of found patterns
  template <PatternSeedSelector<Bucket_t, nGroups> SeedSelector_t>
  std::vector<OutputPattern> findPatterns(
      const SearchTree_t& treeData,
      const SeedSelector_t& seedSelector,
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
  ///  @param orderedSpacepoints: Search tree with spacepoints ordered by their corresponding coordinates
  ///  @param seedSelector: The seed selector to use for selecting the initial seed
  ///  @param beamspotInfo: Information about the beam spot
  ///  @return: resulting vector of PatternStateAuxs successfully built
  template <PatternSeedSelector<Bucket_t, nGroups> SeedSelector_t>
  PatternStateAuxVec findPatternsInEta(const SearchTree_t& orderedSpacepoints,
                                       const SeedSelector_t& seedSelector,
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
  /// @brief Method to sort the hits in a pattern by their measurement layer
  ///  @param a: first hit
  ///  @param b: second hit
  ///  @return: weak ordering of the two hits by their measurement layer
  std::weak_ordering compareLayers(const Hit_t& a, 
                                   const Hit_t& b) const;
  /// @brief Default hit overlap: same space point, or same layer of the same local frame
  ///        with bending coordinates compatible within the measurement uncertainty
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
  ///  @param seedCoords: The coordinates of the seed hit
  ///  @return: The number of patterns that contain the hit
  unsigned countPatterns(const PatternStateAuxVec& patterns, const Hit_t& seed,
                         const SearchTree_t::coordinate_t& seedCoords) const;
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
  /// @brief Type alias for the table caching the ordering of the groups along the trajectory
  using GroupOrderTable_t = std::array<std::array<std::optional<bool>, nGroups>, nGroups>;
  /// @brief Cached Config::groupOrder: m_groupOrder[a][b] is true if group a comes before
  ///        group b, false if after, and empty if the order is not fixed
  GroupOrderTable_t m_groupOrder{};
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
