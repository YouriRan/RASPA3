module;

export module mc_moves_move_types;

import std;

export struct Move
{
  enum class Types : std::size_t
  {
    Translation = 0,
    RandomTranslation = 1,
    Rotation = 2,
    RandomRotation = 3,
    VolumeChange = 4,
    ReinsertionCBMC = 5,
    PartialReinsertionCBMC = 6,
    IdentityChangeCBMC = 7,
    Swap = 8,
    SwapCBMC = 9,
    SwapCFCMC = 10,
    SwapCBCFCMC = 11,
    GibbsVolume = 12,
    GibbsSwapCBMC = 13,
    GibbsSwapCFCMC = 14,
    Widom = 15,
    WidomCFCMC = 16,
    WidomCBCFCMC = 17,
    ParallelTempering = 18,
    HybridMC = 19,
    VolumeNCMC = 20,
    SwapNCMC = 21,
    SwapNonEqCBMC = 22,
    Count = 23
  };

  static std::array<std::string, std::to_underlying(Move::Types::Count)> moveNames;
};

export inline std::unordered_set<Move::Types> componentMoves = {
    Move::Types::Translation,        Move::Types::RandomTranslation,
    Move::Types::Rotation,           Move::Types::RandomRotation,
    Move::Types::ReinsertionCBMC,    Move::Types::PartialReinsertionCBMC,
    Move::Types::IdentityChangeCBMC, Move::Types::Swap,
    Move::Types::SwapCBMC,           Move::Types::SwapCFCMC,
    Move::Types::SwapCBCFCMC,        Move::Types::Widom,
    Move::Types::WidomCFCMC,         Move::Types::WidomCBCFCMC,
    Move::Types::SwapNCMC,           Move::Types::SwapNonEqCBMC};

export inline std::unordered_set<Move::Types> systemMoves = {Move::Types::VolumeChange, Move::Types::HybridMC,
                                                             Move::Types::VolumeNCMC};

export inline std::unordered_set<Move::Types> crossSystemMoves = {
    Move::Types::GibbsVolume, Move::Types::GibbsSwapCBMC, Move::Types::GibbsSwapCFCMC, Move::Types::ParallelTempering};

export inline std::unordered_set<Move::Types> groupMoves = {};
