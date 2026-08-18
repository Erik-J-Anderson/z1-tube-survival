#pragma once


#include "Tube_Survival.hpp"

void InitializeAccumulator(
    SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result
);

void AccumulateSegmentSurvival(SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result);

void FinalizeAccumulator(
    SegmentSurvivalFunction& accumulator
);