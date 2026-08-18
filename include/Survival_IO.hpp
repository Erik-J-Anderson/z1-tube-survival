#pragma once

#include "Tube_Survival.hpp"

#include <string>


void WriteSegmentSurvivalFunction(
    const SegmentSurvivalFunction& result,
    const std::string& filename
);


void WriteTubeSurvivalFunction(
    const TubeSurvivalFunction& result,
    const std::string& filename
);

void WriteEndRetractionFunction(
    const EndRetractionFunction& result,
    const std::string& filename
);