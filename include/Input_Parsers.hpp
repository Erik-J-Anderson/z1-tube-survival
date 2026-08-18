#pragma once

#include "Primitive_Path_Trajectory.hpp"

#include <string>
#include <vector>


// Input Parser Helpers

std::vector<double> ParseDoubleList(const std::string& input); // Parse "3.0,5.0,7.0" -> vector<double> For parsing tube diameters and lags

Vec3 ParseVec3(const std::string& input); // Parse "20.6,20.6,39.2" -> Vec3 For parsing box center coordinates

void SetFixedBoxCenter(std::vector<Box>& frame_boxes, const Vec3& center); // Reconstruct each box origin from a fixed physical box center


