#pragma once

#include <cstdint>
#include <fstream>


#include "common.h"


#define QUADRANT_CHROMOSOME_SUCCESS				0
#define QUADRANT_CHROMOSOME_BAD_PATH			1
#define QUADRANT_CHROMOSOME_FILE_SHORT			2
#define QUADRANT_CHROMOSOME_MISSING_QUADRANTS	3


struct QuadrantChromosomeParseResult {
	uint32_t** SubnetworkIndices = nullptr;
	uint32_t QuadrantCount = 0;
	uint32_t SubnetworksPerQuadrant = 0;
	uint32_t ReturnCode = QUADRANT_CHROMOSOME_SUCCESS;
};

DLL_EXPORT QuadrantChromosomeParseResult ParseQuadrantChromosome(const char* Filepath);
DLL_EXPORT void FreeQuadrantParseResult(QuadrantChromosomeParseResult Result);