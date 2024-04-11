#pragma once

#include "subnetwork.hpp"
#include "quadrant.hpp"
#include "connections.hpp"
#include "common.h"


#define GENERIC_PARSE_SUBNETWORKS		0
#define GENERIC_PARSE_QUADRANTS			1
#define GENERIC_PARSE_CONNECTIONS		2
#define GENERIC_PARSE_BAD_PATH			3
#define GENERIC_PARSE_UNRECOGNISED		4
#define GENERIC_PARSE_SHORT				5

struct GenericChromosomeParseResult {
	uint32_t ReturnCode;
	union {
		SubnetworkChromosomeParseResult SCPR;
		QuadrantChromosomeParseResult QCPR;
		ConnectionsChromosomeParseResult  CCPR;
	} ParseResult;
};

DLL_EXPORT GenericChromosomeParseResult ParseGenericChromosome(const char* Filepath);