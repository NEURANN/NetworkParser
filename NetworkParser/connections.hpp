#pragma once

#include <cstdint>
#include <fstream>

#include "common.h"


struct ConnectionGene {
	float Weight = 0.0f;
	uint32_t SourceSubnetworkIndex = 0;
	uint8_t SourceOutputIndex = 0;
	uint32_t TargetSubnetworkIndex = 0;
	uint8_t TargetInputIndex = 0;
};

struct QuadrantConnections {
	ConnectionGene* ConnectionGenes = nullptr;
	uint32_t ConnectionGeneCount = 0;
	uint32_t SourceQuadrantIndex = 0;
	uint32_t TargetQuadrantIndex = 0;
};


#define CONNECTIONS_CHROMOSOME_SUCCESS				0
#define CONNECTIONS_CHROMOSOME_BAD_PATH				1
#define CONNECTIONS_CHROMOSOME_SHORT				2

struct ConnectionsChromosomeParseResult {
	QuadrantConnections* QuadrantConnectionsArray = nullptr;
	uint32_t QuadrantConnectionsCount = 0;
	uint32_t ReturnCode = CONNECTIONS_CHROMOSOME_SUCCESS;
	uint32_t AdditionalInfo = -1;
};

DLL_EXPORT ConnectionsChromosomeParseResult ParseConnectionsChromosome(const char* Filepath);
DLL_EXPORT void FreeConnectionsParseResult(ConnectionsChromosomeParseResult Result);


#define CONNECTIONS_CHROMOSOME_QUADRANT_SHORT		3

uint32_t ParseQuadrantConnections(
	std::ifstream& ChromosomeFile,
	QuadrantConnections* ResultLocation);


#define CONNECTIONS_GENE_SUCCESS					4
#define CONNECTIONS_GENE_SHORT						5

uint32_t ParseConnectionGene(
	std::ifstream& ChromosomeFile,
	ConnectionGene* ResultLocation);