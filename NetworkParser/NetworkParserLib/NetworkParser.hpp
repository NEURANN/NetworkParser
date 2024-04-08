#pragma once

#include <cstdint>


#define DLL_IMPORT extern "C" __declspec(dllimport)


/*-----SUBNETWORK CHROMOSOME-----*/
#define SOURCE_TYPE_INPUT		0b00000001
#define SOURCE_TYPE_HIDDEN		0b00000010
#define TARGET_TYPE_HIDDEN		0b00000100
#define TARGET_TYPE_OUTPUT		0b00001000

#define SUBNETWORK_CHROMOSOME_SUCCESS			0
#define SUBNETWORK_CHROMOSOME_BAD_PATH			1
#define SUBNETWORK_CHROMOSOME_FILE_SHORT		2
#define SUBNETWORK_CHROMOSOME_MISSING_GENES		3
#define SUBNETWORK_CHROMOSOME_BAD_GENE			4

#define SUBNETWORK_GENE_SUCCESS					0
#define SUBNETWORK_GENE_SHORT					1
#define SUBNETWORK_GENE_EOF						2
#define SUBNETWORK_GENE_MISSING_CONNECTIONS		3


struct InternalConnectionCodon {
	uint32_t CodonIndex = 0;
	float Weight = 0.0f;
	uint8_t Source = 0;
	uint8_t Target = 0;
	uint8_t Types = SOURCE_TYPE_INPUT | SOURCE_TYPE_HIDDEN;
};

struct SubnetworkGene {
	uint64_t GeneIndex = 0;
	InternalConnectionCodon* Codons = nullptr;
	uint32_t CodonCount = 0;
};

struct SubnetworkChromosomeParseResult {
	SubnetworkGene* Genes = nullptr;
	uint32_t GeneCount = 0;
	unsigned int ReturnCode = SUBNETWORK_CHROMOSOME_SUCCESS;
	signed   int AdditionalInfo = 0;
};	

DLL_IMPORT SubnetworkChromosomeParseResult ParseSubnetworkChromosome(const char* Filepath);
DLL_IMPORT void FreeSubnetworkParseResult(SubnetworkChromosomeParseResult Result);



/*-----QUADRANT CHROMOSOME-----*/
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

DLL_IMPORT QuadrantChromosomeParseResult ParseQuadrantChromosome(const char* Filepath);
DLL_IMPORT void FreeQuadrantParseResult(QuadrantChromosomeParseResult Result);