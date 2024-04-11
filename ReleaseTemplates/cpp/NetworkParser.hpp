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



/*-----CONNECTIONS CHROMOSOME-----*/
#define CONNECTIONS_CHROMOSOME_SUCCESS				0
#define CONNECTIONS_CHROMOSOME_BAD_PATH				1
#define CONNECTIONS_CHROMOSOME_SHORT				2
#define CONNECTIONS_CHROMOSOME_QUADRANT_SHORT		3


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

struct ConnectionsChromosomeParseResult {
	QuadrantConnections* QuadrantConnectionsArray = nullptr;
	uint32_t QuadrantConnectionsCount = 0;
	unsigned int ReturnCode = CONNECTIONS_CHROMOSOME_SUCCESS;
	signed   int AdditionalInfo = -1;
};

DLL_IMPORT ConnectionsChromosomeParseResult ParseConnectionsChromosome(const char* Filepath);
DLL_IMPORT void FreeConnectionsParseResult(ConnectionsChromosomeParseResult Result);



/*-----GENERIC CHROMOSOME-----*/
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

DLL_IMPORT GenericChromosomeParseResult ParseGenericChromosome(const char* Filepath);