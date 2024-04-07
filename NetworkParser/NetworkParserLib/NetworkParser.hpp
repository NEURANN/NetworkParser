#pragma once

#include <cstdint>


#define DLL_IMPORT extern "C" __declspec(dllimport)


#define SOURCE_TYPE_INPUT		0b00000001
#define SOURCE_TYPE_HIDDEN		0b00000010
#define TARGET_TYPE_HIDDEN		0b00000100
#define TARGET_TYPE_OUTPUT		0b00001000

#define CHROMOSOME_SUCCESS			0
#define CHROMOSOME_BAD_PATH			1
#define CHROMOSOME_FILE_SHORT		2
#define CHROMOSOME_MISSING_GENES	3
#define CHROMOSOME_BAD_GENE			4

#define GENE_SUCCESS				0
#define GENE_SHORT					1
#define GENE_EOF					2
#define GENE_MISSING_CONNECTIONS	3


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
	unsigned int ReturnCode = CHROMOSOME_SUCCESS;
	signed   int AdditionalInfo = 0;
};	

DLL_IMPORT SubnetworkChromosomeParseResult ParseSubnetworkChromosome(const char* Filepath);
DLL_IMPORT void FreeParseResult(SubnetworkChromosomeParseResult Result);