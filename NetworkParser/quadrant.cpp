#include "pch.h"
#include "quadrant.hpp"

QuadrantChromosomeParseResult ParseQuadrantChromosome(const char* Filepath) {
	QuadrantChromosomeParseResult result = { };


	//attempt to open the chromosome file
	std::ifstream chromosomeFile;
	chromosomeFile.open(Filepath, std::ifstream::in | std::ifstream::binary);
	if (!chromosomeFile) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_BAD_PATH;
		return result;
	}

	//validate the chromosome is a quadrants chromosome
	char* buf = new char[4];
	chromosomeFile.read(buf, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_FILE_SHORT;
		return result;
	}
	if (strcmp(buf, "QUAD") != 0) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_BAD_PATH;
		return result;
	}

	//stream successfully opened, first 4 bytes are the quadrant count
	chromosomeFile.read((char*)&result.QuadrantCount, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_FILE_SHORT;
		return result;
	}

	//next 4 bytes are the subnetworks per quadrant
	//seems redundant, but necessary so that a chromosome is intelligible on its own
	chromosomeFile.read((char*)&result.SubnetworksPerQuadrant, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_FILE_SHORT;
		return result;
	}

	//if the chromosome is as it has described itself, we can allocate the entire
	//set of quadrants into memory, then section it out into a 2D array from there
	uint32_t nIndices = result.SubnetworksPerQuadrant* result.QuadrantCount;
	uint32_t* subnetworkIndices = new uint32_t[nIndices];
	chromosomeFile.read((char*)subnetworkIndices, nIndices * sizeof(uint32_t));
	if (chromosomeFile.gcount() != nIndices * sizeof(uint32_t)) {
		chromosomeFile.close();
		result.ReturnCode = QUADRANT_CHROMOSOME_MISSING_QUADRANTS;
		return result;
	}

	//carve up the 1D array into a 2D array for convenient processing
	result.SubnetworkIndices = new uint32_t*[result.QuadrantCount];
	for (uint32_t i = 0; i < result.QuadrantCount; i++) {
		result.SubnetworkIndices[i] = subnetworkIndices + (i * result.SubnetworksPerQuadrant);
	}

	chromosomeFile.close();
	result.ReturnCode = QUADRANT_CHROMOSOME_SUCCESS;
	return result;
}

void FreeQuadrantParseResult(QuadrantChromosomeParseResult Result) {
	if (Result.SubnetworkIndices == nullptr) {
		return;
	}
	if (Result.SubnetworkIndices[0] == nullptr) {
		return;
	}

	//the subnetwork indices 2D array was initially allocated as a 1D array.
	//therefore, only the first pointer to an array, and the array itself, must be freed.
	delete Result.SubnetworkIndices[0];
	delete Result.SubnetworkIndices;
}