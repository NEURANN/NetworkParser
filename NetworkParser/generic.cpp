#include "pch.h"
#include "generic.hpp"

GenericChromosomeParseResult ParseGenericChromosome(const char* Filepath) {
	GenericChromosomeParseResult result = {  };
	
	//attempt to open the chromosome file
	std::ifstream chromosomeFile;
	chromosomeFile.open(Filepath, std::ifstream::in | std::ifstream::binary);
	if (!chromosomeFile) {
		chromosomeFile.close();
		result.ReturnCode = GENERIC_PARSE_BAD_PATH;
		return result;
	}

	//read the first 4 bytes as a string
	char *buf = new char[4];
	chromosomeFile.read(buf, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = GENERIC_PARSE_SHORT;
		return result;
	}

	//parse the chromosome based on the determined type
	chromosomeFile.close();
	if (strcmp(buf, "SUBN")) {
		result.ParseResult.SCPR = ParseSubnetworkChromosome(Filepath);
		result.ReturnCode = GENERIC_PARSE_SUBNETWORKS;
	}
	else if (strcmp(buf, "QUAD")) {
		result.ParseResult.QCPR = ParseQuadrantChromosome(Filepath);
		result.ReturnCode = GENERIC_PARSE_QUADRANTS;
	}
	else if (strcmp(buf, "CONN")) {
		result.ParseResult.CCPR = ParseConnectionsChromosome(Filepath);
		result.ReturnCode = GENERIC_PARSE_CONNECTIONS;
	}
	else {
		result.ReturnCode = GENERIC_PARSE_UNRECOGNISED;
	}

	return result;
}