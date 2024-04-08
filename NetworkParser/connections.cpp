#include "pch.h"
#include "connections.hpp"


ConnectionsChromosomeParseResult ParseConnectionsChromosome(const char* Filepath) {
	ConnectionsChromosomeParseResult result = { };
	
	//attempt to open the chromosome file
	std::ifstream chromosomeFile;
	chromosomeFile.open(Filepath, std::ifstream::in | std::ifstream::binary);
	if (!chromosomeFile) {
		chromosomeFile.close();
		result.ReturnCode = CONNECTIONS_CHROMOSOME_BAD_PATH;
		return result;
	}

	//stream successfully opened, first 4 bytes are the quadrant connections count
	chromosomeFile.read((char*)&result.QuadrantConnectionsCount, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = CONNECTIONS_CHROMOSOME_SHORT;
		return result;
	}

	//attempt to read the quadrant connections from disk into the newly created array
	result.QuadrantConnectionsArray = new QuadrantConnections[result.QuadrantConnectionsCount];
	for (uint32_t i = 0; i < result.QuadrantConnectionsCount; i++) {
		int status = ParseQuadrantConnections(chromosomeFile, result.QuadrantConnectionsArray + i);

		if (status != CONNECTIONS_GENE_SUCCESS) {
			chromosomeFile.close();
			result.AdditionalInfo = i;
			result.ReturnCode = CONNECTIONS_CHROMOSOME_QUADRANT_SHORT;
			return result;
		}
	}


	chromosomeFile.close();
	result.ReturnCode = CONNECTIONS_CHROMOSOME_SUCCESS;
	return result;
}

void FreeConnectionsParseResult(ConnectionsChromosomeParseResult Result) {
	if (Result.QuadrantConnectionsArray == nullptr) {
		return;
	}

	for (uint32_t i = 0; i < Result.QuadrantConnectionsCount; i++) {
		if (Result.QuadrantConnectionsArray[i].ConnectionGenes != nullptr) {
			delete Result.QuadrantConnectionsArray[i].ConnectionGenes;
		}
	}

	delete Result.QuadrantConnectionsArray;
}


uint32_t ParseQuadrantConnections(
		std::ifstream& ChromosomeFile,
		QuadrantConnections* ResultLocation) {
	ChromosomeFile.read((char*)&(ResultLocation->SourceQuadrantIndex), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_CHROMOSOME_QUADRANT_SHORT;
	}
	ChromosomeFile.read((char*)&(ResultLocation->TargetQuadrantIndex), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_CHROMOSOME_QUADRANT_SHORT;
	}
	ChromosomeFile.read((char*)&(ResultLocation->ConnectionGeneCount), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_CHROMOSOME_QUADRANT_SHORT;
	}

	ResultLocation->ConnectionGenes = new ConnectionGene[ResultLocation->ConnectionGeneCount];
	for (uint32_t i = 0; i < ResultLocation->ConnectionGeneCount; i++) {
		uint32_t result = ParseConnectionGene(ChromosomeFile, ResultLocation->ConnectionGenes + i);

		if (result != CONNECTIONS_GENE_SUCCESS) {
			return result;
		}
	}

	return CONNECTIONS_GENE_SUCCESS;
}

uint32_t ParseConnectionGene(
		std::ifstream& ChromosomeFile,
		ConnectionGene* ResultLocation) {
	ChromosomeFile.read((char*)&(ResultLocation->Weight), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_GENE_SHORT;
	}

	ChromosomeFile.read((char*)&(ResultLocation->SourceSubnetworkIndex), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_GENE_SHORT;
	}

	ChromosomeFile.read((char*)&(ResultLocation->SourceOutputIndex), 1);
	if (ChromosomeFile.gcount() != 1) {
		return CONNECTIONS_GENE_SHORT;
	}

	ChromosomeFile.read((char*)&(ResultLocation->TargetSubnetworkIndex), 4);
	if (ChromosomeFile.gcount() != 4) {
		return CONNECTIONS_GENE_SHORT;
	}

	ChromosomeFile.read((char*)&(ResultLocation->TargetInputIndex), 1);
	if (ChromosomeFile.gcount() != 1) {
		return CONNECTIONS_GENE_SHORT;
	}

	return CONNECTIONS_GENE_SUCCESS;
}