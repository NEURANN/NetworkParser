#include "pch.h"
#include "subnetwork.hpp"


SubnetworkChromosomeParseResult ParseSubnetworkChromosome(const char* Filepath) {
	SubnetworkChromosomeParseResult result = { };


	//attempt to open the chromosome file
	std::ifstream chromosomeFile;
	chromosomeFile.open(Filepath, std::ifstream::in | std::ifstream::binary);
	if (!chromosomeFile) {
		chromosomeFile.close();
		result.ReturnCode = SUBNETWORK_CHROMOSOME_BAD_PATH;
		result.AdditionalInfo = -1;
		return result;
	}

	//validate the chromosome is a subnetworks chromosome
	char* buf = new char[4];
	chromosomeFile.read(buf, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = SUBNETWORK_CHROMOSOME_FILE_SHORT;
		return result;
	}
	if (strncmp(buf, "SUBN", 4) != 0) {
		chromosomeFile.close();
		result.ReturnCode = SUBNETWORK_CHROMOSOME_BAD_PATH;
		return result;
	}

	//stream successfully opened, first 4 bytes are the subnetwork gene count,
	//i.e. how many times we try to read an individual gene
	chromosomeFile.read((char*)&result.GeneCount, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = SUBNETWORK_CHROMOSOME_FILE_SHORT;
		result.AdditionalInfo = -1;
		return result;
	}


	//iterate through each gene until the end of the genome
	result.Genes = new SubnetworkGene[result.GeneCount];
	for (uint32_t i = 0; i < result.GeneCount; i++) {
		int status = ParseSubnetworkGene(chromosomeFile, result.Genes + i);

		switch (status) {
		case SUBNETWORK_GENE_SUCCESS:
			continue;

		case SUBNETWORK_GENE_EOF:
			chromosomeFile.close();
			result.ReturnCode = SUBNETWORK_CHROMOSOME_MISSING_GENES;
			result.AdditionalInfo = i;
			return result;

		case SUBNETWORK_GENE_SHORT:
		case SUBNETWORK_GENE_MISSING_CONNECTIONS:
			chromosomeFile.close();
			result.ReturnCode = SUBNETWORK_CHROMOSOME_BAD_GENE;
			result.AdditionalInfo = i;
			return result;
		}
	}

	//successful chromosome parse
	chromosomeFile.close();
	result.ReturnCode = SUBNETWORK_CHROMOSOME_SUCCESS;
	result.AdditionalInfo = -1;
	return result;
}

void FreeSubnetworkParseResult(SubnetworkChromosomeParseResult Result) {
	if (Result.Genes == nullptr) {
		return;
	}

	for (uint32_t i = 0; i < Result.GeneCount; i++) {
		delete Result.Genes[i].Codons;
		Result.Genes[i].Codons = nullptr;
	}
	delete Result.Genes;
	Result.Genes = nullptr;
}


int ParseSubnetworkGene(std::ifstream& ChromosomeFile, SubnetworkGene* ResultLocation) {
	//read in the gene parameters
	ChromosomeFile.read((char*)&ResultLocation->GeneIndex, 8);
	if (ChromosomeFile.gcount() != 8) {
		return SUBNETWORK_GENE_SHORT;
	}
	ChromosomeFile.read((char*)&ResultLocation->CodonCount, 4);
	if (ChromosomeFile.gcount() != 4) {
		return SUBNETWORK_GENE_SHORT;
	}

	//iterate through each connection codon until the end of the subnetwork gene
	ResultLocation->Codons = 
		new InternalConnectionCodon[ResultLocation->CodonCount];
	for (uint32_t i = 0; i < ResultLocation->CodonCount; i++) {
		ChromosomeFile.peek();
		if (ChromosomeFile.eof()) {
			return SUBNETWORK_GENE_EOF;
		}

		int status = ParseInternalConnectionCodon(ChromosomeFile, 
			ResultLocation->Codons + i);

		if (status != SUBNETWORK_GENE_SUCCESS) {
			return status;
		}
	}
	return SUBNETWORK_GENE_SUCCESS;
}

int ParseInternalConnectionCodon(
		std::ifstream& ChromosomeFile,
		InternalConnectionCodon* ResultLocation) {

	//read in the codon parameters
	ChromosomeFile.read((char*)&ResultLocation->CodonIndex, 4);
	if (ChromosomeFile.gcount() != 4) {
		return SUBNETWORK_GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Source, 1);
	if (ChromosomeFile.gcount() != 1) {
		return SUBNETWORK_GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Target, 1);
	if (ChromosomeFile.gcount() != 1) {
		return SUBNETWORK_GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Types, 1);
	if (ChromosomeFile.gcount() != 1) {
		return SUBNETWORK_GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Weight, 4);
	if (ChromosomeFile.gcount() != 4) {
		return SUBNETWORK_GENE_MISSING_CONNECTIONS;
	}

	return SUBNETWORK_GENE_SUCCESS;
}
