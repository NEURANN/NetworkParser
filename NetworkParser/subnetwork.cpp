#include "pch.h"
#include "subnetwork.hpp"

//TEMPORARY FOR DEBUGGING. Not required for file.
//#include <iostream>

//helper function to destroy a gene array and the genes' codon arrays
auto CleanupGenes = [](SubnetworkGene* Genes, uint32_t Count) {
	if (Genes == nullptr) {
		return;
	}

	for (uint32_t i = 0; i < Count; i++) {
		delete Genes[i].Codons;
	}
	delete Genes;
};

SubnetworkChromosomeParseResult ParseSubnetworkChromosome(const char* Filepath) {
	SubnetworkChromosomeParseResult result = { };


	//attempt to open the chromosome file
	std::ifstream chromosomeFile;
	chromosomeFile.open(Filepath, std::ifstream::in | std::ifstream::binary);
	if (!chromosomeFile) {
		chromosomeFile.close();
		result.ReturnCode = CHROMOSOME_BAD_PATH;
		return result;
	}

	//stream successfully opened, first 4 bytes are the subnetwork gene count,
	//i.e. how many times we try to read an individual gene
	chromosomeFile.read((char*)&result.GeneCount, 4);
	if (chromosomeFile.gcount() != 4) {
		chromosomeFile.close();
		result.ReturnCode = CHROMOSOME_FILE_SHORT;
		return result;
	}


	//iterate through each gene until the end of the genome
	result.Genes = new SubnetworkGene[result.GeneCount];
	for (uint32_t i = 0; i < result.GeneCount; i++) {
		int status = ParseSubnetworkGene(chromosomeFile, result.Genes + i);

		switch (status) {
		case GENE_SUCCESS:
			continue;

		case GENE_EOF:
			chromosomeFile.close();
			CleanupGenes(result.Genes, result.GeneCount);
			result.ReturnCode = CHROMOSOME_MISSING_GENES;
			result.AdditionalInfo = i;
			return result;

		case GENE_SHORT:
		case GENE_MISSING_CONNECTIONS:
			chromosomeFile.close();
			CleanupGenes(result.Genes, result.GeneCount);
			result.ReturnCode = CHROMOSOME_BAD_GENE;
			result.AdditionalInfo = i;
			return result;
		}
	}

	//successful chromosome parse
	result.ReturnCode = CHROMOSOME_SUCCESS;
	return result;
}

void FreeParseResult(SubnetworkChromosomeParseResult Result) {
	//CleanupGenes(Result.Genes, Result.GeneCount);
}


int ParseSubnetworkGene(std::ifstream& ChromosomeFile, SubnetworkGene* ResultLocation) {
	//read in the gene parameters
	ChromosomeFile.read((char*)&ResultLocation->GeneIndex, 8);
	if (ChromosomeFile.gcount() != 8) {
		return GENE_SHORT;
	}
	ChromosomeFile.read((char*)&ResultLocation->CodonCount, 4);
	if (ChromosomeFile.gcount() != 4) {
		return GENE_SHORT;
	}

	//iterate through each connection codon until the end of the subnetwork gene
	ResultLocation->Codons = 
		new InternalConnectionCodon[ResultLocation->CodonCount];
	for (uint32_t i = 0; i < ResultLocation->CodonCount; i++) {
		ChromosomeFile.peek();
		if (ChromosomeFile.eof()) {
			return GENE_EOF;
		}

		int status = ParseInternalConnectionCodon(ChromosomeFile, 
			ResultLocation->Codons + i);

		if (status != GENE_SUCCESS) {
			return status;
		}
	}
	return GENE_SUCCESS;
}

int ParseInternalConnectionCodon(
		std::ifstream& ChromosomeFile,
		InternalConnectionCodon* ResultLocation) {

	//read in the codon parameters
	ChromosomeFile.read((char*)&ResultLocation->CodonIndex, 4);
	if (ChromosomeFile.gcount() != 4) {
		return GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Source, 1);
	if (ChromosomeFile.gcount() != 1) {
		return GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Target, 1);
	if (ChromosomeFile.gcount() != 1) {
		return GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Types, 1);
	if (ChromosomeFile.gcount() != 1) {
		return GENE_MISSING_CONNECTIONS;
	}
	ChromosomeFile.read((char*)&ResultLocation->Weight, 4);
	if (ChromosomeFile.gcount() != 4) {
		return GENE_MISSING_CONNECTIONS;
	}

	return GENE_SUCCESS;
}
