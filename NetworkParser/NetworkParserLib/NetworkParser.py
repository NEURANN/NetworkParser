import ctypes
from enum import Enum

import logging
logger = logging.getLogger("NetworkParser")

SOURCE_TYPE_INPUT =   0b00000001
SOURCE_TYPE_HIDDEN =  0b00000010
TARGET_TYPE_HIDDEN =  0b00000100
TARGET_TYPE_OUTPUT =  0b00001000


#load the dll on import
logger.info("Loading library file")
__network_parser_dll = ctypes.cdll.LoadLibrary("./NetworkParser.dll")
logger.info("Library file loaded.")

class C_InternalConnectionCodon(ctypes.Structure):
    _fields_ = [
        ("CodonIndex", ctypes.c_uint32),
        ("Weight", ctypes.c_float),
        ("Source", ctypes.c_uint8),
        ("Target", ctypes.c_uint8),
        ("Types", ctypes.c_uint8)
    ]   

class InternalConnectionCodon:
    def __init__(self, c_result):
        self.codon_index = int(c_result.CodonIndex)
        self.weight = float(c_result.Weight)
        self.source = int(c_result.Source)
        self.target = int(c_result.Target)
        self.types = int(c_result.Types)
        
        #separates out the SOURCE_ and TARGET_ types.
        #                  bitmask for bits we care about          bits we're testing against
        self.source_type = (SOURCE_TYPE_INPUT | SOURCE_TYPE_HIDDEN) & self.types
        self.target_type = (TARGET_TYPE_HIDDEN | TARGET_TYPE_OUTPUT) & self.types

    def __str__(self):
        if self.source_type == SOURCE_TYPE_INPUT:
            source_str = "Source type = Input"
        else:
            source_str = "Source type = Hidden"

        if self.target_type == TARGET_TYPE_HIDDEN:
            target_str = "Target type = Hidden"
        else:
            target_str = "Target type = Output"

        return f"""InternalConnectionCodon:
Codon index = {self.codon_index}
\t{source_str}, Source index = {self.source}
\t{target_str}, Target index = {self.target}
\tWeight = {self.weight}
\tTypes = {self.types :08b}
"""
    
    def __repr__(self):
        source_char = "i" if self.source_type == SOURCE_TYPE_INPUT else "h"
        target_char = "h" if self.target_type == TARGET_TYPE_HIDDEN else "o"
        return f"{self.codon_index}: {source_char}{self.source} -- {self.weight:.2f}--> {target_char}{self.target}"
    

class C_SubnetworkGene(ctypes.Structure):
    _fields_ = [
        ("GeneIndex", ctypes.c_uint64),
        ("Codons", ctypes.POINTER(C_InternalConnectionCodon)),
        ("CodonCount", ctypes.c_uint32)
    ] 

class SubnetworkGene:
    def __init__(self, c_result):
        self.gene_index = int(c_result.GeneIndex)
        self.codon_count = int(c_result.CodonCount)
        self.codons = []

        #ctypes pointers can be used as arrays, so to get the whole array
        #we take it as a slice to prevent going out of bounds
        for codon in c_result.Codons[:self.codon_count]:
            self.codons.append(InternalConnectionCodon(codon))
    
    def __str__(self):
        return f"""SubnetworkGene:
\tGene index = {self.gene_index}
\tCodon count = {self.codon_count}
\tCodons[:3] = {self.codons[:3]}
"""
    
    def __repr__(self):
        return f"{self.gene_index}, {self.codon_count} codons: {self.codons[:3]}"


class C_SubnetworkChromosomeParseResult(ctypes.Structure):
    _fields_ = [
        ("Genes", ctypes.POINTER(C_SubnetworkGene)),
        ("GeneCount", ctypes.c_uint32),
        ("ReturnCode", ctypes.c_uint),
        ("AdditionalInfo", ctypes.c_int)
    ]

class SubnetworkChromosomeParseResult:
    class Retcodes(Enum):
        CHROMOSOME_SUCCESS =        0
        CHROMOSOME_BAD_PATH =       1
        CHROMOSOME_FILE_SHORT =     2
        CHROMOSOME_MISSING_GENES =  3
        CHROMOSOME_BAD_GENE =       4

    def __init__(self, c_result, path="No path provided"):
        self.genes = []
        self.gene_count = int(c_result.GeneCount)
        self.return_code = int(c_result.ReturnCode)
        self.additional_info = int(c_result.AdditionalInfo)
        
        if c_result.ReturnCode == SubnetworkChromosomeParseResult.Retcodes.CHROMOSOME_SUCCESS.value:
            #ctypes pointers can be used as arrays, so to get the whole array
            #we take it as a slice to prevent going out of bounds
            for c_gene in c_result.Genes[:self.gene_count]:
                self.genes.append(SubnetworkGene(c_gene))

            logger.debug(f"Finished parsing \"{path}\": {self.gene_count} genes found, additional info = {self.additional_info}")
                
        else:
            ret_str = str(SubnetworkChromosomeParseResult.Retcodes(c_result.ReturnCode))
            logger.error(f"Error parsing \"{path}\": {ret_str}. {self.gene_count} genes found, additional info: {self.additional_info}")
         
    def __str__(self):
        return f"""SubnetworkChromosomeParseResult: 
\tFirst gene slice = {self.genes[:1]},
\tGene count = {self.gene_count} ({len(self.genes) = })
\tReturn code = {self.return_code}
\tAdditional info = {self.additional_info}"""


#now that the types have been defined, set up our library functions
logger.info("Preparing library functions")
__parse_subnetwork_chromosome = __network_parser_dll.ParseSubnetworkChromosome
__parse_subnetwork_chromosome.restype = C_SubnetworkChromosomeParseResult
__parse_subnetwork_chromosome.argtypes = [
    ctypes.c_char_p 
]

__free_parse_result = __network_parser_dll.FreeParseResult
__free_parse_result.argtypes = [
    C_SubnetworkChromosomeParseResult
]
logger.info("Prepared library functions")

def parse_subnetwork_chromosome(path):
    logger.debug(f"Attempting to parse\"{path}\"")
    c_result = __parse_subnetwork_chromosome(bytes(path, "ASCII"))
    result = SubnetworkChromosomeParseResult(c_result, path=path)

    logger.debug(f"Freeing C_SubnetworkChromosomeParseResult for \"{path}\"")
    __free_parse_result(c_result)
    logger.debug(f"Freed C_SubnetworkChromosomeParseResult for \"{path}\"")
    return result



if __name__ == "__main__":
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(filename="test_log.log", level=logging.DEBUG, format=FORMAT)

    test_resource_dir = "Tests/"
    test_files = ["Doesn't exist.chr", "Short.chr", "Missing genes.chr", "Bad genes.chr", "Success.chr"]
    for file in test_files:
        print(file)
        print(str(parse_subnetwork_chromosome(test_resource_dir + file)))
        print()