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



#-----SUBNETWORK CHROMOSOME-----
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
        SUCCESS =        0
        BAD_PATH =       1
        FILE_SHORT =     2
        MISSING_GENES =  3
        BAD_GENE =       4

    def __init__(self, c_result, path="No path provided"):
        self.genes = []
        self.gene_count = int(c_result.GeneCount)
        self.return_code = SubnetworkChromosomeParseResult.Retcodes(c_result.ReturnCode)
        self.additional_info = int(c_result.AdditionalInfo)
        
        if self.return_code == SubnetworkChromosomeParseResult.Retcodes.SUCCESS:
            #ctypes pointers can be used as arrays, so to get the whole array
            #we take it as a slice to prevent going out of bounds
            for c_gene in c_result.Genes[:self.gene_count]:
                self.genes.append(SubnetworkGene(c_gene))

            logger.debug(f"Finished parsing subnetwork chromosome \"{path}\": {self.gene_count} genes found, additional info = {self.additional_info}")
                
        else:
            logger.error(f"Error parsing subnetwork chromosome \"{path}\": {self.return_code}. {self.gene_count} genes found, additional info: {self.additional_info}")
         
    def __str__(self):
        return f"""SubnetworkChromosomeParseResult: 
\tReturn code = {self.return_code}
\tAdditional info = {self.additional_info}
\tGene count = {self.gene_count} ({len(self.genes) = })
\tFirst gene slice = {self.genes[:1]}"""


#now that the types have been defined, set up our library functions
logger.info("Preparing subnetwork chromosome library functions")
__parse_subnetwork_chromosome = __network_parser_dll.ParseSubnetworkChromosome
__parse_subnetwork_chromosome.restype = C_SubnetworkChromosomeParseResult
__parse_subnetwork_chromosome.argtypes = [
    ctypes.c_char_p 
]

__free_subnetwork_parse_result = __network_parser_dll.FreeSubnetworkParseResult
__free_subnetwork_parse_result.argtypes = [
    C_SubnetworkChromosomeParseResult
]
logger.info("Prepared subnetwork chromosome library functions")



#-----QUADRANT CHROMOSOME-----
class C_QuadrantChromosomeParseResult(ctypes.Structure):
    _fields_ = [
        ("SubnetworkIndices", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint32))),
        ("QuadrantCount", ctypes.c_uint32),
        ("SubnetworksPerQuadrant", ctypes.c_uint32),
        ("ReturnCode", ctypes.c_uint32)
    ]


class QuadrantChromosomeParseResult:
    class Retcodes(Enum):
        SUCCESS =				0
        BAD_PATH =  			1
        FILE_SHORT =			2
        MISSING_QUADRANTS = 	3

    def __init__(self, c_result, path="No path provided"):
        self.quadrants = []
        self.quadrant_count = int(c_result.QuadrantCount)
        self.subnetworks_per_quadrant = int(c_result.SubnetworksPerQuadrant)
        self.return_code = QuadrantChromosomeParseResult.Retcodes(c_result.ReturnCode)

        if self.return_code == QuadrantChromosomeParseResult.Retcodes.SUCCESS:
            self.quadrants = [
                [c_result.SubnetworkIndices[i][j] for j in range(self.subnetworks_per_quadrant)] 
                for i in range(self.quadrant_count)
            ]
            logger.debug(f"Finished parsing subnetwork chromosome \"{path}\": {self.quadrant_count} quadrants found, {self.subnetworks_per_quadrant} subnetworks per quadrant")
        else:
            logger.error(f"Error parsing quadrant chromosome \"{path}\": {self.return_code}. {self.quadrant_count} quadrants found, {self.subnetworks_per_quadrant} subnetworks per quadrant")


    def __str__(self):
        return f"""QuadrantChromosomeParseResult:
Return code: {self.return_code}
Subnetworks per quadrant: {self.subnetworks_per_quadrant}
Quadrant count: {self.quadrant_count}
First quadrant slice: {self.quadrants[:1]}
"""


logger.info("Preparing quadrant chromosome library functions")
__parse_quadrant_chromosome = __network_parser_dll.ParseQuadrantChromosome
__parse_quadrant_chromosome.restype = C_QuadrantChromosomeParseResult
__parse_quadrant_chromosome.argtypes = [
    ctypes.c_char_p
]

__free_quadrant_parse_result = __network_parser_dll.FreeQuadrantParseResult
__free_quadrant_parse_result.argtypes = [
    C_QuadrantChromosomeParseResult
]
logger.info("Prepared quadrant chromosome library functions")



#TODO should these be staticmethods?
#-----HELPER FUNCTIONS-----
def parse_subnetwork_chromosome(path):
    logger.debug(f"Attempting to parse subnetwork chromosome \"{path}\"")
    c_result = __parse_subnetwork_chromosome(bytes(path, "ASCII"))
    result = SubnetworkChromosomeParseResult(c_result, path=path)

    logger.debug(f"Freeing C_SubnetworkChromosomeParseResult for \"{path}\"")
    __free_subnetwork_parse_result(c_result)
    logger.debug(f"Freed C_SubnetworkChromosomeParseResult for \"{path}\"")

    return result

def parse_quadrant_chromosome(path):
    logger.debug(f"Attempting to parse quadrant chromosome \"{path}\"")
    c_result = __parse_quadrant_chromosome(bytes(path, "ASCII"))
    result = QuadrantChromosomeParseResult(c_result, path=path)

    logger.debug(f"Freeing C_QuadrantChromosomeParseResult for \"{path}\"")
    __free_quadrant_parse_result(c_result)
    logger.debug(f"Freed C_QuadrantChromosomeParseResult for \"{path}\"")

    return result


if __name__ == "__main__":
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(filename="test_log.log", level=logging.DEBUG, format=FORMAT)

    #test subnetwork chromosome parsing
    print("----- Testing subnetwork parsing -----")
    test_resource_dir = "Tests/subnetworks/"
    test_files = ["Bad path.chr", "Short.chr", "Missing genes.chr", "Bad genes.chr", "Success.chr"]
    for file in test_files:
        print(file)
        print(str(parse_subnetwork_chromosome(test_resource_dir + file)))
        print()

    print()

    #test quadrant chromosome parsing
    print("----- Testing quadrant parsing -----")
    test_resource_dir = "Tests/quadrants/"
    test_files = ["Bad path.chr", "Short.chr", "Missing quadrants.chr", "Success.chr"]
    for file in test_files:
        print(file)
        print(str(parse_quadrant_chromosome(test_resource_dir + file)))
        print()