import ctypes
import pathlib

from enum import Enum

import logging
logger = logging.getLogger("NetworkParser")

SOURCE_TYPE_INPUT =   0b00000001
SOURCE_TYPE_HIDDEN =  0b00000010
TARGET_TYPE_HIDDEN =  0b00000100
TARGET_TYPE_OUTPUT =  0b00001000


#load the dll on import
logger.info("Loading library file")
import os.path
libpath = os.path.join(__file__, + "NetworkParser.dll")
__network_parser_dll = ctypes.cdll.LoadLibrary(libpath)
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
    

    @staticmethod
    def from_file(path):
        global _parse_subnetwork_chromosome
        global _free_subnetwork_parse_result


        logger.debug(f"Attempting to parse subnetwork chromosome \"{path}\"")
        c_result = _parse_subnetwork_chromosome(bytes(path, "ASCII"))
        result = SubnetworkChromosomeParseResult(c_result, path=path)

        logger.debug(f"Freeing C_SubnetworkChromosomeParseResult for \"{path}\"")
        _free_subnetwork_parse_result(c_result)
        logger.debug(f"Freed C_SubnetworkChromosomeParseResult for \"{path}\"")

        return result


#now that the types have been defined, set up our library functions
logger.info("Preparing subnetwork chromosome library functions")
_parse_subnetwork_chromosome = __network_parser_dll.ParseSubnetworkChromosome
_parse_subnetwork_chromosome.restype = C_SubnetworkChromosomeParseResult
_parse_subnetwork_chromosome.argtypes = [
    ctypes.c_char_p 
]

_free_subnetwork_parse_result = __network_parser_dll.FreeSubnetworkParseResult
_free_subnetwork_parse_result.argtypes = [
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
            logger.debug(f"Finished parsing quadrant chromosome \"{path}\": {self.quadrant_count} quadrants found, {self.subnetworks_per_quadrant} subnetworks per quadrant")
        
        else:
            logger.error(f"Error parsing quadrant chromosome \"{path}\": {self.return_code}. {self.quadrant_count} quadrants found, {self.subnetworks_per_quadrant} subnetworks per quadrant")


    def __str__(self):
        return f"""QuadrantChromosomeParseResult:
\tReturn code: {self.return_code}
\tSubnetworks per quadrant: {self.subnetworks_per_quadrant}
\tQuadrant count: {self.quadrant_count}
\tFirst quadrant slice: {self.quadrants[:1]}"""
    
    @staticmethod
    def from_file(path):
        global _parse_quadrant_chromosome
        global _free_quadrant_parse_result


        logger.debug(f"Attempting to parse quadrant chromosome \"{path}\"")
        c_result = _parse_quadrant_chromosome(bytes(path, "ASCII"))
        result = QuadrantChromosomeParseResult(c_result, path=path)

        logger.debug(f"Freeing C_QuadrantChromosomeParseResult for \"{path}\"")
        _free_quadrant_parse_result(c_result)
        logger.debug(f"Freed C_QuadrantChromosomeParseResult for \"{path}\"")

        return result


logger.info("Preparing quadrant chromosome library functions")
_parse_quadrant_chromosome = __network_parser_dll.ParseQuadrantChromosome
_parse_quadrant_chromosome.restype = C_QuadrantChromosomeParseResult
_parse_quadrant_chromosome.argtypes = [
    ctypes.c_char_p
]

_free_quadrant_parse_result = __network_parser_dll.FreeQuadrantParseResult
_free_quadrant_parse_result.argtypes = [
    C_QuadrantChromosomeParseResult
]
logger.info("Prepared quadrant chromosome library functions")



#-----CONNECTIONS CHROMOSOME-----
class C_ConnectionGene(ctypes.Structure):
    _fields_ = [
        ("Weight", ctypes.c_float),
        ("SourceSubnetworkIndex", ctypes.c_uint32),
        ("SourceOutputIndex", ctypes.c_uint8),
        ("TargetSubnetworkIndex", ctypes.c_uint32),
        ("TargetInputIndex", ctypes.c_uint8)
    ]

class ConnectionGene:
    def __init__(self, c_result):
        self.weight = float(c_result.Weight)
        self.source_subnetwork_index = int(c_result.SourceSubnetworkIndex)
        self.source_output_index = int(c_result.SourceOutputIndex)
        self.target_subnetwork_index = int(c_result.TargetSubnetworkIndex)
        self.target_input_index = int(c_result.TargetInputIndex)
        
    def __str__(self):
        return f"""ConnectionGene:
\tWeight: {self.weight}
\tSource: {self.source_subnetwork_index}.{self.source_output_index}
\tTarget: {self.target_subnetwork_index}.{self.target_input_index}
"""

    def __repr__(self):
        return f"{self.source_subnetwork_index}[{self.source_output_index}] -- {self.weight}--> {self.target_subnetwork_index}[{self.target_input_index}]"

class C_QuadrantConnections(ctypes.Structure):
    _fields_ = [
        ("ConnectionGenes", ctypes.POINTER(C_ConnectionGene)),
        ("ConnectionGeneCount", ctypes.c_uint32),
        ("SourceQuadrantIndex", ctypes.c_uint32),
        ("TargetQuadrantIndex", ctypes.c_uint32),
    ]

class QuadrantConnections:
    def __init__(self, c_result):
        self.connection_genes = []
        self.connection_gene_count = int(c_result.ConnectionGeneCount)
        self.source_quadrant_index = int(c_result.SourceQuadrantIndex)
        self.target_quadrant_index = int(c_result.TargetQuadrantIndex)
        self.connection_genes = [
            ConnectionGene(c_cg) 
            for c_cg in c_result.ConnectionGenes[:self.connection_gene_count]
        ]

    def __str__(self):
        return f"""QuadrantConnections:
\tSource quadrant: {self.source_quadrant_index}
\tTarget quadrant: {self.target_quadrant_index}
\tConnection gene count: {self.connection_gene_count}
\tConnection genes slice: {self.connection_genes[:3]}
"""
    
    def __repr__(self):
        return f"{self.source_quadrant_index}->{self.target_quadrant_index}: {self.connection_gene_count}@{self.connection_genes[:3]}"

class C_ConnectionsChromosomeParseResult(ctypes.Structure):
    _fields_ = [
        ("QuadrantConnectionsArray", ctypes.POINTER(C_QuadrantConnections)),
        ("QuadrantConnectionsCount", ctypes.c_uint32),
        ("ReturnCode", ctypes.c_uint),
        ("AdditionalInfo", ctypes.c_int),
    ]

class ConnectionsChromosomeParseResult:
    class Retcodes(Enum):
        SUCCESS =           0
        BAD_PATH =          1
        SHORT =             2
        QUADRANT_SHORT =    3

    def __init__(self, c_result, path="No path provided"):
        self.quadrant_connections_array = []
        self.quadrant_connections_count = int(c_result.QuadrantConnectionsCount)
        self.return_code = ConnectionsChromosomeParseResult.Retcodes(c_result.ReturnCode)
        self.additional_info = int(c_result.AdditionalInfo)

        if self.return_code == ConnectionsChromosomeParseResult.Retcodes.SUCCESS:
            self.quadrant_connections_array = [
                QuadrantConnections(c_cr)
                for c_cr in c_result.QuadrantConnectionsArray[:self.quadrant_connections_count]
            ]
            logger.debug(f"Finished parsing connections chromosome \"{path}\": {self.quadrant_connections_count} Quadrant connections found")

        else:
            logger.error(f"Error parsing connections chromosome \"{path}\": {self.return_code}, Additional info = {self.additional_info}")

    def __str__(self):
        return f"""
\tReturn code: {self.return_code}
\tAdditional info: {self.additional_info}
\tQuadrant connections count: {self.quadrant_connections_count}
\tFirst quadrant connections slice: {self.quadrant_connections_array[:1]}
"""

    @staticmethod
    def from_file(path):
        global _parse_connections_chromosome
        global _free_connections_parse_result


        logger.debug(f"Attempting to parse connections chromosome \"{path}\"")
        c_result = _parse_connections_chromosome(bytes(path, "ASCII"))
        result = ConnectionsChromosomeParseResult(c_result, path=path)

        logger.debug(f"Freeing C_ConnectionsChromosomeParseResult for \"{path}\"")
        _free_connections_parse_result(c_result)
        logger.debug(f"Freed C_ConnectionsChromosomeParseResult for \"{path}\"")

        return result


logger.info("Preparing quadrant chromosome library functions")
_parse_connections_chromosome = __network_parser_dll.ParseConnectionsChromosome
_parse_connections_chromosome.restype = C_ConnectionsChromosomeParseResult
_parse_connections_chromosome.argtypes = [
    ctypes.c_char_p
]

_free_connections_parse_result = __network_parser_dll.FreeConnectionsParseResult
_free_connections_parse_result.argtypes = [
    C_ConnectionsChromosomeParseResult
]
logger.info("Prepared quadrant chromosome library functions")



#helper class that stores all aspects of a network's genome
class NetworkGenome:
    class ChromosomeParseException(Exception):
        pass

    def __init__(self, path):
        self.__parse_subnetworks(path)
        self.__parse_quadrants(path)
        self.__parse_connections(path)

    def __parse_subnetworks(self, path):
        subnetworks_path = pathlib.Path.joinpath(path, "subnetworks.chr")
        subnetworks_result = SubnetworkChromosomeParseResult.from_file(subnetworks_path)

        if subnetworks_result.return_code == SubnetworkChromosomeParseResult.Retcodes.SUCCESS:
            self.subnetwork_genes = subnetworks_result.genes
        else:
            exception_message = f"{subnetworks_result.return_code}: Additional info = {subnetworks_result.additional_info}"
            raise NetworkGenome.ChromosomeParseException(exception_message)
        
    def __parse_quadrants(self, path):
        quadrants_path = pathlib.Path.joinpath(path, "quadrants.chr")
        quadrants_result = QuadrantChromosomeParseResult.from_file(quadrants_path)

        if quadrants_result.return_code == QuadrantChromosomeParseResult.Retcodes.SUCCESS:
            self.quadrant_definitions = quadrants_result.quadrants
            self.subnetworks_per_quadrant = quadrants_result.subnetworks_per_quadrant
        else:
            exception_message = f"{quadrants_result.return_code}"
            raise NetworkGenome.ChromosomeParseException(exception_message)

    def __parse_connections(self, path):
        connections_path = pathlib.Path.joinpath(path, "connections.chr")
        connections_result = ConnectionsChromosomeParseResult.from_file(connections_path)

        if connections_result.return_code == ConnectionsChromosomeParseResult.Retcodes.SUCCESS:
            self.quadrant_connections = connections_result.quadrant_connections_array
        else:
            exception_message = f"{connections_result.return_code}: Additional info = {connections_result.additional_info}"
            raise NetworkGenome.ChromosomeParseException(exception_message)
