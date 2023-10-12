'''For generating chemically-specified SDF files from PDB and monomer JSON structure files'''

import warnings
warnings.catch_warnings(record=True)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import argparse
import logging
from polysaccharide2.genutils.logutils.IOHandlers import LOG_FORMATTER
logging.basicConfig(
    level=logging.INFO,
    format =LOG_FORMATTER._fmt,
    datefmt=LOG_FORMATTER.datefmt,
    force=True
)
LOGGER = logging.getLogger(__name__)

from pathlib import Path
from openff.toolkit import Topology

from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.openfftools import TKREGS
from polysaccharide2.openfftools import topology
from polysaccharide2.monomers import MonomerGroup
from polysaccharide2.residues.partition import partition


def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved', type=Path, required=True)
    parser.add_argument('-ldir', '--logging_directory'  , help='Directory into which log files should be output', type=Path)
    parser.add_argument('-pdb' , '--pdb_path'         , help='Path to PDB file from which to read molecule structure', type=Path)
    parser.add_argument('-mono', '--monomer_path'     , help='Path to JSON file from which to read monomer chemistry', type=Path)
    parser.add_argument('-n'   , '--mol_name'         , help='Name to assign to chemically-specified molecule'       , type=str)

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    if not args.logging_directory:
        args.logging_directory = args.working_directory
    else:
        args.logging_directory.mkdir(exist_ok=True)
    
    assert(args.pdb_path.suffix == '.pdb')
    assert(args.monomer_path.suffix == '.json')

    return args

def main() -> None:
    '''Define main code body'''
    args = parse_args()
    with MSFHandlerFlex(args.logging_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        LOGGER.info(f'Assigning chemistry for molecule "{args.mol_name}"')
        # generate substructure cover and partition
        monogrp = MonomerGroup.from_file(args.monomer_path)
        offtop = Topology.from_pdb(args.pdb_path, _custom_substructures=monogrp.monomers, toolkit_registry=TKREGS['The RDKit'])
        was_partitioned = partition(offtop)
        assert(was_partitioned)

        # assign properties to Molecule
        offmol = topology.get_largest_offmol(offtop)
        offmol.name = args.mol_name
        offmol.properties['solvent'] = None
        offmol.properties['charge_method'] = None

        # save partitioned Topology
        mol_path = assemble_path(args.working_directory, args.mol_name, extension='sdf')
        print(str(mol_path))
        topology.topology_to_sdf(mol_path, offtop=offtop, toolkit_registry=TKREGS['The RDKit'])

# main code
if __name__ == '__main__':
    main()