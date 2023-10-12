'''Add solvent to a Topology'''

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

import numpy as np
from pathlib import Path

import openmm.unit
from openmm.unit import gram, centimeter

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.genutils.fileutils.pathutils import assemble_path

from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.solvation import solvents, boxvectors
from polysaccharide2.openfftools.solvation.packing import pack_topology_with_solvent


def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved', type=Path, required=True)
    parser.add_argument('-ldir', '--logging_directory'  , help='Directory into which log files should be output', type=Path)
    parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
    parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)
    
    parser.add_argument('-solv', '--solvent'            , help='The name of the solvent to use to pack the box', type=str, default='water_TIP3P')
    parser.add_argument('-rho' , '--density'            , help='Density (in grams/cm**3) up to which the desired box should be packed', type=float, default=1.0)

    parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', nargs=3, type=float)
    parser.add_argument('-bdu' , '--box_dimension_unit' , help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
    parser.add_argument('-exc' , '--exclusion'          , help='Name of the ForceField XML to use for MD', type=float, default=1.0)
    parser.add_argument('-excu', '--exclusion_unit'     , help='Unit to use when assigning box exlusion (default nanometer)', type=str, default='nanometer')

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    if not args.logging_directory:
        args.logging_directory = args.working_directory
    else:
        args.logging_directory.mkdir(exist_ok=True)
    
    assert(args.sdf_path.suffix == '.sdf')
    args.density = args.density * (gram/centimeter**3)
    args.solvent = getattr(solvents, args.solvent)

    box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
    box_dims = np.array(args.box_dimensions) * box_dim_unit
    args.box_vecs = boxvectors.box_vectors_flexible(box_dims)

    exclusion_unit = getattr(openmm.unit, args.exclusion_unit)
    args.exclusion = args.exclusion * exclusion_unit

    return args

def main() -> None:
    '''Define main code body'''
    args = parse_args()
    with MSFHandlerFlex(args.logging_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        offmol = topology.get_largest_offmol(offtop)
        mol_name = offmol.name

        solv_top = pack_topology_with_solvent(offtop, args.solvent, box_vecs=args.box_vecs, density=args.density, exclusion=args.exclusion)
        solv_top_path = assemble_path(args.working_directory, mol_name, extension='sdf', postfix=f'solv_{args.solvent.name}')
        print(str(solv_top_path)) # for shell return-catching
        topology.topology_to_sdf(solv_top_path, solv_top)

if __name__ == '__main__':
    main()