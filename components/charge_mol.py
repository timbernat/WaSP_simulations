'''Assign charges to a chemically-defined Topology by a method of choice'''

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

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.genutils.fileutils.pathutils import assemble_path
from polysaccharide2.openfftools import topology
from polysaccharide2.openfftools.pcharge import MolCharger
from polysaccharide2.residues.rescharge.rctypes import ChargesByResidue


def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved', type=Path, required=True)
    parser.add_argument('-ldir', '--logging_directory'  , help='Directory into which log files should be output', type=Path)
    parser.add_argument('-sdf' , '--sdf_path'         , help='Path to SDF file from which to read parameterized molecule structure', type=Path)

    parser.add_argument('-cmet', '--charging_method'      , help='Charge assignment method to use to compute charges', type=str)
    parser.add_argument('-lc'  , '--library_charge_path'  , help='Path to a JSON file containing library charges for the desired molecule', type=Path)

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    if not args.logging_directory:
        args.logging_directory = args.working_directory
    else:
        args.logging_directory.mkdir(exist_ok=True)

    assert(args.sdf_path.suffix == '.sdf')
    assert(args.charging_method in MolCharger.subclass_registry)
    assert(args.library_charge_path.suffix == '.json')

    return args

def main() -> None:
    '''Define main code body'''
    args = parse_args()
    with MSFHandlerFlex(args.logging_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        offmol = topology.get_largest_offmol(offtop)
        mol_name = offmol.name

        ChargerType = MolCharger.subclass_registry[args.charging_method]
        if args.charging_method == 'RCT':
            if (args.library_charge_path is not None) and (args.library_charge_path.exists()):
                charger = ChargerType(ChargesByResidue.from_file(args.library_charge_path))
            else:
                raise FileExistsError(f'Invalid library charge file "{args.library_charge_file}" provided')
        else:
            charger = ChargerType()

        cmol = charger.charge_molecule(offmol, in_place=False)
        charged_top_path = assemble_path(args.working_directory, mol_name, extension='sdf', postfix=args.charging_method)
        print(str(charged_top_path)) # for shell return-catching
        topology.topology_to_sdf(charged_top_path, cmol.to_topology())

if __name__ == '__main__':
    main()