'''For generating library charges according to the RCT strategy for a given monomer chemistry'''

import warnings
warnings.catch_warnings(record=True)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import argparse
import logging
from polymerist.genutils.logutils.IOHandlers import LOG_FORMATTER
logging.basicConfig(
    level=logging.INFO,
    format =LOG_FORMATTER._fmt,
    datefmt=LOG_FORMATTER.datefmt,
    force=True
)
LOGGER = logging.getLogger(__name__)

from pathlib import Path
from openff.toolkit import Molecule, Topology

from polymerist.genutils.decorators.functional import allow_string_paths
from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.genutils.logutils.IOHandlers import MSFHandlerFlex

from polymerist.monomers import MonomerGroup
from polymerist.polymers import estimation, building
from polymerist.polymers.exceptions import MorphologyError

from polymerist.residues.rescharge.calculation import compute_residue_charges
from polymerist.residues.rescharge.rctypes import ChargesByResidue

from polymerist.openfftools.pcharge import MolCharger
from polymerist.openfftools import topology, TKREGS
from polymerist.residues.partition import partition


def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved', type=Path, required=True)
    parser.add_argument('-ldir', '--logging_directory'  , help='Directory into which log files should be output', type=Path)
    parser.add_argument('-mono', '--monomer_path'     , help='Path to JSON file from which to read monomer chemistry', type=Path)
    parser.add_argument('-n'   , '--mol_name'         , help='Name to assign to chemically-specified molecule'       , type=str)
    parser.add_argument('-a'   , '--affix'            , help='Additional text to attach to the end of the molecule name provided when naming ', type=str, default='reduced')

    parser.add_argument('-cmet', '--charging_method'  , help='Base charge assignment method from which to extract library charges', type=str)
    parser.add_argument('-N'   , '--chain_length'     , help='Maximum number of atoms to include in reduced chain'   , type=int, default=150)
    parser.add_argument('-k'   , '--keep_pdb'         , help='If set, will preserve the intermediary PDB file for the reduction', action='store_true')

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    if not args.logging_directory:
        args.logging_directory = args.working_directory
    else:
        args.logging_directory.mkdir(exist_ok=True)
    
    assert(args.monomer_path.suffix == '.json')
    assert(args.charging_method in MolCharger.subclass_registry)
    args.delete_pdb = not args.keep_pdb

    return args

@allow_string_paths
def rct_protocol(working_dir : Path, mol_name : str, monogroup : MonomerGroup, N : int, charger : MolCharger, delete_pdb : bool=True) -> tuple[ChargesByResidue, Molecule]:
    '''
    Generates library charges for a monomer group given terminal group head-tail orientations and a maximum chain length
    If delete_pdb=True, will remove the working pdb path after charges are generated
    as currently implemented, only supports monomer groups which constitute a linear homopolymer
    '''
    if not monogroup.is_linear:
        raise MorphologyError('RCT currently only supports linear homopolymers')
    
    DOP = estimation.estimate_DOP_lower(monogroup, max_chain_len=N)
    chain = building.build_linear_polymer(monogroup, DOP)

    pdb_path = assemble_path(working_dir, mol_name, extension='pdb')
    building.mbmol_to_openmm_pdb(pdb_path, chain) # output PDB file to disc 
    offtop = Topology.from_pdb(pdb_path, _custom_substructures=monogroup.monomers, toolkit_registry=TKREGS['The RDKit']) # load custom substructures - raises error if PDB has issues
    if delete_pdb:
        pdb_path.unlink()
        
    was_partitioned = partition(offtop) 
    assert(was_partitioned)
    offmol = topology.get_largest_offmol(offtop)
    offmol.name = mol_name

    cmol = charger.charge_molecule(offmol, in_place=False)
    res_chgs = compute_residue_charges(cmol, monogroup)

    return res_chgs, cmol

def main() -> None:
    '''Define main code body'''
    args=parse_args()
    with MSFHandlerFlex(args.logging_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        chgr = MolCharger.subclass_registry[args.charging_method]()
        monogroup = MonomerGroup.from_file(args.monomer_path)
        lib_chgs, cmol_redux = rct_protocol(args.working_directory, args.mol_name, monogroup, args.chain_length, charger=chgr, delete_pdb=args.delete_pdb)
        mol_name_redux = f'{cmol_redux.name}{"_" if args.affix else ""}{args.affix}'
        cmol_redux.name = mol_name_redux

        lib_chg_path = assemble_path(args.working_directory, args.mol_name, extension='json', postfix='residue_charges')
        lib_chgs.to_file(lib_chg_path)

        # save small molecule with explicit charges for benchmark
        sdf_path_exact = assemble_path(args.working_directory, mol_name_redux, extension='sdf', postfix=args.charging_method)
        topology.topology_to_sdf(sdf_path_exact, cmol_redux.to_topology())
        print(f'{lib_chg_path} {sdf_path_exact}')

# main code
if __name__ == '__main__':
    main()