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

from pathlib import Path
import numpy as np

import openmm.unit
from openmm.unit import nanometer, Unit

import mdtraj
from polysaccharide2.analysis import mdtrajutils

from polysaccharide2.genutils.logutils.IOHandlers import MSFHandlerFlex
from polysaccharide2.genutils.fileutils.pathutils import assemble_path

from polysaccharide2.openmmtools.serialization import SimulationPaths
from polysaccharide2.openmmtools.parameters import SimulationParameters

from polysaccharide2.openfftools import topology
from polysaccharide2.openmmtools import parameters, execution
from polysaccharide2.openfftools.solvation import boxvectors
from polysaccharide2.openfftools.omminter import openff_topology_to_openmm


class ParseKwargs(argparse.Action):
    '''Parser action for dict-like key-value pairs'''
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

def parse_args() -> argparse.Namespace:
    '''Accept and pre-process arguments from command line'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-wdir', '--working_directory'  , help='Directory into which files generated should be saved', type=Path, required=True)
    parser.add_argument('-ldir', '--logging_directory'  , help='Directory into which log files should be output', type=Path)
    parser.add_argument('-sdf' , '--sdf_path'           , help='Path to SDF file from which to read parameterized molecule structure', type=Path)
    parser.add_argument('-a'   , '--affix'              , help='Additional text to attach to the end of the molecule name provided when naming ', type=str)

    parser.add_argument('-bd'  , '--box_dimensions'     , help='XYZ dimensions of the desired periodic box for the system', nargs=3, type=float)
    parser.add_argument('-bdu' , '--box_dimension_unit' , help='Unit to use when assigning box dimensions (default nanometer)', type=str, default='nanometer')
    
    parser.add_argument('-sp'  , '--sim_param_paths'    , help='Path to the file storing the simulation parameters for the anneal', nargs='*', action=ParseKwargs)
    parser.add_argument('-ff'  , '--forcefield'         , help='Name of the ForceField XML to use for MD', type=str)
    parser.add_argument('-az'  , '--analyze'            , help='List of production steps whose trajectories should be analyzed', nargs='*')

    parser.add_argument('-rmin', '--min_radius'         , help='Lower bound of RDF radii, in specified units (default = nanometer)', type=float, default=0.0)
    parser.add_argument('-rmax', '--max_radius'         , help='Upper bound of RDF radii, in specified units (default = nanometer)', type=float, default=0.0)
    parser.add_argument('-ru'  , '--radius_unit'        , help='The unit to ascribe to the radii used to calculate the RDF', type=str)

    # post-process args as needed
    args = parser.parse_args()
    args.working_directory.mkdir(exist_ok=True)
    if not args.logging_directory:
        args.logging_directory = args.working_directory
    else:
        args.logging_directory.mkdir(exist_ok=True)
    
    assert(args.sdf_path.suffix == '.sdf')
    box_dim_unit = getattr(openmm.unit, args.box_dimension_unit)
    box_dims = np.array(args.box_dimensions) * box_dim_unit
    args.box_vecs = boxvectors.box_vectors_flexible(box_dims)

    args.schedule = {
        step_name : parameters.SimulationParameters.from_file(path_str)
            for step_name, path_str in args.sim_param_paths.items()
    }
    args.radius_unit = getattr(openmm.unit, args.radius_unit)

    return args

def analyze_trajectory(prefix : str, sim_paths : SimulationPaths, min_rad : float=0.0, max_rad : float=2.0, rad_unit : Unit=nanometer) -> None:
    '''Generate property time series' and RDF data from Simulation files'''
    working_dir = sim_paths.paths_path.parent # output files into the same folder the Simulation paths live in
    sim_params = SimulationParameters.from_file(sim_paths.parameters_path)

    traj = mdtraj.load(sim_paths.trajectory_path, top=sim_paths.topology_path)
    traj_no_solv = traj.remove_solvent(inplace=False)

    # computing and saving shape property time series'
    prop_data = mdtrajutils.acquire_time_props(traj_no_solv, time_points=sim_params.integ_params.time_points) 
    LOGGER.info('Computed time series\' from trajectory')
    time_data_path = assemble_path(working_dir, prefix, extension='csv', postfix='time_series')
    prop_data.to_csv(time_data_path, index=False)
    sim_paths.time_data_path = time_data_path

    # RDFs
    ## determine IDs of relevant atom pairs
    pair_dict = {
        'chain O - water O' : traj.top.select_pairs('not water and element O', 'water and element O')
    }
    if 'N' in mdtrajutils.unique_elem_types(traj):
        pair_dict['chain N - water O'] = traj.top.select_pairs('not water and element N', 'water and element O')

    ## computing and saving RDFs 
    rdf_data = mdtrajutils.acquire_rdfs(traj, pair_dict, min_rad=min_rad, max_rad=max_rad, rad_unit=rad_unit)
    LOGGER.info('Computed radial distribution functions (RDFs) from trajectory')
    rdf_data_path = assemble_path(working_dir, prefix, extension='csv', postfix='rdfs')
    rdf_data.to_csv(rdf_data_path, index=False)
    sim_paths.spatial_data_path = rdf_data_path

    # updating path references
    sim_paths.to_file(sim_paths.paths_path) # update simulation paths

def main() -> None:
    '''Define main code body'''
    args = parse_args()
    with MSFHandlerFlex(args.logging_directory, proc_name=Path(__file__).stem, loggers='all') as log_handler:
        offtop = topology.topology_from_sdf(args.sdf_path)
        ommtop, ommsys, ommpos = openff_topology_to_openmm(offtop, forcefield=args.forcefield, box_vecs=args.box_vecs)
        history = execution.run_simulation_schedule(args.working_directory, args.schedule, ommtop, ommsys, ommpos, return_history=True)

        for step_to_analyze in args.analyze:
            sim_paths = history[step_to_analyze]['paths']
            analyze_trajectory(step_to_analyze, sim_paths, min_rad=args.min_radius, max_rad=args.max_radius, rad_unit=args.radius_unit)

if __name__ == '__main__':
    main()