#!/bin/bash

echo $PATH > fake_path.txt

# ============================================================================================================================
# DESC : orchestrates parallel SBATCH dispatch of entire chemistry- charge- and replicate-matrix simulation battery
# Accepts no parameters, but rather defined them locally for permanence
# ============================================================================================================================

# ============================================================================================================================
# parameter definitions
# ============================================================================================================================
slurm_log_dir='slurm_logs'
root_dir='sbatch_test'
clear_existing=true

struct_src='polymer_structures'
pdb_src="$struct_src/pdb"
mono_src="$struct_src/monomers"

mol_names=(
    'paam_modified'
    # 'peg_modified'
    # 'pnipam_modified'
)

# charging methods
N_max=150
rct_charge_method='Espaloma-AM1-BCC'
# rct_charge_method='AM1-BCC-ELF10'
charge_methods=(
    # 'AM1-BCC-ELF10'
    # 'Espaloma-AM1-BCC'
    'RCT'
)

# periodic box size
box_x=6.0
box_y=6.0
box_z=6.0
box_unit='nanometer'

# simulation parameter sets
num_confs=1
forcefield='openff-2.0.0'

param_dir='sim_param_sets'
anneal_param_path=$param_dir/anneal_params.json
equil_param_path=$param_dir/equilibration_params.json
prod_param_path=$param_dir/production_lite_params.json

# solvation 
solvent='water_TIP3P'
density_gcm3=0.997
exclusion=1.0
exclusion_unit='nanometer'

# RDF radii
rmin=0.0
rmax=2.0
rad_unit='nanometer'

# ============================================================================================================================
# setting up directories and files
# ============================================================================================================================
if [ -e $root_dir ] && [ $clear_existing == true ]; then
    rm -r $root_dir
    echo "Removed existing directory $root_dir"
fi

mkdir -p $root_dir # make root directory (exists ok)
echo "Touched directory $root_dir"

mkdir -p $slurm_log_dir

# ============================================================================================================================
# code execution
# ============================================================================================================================
for mol_name in "${mol_names[@]}"; do # iterate over list of names
    mol_dir=$root_dir/$mol_name
    root_log_dir="$mol_dir/Logs"

    param_args="$mol_name $mol_dir $root_log_dir $pdb_src/$mol_name.pdb $mono_src/$mol_name.json $N_max $rct_charge_method"
    echo $param_args
    job_id_param=$(sbatch --parsable --job-name par_$mol_name --output $slurm_log_dir/param_$mol_name.log param.job $param_args )
    echo $job_id_param

    # for charge_method in "${charge_methods[@]}"; do
    #     charge_dir=$mol_dir/$charge_method
    #     charge_log_dir="$charge_dir/Logs"

    #     charge_path=$(python -m components.charge_mol -wdir $charge_dir -ldir $charge_log_dir  -sdf $sdf_path -cmet $charge_method -lc $lib_chg_path)

    #     for ((i=1; i<=$num_confs; i++)); do # can't use {1..$num_confs} (variable upper bound returns curly braces as literal)
    #         conf_name="conformer_$i"
    #         conf_dir=$charge_dir/$conf_name
    #         conf_log_dir="$conf_dir/Logs"

    #         conf_path=$(python -m components.anneal  -wdir $conf_dir -ldir $conf_log_dir -sdf $charge_path -a $conf_name -bd $box_x $box_y $box_z -bdu $box_unit -sp $anneal_param_path -ff $forcefield)
    #         solv_path=$(python -m components.solvate -wdir $conf_dir -ldir $conf_log_dir -sdf $conf_path   -a $conf_name -solv $solvent -rho $density_gcm3 -bd $box_x $box_y $box_z -bdu $box_unit -exc $exclusion -excu $exclusion_unit)
    #         python -m components.simulate -wdir $conf_dir -ldir $conf_log_dir -sdf $solv_path -a $conf_name -bd $box_x $box_y $box_z -bdu $box_unit -sp "equilibration=$equil_param_path" "production=$prod_param_path" -ff $forcefield -az 'production' -rmin $rmin -rmax $rmax -ru $rad_unit
    #     done
    # done
done
