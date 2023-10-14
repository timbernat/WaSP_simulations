#!/bin/sh

# ============================================================================================================================
# DESC : orchestrates parallel SBATCH dispatch of entire chemistry- charge- and replicate-matrix simulation battery
# Accepts no parameters, but rather defined them locally for permanence
# ============================================================================================================================

# ============================================================================================================================
# parameter definitions
# ============================================================================================================================
slurm_log_dir='slurm_logs'
root_dir='wasp_sims'
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
    'Espaloma-AM1-BCC'
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
    ioparam=$mol_dir/argio.txt # temporary text file for passing data between sbatch calls

    param_args="$mol_dir $root_log_dir $ioparam $mol_name $pdb_src/$mol_name.pdb $mono_src/$mol_name.json $N_max $rct_charge_method"
    job_id_param=$(sbatch --parsable --job-name par_$mol_name --output ${slurm_log_dir}/param_${mol_name}.log param.job ${param_args} )

    for charge_method in "${charge_methods[@]}"; do
        charge_dir=$mol_dir/$charge_method
        charge_log_dir="$charge_dir/Logs"
        iocharge=$charge_dir/argio.txt

        charge_args="$charge_dir $charge_log_dir $ioparam $iocharge $charge_method"
        job_id_charge=$(sbatch --parsable --dependency "afterok:$job_id_param" --job-name $charge_method --output ${slurm_log_dir}/charge_${mol_name}_${charge_method}.log charge.job $charge_args)

        for ((i=1; i<=$num_confs; i++)); do # can't use {1..$num_confs} (variable upper bound returns curly braces as literal)
            conf_name="conf${i}"
            conf_dir=$charge_dir/$conf_name
            conf_log_dir="$conf_dir/Logs"

            conf_args="$conf_dir $conf_log_dir $iocharge $conf_name $box_x $box_y $box_z $box_unit $density_gcm3 $solvent $exclusion $exclusion_unit $forcefield $anneal_param_path $equil_param_path $prod_param_path $rmin $rmax $rad_unit"
            echo $conf_args
            job_id_conf=$(sbatch --dependency "afterok:$job_id_charge" --job-name ${conf_name}${mol_name} --output ${slurm_log_dir}/${mol_name}_${charge_method}_${conf_name}.log sims.job $conf_args) # not actually necessary to catch job ID (no subsequent jobs), just keeps logging clean
        done
    done
done
