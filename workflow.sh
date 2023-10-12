#!/bin/bash
# ===================================================
# Shell-based wrapper for polymer simulation workflow
# ===================================================


# parameter definitions
# ===================================================
root_dir='workflow_test'
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
N=150
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


# setting up directories and files
# ===================================================
if [ -e $root_dir ] && [ $clear_existing == true ]; then
    rm -r $root_dir
    echo "Removed existing directory $root_dir"
fi

mkdir -p $root_dir # make root directory (exists ok)
echo "Touched directory $root_dir"


# actual code
# ===================================================
for mol_name in "${mol_names[@]}"; do
    mol_dir=$root_dir/$mol_name
    root_log_dir="$mol_dir/Logs"
    sdf_path=$(python -m components.assign_chem   -wdir $mol_dir -ldir $root_log_dir -pdb $pdb_src/$mol_name.pdb -mono $mono_src/$mol_name.json -n $mol_name)
    rct_paths=$(python -m components.rct_protocol -wdir $mol_dir -ldir $root_log_dir -mono $mono_src/$mol_name.json -n $mol_name -cmet $rct_charge_method -N $N)
    
    echo ${rct_paths[0]}
    # echo ${rct_paths[1]}

    for charge_method in "${charge_methods[@]}"; do
        charge_dir=$mol_dir/$charge_method
        charge_log_dir="$charge_dir/Logs"
        # TODO : also add RCT charge of reduction
        charge_path=$(python -m components.charge_mol -wdir $charge_dir -ldir $charge_log_dir  -sdf $sdf_path -cmet $charge_method -lc $mol_dir/${mol_name}_residue_charges.json)

        for ((i=1; i<=$num_confs; i++)); do
            conf_name="conformer_$i"
            conf_dir=$charge_dir/$conf_name
            conf_log_dir="$conf_dir/Logs"

            conf_path=$(python -m components.anneal  -wdir $conf_dir -ldir $conf_log_dir -sdf $charge_path -a $conf_name -bd $box_x $box_y $box_z -bdu $box_unit -sp $anneal_param_path -ff $forcefield)
            solv_path=$(python -m components.solvate -wdir $conf_dir -ldir $conf_log_dir -sdf $conf_path   -a $conf_name -solv $solvent -rho $density_gcm3 -bd $box_x $box_y $box_z -bdu $box_unit -exc $exclusion -excu $exclusion_unit)
            python -m components.simulate -wdir $conf_dir -ldir $conf_log_dir -sdf $solv_path -a $conf_name -bd $box_x $box_y $box_z -bdu $box_unit -sp "equilibration=$equil_param_path" "production=$prod_param_path" -ff $forcefield -az 'production' -rmin $rmin -rmax $rmax -ru $rad_unit
        done
    done
done
