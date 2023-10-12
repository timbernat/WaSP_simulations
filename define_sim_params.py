'''Record script for regenerating simulation parameter files from scratch'''

import warnings
warnings.catch_warnings(record=True)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

from pathlib import Path
from openmm.unit import nanosecond, picosecond, femtosecond
from openmm.unit import kelvin, atmosphere
from polysaccharide2.openmmtools.parameters import IntegratorParameters, ThermoParameters, ReporterParameters, SimulationParameters


# DEFINE OUTPUT DIRECTORY HERE
params_dir = Path('sim_param_sets')

# DEFINE PARAMETERS HERE
all_params = {
    'anneal' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtosecond,
            # total_time=1*nanosecond,
            total_time=50*picosecond,
            num_samples=100
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=600*kelvin,
            pressure=1*atmosphere,
            friction_coeff=1/picosecond,
            barostat_freq=25
        ),
        reporter_params=ReporterParameters(
            report_state=True,
            report_checkpoint=True,
            report_state_data=True,
            report_trajectory=True,
        ),
    ),
    'equilibration' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtosecond,
            total_time=50*picosecond,
            num_samples=100
        ),
        thermo_params=ThermoParameters(
            ensemble='NPT',
            temperature=300*kelvin,
            pressure=1*atmosphere,
            friction_coeff=1/picosecond,
            barostat_freq=100
        ),
        reporter_params=ReporterParameters(
            report_state=True,
            report_checkpoint=True,
            report_state_data=True,
            report_trajectory=True,
        ),
    ),
    'production_lite' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtosecond,
            total_time=100*picosecond,
            num_samples=50
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=300*kelvin,
            pressure=1*atmosphere, # doesn't actually matter for NVT
            friction_coeff=1/picosecond,
            barostat_freq=100      # doesn't actually matter for NVT
        ),
        reporter_params=ReporterParameters(
            report_state=True,
            report_checkpoint=True,
            report_state_data=True,
            report_trajectory=True,
        ),
    ),
    'production' : SimulationParameters(
        integ_params=IntegratorParameters(
            time_step=2*femtosecond,
            total_time=100*nanosecond,
            num_samples=5000
        ),
        thermo_params=ThermoParameters(
            ensemble='NVT',
            temperature=300*kelvin,
            pressure=1*atmosphere, # doesn't actually matter for NVT
            friction_coeff=1/picosecond,
            barostat_freq=100      # doesn't actually matter for NVT
        ),
        reporter_params=ReporterParameters(
            report_state=True,
            report_checkpoint=True,
            report_state_data=True,
            report_trajectory=True,
        ),
    ),
}

if __name__ == '__main__':
    params_dir.mkdir(exist_ok=True)
    for name, param_set in all_params.items():
        param_set.to_file(params_dir / f'{name}_params.json')