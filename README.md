# WaSP_simulations

A repository for reproducing the **Wa**ter **S**oluble **P**olymer (WaSP) molecular dynamics simulations implemented with the OpenFf-based parameterization approaches described in our recent paper(https://chemrxiv.org/engage/chemrxiv/article-details/652f6614bda59ceb9acf6393)

## Installation from GitHub with Conda
First create a local clone of this repository and, within it, a clone of the in-house toolkit this repo depends on:
the toolkit, first download and install . Then, run:
```sh

git clone https://github.com/shirtsgroup/WaSP_simulations
cd WaSP_simulations
git clone https://github.com/shirtsgroup/polymerist
```

Next, recreate the conda environment (after installing the [Anaconda Distribution](https://www.anaconda.com/download)) from the packaged requirements:
```sh
conda env create -n polymer_env -f reqs.yml
conda activate polymer_env
```

## Configuring Scripts
