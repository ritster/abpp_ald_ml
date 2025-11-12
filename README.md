# abpp_ald_ml
Pre-processing and ML workflow code for chemoproteomics reactivity data.

TODO: choose a better name for the repository

## Setup Instructions

### 1. Clone the Repository

```bash
git clone git@github.com:ritster/abpp_ald_ml.git
cd abpp_ald_ml
```

For detailed instructions on setting up ssh keys for GitHub use, see [these docs](https://docs.github.com/en/authentication/connecting-to-github-with-ssh).

If you simply want to view and run repository contents (not propose changes), you can clone via the web URL:

```bash
git clone https://github.com/ritster/abpp_ald_ml.git
cd abpp_ald_ml
```

### 2. Create and Activate the Conda Environment

This project uses Conda for environment management to ensure reproducibility across systems.

```bash
conda env create -f env_conda_specific.yaml
conda activate chemoML
```

If you don’t have Conda installed, you can get it via [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install).

### 3. Install the Package in Editable Mode

```bash
pip install -e .
```
