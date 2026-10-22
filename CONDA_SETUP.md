# Environment setup

The manuscript workflows use Python 3.9.15, NumPy 1.21.6, Biopython 1.79, SciPy 1.8.1 and Matplotlib 3.6.2. Conda installs these dependencies; pip installs the separate Sinabro 3.0.0 checkout.

## Get both repositories

```bash
git clone https://github.com/StudyingAnt/sinabro.git
git clone https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds.git
cd Simple_robustness_analysis_of_cds
```

For an exact code snapshot, use the package commit in [software_manifest.json](software_manifest.json) and record this repository's commit. The publication-oriented package checkout retains the same scientific source fingerprint as the validated calculations. Existing working-branch clones can continue with `git pull --ff-only`.

## Create a named environment

```bash
conda env create --file environment.yml
conda activate sinabro-revision
python -m pip install --no-deps --no-build-isolation --force-reinstall ../sinabro
```

The environment identifier `sinabro-revision` is retained for existing installations and scheduler defaults. It is a location/name, not an analysis mode. A new name may be selected with `conda env create --name YOUR_ENV_NAME --file environment.yml`, followed by activation of that name.

## Install at a chosen absolute path

The named installation above remains valid. When storage policy requires a specific environment directory, use a prefix instead:

```bash
conda env create --prefix /PATH/TO/ENVS/sinabro-revision --file environment.yml
conda activate /PATH/TO/ENVS/sinabro-revision
python -m pip install --no-deps --no-build-isolation --force-reinstall ../sinabro
```

Replace `/PATH/TO/ENVS/` with the desired absolute parent directory. Activate the environment by that same full path. Existing environments do not need recreation after documentation updates.

Use a regular package installation as shown above when the package and analysis checkouts share a parent directory. With some editable installations, the outer `sinabro/` checkout can shadow the installed Python package. Reinstalling with the regular command fixes that import precedence; the scientific dependencies can remain unchanged. From the shared parent directory the equivalent package path is `./sinabro`.

## Verify the installation

```bash
python - <<'PYCODE'
import platform, numpy, Bio, scipy, matplotlib, sinabro
from sinabro import Trajectory, RobustnessComputer
versions = {
    'Python': ('3.9.15', platform.python_version()),
    'NumPy': ('1.21.6', numpy.__version__),
    'Biopython': ('1.79', Bio.__version__),
    'SciPy': ('1.8.1', scipy.__version__),
    'Matplotlib': ('3.6.2', matplotlib.__version__),
    'Sinabro': ('3.0.0', sinabro.__version__),
}
for name, (expected, actual) in versions.items():
    print(name, actual)
    assert actual == expected, (name, expected, actual)
assert sinabro.__file__ is not None
print('Package source:', sinabro.__file__)
PYCODE

python -m unittest discover -s tests -v
```

Tests use bounded scientific/software cases; they do not submit full genome calculations. Package tests are documented in the separate Sinabro repository. CI configuration is under `.github/workflows/`.

## SeaWulf environment activation

Where required by the site's installed module configuration:

```bash
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ABSOLUTE/PATH/TO/YOUR/ENVIRONMENT
export SINABRO_ENV="$CONDA_PREFIX"
export SINABRO_PYTHON="$CONDA_PREFIX/bin/python"
export SINABRO_SOURCE="$(cd ../sinabro && pwd)"
```

Export the environment prefix before submitting the Slurm wrappers, particularly when using a custom installation path. Reference, package and output directories must be visible from compute nodes. The wrappers' standard project layout is described in [the reproduction guide](MANUSCRIPT_REPRODUCTION.md).

## Preserve the resolved environment

```bash
mkdir -p results/environment_records
conda env export > results/environment_records/conda-resolved.yml
conda list --explicit > results/environment_records/conda-explicit.txt
python -m pip freeze > results/environment_records/pip-freeze.txt
git rev-parse HEAD > results/environment_records/analysis-commit.txt
git -C ../sinabro rev-parse HEAD > results/environment_records/sinabro-commit.txt
```

The YAML pins scientific-library versions rather than every binary build/transitive dependency. Preserve these exports, manifests and input hashes with each campaign. Keep records from distinct environments in separate directories.
