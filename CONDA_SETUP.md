# Conda setup for the manuscript revision

Use this environment for GENCODE preprocessing, the timing pilots and their evaluation, and the SeaWulf Slurm simulation runner. It installs Sinabro 3.0.0 from its separate Git repository. Later historical analysis and R scripts still require separate review and dependencies.

## 1. Make Conda available on the HPC login node

If you already use your own working Conda installation, use it. Otherwise, the [SeaWulf package-management guide](https://rci.stonybrook.edu/docs/python-package-management) documents:

```bash
module load anaconda/3
```

Check available module names with `module avail anaconda` if that name is unavailable; use the name actually provided by the cluster. Then check:

```bash
conda --version
```

If activation later asks you to initialize Conda, initialize only the current shell (Bash or Zsh):

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
```

## 2. Get the two repositories

For a fresh setup, clone them into the same parent directory on storage accessible to both login and compute nodes:

```bash
mkdir -p "$HOME/sinabro_revision"
cd "$HOME/sinabro_revision"
git clone --branch work/package-refactor-2026-10-14 https://github.com/StudyingAnt/sinabro.git
git clone --branch work/reproduce-manuscript-2026-10-14 https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds.git
cd Simple_robustness_analysis_of_cds
```

If you cloned these branches already, enter each checkout and use `git pull --ff-only` instead of cloning again. Keep local modifications before resolving any pull error. Run the remaining commands from the **analysis repository root**, with the package checkout at `../sinabro`.

## 3. Create the environment and install Sinabro

```bash
conda env create --file environment.yml
conda activate sinabro-revision
python -m pip install --no-deps --no-build-isolation --force-reinstall ../sinabro
```

Create this environment once. On later logins, load the same Conda module if applicable and run `conda activate sinabro-revision`. If the name already exists, check it before proceeding; choose a new name with `conda env create --name YOUR_NEW_NAME --file environment.yml` if needed, then activate that name.

The [Conda environment guide](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html) supports creating environments from YAML. The `conda env create` form is retained here for compatibility with older cluster Conda versions.

| Component | Pinned version | Purpose |
|---|---|---|
| Python | 3.9.15 | Recorded historical environment and checked revision pilots |
| NumPy | 1.21.6 | Numerical calculations and random generators |
| Biopython | 1.79 | Sequence handling and translation |
| SciPy | 1.8.1 | Optional Sinabro analysis functions |
| Matplotlib | 3.6.2 | Pilot evaluation plots |
| Sinabro | 3.0.0 | Installed from the separate package checkout |

The recorded full historical environment has Python 3.9.15; an older preprocessing README stated 3.9.16. The checked revision workflow uses 3.9.15. The YAML pins these scientific-library versions, not every transitive dependency or binary build; export the resolved HPC environment as shown below.

Conda installs the dependencies first. Pip installs only Sinabro from the local checkout: `--no-deps` prevents dependency replacement, and `--no-build-isolation` uses the environment's installed build tools. `--force-reinstall` replaces any previously installed Sinabro, including an editable installation or another checkout with the same version number. This regular installation copies the package into the environment's `site-packages`; source edits require running the installation command again.

The default PEP 660 editable installation (`pip install -e ...`) can import the outer repository directory as a namespace package when Python runs from the repositories' shared parent directory. The symptom is `sinabro.__file__ is None` and a missing `sinabro.__version__`, even though installation succeeded. Regular installation avoids this reproduced conflict; see the [Setuptools editable-install limitations](https://setuptools.pypa.io/en/latest/userguide/development_mode.html#limitations). No `conda develop` or `import sinabro.sinabro` is needed with this version.

Record the package commit before each run and keep the checkout fixed while jobs run. The pilot and Slurm runners explicitly load the checkout selected by `--sinabro-source` or `SINABRO_SOURCE`, so keep that checkout consistent with the installed package and reinstall after changing it.

### Optional: install the environment at a specific absolute path

If you need to store the environment outside Conda's default environment directory, use `--prefix` instead of the default creation command above. Replace `/PATH/TO/ENVS` with your actual absolute directory on storage accessible to both login and compute nodes. Run these commands from the analysis repository root:

```bash
conda env create --prefix /PATH/TO/ENVS/sinabro-revision --file environment.yml
conda activate /PATH/TO/ENVS/sinabro-revision
python -m pip install --no-deps --no-build-isolation --force-reinstall ../sinabro
```

Choose either the default named environment or this path-based installation; you do not need to create both. The explicit `--prefix` takes precedence over `name: sinabro-revision` in `environment.yml`. On later logins, activate this environment by its absolute path:

```bash
conda activate /PATH/TO/ENVS/sinabro-revision
```

The checks and exports in the following sections apply to either installation method. Before job submission, use the same command shown in section 6:

```bash
export SINABRO_PYTHON="$CONDA_PREFIX/bin/python"
```

For this path-based installation, it resolves to `/PATH/TO/ENVS/sinabro-revision/bin/python` after replacing the placeholder. Export `SINABRO_ENV="$CONDA_PREFIX"` before submitting the Slurm wrapper so it activates the same custom environment inside the batch job.

## 4. Check the environment without running simulations

```bash
python - <<'PY'
import platform
import numpy
import Bio
import scipy
import matplotlib
import sinabro
from sinabro import Trajectory, RobustnessComputer
from sinabro.core import signature_path

expected = {
    'Python': ('3.9.15', platform.python_version()),
    'NumPy': ('1.21.6', numpy.__version__),
    'Biopython': ('1.79', Bio.__version__),
    'SciPy': ('1.8.1', scipy.__version__),
    'Matplotlib': ('3.6.2', matplotlib.__version__),
    'Sinabro': ('3.0.0', sinabro.__version__),
}
for name, (wanted, actual) in expected.items():
    print(f'{name}: {actual}')
    assert actual == wanted, f'{name}: expected {wanted}, got {actual}'
assert signature_path('3.3.1', 'GRCh37').is_file()
print('Sinabro source:', sinabro.__file__)
print('Environment check passed; no trajectories simulated.')
PY

python -m unittest discover -s tests -v
```

The last command runs the preprocessing and HPC-runner checkpoint/argument tests. These commands do not execute the full simulation or submit jobs. With the regular installation above, `Sinabro source` should point into the activated environment's `lib/python3.9/site-packages/sinabro/__init__.py`. It should not be `None`. A successful local check does not establish that the same environment has been created on HPC; run it there too.

Also verify the import from the repositories' shared parent directory. The subshell returns you to the analysis directory afterward:

```bash
(
  cd ..
  python -c "import sinabro; print(sinabro.__file__); assert sinabro.__file__ is not None; assert sinabro.__version__ == '3.0.0'; from sinabro import Trajectory, RobustnessComputer; print('Parent-directory import: OK')"
)
```

### Repair an existing editable installation

Keep the current Conda environment. From the analysis repository root, replace the editable installation and repeat the checks above:

```bash
python -m pip install --no-deps --no-build-isolation --force-reinstall ../sinabro
```

If you are instead in the shared parent directory (for example, `sinabro_pub`), the package path is `./sinabro`:

```bash
python -m pip install --no-deps --no-build-isolation --force-reinstall ./sinabro
```

These are the same repair command with different relative paths; run the one matching your current directory. Repeating the old `-e` command does not fix this import-precedence conflict. No environment recreation or scientific-library upgrade is required.

## 5. Record the actual HPC environment

```bash
mkdir -p results/environment_records
conda env export > results/environment_records/conda-resolved.yml
conda list --explicit > results/environment_records/conda-linux-64-explicit.txt
python -m pip freeze > results/environment_records/pip-freeze.txt
git rev-parse HEAD > results/environment_records/analysis-commit.txt
git -C ../sinabro rev-parse HEAD > results/environment_records/sinabro-commit.txt
```

The explicit Conda export records resolved package builds for the same platform. The separate package commit is required to reconstruct the local source used for installation and by the simulation runners. These files are generated provenance under Git-ignored `results/`; retain them with the run. Use a new recording directory when preserving records from multiple environments.

## 6. Prepare the input, then submit manually

Git does not include FASTA inputs or generated results. Transfer the checked CDS-plus FASTA to HPC, or follow the [download and preprocessing instructions](01_data_downloading_and_preprocessing/README.md) using this environment. A Conda installation alone does not run either stage.

Before following [HPC_SUBMISSION.md](02_run_sinabro/HPC_SUBMISSION.md), activate the environment in the same shell used for submission. The following only prepares variables and the log directory; it does **not** submit a job:

```bash
export SINABRO_PYTHON="$CONDA_PREFIX/bin/python"
export SINABRO_SOURCE="$(cd ../sinabro && pwd)"
export SINABRO_ENV="$CONDA_PREFIX"
export SINABRO_FASTA="$PWD/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa"
export SINABRO_OUTPUT="$PWD/results/hpc_gencode_v40_n1000"
mkdir -p logs
```

The paths above illustrate the analysis repository's local layout. On HPC, use the project-level `data/processed/` and `runs/` paths in the submission guide. The Slurm wrapper activates `SINABRO_ENV` inside each batch job and uses `SINABRO_PYTHON` if specified. The environment, source and input paths must all be visible from compute nodes. Submit from the project root only when ready, as described in the submission guide.
