# Single structure refinement: example of `1wou`

1. Download a PDB file (coordinates) and CIF file (structure factors) from the RCSB database.

2. Convert CIF to MTZ using Phenix and expand it to P1 space group:
```bash
phenix.cif_as_mtz 1WOU-sf.cif --ignore_bad_sigmas --merge --output_file_name=1wou-sf.mtz
phenix.reflection_file_converter --expand_to_p1 1wou-sf.mtz --write_mtz_amplitudes --mtz_root_label="FOBS" --label="FOBS" --generate_r_free_flags --non_anomalous --mtz 1wou.mtz
```

3. Create topology (and coordinate) file from PDB:

   * Put the PDB and P1-MTZ files into `data/input/1wou` directory of the cloned repo.

   * `source` AMBER, since we rely on its Python libraries.

   * Set up Python virtual environment and install dependencies (to avoid potential conflicts, use `amber.python`):
   ```bash
   amber.python -m venv venv
   source venv/bin/activate
   pip install -U pip wheel setuptools
   pip install -r requirements.txt
   ```

   * Add the repo's directory into `PYTHONPATH` variable
   ```bash
   export PYTHONPATH=$PYTHONPATH:$(pwd)
   ```

   * Run the script to prepare a model
   ```bash
   python tools/prepare-structures.py
   ```
   The results will be written to `data/amber-topology/1wou`.

4. Prepare the refinement job:
```bash
python init.py
```
First, as a preparatory step for the x-ray restrained MD run, one needs to convert MTZ into a simple text format (described in the AMBER manual). File conversion is performed by `write_sf_dat_file()` method of `init.py`. Second, one needs to create x-ray-specific topology file. The expansion of the topology file (adding crystallographic parameters) is performed by `prepare_xray_prmtop()` method  of `init.py`. The files are written to `data/output/1wou`. Note that the script automatically reads the number of residues from the topology file and uses it during the refinement (`heating`, `evolution` and `cooling` stages).

If you wish to employ a 2x2x2 supercell model instead of the unit cell model, please, open the file `init.py` and add `.sc` postfix to the name of the folder in the `main()` function (i.e. replace `prepared` with `prepared.sc`).

5. Finally, run the refinement job locally by executing:
```bash
python run_locally.py
```
Alternatively, if you want to use a different machine, you should employ `LocalSlurmWorker` remote runner and execute `python run_remotely.py` instead. NB: Don't forget to adjust paths that are `source`d and `export`ed (in essence, environmental variables should be defined similar to step 3). The results will be written to `data/output/1wou`.


# Water picking, B-factors refinement and MolProbity reports generation

Historically, these tasks were executed separately from the Amber-based pipeline. The related files can be found in the `tools/phenix_tasks` directory. Please, adjust the paths in the `run_phenix.sh` accordingly before running it. The first argument of the script should be the path to the PDB file, the second one - the path to the MTZ file.


# Batch refinement (summary of commands)

```bash
# Copy the structures (PDB and MTZ files) from a given directory to `data/input/`
python tools/copy-search-results.py

# Create Amber coordinate/topology (rst7, parm7 and corresponding pdb) files
python tools/prepare-structures.py

# Generate MD protocols
# Note: change the list of structure you wish to iterate through
python init.py

# Run MD protocols
python run_locally.py
```

# Development of the pipeline

### Set up environment

```bash
python3.9 -m venv venv
source venv/bin/activate
pip install -U pip wheel setuptools
```

### Install dependencies

```bash
pip install pip-tools
pip-sync requirements-dev.txt
pre-commit install
```

### Run linters

```bash
# Pre-commit hooks
pre-commit run --all-files
```

### Update requirements

```bash
./dev-tools/rebuild-requirements.sh
```
