# Protein loop reconstruction by using MD simulations and crystallographic data

_This repository provides results and code for the project "Protein loop reconstruction by using MD simulations and
crystallographic data"_

**Author:**

* Polina Vaganova _(Saint Petersburg State University, Saint Petersburg, Russia)_

**Supervisors:**

* Olga Lebedenko _(BioNMR laboratory, Saint Petersburg State University, Saint Petersburg, Russia)_
* Nikolai Skrynnikov _(BioNMR laboratory, Saint Petersburg State University, Saint Petersburg, Russia)_

**Table of contents:**

- [Background](#sec1) </br>
- [Methods](#sec2) </br>
- [System requirements](#sec3) </br>
- [Repository structure](#sec4) </br>
- [Single structure loop building: example of `1k33`](#sec5) </br>
- [References](#sec6) </br>

<a name="sec1"></a>

### Background:

To date, X-ray crystallography remains a primary experimental method for elucidation of three-dimensional (3D) protein
structure. The diffraction process is fundamentally different from microscopic imaging; crystallography is not an
imaging technique. The crucial difference is that visible light scattered from objects can be focused through refractive
lenses to create a magnified image of the object. This is not the case for X-rays, which are also electromagnetic
radiation but of several orders of magnitude shorter wavelength and correspondingly higher energy; the
refractive index of X-rays in different materials si essentially equal and close to unity, and no refractive lenses can
be constructed for X-rays. Instead, the electron density of the scattering molecular structure must be reconstructed by
special technique (Fourier transform). In itself, the Fourier reconstruction from reciprocal diffraction space back
into direct molecular space poses no difficulties in principle, with the unfortunate qualifier that for this type of
reconstruction two terms are needed as Fourier coefficients: the structure factor amplitudes, readily accessible in the
form of the square root of the measured ($F_{obs}$) and corrected ($F_{calc}$) diffraction spot intensities (Figure 1);
and as a second term for each observed diffraction spot, its relative phase angle.

|    ![xray-experiment.png](/figures/xray-experiment.png)     |
|:-----------------------------------------------------------:|
| *Figure 1. The principle of X-ray structure determination.* |

These phase angles are not directly accessible and must be adjusted during special procedure, called refinement. In
refinement procedure we are adjusting phase angles to obtain a best fit between the observed structure factor
amplitudes (($F_{obs}$)) and the computed model structure factor amplitudes ($F_{calc}$). The overall fit between
diffraction data and model is numerically quantified by a global linear residual (the R-value) between the scaled
structure factor amplitudes $F_{obs}$ and $F_{calc}$:
$R = \frac{{\sum\limits_{hkl} |F_{\text{obs}} - F_{\text{calc}}|}}{{\sum\limits_{hkl} |F_{\text{obs}}|}}$, where $hkl$
Miller indexes for each of reflection spots.

The X-ray crystallography technique operates with static images of molecular structures. The 3D structure is formed by
secondary struture elements (e.g. α-helices and β-strands) as well as loosely
structured elements such as loops. Protein loops are often functionally important. Their conformational plasticity is
key to molecular recognition, allosteric control, ligand binding, and signaling.

However, the structural variability of protein loops presents a challenge for the X-ray crystallography
which is normally limited to static structural models. As a consequence, mobile loops are
often absent from the crystallographic structures deposited in the Protein Data Bank (PDB).
Nevertheless, residual electron density associated with such loops can be used to rebuild them
in a form of conformational ensemble. To this end, we have used Molecular Dynamics (MD)
simulations additionally guided by the experimental diffraction data. The appropriate
crystallographic refinement procedure has been developed in our laboratory using the
biomolecular MD simulation platform Amber [[1]](#arx).

**Goal:** The main goal of this project is to reconstruct the protein loop ensemble using Molecular Dynamics (MD)
simulation guided by the experimental diffraction data.

<a name="sec2"></a>

### Methods:

Protein structures with missing loop regions and different levels of residual electron density at the disordered loop
sites have been identified by the automated parsing of [RCSB database](https://www.rcsb.org/). An example of one such
structure (PDB: `1k33`) is illustrated in Figure 2.

$R_{work}$ / $R_{free}$ for this deposited structure are 0.206 / 0.226.

Initial loop conformations were generated using programs such as Modeller or Rapper (Figure 3). This step doesn't
operate with
electron density data and thus the resulting structure has worse scores (higher is worse).

$R_{work}$ / $R_{free}$ for Modeller resulting structure are 0.316 / 0.317.

Subsequently, MD-based crystallographic refinement procedure has been performed on
structural models in a form of crystal unit cells or supercells (Figure 4).

The refined models involving multiple loop conformations have been validated against the available electron density
maps using [phenix.molprobity](https://phenix-online.org/documentation/reference/molprobity_tool.html) module.

| ![1k33_without_loop.png](/figures/1k33_without_loop.png) |     ![1k33_with_loop.png](/figures/1k33_loop_modeller.png)      | ![1k33_ucell.png](/figures/1k33_ucell.png) |   ![1k33_loop_ensemble.png](/figures/1k33_loop_ensemble.png)   |
|:--------------------------------------------------------:|:---------------------------------------------------------------:|:------------------------------------------:|:--------------------------------------------------------------:|
|  *Figure 2. Structure of `1k33` without loop residues.*  | *Figure 3. Structure of `1k33` with built loop residues (red).* | *Figure 4. Structure of `1k33` unit cell.* | *Figure 5. Loop ensemble of `1k33` with difference Fo-Fc map.* |

After refinement procedure we observed the dramatic improvement in $R_{work}$ / $R_{free}$ from 0.206 / 0.226 (RCSB) to
0.181 / 0.202 in our new approach.

<a name="sec3"></a>

### System requirements:

**Key packages and programs:**

- [Python](https://www.python.org/downloads/) (>= 3.9)
- [Amber22](https://ambermd.org/index.php) (a build with MPI is employed)
- [Modeller](https://salilab.org/modeller/download_installation.html) (10.4) package for loop modelling
- [phenix](https://phenix-online.org) software package for macromolecular structure determination
- [pyxmolpp2](https://github.com/sizmailov/pyxmolpp2) (1.6.0) in-house python library for processing molecular
  structures and MD trajectories
- `slurm` (20.11.8) cluster management and job scheduling system
- other python libraries used are listed in `requirements.txt`

<a name="sec4"></a>

### Repository structure:

- [0_prepare_annotation](0_prepare_annotation) contains rcsb annotation results
- [1_annotated_rcsb](1_annotated_rcsb) contains annotation only for monomers with different symmetry operations
- [2_rcsb_data](2_rcsb_data) contains examples for input files
- [3_loop_building](3_loop_building) contains examples for building initial loop conformation
- [4_protocol_run](4_protocol_run) contains examples for running MD-based refinement protocol
- [arx](arx) contains special module for AMBER
- [utils](utils) contains python scripts for annotations, filtering, loop building and input preparation

<a name="sec5"></a>

### Single structure loop building: example of `1k33`

1. Download a PDB file (coordinates), CIF file (structure factors) and FASTA file (aminoacid sequence) from
   the [RCSB database](https://www.rcsb.org/). Place files into corresponding folders in `2_rcsb_data` dir.

2. Go to `2_rcsb_data/mtz` dir and convert CIF to MTZ using Phenix. Then expand it to P1 space group:

```bash
cd 2_rcsb_data/mtz
phenix.cif_as_mtz ../cif/1k33-sf.cif --ignore_bad_sigmas --merge --output_file_name=1k33-sf.mtz
phenix.reflection_file_converter --expand_to_p1 1k33-sf.mtz --write_mtz_amplitudes --mtz_root_label="FOBS" --label="FOBS" --generate_r_free_flags --non_anomalous --mtz 1k33.mtz
```

3. Build initial loop using Modeller. Please, adjust the `pdb_ids` variable in the `run_phenix.sh` accordingly before
   running it.

```bash
cd ../..
python utils/modeller/model.py
```

4. Create topology (and coordinate) file from PDB:

    * Copy pdn and P1-MTZ files into input dir:
    ```bash
   python utils/prepare_inputs/1_copy_inputs.py
    ```

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
   python utils/prepare_inputs/2_prepare-structures.py
   ```
   The results will be written to `4_protocol_run/amber-topology/1k33/`.

4. Prepare the refinement job:

```bash
python init.py
```

First, as a preparatory step for the x-ray restrained MD run, one needs to convert MTZ into a simple text format (
described in the AMBER manual). File conversion is performed by `write_sf_dat_file()` method of `init.py`. Second, one
needs to create x-ray-specific topology file. The expansion of the topology file (adding crystallographic parameters) is
performed by `prepare_xray_prmtop()` method of `init.py`. The files are written to `4_protocol_run/output/1k33/`. Note
that the
script automatically reads the number of residues from the topology file and uses it during the
refinement (`heating`, `evolution` and `cooling` stages).

If you wish to employ a 2x2x2 supercell model instead of the unit cell model, please, open the file `init.py` and
add `.sc` postfix to the name of the folder in the `main()` function (i.e. replace `prepared` with `prepared.sc`).

5. Finally, run the refinement job locally by executing:

```bash
python run_locally.py
```

Alternatively, if you want to use a different machine, you should employ `LocalSlurmWorker` remote runner and
execute `python run_remotely.py` instead. NB: Don't forget to adjust paths that are `source`d and `export`ed (in
essence, environmental variables should be defined similar to step 3). The results will be written
to `4_protocol_run/output/1k33`.

6. Water picking, B-factors refinement and MolProbity reports generation

Historically, these tasks were executed separately from the AMBER-based pipeline. The related files can be found in
the `4_protocol_run/phenix_refinment` directory. Please, adjust the paths in the `run_phenix_remotely.py` accordingly
before running it (change `{your_path_to_repo}` with your actually path).

If you want to run locally without slurm then launch `residue_renamer.py` on your pdb, adjust the paths in
the `phenix.sh` and launch it. The first argument of the script should be the path to the PDB file, the second one - the
path to the MTZ file.

<a name="sec6"></a>

### References

<a name="arx">[1]</a>
Mikhailovskii, O., Xue, Y., Skrynnikov, N. R. Modeling a Unit Cell: Crystallographic Refinement Procedure Using the
Biomolecular MD Simulation Platform Amber. 2022. IUCrJ, 9 (1): 114–133. https://doi.org/10.1107/S2052252521011891.
