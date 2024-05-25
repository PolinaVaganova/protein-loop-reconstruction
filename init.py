from __future__ import annotations
import os
import shutil
from pathlib import Path
from typing import List, Union
import pandas as pd
import gemmi
from amber_runner.MD import MdProtocol, PmemdCommand, SingleSanderCall, Step
from remote_runner import Task


class Prepare(Step):
    @classmethod
    def write_sf_dat_file(
            cls, mtz: "gemmi.Mtz", output_filename: Union[os.PathLike, str], cell_size: int
    ):
        """Write .tab file for pmemd.arx
        :param mtz: mtz file with P1 symmetry
        :param output_filename: output .dat filename
        :param cell_size: size of supercell
        """
        import re

        assert mtz.spacegroup.hm == "P 1"

        def missing_column(label):
            raise RuntimeError(f"MTZ file missing {label} column")

        R_FREE_FLAG = None
        r_free_pattern_string = r"(r.*free)|(free.*r)"
        r_free_pattern = re.compile(r_free_pattern_string, flags=re.IGNORECASE)

        for column in mtz.columns:  # type: gemmi.Mtz.Column
            if r_free_pattern.search(column.label):
                R_FREE_FLAG = column

        if R_FREE_FLAG is None:
            raise RuntimeError(
                f"MTZ file missing R-FREE-FLAG column "
                f"(pattern: `{r_free_pattern_string}`)"
                f"\nPresent columns: {[column.label for column in mtz.columns]}"
            )

        H, K, L, FOBS, SIGMA_FOBS = [
            mtz.column_with_label(label)  # or missing_column(label)
            for label in ("H", "K", "L", "FOBS", "SIGFOBS")
        ]
        # if SIGMA_FOBS was not found, fill it with ones
        # currently it's not used anyway
        if not SIGMA_FOBS:
            import numpy as np

            SIGMA_FOBS = np.ones(len(FOBS))

        n_positive_r_flags = sum(R_FREE_FLAG)
        flag_is_one = n_positive_r_flags > len(R_FREE_FLAG) / 2

        with open(output_filename, "w") as output:
            output.write(f"{len(H)} 0\n")

            for h, k, l, fobs, sigma, r_flag in zip(
                    H, K, L, FOBS, SIGMA_FOBS, R_FREE_FLAG
            ):
                r = r_flag if flag_is_one else 1 - r_flag
                output.write(
                    f"{h * cell_size:3.0f} {k * cell_size:3.0f} {l * cell_size:3.0f} "
                    f"{fobs:15.8e} {sigma:15.8e} {r:1.0f}\n"
                )

    def __init__(
            self,
            name: str,
            parm7: Path,
            rst7: Path,
            pdb: Path,
            mtz: Path,
            xray_weight_target: float,
            cell_size: int,
    ):
        Step.__init__(self, name)

        self.charge = None
        self.input_pdb_path: Path = pdb
        self.input_mtz_path: Path = mtz
        self.wbox_prmtop = parm7
        self.wbox_inpcrd_path = rst7
        self.wbox_pdb = pdb
        self.xray_weight_target = xray_weight_target
        self.cell_size = cell_size

    @property
    def wbox_xray_prmtop_path(self):
        return self.step_dir / "wbox.xray.parm7"

    @property
    def structure_factors_dat(self):
        return str(self.step_dir / "sf.dat")

    def prepare_structure_factors(self):
        import gemmi

        mtz = gemmi.read_mtz_file(str(self.input_mtz_path))
        self.write_sf_dat_file(
            mtz=mtz,
            output_filename=self.structure_factors_dat,
            cell_size=self.cell_size,
        )

    def prepare_xray_prmtop(self):
        from arx.utils import check_call

        tmp_parm = self.step_dir / "tmp.parm7"
        parmed_add_xray_parameters_in = self.step_dir / "parmed.add-xray-parameters.in"
        parmed_in = f"""
addPdb {self.wbox_pdb} elem strict allicodes
lmod
parmout {tmp_parm}
go
"""

        with open(parmed_add_xray_parameters_in, "w") as f:
            f.write(parmed_in)
        check_call(
            [
                "parmed",
                "--overwrite",
                str(self.wbox_prmtop),
                str(parmed_add_xray_parameters_in),
            ]
        )

        check_call(
            [
                "add_xray",
                "-i",
                str(tmp_parm),
                "-o",
                str(self.wbox_xray_prmtop_path),
                "-scattering",
                "xray",
            ]
        )

    def prepare_files_for_next_stages(self, md: RefinementProtocol):
        import gemmi
        from amber_runner.inputs import AmberInput

        def count_polymer_residues(
                st: gemmi.Structure,
        ) -> int:
            non_polymer_residue_names = ["WAT", "Cl-", "Na+"]
            count = 0
            for model in st:
                for chain in model:
                    for residue in chain:  # type: gemmi.Residue
                        if residue.name in non_polymer_residue_names:
                            return count
                        else:
                            count += 1
            return count

        # Set global attributes
        md.sander.prmtop = str(self.wbox_xray_prmtop_path)
        md.sander.inpcrd = str(self.wbox_inpcrd_path)
        md.sander.refc = str(self.wbox_inpcrd_path)

        n_polymer_residues = count_polymer_residues(gemmi.read_pdb(str(self.wbox_pdb)))

        # Configure minimize
        md.minimize.input.cntrl(
            imin=1,
            maxcyc=500,
            ncyc=150,
            ntb=1,
            ntr=0,
            cut=8.0,
        )
        # Configure heating
        md.heating.input.cntrl(
            cut=8.0,
            dt=0.002,
            imin=0,
            ioutfm=1,
            irest=0,
            nstlim=10000,
            ntb=1,
            ntc=2,
            ntf=2,
            ntpr=50,
            ntr=1,
            ntt=1,
            ntwr=1000,
            ntwx=200,
            ntx=1,
            temp0=298.0,
            tempi=0.0,
        )
        md.heating.input.pin(
            AmberInput.GroupSelection(
                title="Keep protein fixed with weak restraints",
                weight=10.0,
                residue_id_ranges=[(1, n_polymer_residues)],
            )
        )

        # Configure evolution
        nstlim = 5000
        md.evolution.input.cntrl(
            imin=0,
            irest=1,
            ntx=5,
            iwrap=1,
            ntb=1,
            ntt=3,
            gamma_ln=3.0,
            ig=-1,
            tempi=298.0,
            temp0=298.0,
            ntp=0,
            pres0=1.0,
            taup=2.0,
            cut=8.0,
            ntr=0,
            ntc=2,
            ntf=2,
            nstlim=nstlim,
            nscm=100,
            dt=0.002,
            ntpr=100,
            ntwx=100,
            ntwr=5000,
            ioutfm=1,
        )

        md.evolution.input._get("xray")(
            spacegroup_name="P1",
            pdb_infile=str(self.wbox_pdb),
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            atom_selection_mask=f":1-{n_polymer_residues}",
            xray_weight_initial=0.0,
            xray_weight_final=self.xray_weight_target,
            target="ml",
            bulk_solvent_model="afonine-2013",
        )

        # Configure cool
        md.cooling.input.cntrl(
            imin=0,
            ntx=5,
            irest=1,
            iwrap=1,
            nstlim=nstlim,
            dt=0.002,
            ntf=2,
            ntc=2,
            tempi=298.0,
            temp0=0.0,
            ntpr=100,
            ntwx=100,
            cut=8.0,
            ntb=1,
            ntp=0,
            ntt=3,
            gamma_ln=2.0,
            nscm=200,
            nmropt=1,
        )

        md.cooling.input._get("xray")(
            spacegroup_name="P1",
            pdb_infile=str(self.wbox_pdb),
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            atom_selection_mask=f":1-{n_polymer_residues}",
            xray_weight_initial=self.xray_weight_target,
            xray_weight_final=self.xray_weight_target,
            target="ml",
            bulk_solvent_model="afonine-2013",
        )

        cool_steps = []
        start = 0
        steps_inc = 500
        steps_inc_steady = 125
        temp = 300.0
        steps = nstlim // (steps_inc + steps_inc_steady)
        temp_inc = -temp / steps

        for i in range(steps):
            if i == 0:
                cool_steps.append((start, start + steps_inc, 298.0, temp + temp_inc))
                cool_steps.append(
                    (
                        start + steps_inc + 1,
                        start + steps_inc + steps_inc_steady,
                        temp + temp_inc,
                        temp + temp_inc,
                    )
                )
            else:
                cool_steps.append((start + 1, start + steps_inc, temp, temp + temp_inc))
                cool_steps.append(
                    (
                        start + steps_inc + 1,
                        start + steps_inc + steps_inc_steady,
                        temp + temp_inc,
                        temp + temp_inc,
                    )
                )
            start += steps_inc + steps_inc_steady
            temp += temp_inc
        for istep1, istep2, value1, value2 in cool_steps:
            md.cooling.input.varying_conditions.add(
                type="TEMP0", istep1=istep1, istep2=istep2, value1=value1, value2=value2
            )

    def run(self, md: RefinementProtocol):
        self.prepare_structure_factors()
        self.prepare_xray_prmtop()
        self.prepare_files_for_next_stages(md)


class ConvertToPdb(Step):
    def run(self, md: "RefinementProtocol"):
        import tempfile

        from arx.prepare import (
            copy_coordinates,
            read_pdb,
            remove_ligands_and_water,
            write_pdb,
        )
        from arx.utils import check_call

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_pdb = Path(tmp_dir) / "tmp.pdb"
            # Don't change directory to avoid problems
            # with relative paths in md.sander.*
            with open(tmp_pdb, "wb") as f:
                check_call(
                    [
                        "ambpdb",
                        # ambpdb doesn't work well with xray-modified prmtop
                        # use original parm7
                        "-p",
                        str(md.prepare.wbox_prmtop),
                        "-c",
                        str(md.sander.inpcrd),
                    ],
                    stdout=f,
                )
            final = read_pdb(tmp_pdb)

        initial = read_pdb(md.prepare.wbox_pdb)
        result = copy_coordinates(initial, reference=final)
        write_pdb(result, self.step_dir / "final.pdb")
        dry_result = remove_ligands_and_water(result)
        write_pdb(dry_result, self.step_dir / "final_no_water_no_ligands.pdb")


class RefinementProtocol(MdProtocol):
    def __init__(
            self,
            pdb: Path,
            mtz: Path,
            parm7: Path,
            rst7: Path,
            output_dir: Path,
            xray_weight_target: float,
            cell_size: int,
    ):
        wd = output_dir
        wd.mkdir(mode=0o755, exist_ok=True, parents=True)
        MdProtocol.__init__(self, name=str(wd).split("/")[-2] + "_amb", wd=wd)

        self.sander = PmemdCommand()
        self.sander.executable = ["pmemd.cuda"]
        self.sander.allow_small_box = True

        self.prepare = Prepare(
            name="prepare",
            pdb=pdb,
            mtz=mtz,
            parm7=parm7,
            rst7=rst7,
            xray_weight_target=xray_weight_target,
            cell_size=cell_size,
        )

        self.minimize = SingleSanderCall("minimize")
        self.heating = SingleSanderCall("heating")
        self.evolution = SingleSanderCall("evolution")
        self.cooling = SingleSanderCall("cooling")
        self.convert_to_pdb = ConvertToPdb("convert_to_pdb")


def create_tasks(topology_dirs: List[Path],
                 small_molecules_dir: Path,
                 total_output_dir: Path,
                 path_to_mtz_file: Path
                 ):
    tasks = []
    # collect structures with missing files
    bad_tasks_parm = []
    bad_tasks_pdb = []
    bad_tasks_mtz = []

    cell_size = 1
    xray_weight = 1.0

    for topology_dir in topology_dirs:
        output_dir = total_output_dir / os.path.basename(topology_dir)
        input_copy = output_dir / "inputs"
        pdb = input_copy / "wbox.pdb"
        mtz = input_copy / os.path.basename(path_to_mtz_file)
        rst7 = input_copy / "wbox.rst7"
        parm7 = input_copy / "wbox.parm7"

        # skip execution if ferined structure already exist
        if (
                not (topology_dir / "wbox.rst7").exists()
                or not (topology_dir / "wbox.parm7").exists()
        ):
            bad_tasks_parm.append(os.path.basename(topology_dir))
        elif not (topology_dir / "wbox.pdb").exists():
            bad_tasks_pdb.append(os.path.basename(topology_dir))
        elif not (path_to_mtz_file).exists():
            bad_tasks_mtz.append(os.path.basename(topology_dir))
        else:
            input_copy.mkdir(exist_ok=True, parents=True)
            shutil.copy(topology_dir / "wbox.pdb", pdb)
            shutil.copy(path_to_mtz_file, mtz)
            shutil.copy(topology_dir / "wbox.rst7", rst7)
            shutil.copy(topology_dir / "wbox.parm7", parm7)

            md = RefinementProtocol(
                pdb=pdb.relative_to(output_dir),
                mtz=mtz.relative_to(output_dir),
                rst7=rst7.relative_to(output_dir),
                parm7=parm7.relative_to(output_dir),
                output_dir=output_dir,
                xray_weight_target=xray_weight,
                cell_size=cell_size,
            )
            md.save(md.wd / "state.dill")
            tasks.append(md)
    print(f"Structures with missing parm file: {', '.join(bad_tasks_parm)}")
    print(f"Structures with missing pdb file: {', '.join(bad_tasks_pdb)}")
    print(f"Structures with missing mtz file: {', '.join(bad_tasks_mtz)}")
    print(f"Prepared {len(tasks)} jobs")

    return tasks


def run_all_locally(tasks: List[Task]):
    import remote_runner
    from remote_runner import LocalWorker, Pool

    remote_runner.log_to(".remote-runner.log", level="DEBUG")
    workers = [LocalWorker()]
    Pool(workers).run(tasks)


def run_sequentially_inplace(tasks: List[Task]):
    from remote_runner.utility import ChangeDirectory

    for md in tasks:
        with ChangeDirectory(md.wd):
            md.run()


def main():
    # settings
    path_to_annotation_csv = "/home/polina/xray-refinement/1_annotated_rcsb/1u3y_1k34.csv"
    run_dir = "4_protocol_run"
    small_molecules_dir = os.path.join(run_dir, "small-molecules")
    pdb_df = pd.read_csv(path_to_annotation_csv, sep=",")
    pdb_ids = pdb_df["pdb_id"]
    postfix = "_rapper_saved_split"

    # iterate over available pdb codes and initialize protocol
    for pdb_id in pdb_ids:
        path_to_main_dir = os.path.join(run_dir, "amber-topology", f"{pdb_id}{postfix}")

        topology_dirs = [Path(os.path.join(path_to_main_dir, dirname))
                         for dirname in os.listdir(path_to_main_dir)]

        output_dir = os.path.join(run_dir, "output", f"{pdb_id}{postfix}")
        path_to_mtz_file = os.path.join(run_dir, "input", f"{pdb_id}{postfix}", f"{pdb_id}.mtz")

        tasks = create_tasks(topology_dirs=topology_dirs,
                             small_molecules_dir=Path(small_molecules_dir),
                             total_output_dir=Path(output_dir),
                             path_to_mtz_file=Path(path_to_mtz_file))

        assert tasks


if __name__ == "__main__":
    main()
