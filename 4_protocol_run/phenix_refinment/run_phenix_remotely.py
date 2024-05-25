import os
from cmd_runner.runner import Pool, LocalSlurmWorker, Cmd


def main():
    # Define tasks
    structures_names = ["1K33_1|Chain.B99990034/"
                        ]

    path_to_structures = (
        "{your_path_to_repo}/4_protocol_run/output/1k33/"
    )

    path_to_data = "{your_path_to_repo}/2_rcsb_data"
    cmds = []

    for structures_name in structures_names:

        os.makedirs(structures_name, exist_ok=True)

        pdb_id = structures_name[:4].lower()

        path_to_amber_final_pdb = os.path.join(
            path_to_structures,
            structures_name,
            "5_convert_to_pdb",
            "final_no_water_no_ligands.pdb",
        )

        path_to_mtz = os.path.join(path_to_data, "mtz", f"{pdb_id}.mtz")

        if not os.path.exists(os.path.join(structures_name, "phenix.sh")):
            os.symlink(
                f"{os.getcwd()}/phenix.sh",
                os.path.join(structures_name, "phenix.sh"),
            )

        cmds.append(
            Cmd(
                cmd=f"python ../residue_renamer.py {path_to_amber_final_pdb} input.pdb",
                workdir=f"{structures_name}",
            )
        )

        cmds.append(
            Cmd(
                cmd=f"bash phenix.sh input.pdb {path_to_mtz}",
                workdir=f"{structures_name}",
            )
        )

    # Configure workers
    workers = [
        LocalSlurmWorker(
            remote_user_rc="""
source /home/olebedenko/xray-refinment/venv_amber22_with_patches/bin/activate
""",
            sbatch_args=[
                f"--job-name={structures_name}_phenix.refine",
                "--ntasks=1",
                f"--output={structures_name}.out",
                f"--error={structures_name}.err"
            ],
        )
        for _ in range(len(structures_names))
    ]

    Pool(workers).run(cmds)


if __name__ == "__main__":
    main()
