from pathlib import Path
import remote_runner
from remote_runner import Pool, LocalSlurmWorker, Task
import pandas as pd

def main():
    remote_runner.log_to(".remote-runner.log")
    
    # path_to_csv = "/home/polina/xray-refinement/1_annotated_rcsb/monomers_to_ucell.csv"
    # pdb_df = pd.read_csv(path_to_csv, sep=",")
    # pdb_ids = pdb_df["pdb_id"]
    pdb_ids = ["6vo7"]
    paths_to_run_dirs = []
    postfix = "modeller_saved_split"
    # pdb_ids = ["1u3y_rapper_saved_split"]

    for pdb_id in pdb_ids:
        paths_to_run_dirs += sorted(Path("./").glob(f"4_protocol_run/output/{pdb_id}_{postfix}/6VO7_1|Chain.B99990095/state.dill"))

    workers = [LocalSlurmWorker(
        remote_user_rc="""
export PYTHONPATH=$PYTHONPATH:$(pwd)
source /opt/amber22_with_patches/amber.sh
source /home/olebedenko/xray-refinment/venv_amber22_with_patches/bin/activate
""",
        sbatch_args=[
            "--cpus-per-task=2",  # cpu here means threads
            "--partition=full",
            "--gpus=1",
            # "--nodelis=bionmr-mom-001,bionmr-mom-002,bionmr-mom-003,bionmr-mom-004,bionmr-mom-005",
        ]
    ) for _ in range(len(pdb_ids))]

    tasks = [
        Task.load(path)
        for path in paths_to_run_dirs
    ]
    print(len(tasks))
    Pool(workers).run(tasks)


if __name__ == "__main__":
    main()
