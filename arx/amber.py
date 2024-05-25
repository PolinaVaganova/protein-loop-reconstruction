import tempfile
from pathlib import Path
from typing import Iterator, Tuple

import gemmi

from arx.prepare import read_pdb, write_pdb
from arx.utils import chdir, check_call

from .prepare import add_missing_b_factors


def find_ss_bond_pairs(st: gemmi.Structure) -> Iterator[Tuple[int, int]]:
    cyx_residues = [
        residue for chain in st[0] for residue in chain if residue.name == "CYX"
    ]

    for i in range(len(cyx_residues)):
        ri = cyx_residues[i]
        sgi = ri["SG"][0]
        for j in range(i + 1, len(cyx_residues)):
            rj = cyx_residues[j]
            sgj = rj["SG"][0]
            if sgi.pos.dist(sgj.pos) < 2.3:
                yield ri.seqid.num, rj.seqid.num


def get_ss_bond_commands(st: gemmi.Structure) -> str:
    return "\n".join(f"bond wbox.{i}.SG wbox.{j}.SG" for i, j in find_ss_bond_pairs(st))


def create_topology_and_input(
    st: gemmi.Structure, parm7_path: Path, rst7_path: Path
) -> gemmi.Structure:
    bond_commands = get_ss_bond_commands(st)
    config = f"""
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p

wbox = loadpdb input.pdb

{bond_commands}

set default nocenter on
setBox wbox vdw 1.0
saveamberparm wbox {parm7_path} {rst7_path}
savepdb wbox wbox.pdb
quit
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            with open("tleap.in", "w") as f:
                f.write(config)

            write_pdb(st, "input.pdb")

            check_call(["tleap", "-s", "-f", "tleap.in"])

            check_call(
                [
                    "ChBox",
                    "-c",
                    str(rst7_path),
                    "-o",
                    str(rst7_path),
                    "-X",
                    str(st.cell.a),
                    "-Y",
                    str(st.cell.b),
                    "-Z",
                    str(st.cell.c),
                    "-al",
                    str(st.cell.alpha),
                    "-bt",
                    str(st.cell.beta),
                    "-gm",
                    str(st.cell.gamma),
                ]
            )

            result = read_pdb("wbox.pdb")
    result.spacegroup_hm = st.spacegroup_hm
    result.cell = st.cell
    return add_missing_b_factors(result, reference=st)
