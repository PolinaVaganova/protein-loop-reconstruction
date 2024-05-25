import os
import re
import subprocess
import tempfile
import typing
from pathlib import Path

import gemmi
import numpy as np
import parmed
import propka.run
import sander
from pdb4amber.residue import AMBER_SUPPORTED_RESNAMES, RESPROT
from scipy import optimize

from hybrid_36 import hy36encode

from .utils import chdir, check_call

AnyPath = typing.Union[str, bytes, os.PathLike]

tleap_in = """
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3

HOH = SPC
WAT = SPC

loadAmberParams frcmod.ionsjc_spce
loadAmberParams frcmod.ions1lm_1264_spce
loadAmberParams frcmod.spce

wbox = loadpdb input.pdb
set default nocenter on
setBox wbox vdw 1.0
saveamberparm wbox wbox.prmtop wbox.inpcrd
savepdb wbox wbox.dry.pdb
quit
"""


def add_ions(
    st: gemmi.Structure, ion: gemmi.Structure, count: int, *, n_protein_atoms: int
) -> gemmi.Structure:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            write_pdb(ion, "ion.pdb")
            check_call(
                [
                    "AddToBox",
                    "-c",
                    "input.pdb",
                    "-a",
                    "ion.pdb",
                    "-na",
                    str(count),
                    "-o",
                    "result.pdb",
                    "-P",
                    str(n_protein_atoms),
                    "-RP",
                    str(3.0),
                    "-RW",
                    str(6.0),
                    "-G",
                    str(0.2),
                    "-V",
                    str(1),
                ]
            )
            return read_pdb("result.pdb")


def add_water(st: gemmi.Structure, water: gemmi.Structure) -> gemmi.Structure:
    n_non_water_atoms = st[0].count_atom_sites()
    count = estimate_water_molecules(st)
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            write_pdb(water, "water.pdb")
            check_call(
                [
                    "AddToBox",
                    "-c",
                    "input.pdb",
                    "-a",
                    "water.pdb",
                    "-na",
                    str(count),
                    "-o",
                    "result.pdb",
                    "-P",
                    str(n_non_water_atoms),
                    "-RP",
                    str(3.0),
                    "-RW",
                    str(3.0),
                    "-G",
                    str(0.2),
                    "-V",
                    str(1),
                ]
            )
            return read_pdb("result.pdb")


def calc_total_charge(st: gemmi.Structure) -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            with open("tleap.in", "w") as f:
                f.write(tleap_in)

            with open("parmed.in", "w") as f:
                f.write("summary\n")
                f.write("quit\n")

            check_call(["tleap", "-s", "-f", "tleap.in"])

            parmed_output = subprocess.check_output(
                [
                    "parmed",
                    "wbox.prmtop",
                    "parmed.in",
                ]
            ).decode("utf-8")

    for line in parmed_output.split("\n"):
        if "Total charge (e-)" in line:
            return int(float(line.split(":")[-1].strip().split()[0]))

    raise RuntimeError("Can't find Total charge in parmed output")


def estimate_water_molecules(st: gemmi.Structure) -> int:
    psv = 0.74
    MW = estimate_weight(st)
    nau = 1
    # 100K->RT volume expansion coefficient
    expand = 1.0  # no expansion
    # ==============================================================================
    # Volume of protein
    VP = MW * nau * psv * 1e24 / 6.02e23
    V = st.cell.volume
    # print("Water content: %.3f" % (1 - VP / V))
    V *= expand
    # estimated water molecule numbers
    nwat = (V - VP) * 1e-24 * 1.00 / 18.0 * 6.02e23
    assert nwat > 10, f"Too low number of water molecules: {nwat}"
    return int(nwat)


def estimate_weight(st: gemmi.Structure) -> float:
    mass = 0.0
    aname_map = dict([(n, gemmi.Element(n).weight) for n in ["Na+", "Cl-"]])
    element_map = dict(
        [(n, gemmi.Element(n).weight) for n in ["P", "H", "C", "S", "O", "N"]]
    )
    for m in st:
        for c in m:
            for r in c:
                for a in r:
                    if a.name in aname_map:
                        mass += aname_map[a.name]
                    else:
                        mass += element_map[a.name[0]]

    return mass


def expand_non_crystallographic_symmetries(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.expand_ncs(how=gemmi.HowToNameCopiedChain.AddNumber)
    return result


def expand_crystallographic_symmetries(st: gemmi.Structure) -> gemmi.Structure:
    g: gemmi.SpaceGroup = gemmi.find_spacegroup_by_name(st.spacegroup_hm)

    def apply_to_chain(chain, op, cell):
        for residue in chain:
            for atom in residue:
                fractional = cell.fractionalize(atom.pos)
                transformed = gemmi.Fractional(*op.apply_to_xyz(fractional.tolist()))
                cartesian = cell.orthogonalize(transformed)
                atom.pos = cartesian

    result = gemmi.Structure()

    for model in st:
        new_model = gemmi.Model(model.name)
        for op in g.operations():
            for chain in model:
                new_model.add_chain(chain, pos=-1)
                new_chain = new_model[-1]
                apply_to_chain(new_chain, op, cell=st.cell)
        result.setup_entities()
        result.add_model(new_model, pos=-1)
    result.cell = st.cell
    result.spacegroup_hm = "P1"  # st.spacegroup_hm
    return result


def expand_supercell(st: gemmi.Structure, supercell_size: int) -> gemmi.Structure:
    if supercell_size < 2:
        return st

    assert st.spacegroup_hm == "P1"

    count = 0
    for model in st:
        for _ in model:
            count += 1
    # check if there is a single model in the structure to be expanded
    if count > 2:
        return st

    def apply_to_chain(chain, cell, dim, shift=1):
        for residue in chain:
            for atom in residue:
                fractional = cell.fractionalize(atom.pos)
                transformed = fractional
                transformed[dim] += shift
                cartesian = cell.orthogonalize(transformed)
                atom.pos = cartesian

    result = gemmi.Structure()

    for model in st:
        for dim in [0, 1, 2]:
            new_model = gemmi.Model(model.name)
            for chain in model:
                new_model.add_chain(chain, pos=-1)
            for chain in model:
                for shift in range(1, supercell_size):
                    new_model.add_chain(chain, pos=-1)
                    new_chain = new_model[-1]
                    apply_to_chain(new_chain, cell=st.cell, dim=dim, shift=shift)
            model = new_model
        result.add_model(new_model, pos=-1)
    result.cell = st.cell
    result.cell.set(
        result.cell.a * supercell_size,
        result.cell.b * supercell_size,
        result.cell.c * supercell_size,
        result.cell.alpha,
        result.cell.beta,
        result.cell.gamma,
    )
    return result


def get_target_ph(st: gemmi.Structure) -> float:
    def extract200(line):
        print(line)
        m = re.findall(r"\d*\.\d+|\d+", line.split(":")[-1].upper().strip())
        m = [float(m_) for m_ in m]
        letter_presence = (
            len(re.findall(r"[a-zA-z]", line.split(":")[-1].upper().strip())) > 0
            and len(m) > 0
        )
        if len(m) == 1 or letter_presence:
            return m[0]
        elif len(m) > 1:
            return sum(m) / len(m)
        return None

    def extract280(line):
        print(line)
        m = re.findall(r"PH (?:\d*\.\d+|\d+)", line.upper().strip())
        m = [float(m_[2:]) for m_ in m]
        if len(m) > 0:
            return sum(m) / len(m)
        return None

    header = st.make_pdb_headers().splitlines()
    # print(header)
    ph = None
    for header_line in header:
        if "REMARK 200  PH " in header_line.upper():
            ph = extract200(header_line.strip())
            if ph is not None:
                break
    if ph is None:
        l_accumulate = ""
        accumulate = False
        for header_line in header:
            if "REMARK 280 CRYSTALLIZATION" in header_line.upper():
                accumulate = True
            if accumulate:
                l_accumulate += header_line.upper()[11:-1]
                if "280" not in header_line:
                    accumulate = False
        ph = extract280(l_accumulate.strip())
    if ph is not None:
        result = ph
    else:
        result = 7.5
    return result


def neutralize_with_ions(
    st: gemmi.Structure, negative_ion: gemmi.Structure, positive_ion: gemmi.Structure
) -> gemmi.Structure:
    total_charge = calc_total_charge(st)
    if total_charge > 0:
        ion = negative_ion
    else:
        ion = positive_ion
    count = abs(total_charge)
    n_protein_atoms = st[0].count_atom_sites()
    return add_ions(st, ion, count=count, n_protein_atoms=n_protein_atoms)


def find_gaps(st: gemmi.Structure) -> (gemmi.Structure, typing.List[int]):
    # N.B.: following only finds gaps in protein chains!
    # H.N: Assume that residue has all 3 atoms: CA, C, and N
    respro_nocap = set(RESPROT) - {"ACE", "NME"}
    result = gemmi.Structure()
    is_ter = False
    gaplist = []
    for model in st:
        new_model = gemmi.Model(model.name)
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            # at this point there are no empty chains
            # N.B. the procedure will fail for free-floating residues!
            for i in range(len(chain) - 1):
                residue = chain[i]
                if residue.name in respro_nocap:
                    C_atom = residue.find_atom("C", " ", gemmi.Element("C"))
                    N_atom = chain[i + 1].find_atom("N", " ", gemmi.Element("N"))
                    gap = C_atom.pos.dist(N_atom.pos)
                    if gap > 2.5:
                        gaplist.append(residue.seqid.num)
                        print(residue.seqid.num, residue.name, gap)
                        is_ter = True
                    else:
                        is_ter = False
                new_chain.add_residue(residue)
                if is_ter:
                    new_model.add_chain(new_chain)
                    new_chain = gemmi.Chain(chain.name)
            new_chain.add_residue(chain[-1])
            new_model.add_chain(new_chain)
        result.add_model(new_model)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result, gaplist


def apply_additional_ters(
    st: gemmi.Structure, gaplist: typing.List[int]
) -> gemmi.Structure:
    result = gemmi.Structure()
    for model in st:
        new_model = gemmi.Model(model.name)
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            for i in range(len(chain) - 1):
                residue = chain[i]
                if residue.seqid.num in gaplist:
                    is_ter = True
                else:
                    is_ter = False
                new_chain.add_residue(residue)
                if is_ter:
                    print(residue.seqid.num)
                    new_model.add_chain(new_chain)
                    new_chain = gemmi.Chain(chain.name)
            new_chain.add_residue(chain[-1])
            new_model.add_chain(new_chain)
        result.add_model(new_model)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def extract_pkas(protein, conformation, parameters, target_ph):
    rename_residues = ["GLU", "ASP", "HIS", "CYS", "LYS"]
    rename_map_acids = {
        "GLU": "GLH",
        "ASP": "ASH",
    }
    protonation = {}
    str_ = ""
    # printing pKa summary
    for residue_type in parameters.write_out_order:
        for group in protein.conformations[conformation].groups:
            if group.residue_type == residue_type:
                # print(group.coupled_titrating_group)
                str__ = (
                    f"{group.residue_type:>9s} {group.atom.res_num:>4d} "
                    f"{group.atom.chain_id:>9s} {group.pka_value:8.2f}\n"
                )
                str_ += str__
                if group.residue_type in rename_residues:
                    if group.pka_value < target_ph:  # deprotonated
                        # protonation[
                        #     (
                        #         group.residue_type,
                        #         group.atom.chain_id,
                        #         group.atom.res_num,
                        #     )
                        # ] = False
                        if group.residue_type == "CYS":
                            protonation[
                                (
                                    group.residue_type,
                                    group.atom.chain_id,
                                    group.atom.res_num,
                                )
                            ] = "CYM"
                        elif group.residue_type == "LYS":
                            protonation[
                                (
                                    group.residue_type,
                                    group.atom.chain_id,
                                    group.atom.res_num,
                                )
                            ] = "LYN"
                    else:
                        # protonation[
                        #     (
                        #         group.residue_type,
                        #         group.atom.chain_id,
                        #         group.atom.res_num,
                        #     )
                        # ] = True
                        if group.residue_type in rename_map_acids:  # protonated
                            protonation[
                                (
                                    group.residue_type,
                                    group.atom.chain_id,
                                    group.atom.res_num,
                                )
                            ] = rename_map_acids[group.residue_type]
                        if group.residue_type == "HIS":
                            protonation[
                                (
                                    group.residue_type,
                                    group.atom.chain_id,
                                    group.atom.res_num,
                                )
                            ] = "HIP"
                # else:
                #     print(str__.strip())
    return protonation


def assign_protonation_states(
    st: gemmi.Structure, target_ph: float = 7.5
) -> gemmi.Structure:
    result = st.clone()
    cys_possible_names = ["CYS", "CYM", "CYX"]
    index: gemmi.NeighborSearch = gemmi.NeighborSearch(st, 3)
    index.populate(include_h=False)
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            tmp_pdb = "tmp.pdb"
            write_pdb(st, tmp_pdb)
            # propka step
            mol = propka.run.single(tmp_pdb, write_pka=False)
            protonation_states = extract_pkas(
                mol, "AVR", mol.version.parameters, target_ph
            )
            chains = set([k[1] for k in protonation_states.keys()])
            for chain in result[0]:
                chain_name = str(chain.name).strip()
                if chain_name in chains:
                    for residue in chain:
                        key = (str(residue.name), chain_name, residue.seqid.num)
                        if key in protonation_states.keys():
                            residue.name = protonation_states[key]
                for residue in chain:
                    # histidines
                    if residue.name == "HIS":
                        atom_name_set = sorted(
                            set(atom.name for atom in residue if atom.is_hydrogen())
                        )
                        # this should be taken care of by now
                        # if set(['HD1', 'HE2']).issubset(atom_name_set):
                        #     residue.name = 'HIP'
                        if "HD1" in atom_name_set:
                            residue.name = "HID"
                        elif "HE2" in atom_name_set:
                            residue.name = "HIE"
                    # disulfides
                    if residue.name in cys_possible_names:
                        atom = residue.find_atom("SG", " ", gemmi.Element("S"))
                        closest = index.find_neighbors(atom, 1.8, 2.3)
                        for mark in closest:
                            ref_at = mark.to_cra(st[0]).atom
                            ref_res = mark.to_cra(st[0]).residue
                            if (
                                ref_res.name in cys_possible_names
                                and ref_at.name == "SG"
                            ):
                                residue.name = "CYX"
                                break
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def add_missing_atoms(st: gemmi.Structure) -> gemmi.Structure:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            input_pdb = "input.pdb"
            result_pdb = "wbox.dry.pdb"
            result = st.clone()
            result.remove_hydrogens()
            write_pdb(result, input_pdb)
            with open("tleap.in", "w") as f:
                f.write(tleap_in)
            check_call(["tleap", "-f", "tleap.in"])
            result = read_pdb(result_pdb)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def prepare_5_terminal_nucleic_acids(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    for model in result:
        for chain in model:
            first_res = True
            for residue in chain:
                if len(residue.name) < 3 and first_res:
                    del_P = False
                    del_OP1 = False
                    del_OP2 = False
                    for atom in residue:
                        if atom.name == "P":
                            del_P = True
                        elif atom.name == "OP1":
                            del_OP1 = True
                        elif atom.name == "OP2":
                            del_OP2 = True
                    if del_P:
                        residue.remove_atom("P", " ", gemmi.Element("P"))
                    if del_OP1:
                        residue.remove_atom("OP1", " ", gemmi.Element("O"))
                    if del_OP2:
                        residue.remove_atom("OP2", " ", gemmi.Element("O"))

                first_res = False
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def minimize_in_python(st: gemmi.Structure, bellymask: str) -> gemmi.Structure:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            tmp_pdb = "tmp.pdb"

            result = st.clone()
            result.cell = st.cell
            result.spacegroup_hm = st.spacegroup_hm

            from arx.amber import create_topology_and_input

            create_topology_and_input(
                st, Path.cwd() / "wbox.prmtop", Path.cwd() / "wbox.inpcrd"
            )
            prmtop = parmed.load_file(
                str(Path.cwd() / "wbox.prmtop"), xyz=str(Path.cwd() / "wbox.inpcrd")
            )

            # sander energy minimization
            inp = sander.gas_input(6)
            # inp.ntc = 2
            # inp.ntf = 2
            # inp.ntr = 1
            # inp.restraint_wt = 10
            # inp.restraintmask = f"!@{bellymask} & !@H="
            # inp.ibelly = 1
            # inp.bellymask=f"@{bellymask} | @H="

            def energy_function(var):
                sander.set_positions(var)
                ene, frc = sander.energy_forces()
                return ene.tot, -np.array(frc)

            with sander.setup(
                str(Path.cwd() / "wbox.prmtop"), prmtop.coordinates, prmtop.box, inp
            ) as sander_context:
                print(sander_context)
                # min_method = 'L-BFGS-B'
                min_method = "CG"
                min_result = optimize.minimize(
                    energy_function,
                    prmtop.coordinates,
                    method=min_method,
                    jac=True,
                    options=dict(maxiter=50, disp=True, gtol=0.01),
                )
                sander.set_positions(min_result.x)
                optimized_coordinates = sander.get_positions(as_numpy=True).reshape(
                    (sander_context.natom, 3)
                )

            prmtop.save(
                tmp_pdb, coordinates=optimized_coordinates, standard_resnames=False
            )

            coords = read_pdb(tmp_pdb)
            result = copy_coordinates(result, coords)

    return result


def apply_to_atoms(
    st: gemmi.Structure, atom_mutator: typing.Callable[[gemmi.Atom], None]
) -> gemmi.Structure:
    result = st.clone()
    for model in result:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_mutator(atom)
    return result


def set_b_factors_to(st: gemmi.Structure, value: float = 0) -> gemmi.Structure:
    def setter(atom: gemmi.Atom):
        atom.b_iso = value

    return apply_to_atoms(st, atom_mutator=setter)


def add_missing_b_factors(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    index: gemmi.NeighborSearch = gemmi.NeighborSearch(reference, 3.0)
    index.populate()

    result = set_b_factors_to(st, -1)
    total_assigned = 0
    shell_step = 0.5
    min_distance = 0
    max_distance = shell_step
    is_first_pass = True
    while True:
        n_assigned = 0
        unassigned_exist = False
        for model in result:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.b_iso > 0:
                            continue
                        unassigned_exist = True
                        closest = index.find_neighbors(atom, min_distance, max_distance)
                        max_b_factor = -1
                        for mark in closest:
                            if is_first_pass:
                                ref_at = mark.to_cra(reference[0]).atom
                            else:
                                ref_at = mark.to_cra(model).atom
                            dist = ref_at.pos.dist(atom.pos)
                            if dist < max_distance:
                                max_b_factor = max(max_b_factor, ref_at.b_iso)
                        if max_b_factor > 0:
                            atom.b_iso = max_b_factor
                            n_assigned += 1
        if is_first_pass:
            index: gemmi.NeighborSearch = gemmi.NeighborSearch(result, 3.0)
            index.populate()
            is_first_pass = False
        total_assigned += n_assigned
        min_distance = max_distance
        max_distance += shell_step
        shell_step = 2.0

        if n_assigned == 0 and not unassigned_exist:
            break

    return result


def add_missing_occupancies(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    def setter(atom: gemmi.Atom):
        atom.occ = 1.0

    return apply_to_atoms(st, setter)


def remove_alternative_conformations(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_alternative_conformations()
    return result


def remove_ligands_and_water(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_ligands_and_waters()
    return result


def remove_water(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_waters()
    return result


def remove_empty_chains(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_empty_chains()
    return result


def remove_hydrogens(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r), (_, _, ref_r) in zip(
        iterate_over_residues(result), iterate_over_residues(reference)
    ):
        # condition for NAs
        if len(ref_r.name) < 3:
            unter_na = st_r.name
            unter_na = unter_na.replace("3", "").replace("5", "")
            condition = unter_na == ref_r.name
        # condition for histidines
        elif ref_r.name.startswith("HI"):
            condition = st_r.name.startswith("HI") and ref_r.name.startswith("HI")
        else:
            condition = st_r.name == ref_r.name
        assert condition
        to_del = []
        for h_atom in st_r:
            if h_atom.is_hydrogen():
                if h_atom.name not in [atom.name for atom in ref_r]:
                    to_del.append(h_atom)
        for h_atom in to_del[::-1]:
            st_r.remove_atom(h_atom.name, " ", gemmi.Element("H"))
        if st_r.name.startswith("HI"):
            st_r.name = ref_r.name
        if len(ref_r.name) < 3:
            st_r.name = ref_r.name
    return result


def read_pdb(path: AnyPath) -> gemmi.Structure:
    return gemmi.read_pdb(str(path), split_chain_on_ter=True)


def write_pdb(st: gemmi.Structure, path: AnyPath):
    return st.write_pdb(str(path), numbered_ter=False, ter_ignores_type=True)


def copy_coordinates(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r, st_a), (_, _, ref_r, ref_a) in zip(
        iterate_over_atoms(result), iterate_over_atoms(reference)
    ):
        assert st_r.name == ref_r.name
        assert st_a.name == ref_a.name
        st_a.pos = ref_a.pos
    return result


def iterate_over_atoms(
    st: gemmi.Structure,
) -> typing.Iterator[typing.Tuple[gemmi.Model, gemmi.Chain, gemmi.Residue, gemmi.Atom]]:
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    yield model, chain, residue, atom


def iterate_over_residues(
    st: gemmi.Structure,
) -> typing.Iterator[typing.Tuple[gemmi.Model, gemmi.Chain, gemmi.Residue]]:
    for model in st:
        for chain in model:
            for residue in chain:
                yield model, chain, residue


def copy_residue_names(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r), (_, _, ref_r) in zip(
        iterate_over_residues(result), iterate_over_residues(reference)
    ):
        if st_r.name != ref_r.name:
            # very crude check that there was no accidental shift so far
            assert st_r.name[:2] == ref_r.name[:2]
            st_r.name = ref_r.name
    return result


def check_chain_and_residue_numbering(
    st: gemmi.Structure, ref: gemmi.Structure, strict: bool = True
) -> bool:
    result = True
    for (_, st_c, st_r), (_, ref_c, ref_r) in zip(
        iterate_over_residues(st), iterate_over_residues(ref)
    ):
        if strict:
            res_names_equal = st_r.name == ref_r.name
        else:
            res_names_equal = st_r.name[:2] == ref_r.name[:2]
        result = (
            result
            and st_c.name == ref_c.name
            and st_r.seqid.num == ref_r.seqid.num
            and res_names_equal
        )
    return result


def retain_only_standard_resnames(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    for _, _, residue in iterate_over_residues(result):
        assert residue.name in AMBER_SUPPORTED_RESNAMES
    return result


def renumber_residues(st: gemmi.Structure) -> gemmi.Structure:
    res_renum = 1
    chain_renum = 1

    result = gemmi.Structure()
    for model in st:
        new_model = gemmi.Model(model.name)
        for chain in model:
            new_chain = gemmi.Chain(hy36encode(2, chain_renum))
            chain_renum += 1
            for residue in chain:
                residue.seqid.num = res_renum
                new_chain.add_residue(residue)
                res_renum += 1
            new_model.add_chain(new_chain)
        result.add_model(new_model)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result
