import gemmi


def pos_str(pos: gemmi.Position) -> str:
    return f"{{{pos.x:8.3f}, {pos.y:8.3f}, {pos.z:8.3f}}}"


def show_diff(
    ref: gemmi.Structure,
    st: gemmi.Structure,
    check_occupancy: bool,
    occupancy_tolerance: float,
    check_b_factor: bool,
    b_factor_tolerance: float,
    check_coordinate: bool,
    coordinate_tolerance: float,
    ignore_model_name: bool,
    ignore_chain_name: bool,
    ignore_residue_id: bool,
    ignore_unknown_atom: bool,
):
    model: gemmi.Model
    chain: gemmi.Chain
    residue: gemmi.Residue
    atom: gemmi.Atom

    def to_key(model, chain, residue, residue_seq_idx, atom):
        result = ()
        if not ignore_model_name:
            result = (*result, model.name)
        if not ignore_chain_name:
            result = (*result, chain.name)
        if ignore_residue_id:
            result = (*result, residue_seq_idx)
        else:
            result = (*result, residue.seqid.num)
        result = (*result, residue.name, atom.name)
        return result

    ref_atom_map = {
        to_key(model, chain, residue, residue_seq_idx, atom): atom
        for model in ref
        for chain in model
        for residue_seq_idx, residue in enumerate(chain)
        for atom in residue
    }
    n_missing = 0
    n_ok = 0
    n_bad = 0
    n_checked = 0
    n_bad_occupancy = 0
    n_bad_b_factor = 0
    n_bad_coordinates = 0
    for model in st:
        for chain in model:
            for residue_seq_idx, residue in enumerate(chain):
                for atom in residue:
                    key = to_key(model, chain, residue, residue_seq_idx, atom)
                    try:
                        ref_atom: gemmi.Atom = ref_atom_map[key]
                    except KeyError:
                        if not ignore_unknown_atom:
                            print(f"Unknown atom {key}")
                            n_missing += 1
                        continue
                    n_checked += 1
                    error_lines = []
                    if check_occupancy:
                        delta = abs(ref_atom.occ - atom.occ)
                        if delta > occupancy_tolerance:
                            n_bad_occupancy += 1
                            error_lines += [
                                f"Occupancy delta: {delta:7.2f}   "
                                f":: {ref_atom.occ:.7.2f}  -> {atom.occ:7.2f} "
                            ]
                    if check_b_factor:
                        delta = abs(ref_atom.b_iso - atom.b_iso)
                        if delta > b_factor_tolerance:
                            n_bad_b_factor += 1
                            error_lines += [
                                f"B-factor delta : {delta:7.2f}  "
                                f"::  {ref_atom.b_iso:7.2f}  -> {atom.b_iso:7.2f} "
                            ]
                    if check_coordinate:
                        delta = ref_atom.pos.dist(atom.pos)
                        if delta > coordinate_tolerance:
                            n_bad_coordinates += 1
                            error_lines += [
                                f"Coord delta    : {delta:8.3f}  "
                                f"::  {pos_str(ref_atom.pos)}  ->  {pos_str(atom.pos)}"
                            ]
                    if error_lines:
                        error_str = "\t" + "\n\t".join(error_lines)
                        print(f"{key}:\n{error_str}")
                        n_bad += 1
                    else:
                        n_ok += 1

    print(f"Checked atoms : {n_checked:8}")
    print(f"Unknown atoms : {n_missing:8}  # i.e. missed in reference structure")
    print(f"Good atoms    : {n_ok:8}")
    print(f"Bad atoms     : {n_bad:8}")
    if check_occupancy:
        print(f"Bad occupancy : {n_bad_occupancy:8}")
    else:
        print("Occupancy check skipped")
    if check_b_factor:
        print(f"Bad b-factors : {n_bad_b_factor:8}")
    else:
        print("B-factor check skipped")
    if check_coordinate:
        print(f"Bad coords    : {n_bad_coordinates:8}")
    else:
        print("Coordinates check skipped")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Show difference")
    parser.add_argument("--check-occupancy", default=False, action="store_true")
    parser.add_argument("--occupancy-tolerance", default=0.0, type=float)
    parser.add_argument("--check-b-factor", default=False, action="store_true")
    parser.add_argument("--b-factor-tolerance", default=0.0, type=float)
    parser.add_argument("--check-coordinate", default=False, action="store_true")
    parser.add_argument("--coordinate-tolerance", default=0.0, type=float)
    parser.add_argument("--ignore-chain-name", default=False, action="store_true")
    parser.add_argument("--ignore-model-name", default=False, action="store_true")
    parser.add_argument("--ignore-residue-id", default=False, action="store_true")
    parser.add_argument("--ignore-unknown-atom", default=False, action="store_true")

    parser.add_argument(dest="reference_pdb", metavar="REFERENCE_PDB")
    parser.add_argument(dest="contest_pdb", metavar="CONTEST_PDB")

    args = parser.parse_args()

    ref = gemmi.read_pdb(args.reference_pdb)
    st = gemmi.read_pdb(args.contest_pdb)

    show_diff(
        ref,
        st,
        args.check_occupancy,
        args.occupancy_tolerance,
        args.check_b_factor,
        args.b_factor_tolerance,
        args.check_coordinate,
        args.coordinate_tolerance,
        args.ignore_model_name,
        args.ignore_chain_name,
        args.ignore_residue_id,
        args.ignore_unknown_atom,
    )


if __name__ == "__main__":
    main()
