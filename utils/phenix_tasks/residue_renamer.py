import sys

import gemmi


def rename_residues_amber_to_standard(st: gemmi.Structure):
    def get_residues(st: gemmi.Structure):
        assert len(st) == 1, "Structure MUST have one MODEL"
        for chain in st[0]:
            for residue in chain:
                yield residue

    rename_map = {
        "GLH": "GLU",
        "ASH": "ASP",
        "HID": "HIS",
        "HIP": "HIS",
        "HIE": "HIS",
        "CYX": "CYS",
        "CYM": "CYS",
        "LYN": "LYS",
        "G5": "G",  # RNA 5-end
        "A5": "A",
        "C5": "C",
        "U5": "U",
        "G3": "G",  # RNA 3-end
        "A3": "A",
        "C3": "C",
        "U3": "U",
        "DG5": "DG",  # DNA 5-end
        "DA5": "DA",
        "DC5": "DC",
        "DT5": "DT",
        "DG3": "DG",  # DNA 3-end
        "DA3": "DA",
        "DC3": "DC",
        "DT3": "DT",
    }

    for res in get_residues(st):
        res.name = rename_map.get(res.name, res.name)


def rename_residues_in_pdb(filename: str):
    input_file_name = filename
    output_file_name = filename.strip(".pdb") + "_renamed.pdb"
    final = gemmi.read_pdb(input_file_name, split_chain_on_ter=True)
    rename_residues_amber_to_standard(final)
    final.write_pdb(output_file_name, numbered_ter=False, ter_ignores_type=True)


if __name__ == "__main__":
    try:
        rename_residues_in_pdb(sys.argv[1])
    except FileNotFoundError:
        print("File does not exist.")
    except ValueError:
        print("File does not have a proper extension.")
