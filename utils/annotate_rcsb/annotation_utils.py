import os
from typing import Union, Iterator, Tuple
import gemmi
import numpy as np

AnyPath = Union[str, bytes, os.PathLike]


class Monomer:
    """
    Represents a monomer structure loaded from a PDB file.

    Attributes:
        st: The structure loaded from the PDB file using Gemmi.

    Methods:
        __init__(path_to_pdb): Initializes the Monomer instance with the structure loaded from the specified PDB file.
        clean(): Cleans the structure.
        gaps_lengths(): Computes the lengths of gaps between residue IDs in the structure.
    """

    def __init__(self, path_to_pdb):
        self.st = gemmi.read_pdb(path_to_pdb)

    def remove_hetatm(self) -> "Monomer":
        for _, _, residue in self.iterate_over_residues():
            if residue.het_flag == "A":
                residue.flag = "A"
        selection = gemmi.Selection("/1").set_residue_flags("A")
        selection.remove_not_selected(self.st)
        return self

    def clean(self) -> "Monomer":
        self.remove_hetatm()
        self.st.add_entity_types(overwrite=True)
        self.st.remove_alternative_conformations()
        self.st.remove_ligands_and_waters()
        self.st.remove_empty_chains()
        return self

    def count_chains(self) -> int:
        return len(self.st[0])

    def write_pdb(self, path: AnyPath) -> None:
        return self.st.write_pdb(
            str(path), gemmi.PdbWriteOptions(minimal=True, numbered_ter=False)
        )

    def get_chain_id(self, missed_residues: list) -> str:
        for chain in self.st[0]:
            for residue in chain:
                if (
                    residue.seqid.num == missed_residues[0] - 1
                    and (residue.seqid.num + 1) != missed_residues[0]
                ):
                    return chain.name

    def get_gaps_start_end(self) -> list[list[Union[int, int]]]:
        residue_ids = []

        for chain in self.st[0]:
            for residue in chain:
                residue_ids.append(residue.seqid.num)

        residue_ids = np.array(residue_ids)
        distances = residue_ids[1:] - residue_ids[:-1] - 1
        gap_idx = np.where(distances > 0)
        gaps_lengths = distances[gap_idx]
        # print('gaps_lengths', gaps_lengths)
        gap_start_residue_id = residue_ids[gap_idx] + 1
        missed_residues = [
            [start, start + length - 1]
            for start, length in zip(gap_start_residue_id, gaps_lengths)
        ]

        return missed_residues

    def get_ss_rids(self) -> list[list[int]]:
        ss_rids = []

        # get sheets rids
        for sheet in self.st.sheets:
            for strand in sheet.strands:
                start_id = strand.start.res_id.seqid.num
                end_id = strand.end.res_id.seqid.num
                ss_rids += [rid for rid in range(start_id, end_id + 1)]

        # get helices rids
        for helix in self.st.helices:
            start_id = helix.start.res_id.seqid.num
            end_id = helix.end.res_id.seqid.num
            ss_rids += [rid for rid in range(start_id, end_id + 1)]
        return ss_rids

    def count_hetatm(self) -> int:
        hetatm_count = 0
        for chain in self.st[0]:
            for residue in chain:
                if residue.het_flag == "H" and residue.name != "HOH":
                    hetatm_count += 1
        return hetatm_count

    def iterate_over_residues(
        self,
    ) -> Iterator[Tuple[gemmi.Model, gemmi.Chain, gemmi.Residue]]:
        for model in self.st:
            for chain in model:
                for residue in chain:
                    yield model, chain, residue

    def count_unk(self) -> int:
        unk_count = 0
        for _, _, residue in self.iterate_over_residues():
            if residue.name == "UNK":
                unk_count += 1
        return unk_count
