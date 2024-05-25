from __future__ import annotations
import gemmi
import numpy as np
import pandas as pd


def get_rfree_flags(
        mtz: "gemmi.Mtz"
):
    """Write .tab file for pmemd.arx
    :param mtz: mtz file with P1 symmetry
    """
    import re

    R_FREE_FLAG = None
    r_free_pattern_string = r"(r.*free)|(free.*r)"
    r_free_pattern = re.compile(r_free_pattern_string, flags=re.IGNORECASE)

    for column in mtz.columns:  # type: gemmi.Mtz.Column
        if r_free_pattern.search(column.label):
            R_FREE_FLAG = column

    n_positive_r_flags = sum(R_FREE_FLAG)
    flag_is_one = n_positive_r_flags > len(R_FREE_FLAG) / 2

    r_free_flags = []

    for r_flag in R_FREE_FLAG:
        r = r_flag if flag_is_one else 1 - r_flag
        r_free_flags.append(r)

    return r_free_flags


def analyze_mtz(
        mtz: "gemmi.Mtz"
):

    r_free_flags_rcsb = np.array(get_rfree_flags(mtz))
    num_zeros = len(r_free_flags_rcsb) - np.count_nonzero(r_free_flags_rcsb)
    num_ones = np.count_nonzero(r_free_flags_rcsb)

    print("Размер rcsb:", len(r_free_flags_rcsb))
    print("Количество нулей в rcsb:", num_zeros)
    print("Количество единиц в rcsb:", num_ones)


def mtz_to_df(
        mtz: "gemmi.Mtz"
):
    mtz_data = np.array(mtz, copy=False)
    mtz_df = pd.DataFrame(data=mtz_data, columns=mtz.column_labels())

    return mtz_df


if __name__ == "__main__":
    # settings
    # add path to annotated .csv file if you need to process multiple pdb files
    # path_to_annotation_csv = "../1_annotated_rcsb/monomers_to_ucell.csv"
    # pdb_df = pd.read_csv(path_to_annotation_csv, sep=",")
    # pdb_ids = pdb_df["pdb_id"]

    # hardcoded. but there is id for my example structure
    pdb_ids = ["1k33"]

    # iterate over available pdb codes and initialize protocol
    for pdb_id in pdb_ids:
        mtz_rcsb = gemmi.read_mtz_file(f'mtz/{pdb_id}-sf.mtz')
        # expand original mtz to P1 symmetry group mtz
        mtz_rcsb.expand_to_p1()
        mtz_rcsb.ensure_asu()
        analyze_mtz(mtz_rcsb)
        mtz_rcsb.write_to_file(f'mtz/{pdb_id}.mtz')
