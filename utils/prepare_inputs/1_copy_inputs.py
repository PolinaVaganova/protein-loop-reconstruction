import os
import gemmi
from typing import Union
from hybrid_36 import hy36encode
from glob import glob
import pandas as pd


AnyPath = Union[str, bytes, os.PathLike]


def renumber_residues(st: gemmi.Structure) -> gemmi.Structure:
    chain_renum = 1
    first_res_num = None

    result = gemmi.Structure()
    for model in st:
        new_model = gemmi.Model(model.name)
        for chain in model:
            new_chain = gemmi.Chain(hy36encode(2, chain_renum))
            chain_renum += 1
            for residue in chain:
                if first_res_num is None:
                    first_res_num = residue.seqid.num
                residue.seqid.num = residue.seqid.num - first_res_num + 1
                new_chain.add_residue(residue)

            new_model.add_chain(new_chain)
        result.add_model(new_model)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def read_pdb(path: AnyPath) -> gemmi.Structure:
    return gemmi.read_pdb(str(path), split_chain_on_ter=True)


def copy_b_factor(structure,
                  original_structure):
    original_b_factor_dict = {}

    for model in original_structure:
        for chain in model:
            for residue in chain:
                original_b_factor_dict[residue.seqid.num] = {}
                for atom in residue:
                    original_b_factor_dict[residue.seqid.num].update({atom.name: atom.b_iso})

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if original_b_factor_dict.get(residue.seqid.num):
                        if original_b_factor_dict[residue.seqid.num].get(atom.name):
                            atom.b_iso = original_b_factor_dict[residue.seqid.num][atom.name]
                        else:
                            atom.b_iso = -1
                    else:
                        atom.b_iso = -1
    return structure


if __name__ == '__main__':

    # settings paths
    option = 'modeller'
    path_to_annotation_csv = '/home/polina/xray-refinement/1_annotated_rcsb/1u3y.csv'
    path_to_rcsb_data = '2_rcsb_data'
    path_to_original_pdbs = os.path.join(path_to_rcsb_data, 'pdb')
    path_to_original_mtz = os.path.join(path_to_rcsb_data, 'mtz')
    data_input_dir = '4_protocol_run/input/'
    path_to_loops = os.path.join('3_loop_modelling/', option)

    # read annotation df
    annotation_df = pd.read_csv(path_to_annotation_csv, sep=',')
    pdb_ids = annotation_df['pdb_id']

    for pdb_id in pdb_ids:
        path_to_out_dir = os.path.join(data_input_dir, f'{pdb_id}')

        if not os.path.exists(path_to_out_dir):
            os.makedirs(path_to_out_dir)

        chain_id = annotation_df[annotation_df['pdb_id'] == pdb_id]['chain_ids']

        # find modelled pdbs
        if option == 'rapper':
            path_to_modelled_loops = os.path.join(path_to_loops, pdb_id, f'grafted_loop')
            path_to_modelled_pdbs = glob(os.path.join(path_to_modelled_loops, f'{pdb_id}_rapper_*.pdb'))

        else:
            path_to_modelled_loops = os.path.join(path_to_loops, pdb_id)
            path_to_modelled_pdbs = glob(os.path.join(path_to_modelled_loops, '*Chain*.pdb'))
            print(len(path_to_modelled_pdbs))

        path_to_modelled_pdbs.sort()

        # open original pdb structure
        path_to_original_pdb = os.path.join(path_to_original_pdbs, f'{pdb_id}.pdb')
        st_original = read_pdb(path_to_original_pdb)

        if option == 'modeller':
            st_original = renumber_residues(st_original)

        # iterate over all modelled structures to prepare it
        for modelled_pdb in path_to_modelled_pdbs:
            st_modelled = read_pdb(modelled_pdb)

            # add b-factors
            st_modelled_b_factors = copy_b_factor(st_modelled, st_original)

            # add headers
            original_header = st_original.make_pdb_headers()

            # write output pdb file
            path_to_fout = os.path.join(path_to_out_dir, os.path.basename(modelled_pdb))
            st_modelled_b_factors.setup_entities()
            st_modelled_b_factors.assign_label_seq_id()
            st_modelled_b_factors.entities = st_original.entities

            st_modelled_b_factors_str = st_modelled_b_factors.make_pdb_string(
                gemmi.PdbWriteOptions(cryst1_record=False, end_record=True, ter_records=True))

            with open(path_to_fout, "w") as fout:
                fout.writelines(original_header)
                fout.writelines(st_modelled_b_factors_str)

            st_test = read_pdb(path_to_fout)

            # make symlink for mtz file
            if not os.path.islink(os.path.join(path_to_out_dir, f'{pdb_id}.mtz')):
                os.symlink(os.path.join(os.getcwd(), path_to_original_mtz, f'gemmi_{pdb_id}_p1.mtz'),
                           os.path.join(path_to_out_dir, f'{pdb_id}.mtz'))
