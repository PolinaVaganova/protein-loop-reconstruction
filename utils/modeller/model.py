import os
import re
import pandas as pd
from modeller import *
from modeller.automodel import *


def get_seqs() -> str:
    with open('alignment.ali', 'r') as ali:
        pattern = r':\s*-?\d+(\.\d+)?'
        start_reading = False
        seqs = []
        seq = ''
        for line in ali:
            if re.search(pattern, line):
                start_reading = True
                continue
            elif start_reading:
                seq += line.strip()
                if '*' in line:
                    seqs.append(seq)
                    seq = ''
                    start_reading = False
        return seqs


def remove_tails() -> None:
    pattern = r':\s*-?\d+(\.\d+)?'
    seq_region = False

    with open('alignment.ali', 'r') as ali:
        lines = ali.readlines()

    seq_num = 0
    first_line_in_fasta_seq = True
    len_start_gap = None
    len_end_gap = None
    for i, line in enumerate(lines):
        if re.search(pattern, line):
            seq_region = True
            seq_num += 1
            continue

        elif seq_region and line.startswith('-') and seq_num == 1:
            new_line = line.strip('-')
            len_start_gap = len(line) - len(new_line)
            lines[i] = new_line

        elif seq_region and line.strip().endswith('-*') and seq_num == 1:
            new_line = line.strip().strip('-*') + '*'
            len_end_gap = len(line) - len(new_line) - 1
            lines[i] = new_line
            seq_region = False

        elif seq_region and seq_num == 2:
            if first_line_in_fasta_seq and len_start_gap:
                lines[i] = line[len_start_gap:]
                first_line_in_fasta_seq = False

            if line.strip().endswith('*') and len_end_gap:
                lines[i] = line.strip()[:-len_end_gap - 1] + '*' + '\n'
                seq_region = False

    with open('alignment.ali', 'w') as file:
        file.writelines(lines)


def get_gap_indexes(seq: str) -> list[int]:
    start_sense_part = False
    start_gap = False
    gap_indexes = []
    for idx, chr in enumerate(seq):
        if chr != '-' and not start_sense_part:
            start_sense_part = True
        if chr == '-' and start_sense_part:
            gap_indexes.append(idx + 1)
            start_gap = True
            start_sense_part = False
        if chr != '-' and start_gap:
            gap_indexes.append(idx)
            break
    return gap_indexes


def extract_name_from_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        first_line = file.readline().strip()
        if first_line.startswith('>'):
            name = first_line.split('|')[0][1:]
            return f'{name}|Chain'
        else:
            print('Invalid FASTA format: Missing header line.')
            return None


if __name__ == "__main__":
    # settings
    cwd = os.getcwd()
    path_to_pdb_dir = '../../../2_rcsb_data/pdb/'
    path_to_fasta_dir = '../../../2_rcsb_data/fasta/'

    # add path to annotated .csv file if you need to process multiple pdb files
    # path_to_csv = ''
    # pdb_df = pd.read_csv(path_to_csv, sep=',')
    # pdb_ids = pdb_df['pdb_id']

    pdb_ids = ['1k33']
    
    for pdb_id in pdb_ids:
        if not os.path.exists(pdb_id):
            os.chdir(cwd)
            path_to_pdb_modelling = os.path.join(cwd, pdb_id)
            os.makedirs(path_to_pdb_modelling, exist_ok=True)
            os.chdir(path_to_pdb_modelling)
            path_to_pdb = os.path.join(path_to_pdb_dir, f'{pdb_id}.pdb')

            if not os.path.exists(f'{pdb_id}.pdb'):
                os.symlink(path_to_pdb, f'{pdb_id}.pdb')

            env = Environ()
            mdl = Model(env)
            mdl.read(file=pdb_id, model_segment=('FIRST:@', 'END:'))
            aln = Alignment(env)
            aln.append_model(mdl, align_codes=pdb_id, atom_files=pdb_id)
            aln.write(file=pdb_id + '.seq')
            path_to_fasta_dir = os.path.join(path_to_fasta_dir, pdb_id)

            if not os.path.exists(f'{pdb_id}.fasta'):
                os.symlink(path_to_fasta_dir, f'{pdb_id}.fasta')

            aln.append(file=f'{pdb_id}.fasta', align_codes='all', alignment_format='FASTA')
            aln.align(gap_penalties_1d=(-600, -400))
            aln.write(file='alignment.ali', align_alignment=False)
            remove_tails()

            seqs = get_seqs()
            gap_idxs = get_gap_indexes(seqs[0])

            log.verbose()

            # directories for input atom files
            env.io.atom_files_directory = ['.', '.']

            seq_name = extract_name_from_fasta(f'{pdb_id}.fasta')


            class MyModel(AutoModel):
                def select_atoms(self):
                    return Selection(self.residue_range(f'{gap_idxs[0]}:A',
                                                        f'{gap_idxs[1]}:A'),
                                     )


            a = MyModel(env, alnfile='alignment.ali',
                        knowns=pdb_id, sequence=seq_name)
            a.starting_model = 1
            a.ending_model = 100
            a.make()
