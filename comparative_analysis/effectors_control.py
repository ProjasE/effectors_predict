import re
import glob
import os
from contextlib import contextmanager
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO

# Secretome Ids
secretome_filepath = 'secretome_all.fasta'
secretome_proteinIds = []

for record in SeqIO.parse(secretome_filepath, "fasta"):
    secretome_proteinIds.append(record.id)

# Orthologs data structure
orthologs_filepath = "09_final_groups/Phytophthora_OrthologousGroups.txt"

orthoGroup2protName = {}

with open(orthologs_filepath) as fp:
    for cnt, line in enumerate(fp):
        ortho_group = re.search('^(\S+):', line).group(1)
        protein_names = re.search(': (.+$)', line).group(1).split(' ')

        protein_names_mod = []
        for pn in protein_names:
            pn_mod = re.search('\|(.+$)', pn).group(1)
            if pn_mod in secretome_proteinIds:
                protein_names_mod.append(pn_mod)

        if protein_names_mod:
            orthoGroup2protName[ortho_group] = protein_names_mod

# Singletons information
singletons_filepath = "09_final_groups/Phytophthora_Singletons.txt"

singletons = []

with open(singletons_filepath) as fp:
    for line in fp:
        line = re.search('\|(.+$)', line).group(1)
        singletons.append(line.rstrip('\n'))
        
def find_orthologs(modified_protein_name):
    ret_dict = {}
    for og in orthoGroup2protName:
        for pn in orthoGroup2protName[og]:
            if pn == modified_protein_name:
                ortho_list = orthoGroup2protName[og].copy()
                ortho_list.remove(modified_protein_name)
                ret_dict[og] = ortho_list
    
    if not ret_dict:
        if modified_protein_name not in singletons:
            raise ValueError('Protein_name not included in Orthologs analysis or protein is not include in the secretome.')
                
    return ret_dict

def in_which_species(modified_protein_name):
    res_dic = find_orthologs(modified_protein_name)
    ret_dic = {}
    for og in res_dic:
        species_names = []
        for pn in res_dic[og]:
            species_names.append(re.search('(^.+)_[0-9]+$', pn).group(1))
        species_names =  list(set(species_names))
        ret_dic[og] = species_names
        
    return ret_dic

# Recover all .txt files in current directory
names_files = []
for f in glob.glob('*trCDS_modified.txt'):
    names_files.append(f)
    
# Load .txt files into a panda Dataframe.
ret_df = pd.DataFrame(columns=['Modified_name', 'Original_name', 'Species_prefix'])
for f in names_files:
    df = pd.read_csv(f, sep='\t')
    df['Species_prefix'] =  re.search('(^.+)_trCDS', f).group(1)
    ret_df = pd.concat([ret_df, df])

controls_path = 'control_positiveEffector'

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

with cd(controls_path):
    names_files = []
    for f in glob.glob('*.xls*'):
        names_files.append(f)

controls_df = pd.DataFrame(columns=['Prot_id', 'Prot_seq'])

for f in names_files:
    # print(f)
    df = pd.read_excel(os.path.join(controls_path, f), names = ['Prot_id', 'Prot_seq'])
    controls_df = pd.concat([controls_df, df])

controls_df =  controls_df[controls_df.Prot_seq != 0 ]
controls_df = controls_df.dropna()

def createSeqRecordRow(x):
    # print(x)
    return SeqRecord(Seq(x.Prot_seq, alphabet=), id =x.Prot_id)

controls_df['Seq_record'] = controls_df.apply(createSeqRecordRow, axis=1)

def blast_protSeq(seqRecord):
    print(seqRecord)
    blastp_cl =  NcbiblastpCommandline(db='Phytophthora_species_proteomes', outfmt="6", num_threads=4, num_alignments=2)
    stdout, stderr = blastp_cl(stdin=str(seqRecord.seq))
    stdout_io = StringIO(stdout)
    
    columns = ['qaccver', 'saccver', 'pident ',
           'length ','mismatch ','gapopen ',
           'qstart', 'qend ','sstart', 'send',
           'evalue', 'bitscore']
    
    prot_df = pd.read_csv(stdout_io, sep='\t', header=None, names=columns)
    prot_df['qaccver'] = seqRecord.id
    
    return prot_df

results_series = controls_df.Seq_record.apply(blast_protSeq)

df_total = pd.concat(results_series.tolist())
print(df_total)

df_best_hits = df_total.loc[0]
df_best_100 = df_best_hits[df_best_hits['pident '] == 100]
print(df_best_100)
df_total.to_csv('balst_rxlr_all.csv', header=True, index=False)
