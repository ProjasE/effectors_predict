import re
import glob
import os
from contextlib import contextmanager
from collections import defaultdict
import csv
import argparse
import yaml

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from upsetplot import from_memberships, plot
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2


class SecretomesLoader(object):

    def __init__(self, secretome_filepath):
        self.secretome_filepath = secretome_filepath
        self.dict = self._build()

    def _build(self):
        secretome_dict = {}

        for record in SeqIO.parse(self.secretome_filepath, "fasta"):
            species_prefix = re.search('(^.+)_[0-9]+', record.id).group(1)
            protein_id = re.search('(^.+_[0-9]+)', record.id).group(1)
            if species_prefix in secretome_dict:
                if protein_id not in secretome_dict[species_prefix]:
                    secretome_dict[species_prefix].append(protein_id)
            else:
                secretome_dict[species_prefix] = [protein_id]

        return secretome_dict


class OrthologsStructure(object):

    def __init__(self, orthologs_filepath, singletones_filepath, secretome_dict, species_clade_dict):
        self.orthologs_filepath = orthologs_filepath
        self.singletons_filepath = singletones_filepath
        self.species_clade_mapping = species_clade_dict

        self.secretome_dict = secretome_dict
        self.orthogroups_dict = self._build_orthogroups_structure()
        self.singletons_list = self._build_singletons_structure()

    def _build_orthogroups_structure(self):
        orthoGroup2protName = {}

        with open(self.orthologs_filepath) as fp:
            for cnt, line in enumerate(fp):
                ortho_group = re.search('^(\S+):', line).group(1)
                protein_names = re.search(': (.+$)', line).group(1).split(' ')

                protein_names_mod = []
                for pn in protein_names:
                    pn_mod = re.search('\|(.+$)', pn).group(1)
                    species_prefix = re.search(
                        '(^.+)_[0-9]+',  pn_mod).group(1)
                    if species_prefix in self.secretome_dict and pn_mod in self.secretome_dict[species_prefix]:
                        protein_names_mod.append(pn_mod)

                if protein_names_mod:
                    orthoGroup2protName[ortho_group] = protein_names_mod

        return orthoGroup2protName

    def _build_singletons_structure(self):
        singletons = []

        with open(self.singletons_filepath) as fp:
            for line in fp:
                line = re.search('\|(.+$)', line).group(1)
                singletons.append(line.rstrip('\n'))

        return singletons

    def find_orthologs(self, modified_protein_name):
        ret_dict = {}

        for og in self.orthogroups_dict:
            for pn in self.orthogroups_dict[og]:
                if pn == modified_protein_name:
                    ortho_list = self.orthogroups_dict[og].copy()
                    ortho_list.remove(modified_protein_name)
                    ret_dict[og] = ortho_list

        if not ret_dict:
            if modified_protein_name not in self.singletons_list:
                raise ValueError('Protein_name not included in Orthologs \
                analysis or protein is not include in the secretome.')

        return ret_dict

    def in_which_species(self, modified_protein_name):
        species = re.search('(^.+)_[0-9]+$', modified_protein_name).group(1)
        res_dict = self.find_orthologs(modified_protein_name)

        ret_dict = {}
        for og in res_dict:
            species_names = []
            for pn in res_dict[og]:
                species_names.append(re.search('(^.+)_[0-9]+$', pn).group(1))
            species_names = list(set(species_names))
            species_names.remove(species)
            ret_dict[og] = species_names

        return ret_dict

    def secretomeProtein_in_which_species(self, secretome_proteinId):

        formated_id = re.search('(^.+_[0-9]+)', secretome_proteinId).group(1)
        return self.in_which_species(formated_id)

    def secretomeProtein_in_which_clade(self, secretome_proteinId):
        res_dict = self.secretomeProtein_in_which_species(secretome_proteinId)
        ret_dict = {}
        for k in res_dict:
            clades = []
            array = res_dict[k]
            for s in array:
                c = self.species_clade_mapping[s]
                clades.append(c)
            clades = list(set(clades))
            ret_dict[k] =  clades

        return ret_dict

    def orthogroups_sets_clades(self):
        clade_groups = defaultdict(lambda: [])
        clades_orthogroups = defaultdict(lambda: [])

        for s in self.secretome_dict:
            try:
                clade = self.species_clade_mapping[s]
            except KeyError:
                pass

            for p_id in self.secretome_dict[s]:
                try:
                    ret_dict = self.secretomeProtein_in_which_clade(p_id)
                    if ret_dict:
                        for k in ret_dict:
                            list_key = ret_dict[k]
                            if list_key:
                                list_key.append(clade)
                                list_key.sort()
                            else:
                                list_key = [clade]
                            list_key = list(set(list_key))
                            tuple_key = tuple(list_key)
                            clade_groups[tuple_key].append(k)
                            clades_orthogroups[clade].append(k)
                    else:
                        tuple_key = tuple([clade])
                        clade_groups[tuple_key].append(p_id)
                        clades_orthogroups[clade].append(p_id)
                except ValueError:
                    pass
        
        return (clade_groups, clades_orthogroups)


    def orthogroups_sets(self):
        species_groups = defaultdict(lambda: [])
        species_orthogroups = defaultdict(lambda: [])

        for s in self.secretome_dict:
            for p_id in self.secretome_dict[s]:
                try:
                    ret_dic = self.secretomeProtein_in_which_species(p_id)
                    if ret_dic:  # there is an orthogroup listed
                        for k in ret_dic:
                            list_key = ret_dic[k]
                            if list_key:  # Orthogroup includes other species
                                list_key.append(s)
                                list_key.sort()
                            else:  # Orthogroup just includes this species
                                list_key = [s]
                            tuple_key = tuple(list_key)
                            species_groups[tuple_key].append(k)
                            species_orthogroups[s].append(k)
                    else:  # No orthogroup listed.  Singleton
                        tuple_key = tuple([s])
                        species_groups[tuple_key].append(p_id)
                        species_orthogroups[s].append(p_id)

                except ValueError:
                    pass

        return (species_groups, species_orthogroups)

    def plot_species_intersections(self, color, ignore_counts=0, orientation='horizontal'):
        memberships = []
        data = []

        species_groups, _ = self.orthogroups_sets()

        for k in species_groups:
            memberships.append(k)
            data.append(len(set(species_groups[k])))

        structured_data = from_memberships(memberships, data=data)

        species_dict = {'P8084_finalAssembly': 'P.betacei',
                        'P_cactorum_10300': 'P.cactorum',
                        'P_infestans_RefSeq': 'P. infestans',
                        'P_palmivora_LILI_trCDS': 'P.palmivora',
                        'P_parasitica_INRA310': 'P.parasitica',
                        'P_ramorum_Pr102': 'P.ramorum',
                        'P_sojae_V3': 'P.sojae'}


        new_names = [species_dict[old_name] for old_name in structured_data.index.names]
        structured_data.index.names = new_names

        structured_data = structured_data[structured_data > ignore_counts].copy()

        p = plot(structured_data,
                 orientation=orientation,
                 show_counts=True,
                 facecolor=color,
                 element_size=40)

        return p

    def plot_venn_diagram_species(self):
        species_groups, _ = self.orthogroups_sets()

        combinations_counts = defaultdict(int)
        for combination in species_groups:
            combinations_counts[combination] = \
                len(set(species_groups[combination]))

        set_names = []
        for combination in combinations_counts:
            if len(combination) == 1:
                set_names.append(combination[0])

        if len(set_names) == 3:
            return venn3(subsets=(combinations_counts[(set_names[0],)],
                               combinations_counts[(set_names[1],)],
                               combinations_counts[(set_names[0], set_names[1])],
                               combinations_counts[(set_names[2],)],
                               combinations_counts[(set_names[0], set_names[2])],
                               combinations_counts[(set_names[1], set_names[2])],
                               combinations_counts[(set_names[0], set_names[1], set_names[2])],
                               ),
                      set_labels=(set_names[0], set_names[1], set_names[2]))

        elif len(set_names) == 2:
            return venn2(subsets=(combinations_counts[(set_names[0],)],
                               combinations_counts[(set_names[1],)],
                               combinations_counts[(set_names[0], set_names[1])]),
                       set_labels=(set_names[0], set_names[1]))
        else:
            return None



    def plot_clades_intersections(self, color):
        memberships = []
        data = []

        clades_groups, _ = self.orthogroups_sets_clades()

        for k in clades_groups:
            memberships.append(k)
            data.append(len(set(clades_groups[k])))

        structured_data = from_memberships(memberships, data=data)

        p = plot(structured_data,
                 orientation='vertical',
                 show_counts=True,
                 facecolor=color,
                 element_size=100)
        return p


    def get_species_intersections_protIds(self):
        species_groups, _ = self.orthogroups_sets()

        species_protIds_dict = {}
        for sg in species_groups:
            p_ids = []
            for og in species_groups[sg]:
                try:
                    p_ids += self.orthogroups_dict[og]
                except KeyError:
                    p_ids += [og]
            species_protIds_dict[sg] = list(set(p_ids))

        return species_protIds_dict

    def get_clade_intersections_protIds(self):
        clades_groups, _ = self.orthogroups_sets_clades()

        clades_protIds_dict = {}
        for cg in clades_groups:
            p_ids = []
            for og in clades_groups[cg]:
                try:
                    p_ids += self.orthogroups_dict[og]
                except KeyError:
                    p_ids += [og]
            clades_protIds_dict[cg] = list(set(p_ids))

        return clades_protIds_dict
                 



def main():

    def write_dict(dict, name):
        w = csv.writer(open(name + ".csv", "w"))
        for key, val in dict.items():
            w.writerow([key, val])

    parser = argparse.ArgumentParser(
        description='Comparative analysis based on Orthologs.')
    parser.add_argument('-p',
                        '--proteins_set_filepath',
                        action="store",
                        dest='proteins_set_filepath',
                        required=True,
                        type=str,
                        help='file with a set of proteins of special interest from the different species to be compared')
    parser.add_argument('-o',
                        '--orthologs_filepath',
                        action="store",
                        dest='orthologs_filepath',
                        required=True,
                        type=str,
                        help='orthologs information file')
    parser.add_argument('-s',
                        '--singletons_filepath',
                        action="store",
                        dest='singletons_filepath',
                        required=True,
                        type=str,
                        help='singletons information file')
    parser.add_argument('-c',
                        '--species2clade',
                        action="store",
                        dest="species2clade",
                        required=True,
                        type=str,
                        help='YAML file mapping species to its respective clade.')
    parser.add_argument('--out_prefix',
                        action="store",
                        dest='out_prefix',
                        required=True,
                        type=str,
                        help='file prefix to use for output files.')
    parser.add_argument('--color',
                        action="store",
                        dest="color",
                        type=str,
                        default='purple',
                        help='color to use in UpSet plots.')
    parser.add_argument('--small_counts',
                        action="store",
                        dest="small_counts",
                        type=int,
                        default=0,
                        help='ignore counts equal or smaller to this parameter'
                             'when rendering UpSet plots.')
    parser.add_argument('--orientation',
                        action="store",
                        dest="orientation",
                        type=str,
                        default='vertical',
                        help='Orientation of UpSet plots.')

    args = parser.parse_args()

    with open(args.species2clade) as f:
        try:
            species2clade = yaml.safe_load(f)
        except yaml.YAMLError as err:
            print(err)

    sl = SecretomesLoader(args.proteins_set_filepath)

    ortho_s = OrthologsStructure(args.orthologs_filepath, 
                                args.singletons_filepath,
                                sl.dict, 
                                species2clade)


    p_s = ortho_s.plot_species_intersections(args.color,
                                             args.small_counts,
                                             args.orientation)
    plt.savefig(args.out_prefix + '_species.png', dpi=300)
    plt.clf()
    p_c = ortho_s.plot_clades_intersections(args.color)
    plt.savefig(args.out_prefix + '_clades.png', dpi=300)
    plt.clf()



    clade_dict = ortho_s.get_clade_intersections_protIds()
    species_dict = ortho_s.get_species_intersections_protIds()

    write_dict(clade_dict, args.out_prefix + '_clade_dict')
    write_dict(species_dict, args.out_prefix + '_species_dict')

    v = ortho_s.plot_venn_diagram_species()
    # if v is not None:
    plt.savefig(args.out_prefix + '_venn.png')
    plt.clf()



if __name__ == "__main__":
    main()
