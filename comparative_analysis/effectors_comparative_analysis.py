import re
import glob
import os
from contextlib import contextmanager
from collections import defaultdict
import csv

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from upsetplot import from_memberships, plot
from matplotlib import pyplot


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

        return (species_groups, species_groups)

    def plot_species_intersections(self, color):
        memberships = []
        data = []

        species_groups, _ = self.orthogroups_sets()

        for k in species_groups:
            memberships.append(k)
            data.append(len(set(species_groups[k])))

        structured_data = from_memberships(memberships, data=data)

        p = plot(structured_data,
                 orientation='vertical',
                 show_counts=True,
                 facecolor=color,
                 element_size=100)

        return p

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

    secretome_filepath = 'CRN_all.fa'
    orthologs_filepath = "09_final_groups/Phytophthora_OrthologousGroups.txt"
    singletons_filepath = "09_final_groups/Phytophthora_Singletons.txt"

    sl = SecretomesLoader(secretome_filepath)


    species2clade = {'P8084_finalAssembly': '1c',
                'P_cactorum_10300': '1a',
                'P_infestans_RefSeq': '1c',
                'P_ramorum_Pr102': '8c',
                'P_palmivora_LILI_trCDS': '4', 
                'P_parasitica_INRA310': '1c',
                'P_sojae_V3': '7b'}

    ortho_s = OrthologsStructure(orthologs_filepath, singletons_filepath, sl.dict, species2clade)

    clade_dict = ortho_s.get_clade_intersections_protIds()
    species_dict = ortho_s.get_species_intersections_protIds()

    write_dict(clade_dict, 'CRN_clade_dict')
    write_dict(species_dict, 'CRN_species_dict')

if __name__ == "__main__":
    main()