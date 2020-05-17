# -*- coding: utf-8 -*-

import os
import pylab
from collections import Counter
import matplotlib.pyplot as plt
from Bio import Phylo
from multiprocessing import Pool
import time

TREE = Phylo.read("insect_tree.newick", format="newick")
# I need to get rid of two branches since I haven't found normal genomes for them.
for branch in TREE.get_terminals():
    if branch.name == "Dpseudoobscura1" or branch.name == "Agambiae1":
        TREE.collapse(branch)
brutal_coditions = True

class Animal(object):
    def __init__(self, name, gene):
        self.name = name
        self.gene = gene
        self.ins = []
        self.dels = []
        self.holes = set()
        self.strange_things = []
        self.ghost_dels = []
    def __str__(self):
        return self.name + "\t" + self.gene

class Indel(object):
    def __init__(self, position, length, name):
        self.position = position
        self.length = length
        self.name = name
        self.corrections = []
        self.total_length = self.length - len(self.corrections)
        self.confidence = {}
        self.key = set([self.position + i for i in range(self.length)]) - set(self.corrections)
        self.ca = None
    def __str__(self):
        if self.key:
            return str(self.name) + "_" + str(min(self.key)) + "_" + str(self.total_length)
        else:
            return str(self.name) + "_" + str(0) + "_" + str(self.total_length)

    def __eq__(self, other):
        if self.name == other.name and self.key == other.key:
            return True
        return False

    def __hash__(self):
        return hash(str(sorted(list(self.key))))

    def correct(self, num, conf):
        self.corrections.append(num)
        self.confidence[num] = conf
        self.total_length = self.length - len(self.corrections)
        self.key = set([self.position + i for i in range(self.length)]) - set(self.corrections)

    def elongate(self, n):
        self.length += n
        self.total_length = self.length - len(self.corrections)
        self.key = set([self.position + i for i in range(self.length)]) - set(self.corrections)

def indel_copy(indel):
    new_indel = Indel(position=indel.position, length=indel.length, name=indel.name)
    if indel.corrections:
        for correction in indel.corrections:
            confidence = indel.confidence[correction]
            new_indel.correct(correction, confidence)
    return new_indel

def nonbci(indel_1, indel_2, animal):
    left = min(indel_1, indel_2, key=lambda indel: indel.position)
    right = max(indel_1, indel_2, key=lambda indel: indel.position)
    primary_length = len(animal.gene[left.position+left.length-1 : right.position-1].replace("-", ""))
    all = animal.ins + animal.dels
    additional_length = 0
    for element in all:
        if left.position < element.position < right.position:
            if element.ca not in TREE.get_path(right.ca):
                if element.name == "In":
                    additional_length -= element.total_length
                elif element.name == "Del":
                    additional_length += element.total_length
                else:
                    raise NameError("Your indel must be either In or Del")
    for element in animal.strange_things:
        if left.position < element.position < right.position:
            additional_length += element.total_length
    for element in animal.ghost_dels:
        if left.position < element.position < right.position and element.ca in TREE.get_path(right.ca):
            additional_length -= element.total_length
    final_length = primary_length + additional_length
    return final_length

def make_list_of_animals(fl):
    """
    :param fl: File with names of animals and their aligned genes (same gene from different animals)
    :return: list of Animal objects (with names and genes given)
    """
    # if computer == "mine":
    # path = "input/test/No_in/"
    # elif computer == "makarich":
    # path = "/mnt/lustre/Blackstone/insect_pipeline_intermediate_data/for_indel_finder/"
    # path = "in/"
    with open(fl, "r") as file:
        lines = list(file.readlines())
        lines = lines[2:]
        animals = []
        for line in lines:
            if "=" in line:
                return None
            ln = line.split()
            name = ln[0].strip()
            gene = ln[1].strip()
            animals.append(Animal(name, gene))
    return animals


def make_holes_list(animals):
    for animal in animals:
        x=0
        start = False
        for ind, letter in enumerate(animal.gene):
            if letter == "-":
                if not start:
                    st = ind + 1
                start = True
                x+=1
            else:
                if start == True:
                    animal.holes.add(Indel(position=st, length=x, name=None))
                    start = False
                x=0
    holes = list()
    for animal in animals:
        for indel in animal.holes:
            if indel not in holes:
                holes.append(indel_copy(indel))
    holes = sorted(list(holes), key=lambda indl: (indl.position, indl.length))
    stop = False
    switch = False
    while not stop:
        if not holes:
            break
        if not switch:
            if holes[0].position == 1:
                holes.remove(holes[0])
            else:
                switch = True
        else:
            while not stop:
                broken = False
                for element in holes:
                    if element.position + element.length - 1 == len(animals[0].gene):
                        holes.remove(element)
                        broken = True
                        break
                if not broken:
                    stop = True
    return set(holes)


def deal_with_holes_second_edition(holes, animals, filename):
    def indel_is_subindel(some_indel, probable_subindel):
        """
        :param some_indel: some indel
        :param probable_subindel: indel, which is tested on being subindel
        :return: True/False
        """
        return True if probable_subindel.key.issubset(some_indel.key) else False

    def diverged_earlier(clade, al):
        """
        :param clade: some clade
        :param al: animal which is tested on earlier divergence
        :return: True/False
        """
        parent = TREE.get_path(clade)[-2]
        leafs_names = {element.name for element in parent.get_terminals()}
        return False if al.name in leafs_names else True
    strange_examples = []
    while len(holes) > 0:
        hole = min(holes, key=lambda element: element.total_length)
        if hole.total_length == 0:
            holes.remove(hole)
            continue
        names_with = set()
        names_without = set()
        for an in animals:
            if hole in an.holes:
                names_with.add(an.name)
            else:
                # I only consider animal not to have hole if al the symbols on this range are letters, not "-"
                no_hole = True
                for z in sorted(list(hole.key)):
                    if an.gene[z-1] == "-":
                        no_hole = False
                if no_hole:
                    names_without.add(an.name)
        if not names_with:
            holes.remove(hole)
            continue
        def my_monophily(set1, set2):
            ca = TREE.common_ancestor(set1)
            leafs_names = {element.name for element in ca.get_terminals()}
            if leafs_names.intersection(set2):
                    return False
            return True

        monophyletic_in = my_monophily(names_without, names_with)
        monophyletic_del = my_monophily(names_with, names_without)
        if monophyletic_del and not monophyletic_in:
            lfs_names = {element.name for element in TREE.common_ancestor(names_with).get_terminals()}
            for animal in animals:
                if animal.name in names_with:
                    new_ind = indel_copy(hole)
                    new_ind.name = "Del"
                    new_ind.ca = TREE.common_ancestor(names_with)
                    animal.dels.append(new_ind)
                elif animal.name in lfs_names:
                    new_ind = indel_copy(hole)
                    new_ind.name = "Del"
                    new_ind.ca = TREE.common_ancestor(names_with)
                    animal.ghost_dels.append(new_ind)
        if monophyletic_in and not monophyletic_del:
            for an in animals:
                if an.name in names_without:
                    new_ind = indel_copy(hole)
                    new_ind.name = "In"
                    new_ind.ca = TREE.common_ancestor(names_without)
                    an.ins.append(new_ind)
                    continue
                if diverged_earlier(TREE.common_ancestor(names_without), an):
                    confident = True
                else:
                    confident = False
                ### Этот кусочек кода не тестировался. Предполагается, что он избавит от редактирования гэпы тех зверей, которые оказались внутри клады
                # зверей с инсерцией. Типа, у них после этого могла произойти делеция, и это не ложные гэпы.
                inside_names_without = TREE.common_ancestor(names_without).get_terminals()
                leafs_names = []
                for lf in inside_names_without:
                    leafs_names.append(lf.name)
                if an.name in inside_names_without:
                    continue
                ### Здесь этот кусочек кода заканчивается
                for h in an.holes:
                    if indel_is_subindel(h, hole):
                        for position in sorted(list(hole.key)):
                            if not confident:
                                with open("output_indel_finder/slow_uncertain_cases.txt", "a") as fl:
                                    text = filename + " " + str(hole) + " " + an.name + " " + str(position) + "\n"
                                    fl.write(text)
                            h.correct(position, confident)
                        h.corrections = sorted(h.corrections)
                        if h not in holes and h.total_length > 0:
                            holes.add(indel_copy(h))
        if (monophyletic_in and monophyletic_del) or (not monophyletic_in and not monophyletic_del):
            strange_examples.append(str(hole))
            for animal in animals:
                if hole in animal.holes and hole not in animal.strange_things:
                    animal.strange_things.append(indel_copy(hole))
        holes.remove(hole)
    with open("output_indel_finder/slow_strange examples.txt", "a") as fl:
        for ex in strange_examples:
            text = filename + " " + str(ex) + "\n"
            fl.write(text)
    return animals



def find_indels_in_animals_second_edition(animals, filename):
    holes = make_holes_list(animals)
    animals = deal_with_holes_second_edition(holes, animals, filename)
    return animals


def make_indel_list(animals):
    indel_list = []
    for animal in animals:
        for indel in animal.ins + animal.dels:
            if indel not in indel_list:
                indel_list.append(indel_copy(indel))
    x = lambda indl: (indl.position, indl.total_length)
    indel_list.sort(key=x)
    return indel_list


def main():
    in_dir = "/mnt/lustre/Blackstone/insect_pipeline_intermediate_data/for_indel_finder/"
    out_dir = "output_indel_finder/"
    # For distributions
    with open("strange examples.txt", "w") as fl:
        fl.write("\n")
    with open("uncertain_cases.txt", "w") as fl:
        fl.write("\n")
    files = os.listdir(in_dir)
    i = 0
    for file in files:
        i += 1
        with open("slow_process.txt", "w") as fl:
            fl.write("%s\n" % file)
            fl.write("%i/%i" % (i, len(files)))
        task = indel_finder(file)
        if task[0]:
            with open(out_dir + "slow_simultaneous.txt", "a") as fl:
                for el in task[0]:
                    fl.write("%20s" % str(el))
                fl.write("\n")
        if task[1]:
            with open(out_dir + "slow_not_simultaneous.txt", "a") as fl:
                for el in task[1]:
                    fl.write("%20s" % str(el))
                fl.write("\n")

def animal_gap_sum(a):
    s = 0
    for delet in a.dels:
        s += delet.total_length
    return s


def stop_in_gene(gene):
    copy = gene[:]
    copy = copy.replace("-", "")
    copy = copy[:-3]
    for q in range(0, len(copy), 3):
        one = True if copy[q:q + 3].upper() == "TAG" else False
        two = True if copy[q:q + 3].upper() == "TAA" else False
        three = True if copy[q:q + 3].upper() == "TGA" else False
        if one or two or three:
            return True
    return False

def flatten(lst):
    """
    Removes all the collections from the list, making their elements parts of the "global" list. Yes, it's an awful
    description.
    :param lst: list - a list you want to flatten. Boy, this description is astonishing.
    :return: list - a list without collections in it.
    """
    result = []
    contains_list = False
    for element in lst:
        if type(element) in (list, tuple):  # To be honest, this looks like cheat since we haven't learned it yet.
            # If there's more elegant way, tell me, please
            result += element
            contains_list = True
        else:
            result.append(element)
    if contains_list:
        return flatten(result)
    else:
        return result

def indel_finder(file):
    frameshifts_simultaneous = set()
    frameshifts_not_simultaneous = set()
    unfiltrated_animals = make_list_of_animals(file)
    # We don't need genes consisting of pure deletion
    unfiltrated_animals = list(filter(lambda an: an.gene != len(an.gene) * "-", unfiltrated_animals))
    ### THE MAIN WORK IS DONE HERE ###
    unfiltrated_animals = list(find_indels_in_animals_second_edition(unfiltrated_animals, file))
    # Filters:
    # We don't need genes with gaps at the start or at the end
    animals = tuple(filter(lambda an: an.gene[0] != "-" and an.gene[-1] != "-", unfiltrated_animals))
    # STOP-codon in the middle of the gene
    animals = tuple(filter(lambda an: not stop_in_gene(an.gene), animals))
    # Nucleotide sum divisible by 3
    animals = tuple(filter(lambda anl: len(anl.gene.replace("-", "")) % 3 == 0, animals))
    # Remove animals without ins, dels, or indels
    animals = tuple(filter(lambda anl: len(anl.dels + anl.ins) != 0, animals))
    file = ".".join(file.split(".")[:-1])
    for animal in animals:
        print(animal.name)
        all = animal.ins + animal.dels
        # Getting rid of indels with length divisible by 3 - we are not interested in those.
        all = list(filter(lambda el: el.total_length % 3 != 0, all))
        # Now we'll get the surroundings of our animal. We need that to get rid of their common indels.
        if animal.name != "Aalbimanus":
            path = TREE.get_path(animal.name)
            parent = path[-2]
            leafs_names = [leaf.name for leaf in parent.get_terminals()]
            neighbours = list(filter(lambda anl: anl.name in leafs_names, unfiltrated_animals))
            indels_of_neighbours = [neighbour.ins + neighbour.dels for neighbour in neighbours]
            indels_of_neighbours = flatten(indels_of_neighbours)
            all = list(filter(lambda indel: indel not in indels_of_neighbours or indel.total_length < 20, all))
        # Now if there are more than 2 indels, we won't proceed - I don't want to sort out which of them compensate \
        # which
        if brutal_coditions and len(all) != 2:
            continue
        for i in range(len(all)):
            for j in range(len(all)):
                left = min(all[i], all[j], key=lambda indel: indel.position)
                right = max(all[i], all[j], key=lambda indel: indel.position)
                if i != j and left.total_length % 3 != 0 and right.total_length % 3 != 0:
                    if (left.name == right.name and (left.total_length + right.total_length) % 3 == 0) \
                            or (left.name != right.name and (left.total_length - right.total_length) % 3 == 0):
                        new_frameshift = (animal.name, file, str(left), str(right), nonbci(all[i], all[j], animal))
                        target_set = frameshifts_simultaneous if left.ca == right.ca else frameshifts_not_simultaneous
                        target_set.add(new_frameshift)
                        for anml in unfiltrated_animals:
                            all_anml = anml.ins + anml.dels
                            if (all[i] in all_anml) and (all[j] in all_anml):
                                new_frameshift = (anml.name, file, str(left), str(right), nonbci(all[i], all[j], anml))
                                target_set.add(new_frameshift)
    return list(frameshifts_simultaneous), list(frameshifts_not_simultaneous)

main()
