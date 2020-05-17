# -*- coding: utf-8 -*-

import os
import pylab
from collections import Counter
import matplotlib.pyplot as plt
from Bio import Phylo
from multiprocessing import Pool
TREE = Phylo.read("tree.newick", format="newick")
brutal_coditions = True

class Animal(object):
    def __init__(self, name, gene):
        self.name = name
        self.gene = gene
        self.ins = []
        self.dels = []
        self.holes = []
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
    path = "/mnt/lustre/Blackstone/pipeline_intermediate_data/for_indel_finder/"
    with open(path + fl, "r") as file:
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
                if start == False:
                    st = ind + 1
                start = True
                x+=1
            else:
                if start == True:
                    animal.holes.append(Indel(position=st, length=x, name=None))
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
    holes.sort(key = lambda indl: (indl.length, indl.position))
    return holes


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
        path = TREE.get_path(clade)
        parent = path[-2]
        leafs = parent.get_terminals()
        leafs_names = []
        for leaf in leafs:
            leafs_names.append(leaf.name)
        return False if al.name in leafs_names else True

    strange_examples = []
    while len(holes) > 0:
        hole = min(holes, key=lambda element: element.total_length)
        the_one = False
        pos = 274
        ln = 8
        if hole.position == pos and hole.total_length == ln:
            the_one = True
            ln = ln + len(hole.corrections)
        # if the_one:
        #     print(hole.length, hole.total_length, hole.corrections, hole.position)
        names_with = []
        names_without = []
        if hole.total_length == 0:
            holes.remove(hole)
            continue
        for an in animals:
            if hole in an.holes:
                names_with.append(an.name)
            else:
                stop = False
                for z in sorted(list(hole.key)):
                    if an.gene[z-1] == "-":
                        stop = True
                if not stop:
                    names_without.append(an.name)
        if not names_with:
            holes.remove(hole)
            continue
        # if the_one:
        #     print("ans_with")
        #     for an in names_with:
        #         print(an)
        #     print("____")
        #     print("ans_without")
        #     for an in names_without:
        #         print(an)
        #     print("________")
        def my_monophily(lst1, lst2):
            ca = TREE.common_ancestor(lst1)
            # if the_one:
            #     ca.color = "blue"
            #     Phylo.draw(TREE)
            leafs = ca.get_terminals()
            leafs_names = []
            for lf in leafs:
                leafs_names.append(lf.name)
            # if the_one:
            #     for name in leafs_names:
            #         for an in animals:
            #             if an.name == name:
            #                 print(an.name, an.gene[pos-2:pos-1+ln+1])
            monophyletic = True
            for n in lst2:
                if n in leafs_names:
                    monophyletic = False
            return monophyletic

        # if the_one:
        #     print("in")
        monophyletic_in = my_monophily(names_without, names_with)
        # if the_one:
        #     print("del")
        monophyletic_del = my_monophily(names_with, names_without)
        # if the_one:
        #     print(monophyletic_in, monophyletic_del)
        if monophyletic_del and not monophyletic_in:
            lfs = TREE.common_ancestor(names_with).get_terminals()
            lfs_names = []
            for lf in lfs:
                lfs_names.append(lf.name)
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
                            holes.append(indel_copy(h))
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
    "done"
    animals = deal_with_holes_second_edition(holes, animals, filename)
    return animals


def filtrate(lst, key):
    temp = []
    for element in lst:
        if key(element):
            temp.append(element)
    lst_copy = lst[:]
    for element in temp:
        lst_copy.remove(element)
    return lst_copy


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
    in_dir = "/mnt/lustre/Blackstone/pipeline_intermediate_data/for_indel_finder/"
    out_dir = "output_indel_finder/"
    # For distributions
    with open("strange examples.txt", "w") as fl:
        fl.write("\n")
    with open("uncertain_cases.txt", "w") as fl:
        fl.write("\n")
    files = os.listdir(in_dir)
    i = 0
    for file in files:
        i+=1
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

def indel_finder(file):
    frameshifts_simultaneous = set()
    frameshifts_not_simultaneous = set()
    discarded_animals = []
    number_of_deletions = []
    sum_length_of_deletions = []
    unfiltrated_animals = make_list_of_animals(file)
    if unfiltrated_animals is None:
        with open("crazy files.txt", "a") as fl:
            fl.write(file + "\n")
        return [], []
    # We don't need genes consisting of pure deletion
    unfiltrated_animals = filtrate(unfiltrated_animals, lambda an: an.gene == len(an.gene) * "-")

    unfiltrated_animals = find_indels_in_animals_second_edition(unfiltrated_animals, file)
    for animal in unfiltrated_animals:
        number_of_deletions.append(len(animal.dels))
        sum_length_of_deletions.append(animal_gap_sum(animal))

    # Filters:
    # We don't need genes with gaps at the start or at the end
    animals = filtrate(unfiltrated_animals, lambda an: an.gene[0] == "-" or an.gene[-1] == "-")

    # STOP-codon in the middle of the gene

    animals = filtrate(animals, lambda an: stop_in_gene(an.gene))
    # Nucleotide sum divisible by 3
    animals = filtrate(animals, lambda anl: len(anl.gene.replace("-", "")) % 3 != 0)
    # We don't want to see genes with to many gaps
    # animals = filtrate(animals, lambda anl: animal_gap_sum(anl) > length_threshold)
    # And, also, genes with too many deletions (this is not a usual event), because they may be sequenced badly or
    # something
    # animals = filtrate(animals, lambda anl: len(anl.dels) > number_threshold)
    # Remove animals without ins, dels, or indels
    animals = filtrate(animals, lambda anl: len(anl.dels + anl.ins) == 0)
    # indel_list = make_indel_list(animals)
    file = ".".join(file.split(".")[:-1])
    # make_good_looking_table(indel_list, animals, fl)

    for animal in animals:
        all = animal.ins + animal.dels
        # Getting rid of indels with length divisible by 3 - we are not interested in those.
        temp = []
        for element in all:
            if element.total_length % 3 != 0:
                temp.append(element)
        all = temp
        # Now we'll get the surroundings of our animal. We need that to get rid of their common indels.
        if animal.name != "petMar2":
            path = TREE.get_path(animal.name)
            parent = path[-2]
            leafs = parent.get_terminals()
            leafs_names = []
            for leaf in leafs:
                leafs_names.append(leaf.name)
            neighbours = []
            for element in leafs_names:
                for animal_neigh in unfiltrated_animals:
                    if animal_neigh.name == element:
                        neighbours.append(animal_neigh)
            temp = []
            for indel in all:
                get_rid = True
                for neighbour in neighbours:
                    if indel not in (neighbour.ins + neighbour.dels):
                        get_rid = False
                if indel.total_length < 20:
                    get_rid = False
                if not get_rid:
                    temp.append(indel)
            all = temp
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
                        if left.ca == right.ca:
                            frameshifts_simultaneous.add(new_frameshift)
                            for anml in unfiltrated_animals:
                                all_anml = anml.ins + anml.dels
                                if (all[i] in all_anml) and (all[j] in all_anml):
                                    new_frameshift = (
                                    anml.name, file, str(left), str(right), nonbci(all[i], all[j], anml))
                                    frameshifts_simultaneous.add(new_frameshift)
                        else:
                            frameshifts_not_simultaneous.add(new_frameshift)
                            for anml in unfiltrated_animals:
                                all_anml = anml.ins + anml.dels
                                if (all[i] in all_anml) and (all[j] in all_anml):
                                    new_frameshift = (anml.name, file, str(left), str(right), nonbci(all[i], all[j], anml))
                                    frameshifts_not_simultaneous.add(new_frameshift)
    return list(frameshifts_simultaneous), list(frameshifts_not_simultaneous)

main()
