# -*- coding: utf-8 -*-
import itertools
import json
import os
from copy import deepcopy
import collections
from pprint import pprint

from ValenceSolver.core.composition_inhouse import CompositionInHouse
from ValenceSolver.core.utils import to_GeneralMat_obj, get_material_valence

from pymatgen.core.periodic_table import Element, Specie


def print_ele_oxi_states(ele):
    tmp_ele = Element(ele)
    print('tmp_ele.icsd_oxidation_states', tmp_ele.icsd_oxidation_states)
    print('Element(ele).oxidation_states', tmp_ele.oxidation_states)

if __name__ == '__main__':
    for ele in Element.__members__:
        print(type(ele), ele)
        print_ele_oxi_states(ele)
        print()

    exit()

    with open('generated/test.json', 'r') as fr:
        solid_state_data = json.load(fr)

    valency_dictionary = {}
    valence_cache = {}
    fails = []

    warning_cases = [191, 823, 1344, 3166, 3685, 4837, 5288, 5502,    ]
    for j, data in enumerate(solid_state_data[55:56] ):
        for chunk in data['synthesis_chunks']:

            for material in chunk['extracted_materials']['targets'] + chunk['extracted_materials']['precursors']:
                try:
                    composition = material['composition']
                    for compound in composition:
                        # if compound['formula'] == 'Ba1.6Ca0.4P2O7':
                        #     print(j)
                        # else:
                        #     continue
                        if compound['formula'] in valency_dictionary:
                            compound['valence'] = valency_dictionary[compound['formula']]
                        else:
                            all_valence = []
                            ele_vars = list(filter(lambda x: len(x[1]) > 0, material['elements_vars'].items()))
                            for sub in itertools.product(*[tmp_var[1] for tmp_var in ele_vars]):
                                element_sub = {
                                    tmp_var[0]: sub[tmp_index]
                                    for tmp_index, tmp_var in enumerate(ele_vars)
                                }
                                # print(compound)
                                # print(material['amounts_vars'])
                                # print(element_sub)

                                element_sub = {}
                                material_obj = to_GeneralMat_obj(
                                    composition=compound,
                                    amounts_vars=material['amounts_vars'],
                                    elements_vars=element_sub,
                                )

                                valence = get_material_valence(material_obj, valence_cache=valence_cache)
                                if valence:
                                    for v in valence:
                                        del v['elements']
                                        v['elements_vars'] = element_sub
                                    all_valence.extend(valence)
                                print(material_obj)
                                print(valence)
                                print()
                            if not len(all_valence):
                                all_valence = None
                            compound['valence'] = all_valence

                            if material['amounts_vars'] == {} and material['elements_vars'] == {} and all_valence:
                                valency_dictionary[compound['formula']] = all_valence

                except:
                    fails.append(dict(
                        idx = j,
                        stage = 'valency'
                    ))