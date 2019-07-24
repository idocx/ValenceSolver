# -*- coding: utf-8 -*-
import itertools
import json
import os
from copy import deepcopy
import collections
from pprint import pprint

from ValenceSolver.core.utils import to_GeneralMat_obj, get_material_valence

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


if __name__ == "__main__":
    with open('rsc/solid_state_data_v11.json', 'r') as fr:
        paras = json.load(fr)
    print('len(paras)', len(paras))

    valence_cache = {}

    for i, para in enumerate(paras):
        for chunk in para['extracted_data']:
            for tmp_mat in chunk['target'] + chunk['precursors'] + chunk['other_materials']:
                for tmp_comp in tmp_mat['composition']:
                    all_valence = []

                    ele_vars = list(filter(lambda x: len(x[1]) > 0, tmp_mat['elements_vars'].items()))
                    for sub in itertools.product(*[tmp_var[1] for tmp_var in ele_vars]):
                        element_sub = {
                            tmp_var[0]: sub[tmp_index]
                            for tmp_index, tmp_var in enumerate(ele_vars)
                        }
                        tmp_mat_obj = to_GeneralMat_obj(
                            composition=tmp_comp,
                            amounts_vars=tmp_mat['amounts_vars'],
                            elements_vars=element_sub,
                        )
                        valence = get_material_valence(tmp_mat_obj, valence_cache=valence_cache)
                        if valence:
                            for v in valence:
                                del v['elements']
                                v['elements_vars'] = element_sub
                            all_valence.extend(valence)
                    if not len(all_valence):
                        all_valence = None
                    tmp_comp['valence'] = all_valence

    if not os.path.exists('generated'):
        os.mkdir('generated')
    with open('generated/solid_state_data_v11_w_valence.json', 'w') as fw:
        json.dump(paras, fw, indent=2)