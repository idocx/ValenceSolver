# -*- coding: utf-8 -*-
import json
import os
import random
from copy import deepcopy
import collections
from pprint import pprint

from ValenceSolver.core.utils import to_GeneralMat_obj, get_material_valence

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


if __name__ == "__main__":
    with open('rsc/upload_v10.json', 'r') as fr:
    # with open('rsc/data_release_v12.json', 'r') as fr:
        reactions = json.load(fr)
    reactions = reactions['reactions']

    print('len(reactions)', len(reactions))
    random.shuffle(reactions)
    valence_cache = {}
    for i, reaction in enumerate(reactions):
    # #     target
    # #     generate Material object from entire composition for general usage
    #     tmp_mat_obj = to_GeneralMat_obj(
    #         composition=reaction['target']['composition'],
    #         amounts_vars=reaction['target']['amounts_vars'],
    #         elements_vars=reaction['reaction']['element_substitution']
    #     )
    #     valence = get_material_valence(tmp_mat_obj, valence_cache=valence_cache)
    #     reaction['target']['valence'] = valence

    #     generate Material object from each dict in composition
    #     only use this when we want to generate valence for each dict in composition
        for tmp_comp in reaction['target']['composition']:
            tmp_mat_obj = to_GeneralMat_obj(
                composition=tmp_comp,
                amounts_vars=reaction['target']['amounts_vars'],
                elements_vars=reaction['reaction']['element_substitution']
            )
            valence = get_material_valence(tmp_mat_obj, valence_cache=valence_cache)
            tmp_comp['valence'] = valence


            if tmp_mat_obj is not None:
                print(valence)

            if valence is None and tmp_mat_obj is not None:
                print(valence)
                error_comps = tmp_mat_obj.get_critical_compositions(
                        skip_wrong_composition=True
                )
                print(reaction['doi'])
                pprint(reaction['targets_string'])
                for c in error_comps:
                    print(c.composition)

    #     precursor
        for tmp_pre in reaction['precursors']:
            for tmp_comp in tmp_pre['composition']:
                tmp_mat_obj = to_GeneralMat_obj(
                    composition=tmp_comp,
                    amounts_vars=tmp_pre['amounts_vars'],
                    elements_vars=reaction['reaction']['element_substitution']
                )
                valence = get_material_valence(tmp_mat_obj, valence_cache=valence_cache)
                tmp_comp['valence'] = valence

    print('len(valence_cache)', len(valence_cache))

    if not os.path.exists('generated'):
        os.mkdir('generated')
    with open('generated/upload_v10_w_valence.json', 'w') as fw:
        json.dump(reactions, fw, indent=2)