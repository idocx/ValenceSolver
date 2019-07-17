# -*- coding: utf-8 -*-
import json
import os
from copy import deepcopy
import collections
from pprint import pprint

from ValenceSolver.core.utils import to_GeneralMat_obj, get_target_valence

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


if __name__ == "__main__":
    with open('rsc/upload_v10.json', 'r') as fr:
        reactions = json.load(fr)
    reactions = reactions['reactions']

    print('len(reactions)', len(reactions))

    valence_cache = {}
    for i in range(len(reactions)):
    #         get reaction
        reaction = reactions[i]

    #     target
    #     generate Material object from entire composition for general usage
    #     tmp_tar_obj = to_GeneralMat_obj(
    #         reaction['target'], 
    #         reaction['reaction']['element_substitution']
    #     )

    #     generate Material object from each dict in composition 
    #     only use this when we want to generate valence for each dict in composition
        for tmp_comp in reaction['target']['composition']:
            tmp_mat = deepcopy(reaction['target'])
            tmp_mat['composition'] = [deepcopy(tmp_comp)]
            tmp_mat['composition'][0]['amount'] = '1.0'
            tmp_mat_obj = to_GeneralMat_obj(
                tmp_mat, 
                reaction['reaction']['element_substitution']
            )
            valence = get_target_valence(tmp_mat_obj, valence_cache=valence_cache)
            if len(valence) == 1: 
                if set(valence[0]['valence'].keys()) == set(tmp_mat_obj.composition.keys()):
                    tmp_comp['valence'] = valence[0]['valence']
                elif isinstance(valence[0]['amounts_vars'], dict):
                    # valence[0]['amounts_vars'] could be either a dict or a list of dict
                    tmp_comp['valence'] = valence
                else:
                    # valence[0]['amounts_vars'] is a list
                    tmp_comp['valence'] = []
                    for i, var in enumerate(valence[0]['amounts_vars']):
                        tmp_comp['valence'].append({
                            'valence': valence[0]['valence'],
                            'amounts_vars': var,
                            'elements': valence[0]['elements'][i],
                        })
            elif len(valence) == 0:
                tmp_comp['valence'] = None
            else:
                tmp_comp['valence'] = valence

    #     precursor
        for tmp_pre in reaction['precursors']:
            for tmp_comp in tmp_pre['composition']:
                tmp_mat = deepcopy(tmp_pre)
                tmp_mat['composition'] = [deepcopy(tmp_comp)]
                tmp_mat['composition'][0]['amount'] = '1.0'
                tmp_mat_obj = to_GeneralMat_obj(
                    tmp_mat, 
                    reaction['reaction']['element_substitution']
                )
                valence = get_target_valence(tmp_mat_obj, valence_cache=valence_cache)
                if len(valence) == 1: 
                    if set(valence[0]['valence'].keys()) == set(tmp_mat_obj.composition.keys()):
                        tmp_comp['valence'] = valence[0]['valence']
                    elif isinstance(valence[0]['amounts_vars'], dict):
                        # valence[0]['amounts_vars'] could be either a dict or a list of dict
                        tmp_comp['valence'] = valence
                    else:
                        # valence[0]['amounts_vars'] is a list
                        tmp_comp['valence'] = []
                        for i, var in enumerate(valence[0]['amounts_vars']):
                            tmp_comp['valence'].append({
                                'valence': valence[0]['valence'],
                                'amounts_vars': var,
                                'elements': valence[0]['elements'][i],
                            })
                elif len(valence) == 0:
                    tmp_comp['valence'] = None
                else:
                    tmp_comp['valence'] = valence  
                    
    if not os.path.exists('generated'):
        os.mkdir('generated')
    with open('generated/upload_v10_w_valence.json', 'w') as fw:
        json.dump(reactions, fw, indent=2)