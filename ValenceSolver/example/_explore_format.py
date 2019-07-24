import itertools
import json

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


def explore_reactions():
    with open('rsc/upload_v10.json', 'r') as fr:
        reactions = json.load(fr)
    reactions = reactions['reactions']

    print('len(reactions)', len(reactions))

    valence_cache = {}
    for i, reaction in enumerate(reactions):
        if reaction['reaction']['element_substitution']:
            print(reaction['reaction']['element_substitution'])
        if i > 1000:
            break

def explore_paragraphs():
    with open('rsc/solid_state_data_v11.json', 'r') as fr:
        paras = json.load(fr)
    print('len(paras)', len(paras))
    for i, para in enumerate(paras):
        for chunk in para['extracted_data']:
            for tmp_mat in chunk['target']:
                if tmp_mat['elements_vars']:
                    print(para['doi'])
                print("tmp_mat['elements_vars']", tmp_mat['elements_vars'])
                ele_vars = list(filter(lambda x: len(x[1]) > 0, tmp_mat['elements_vars'].items()))
                for sub in itertools.product(*[tmp_var[1] for tmp_var in ele_vars]):
                    tmp_var_dict = {
                        tmp_var[0]: sub[tmp_index]
                        for tmp_index, tmp_var in enumerate(ele_vars)
                    }
                    print('tmp_var_dict', tmp_var_dict)
                    print('sub', sub)

if __name__ == "__main__":
    explore_paragraphs()
