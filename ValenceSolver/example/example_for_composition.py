# -*- coding: utf-8 -*-
from pprint import pprint

from pymatgen import Composition

from ValenceSolver.core.composition_inhouse import CompositionInHouse
from ValenceSolver.core.utils import to_GeneralMat_obj, get_material_valence

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


def get_valence_datafile_format():
    composition = {
      "formula": "HCa2Nam-3NbmO3m+1",
      "amount": "1.0",
      "elements": {
        "H": "1.0",
        "Ca": "2.0",
        "Na": "m-3",
        "Nb": "m",
        "O": "3*m+1"
      }
    }
    print(composition)
    tmp_mat_obj = to_GeneralMat_obj(
        composition=composition,
        amounts_vars={
            "m": {
              "values": [],
              "max_value": 6.0,
              "min_value": 3.0
            }
        },
        elements_vars={}
    )
    valence = get_material_valence(tmp_mat_obj)
    pprint(valence)

def get_valence_plain_dict(composition):
    valence_comp = CompositionInHouse(composition)
    valence_comp, inte_factor = valence_comp.get_integer_formula_and_factor()
    valence_comp = CompositionInHouse(valence_comp)
    oxi_state, _, _ = valence_comp.oxi_state_guesses_most_possible()
    print(composition, oxi_state)
    return oxi_state


def get_valence_plain_formula(formula, mode='inhouse'):
    if mode == 'pymatgen':
        valence_comp = Composition(formula)
        valence_comp, inte_factor = valence_comp.get_integer_formula_and_factor()
        valence_comp = Composition(valence_comp)
        oxi_state = valence_comp.oxi_state_guesses()
    else:
        valence_comp = CompositionInHouse(formula)
        valence_comp, inte_factor = valence_comp.get_integer_formula_and_factor()
        valence_comp = CompositionInHouse(valence_comp)
        oxi_state, is_usual, comments = valence_comp.oxi_state_guesses_most_possible()

    print(formula, oxi_state)
    return oxi_state


if __name__ == '__main__':
    composition = {
        "Y": 1.0,
        "Fe": 1.0,
        "O": 3.0,
    }
    get_valence_plain_dict(composition)
    add_zero_valence = True
    get_valence_plain_formula('YFeO3')
    get_valence_plain_formula('Y0.1Fe1.9O3')
    get_valence_plain_formula('Y0.1Fe2O3')
    get_valence_plain_formula('Y0.2Fe2O3')
    get_valence_plain_formula('Y0.5Fe2O3')
    get_valence_plain_formula('Fe1.9Y2.0O3')
    get_valence_plain_formula('YFe2O3')
    get_valence_plain_formula('Eu0.1Y2O3')

