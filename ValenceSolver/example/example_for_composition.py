# -*- coding: utf-8 -*-
from pprint import pprint

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
    oxi_state = valence_comp.oxi_state_guesses_most_possible(all_oxi_states=False)
    print(composition, oxi_state)
    return oxi_state

def get_valence_plain_formula(formula, add_zero_valence=False):
    valence_comp = CompositionInHouse(formula)
    valence_comp, inte_factor = valence_comp.get_integer_formula_and_factor()
    valence_comp = CompositionInHouse(valence_comp)
    oxi_state = valence_comp.oxi_state_guesses_most_possible(
        all_oxi_states=False,
        add_zero_valence=add_zero_valence
    )
    print(formula, oxi_state)
    return oxi_state


if __name__ == '__main__':
    composition = {
        "Y": 1.0,
        "Fe": 1.0,
        "O": 3.0,
    }
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition(composition)
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('YFeO3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Y0.1Fe1.9O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Y0.1Fe2O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Y0.2Fe2O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Y0.5Fe2O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Fe1.9Y2.0O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('YFe2O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Eu0.1Y2O3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Nb(PO4)3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('NaNb0.667Cr0.667(PO4)3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Na2Nb0.333Cr1.333(PO4)3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Na2.5Nb0.1667Cr1.6667(PO4)3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Na3Cr2(PO4)3')
    )
