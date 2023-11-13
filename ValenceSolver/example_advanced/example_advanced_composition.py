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


if __name__ == '__main__':
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition(
            'KClO3',
            all_oxi_states=True,
            # oxi_states_override={
            #     'K': [1,],
            #     'Cl': [-1, 1, 3, 5 7],
            #     'O': [-2,],
            # }
        )
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition(
            'EuBrO2',
            all_oxi_states=True,
            # oxi_states_override={
            #     'Eu': [2,3],
            #     'Br': [-1, 1, 3, 5, 7],
            #     'O': [-2,],
            # }
        )
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('La200Mg100Mn(Ge33O200)3')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Sr63577Mn70641Bi7064O211923')
    )
    print(
        CompositionInHouse.get_most_possible_oxi_state_of_composition('Li118358Mn106522Fe11836P118358O473431')
    )

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
