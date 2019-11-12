# -*- coding: utf-8 -*-
import copy
import unittest
from pprint import pprint

from ValenceSolver.core.composition_inhouse import CompositionInHouse
from ValenceSolver.core.utils import get_valence_single_composition

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


"""
positive cases
10.1016/j.jpowsour.2008.04.046
['Ba0.5Sr0.5Fe1-yCoyO3']
{'Ba': 0.5, 'Sr': 0.5, 'Fe': 1.0, 'O': 3.0}

{'Ba': 2, 'Sr': 2, 'Fe': 4, 'O': -2}

10.1021/cm011292l
SrFeO3

{'Sr': 2, 'Fe': 4, 'O': -2}

SrFeO2.5

{'Sr': 2, 'Fe': 3, 'O': -2}

https://arxiv.org/pdf/cond-mat/0512510.pdf
LaMnO3

{'La': 3, 'Mn': 3, 'O': -2}

10.1016/S0167-2738(01)00978-X
['Sr0.7La0.3FeO3']
{'Sr': 0.7, 'La': 0.3, 'Fe': 1.0, 'O': 3.0}
{'Sr': 2, 'La': 3, 'Fe': 3.7, 'O': -2}

{'Sr': 2, 'La': 3, 'Fe': 3.7, 'O': -2}

10.1039/C7TC00112F
['Cr1.9Mn0.1Al', 'Cr1.8Mn0.2Al', 'Cr1.7Mn0.3Al']
{'Cr': 1.9, 'Mn': 0.1, 'Al': 1.0}
{'Cr': 1.8, 'Mn': 0.2, 'Al': 1.0}
{'Cr': 1.7, 'Mn': 0.3, 'Al': 1.0}

{'Cr': 0, 'Mn': 0, 'Al': 0}

10.1016/S0921-4534(00)00269-0
['Sr2Ca1-xGdxCu2Bi2.4Oy', 'Sr2Ca1-xGdxCu2Bi2.4Oy']
{'Bi': 2.4, 'Sr': 2.0, 'Gd': 1.0, 'Cu': 2.0}
{'Bi': 2.4, 'Sr': 2.0, 'Ca': 1.0, 'Cu': 2.0}

{'Bi': 0, 'Sr': 0, 'Ca': 0, 'Cu': 0}

10.1016/j.matlet.2014.04.029
['Ca2.97Dy0.03Si2O7',
 'Ca2.95Dy0.05Si2O7',
 'Ca2.93Dy0.07Si2O7',
 'Ca2.9Dy0.1Si2O7']
{'Ca': 2.97, 'Si': 2.0, 'O': 7.0, 'Dy': 0.03}
{'Ca': 2.95, 'Si': 2.0, 'O': 7.0, 'Dy': 0.05}
{'Ca': 2.93, 'Si': 2.0, 'O': 7.0, 'Dy': 0.07}
{'Ca': 2.9, 'Si': 2.0, 'O': 7.0, 'Dy': 0.1}

{'Ca': 2, 'Si': 4, 'O': -2, 'Dy': 2}

10.1016/j.jpcs.2008.12.006
['Sr2Ca2-xYxCu3Pb0.3Bi1.7Oy']
{'Bi': 1.7, 'Pb': 0.3, 'Sr': 2.0, 'Ca': 2.0, 'Cu': 3.0}

{'Bi': 0, 'Pb': 0, 'Sr': 0, 'Ca': 0, 'Cu': 0}


10.1016/j.matchemphys.2010.10.027
['Ba0.95La0.033Zr0.09Ti0.91O3']
{'Ba': 0.95, 'Zr': 0.09, 'Ti': 0.91, 'O': 3.0, 'La': 0.033}

{'Ba': 2, 'Zr': 4, 'Ti': 4, 'O': -2, 'La': 3}

10.1007/s10853-008-2897-2
['Ti2NiC']
{'Ni': 1.0, 'Ti': 1.0}
{'Ni': 0, 'Ti': 0}

negative cases
10.1016/S0167-577X(02)00433-0
['BaY2CuO45']
{'Y': 2.0, 'Ba': 1.0, 'Cu': 1.0, 'O': 45.0}

10.1016/j.physc.2013.02.013
['Ba2Ca2Cu3.5O10',
 'Ba2Ca2Ti0.25Cu3.25O10',
 'Ba2Ca2Ti0.5Cu3O10',
 'Ba2Ca2Ti0.75Cu2.75O10']
{'Cu': 3.5, 'Ba': 2.0, 'Ca': 2.0, 'O': 10.0}
{'Cu': 3.25, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.25, 'O': 10.0}
{'Cu': 3.0, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.5, 'O': 10.0}
{'Cu': 2.75, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.75, 'O': 10.0}

questionable
???
10.1016/j.jmmm.2011.06.030
['BaTixFe11.8-2*xCoxO19']
{'Ba': 1.0, 'Fe': 11.8, 'O': 19.0}
x=0.3, (b) x=0.6, (c) x=0.9 and (d) x=1.2


10.1016/j.matlet.2014.04.029
['Ca2.975Si2O7']
{'Ca': 2.975, 'Si': 2.0, 'O': 7.0}
{'Ca': 2, 'Si': 4, 'O': -2}
or {'Ca': 2.xx, 'Si': 4, 'O': -2} ?

10.1021/ic501510u
['BaHo0.1Ce0.9O3']
{'Ba': 1.0, 'Ce': 0.9, 'O': 3.0, 'Ho': 0.1}
{'Ba': 2, 'Ce': 4, 'O': -2, 'Ho': 3 or 4}??? 

10.1016/j.ssi.2004.10.015
['K0.03La1.97Mo2O9',
 'K0.05La1.95Mo2O9',
 'K0.075La1.925Mo2O9',
 'K0.1La1.9Mo2O9',
 'K0.15La1.85Mo2O9']
{'La': 1.97, 'K': 0.03, 'Mo': 2.0, 'O': 9.0}
{'La': 1.95, 'K': 0.05, 'Mo': 2.0, 'O': 9.0}
{'La': 1.925, 'K': 0.075, 'Mo': 2.0, 'O': 9.0}
{'La': 1.9, 'K': 0.1, 'Mo': 2.0, 'O': 9.0}
{'La': 1.85, 'K': 0.15, 'Mo': 2.0, 'O': 9.0}

{'La': 3, 'K': 1 or 3, 'Mo': 6, 'O': -2}

10.1149/1.2988135
['In0.1Sn0.9P2O7']
{'Sn': 0.9, 'In': 0.1, 'P': 2.0, 'O': 7.0}

{'Sn': 4, 'In': 3 or 4, 'P': 5, 'O': -2}

10.1016/j.jpowsour.2017.10.091
['Ba3Ca1.18Nb1.77Ni0.05O9']
{'Ba': 3.0, 'Ca': 1.18, 'Nb': 1.77, 'Ni': 0.05, 'O': 9.0}

{'Ba': 2, 'Ca': 2, 'Nb': 5, 'Ni': ??, 'O': -2}

10.1016/j.ijhydene.2011.10.020
['Ca0.01La0.99Ti0.01Nb0.99O4']
{'La': 0.99, 'Ca': 0.01, 'Nb': 0.99, 'Ti': 0.01, 'O': 4.0}

{'La': 3, 'Ca': 2 or 3, 'Nb': 5, 'Ti': 4 or 5, 'O': -2}

10.1016/j.memsci.2018.01.068
['La5.5Nb0.15W0.45Mo0.4O11.25-δ']
{'La': 5.5, 'W': 0.45, 'Nb': 0.15, 'Mo': 0.4, 'O': 11.25}

{'La': 3, 'W': 6, 'Nb': 5, 'Mo': 6, 'O': -2} deficiency?

"""

class ValenceTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.positive_cases = [
            {
                'doi': '10.1016/j.jpowsour.2008.04.046',
                'composition': {'Ba': 0.5, 'Sr': 0.5, 'Fe': 1.0, 'O': 3.0},
                'valence': {'Ba': 2, 'Sr': 2, 'Fe': 4, 'O': -2},
            },
            {
                'doi': '10.1021/cm011292l',
                'composition': 'SrFeO3',
                'valence': {'Sr': 2, 'Fe': 4, 'O': -2},
            },
            {
                'doi': '10.1021/cm011292l',
                'composition': 'SrFeO2.5',
                'valence': {'Sr': 2, 'Fe': 3, 'O': -2},
            },
            {
                'doi': 'https://arxiv.org/pdf/cond-mat/0512510.pdf',
                'composition': 'LaMnO3',
                'valence': {'La': 3, 'Mn': 3, 'O': -2},
            },
            {
                'doi': '10.1016/S0167-2738(01)00978-X',
                'composition': {'Sr': 0.7, 'La': 0.3, 'Fe': 1.0, 'O': 3.0},
                'valence': {'Sr': 2, 'La': 3, 'Fe': 3.7, 'O': -2},
            },
            {
                'doi': '10.1016/j.matlet.2014.04.029',
                'composition': {'Ca': 2.97, 'Si': 2.0, 'O': 7.0, 'Dy': 0.03},
                'valence': {'Ca': 2, 'Si': 4, 'O': -2, 'Dy': 2},
            },
            {
                'doi': '10.1016/j.matlet.2014.04.029',
                'composition': {'Ca': 2.95, 'Si': 2.0, 'O': 7.0, 'Dy': 0.05},
                'valence': {'Ca': 2, 'Si': 4, 'O': -2, 'Dy': 2},
            },
            {
                'doi': '10.1016/j.matlet.2014.04.029',
                'composition': {'Ca': 2.93, 'Si': 2.0, 'O': 7.0, 'Dy': 0.07},
                'valence': {'Ca': 2, 'Si': 4, 'O': -2, 'Dy': 2},
            },
            {
                'doi': '10.1016/j.matlet.2014.04.029',
                'composition': {'Ca': 2.9, 'Si': 2.0, 'O': 7.0, 'Dy': 0.1},
                'valence': {'Ca': 2, 'Si': 4, 'O': -2, 'Dy': 2},
            },
            {
                'doi': '10.1016/j.matchemphys.2010.10.027',
                'composition': {'Ba': 0.95, 'Zr': 0.09, 'Ti': 0.91, 'O': 3.0, 'La': 0.033},
                'valence': {'Ba': 2, 'Zr': 4, 'Ti': 4, 'O': -2, 'La': 3},
            },
            {
                'doi': '10.1016/j.physc.2013.02.013',
                'composition': {'Cu': 3.5, 'Ba': 2.0, 'Ca': 2.0, 'O': 10.0},
                'valence': {'Cu': 12 / 3.5, 'Ba': 2.0, 'Ca': 2.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.physc.2013.02.013',
                'composition': {'Cu': 3.25, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.25, 'O': 10.0},
                'valence': {'Cu': 11 / 3.25, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.physc.2013.02.013',
                'composition': {'Cu': 3.0, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.5, 'O': 10.0},
                'valence': {'Cu': 10 / 3.0, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.physc.2013.02.013',
                'composition': {'Cu': 2.75, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.75, 'O': 10.0},
                'valence': {'Cu': 9 / 2.75, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.jmmm.2011.06.030',
                'composition': {'Ba': 1.0, 'Fe': 11.8, 'O': 19.0},
                'valence': {'Ba': 2.0, 'Fe': 36 / 11.8, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.jmmm.2011.06.030',
                'composition': {'Ba': 1.0, 'Ti': 0.3, 'Fe': 11.2, 'Co': 0.3, 'O': 19.0},
                'valence': {'Ba': 2.0, 'Ti': 4.0, 'Fe': 3.0, 'Co': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.jmmm.2011.06.030',
                'composition': {'Ba': 1.0, 'Ti': 0.6, 'Fe': 10.6, 'Co': 0.6, 'O': 19.0},
                'valence': {'Ba': 2.0, 'Ti': 4.0, 'Fe': 3.0, 'Co': 3.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.jmmm.2011.06.030',
                'composition': {'Ba': 1.0, 'Ti': 0.9, 'Fe': 10.0, 'Co': 0.9, 'O': 19.0},
                'valence': {'Ba': 2.0, 'Ti': 4.0, 'Fe': 3.0, 'Co': 8 / 3.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.matlet.2014.04.029',
                'composition': {'Ca': 2.975, 'Si': 2.0, 'O': 7.0},
                'valence': {'Ca': 2.0, 'Si': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1021/ic501510u',
                'composition': {'Ba': 1.0, 'Ce': 0.9, 'O': 3.0, 'Ho': 0.1},
                'valence': {'Ba': 2.0, 'Ce': 4.0, 'O': -2.0, 'Ho': 3.0},
            },
            {
                'doi': '10.1016/j.ssi.2004.10.015',
                'composition': {'La': 1.97, 'K': 0.03, 'Mo': 2.0, 'O': 9.0},
                'valence': {'La': 3.0, 'K': 1.0, 'Mo': 6.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.ssi.2004.10.015',
                'composition': {'La': 1.95, 'K': 0.05, 'Mo': 2.0, 'O': 9.0},
                'valence': {'La': 3.0, 'K': 1.0, 'Mo': 6.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.ssi.2004.10.015',
                'composition': {'La': 1.925, 'K': 0.075, 'Mo': 2.0, 'O': 9.0},
                'valence': {'La': 3.0, 'K': 1.0, 'Mo': 6.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.ssi.2004.10.015',
                'composition': {'La': 1.9, 'K': 0.1, 'Mo': 2.0, 'O': 9.0},
                'valence': {'La': 3.0, 'K': 1.0, 'Mo': 6.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.ssi.2004.10.015',
                'composition': {'La': 1.85, 'K': 0.15, 'Mo': 2.0, 'O': 9.0},
                'valence': {'La': 3.0, 'K': 1.0, 'Mo': 6.0, 'O': -2.0},
            },
            {
                'doi': '10.1149/1.2988135',
                'composition': {'Sn': 0.9, 'In': 0.1, 'P': 2.0, 'O': 7.0},
                'valence': {'Sn': 4.0, 'In': 3.0, 'P': 5.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.jpowsour.2017.10.091',
                'composition': {'Ba': 3.0, 'Ca': 1.18, 'Nb': 1.77, 'Ni': 0.05, 'O': 9.0},
                'valence': {'Ba': 2.0, 'Ca': 2.0, 'Nb': 5.0, 'Ni': 4.0, 'O': -2.0},
            },
            {
                'doi': '10.1016/j.ijhydene.2011.10.020',
                'composition': {'La': 0.99, 'Ca': 0.01, 'Nb': 0.99, 'Ti': 0.01, 'O': 4.0},
                'valence': {'La': 3, 'Ca': 2, 'Nb': 5, 'Ti': 4, 'O': -2},
            },
            {
                'doi': '10.1016/j.memsci.2018.01.068',
                'composition': {'La': 5.5, 'W': 0.45, 'Nb': 0.15, 'Mo': 0.4, 'O': 11.25},
                'valence': {'La': 3, 'W': 6, 'Nb': 5, 'Mo': 6, 'O': -2},
            },
        ]

        self.negative_cases = [
            {
                'doi': '10.1016/S0167-577X(02)00433-0',
                'composition': {'Y': 2.0, 'Ba': 1.0, 'Cu': 1.0, 'O': 45.0},
                'valence': None,
            },
            {
                'doi': '10.1039/C7TC00112F',
                'composition': {'Cr': 1.9, 'Mn': 0.1, 'Al': 1.0},
                'valence': {'Cr': 0, 'Mn': 0, 'Al': 0},
            },
            {
                'doi': '10.1039/C7TC00112F',
                'composition': {'Cr': 1.8, 'Mn': 0.2, 'Al': 1.0},
                'valence': {'Cr': 0, 'Mn': 0, 'Al': 0},
            },
            {
                'doi': '10.1039/C7TC00112F',
                'composition': {'Cr': 1.7, 'Mn': 0.3, 'Al': 1.0},
                'valence': {'Cr': 0, 'Mn': 0, 'Al': 0},
            },
            {
                'doi': '10.1016/S0921-4534(00)00269-0',
                'composition': {'Bi': 2.4, 'Sr': 2.0, 'Gd': 1.0,
                                'Cu': 2.0},
                'valence': {'Bi': 0, 'Sr': 0, 'Ca': 0, 'Cu': 0},
            },
            {
                'doi': '10.1016/S0921-4534(00)00269-0',
                'composition': {'Bi': 2.4, 'Sr': 2.0, 'Ca': 1.0,
                                'Cu': 2.0},
                'valence': {'Bi': 0, 'Sr': 0, 'Ca': 0, 'Cu': 0},
            },
            {
                'doi': '10.1016/j.jpcs.2008.12.006',
                'composition': {'Bi': 1.7, 'Pb': 0.3, 'Sr': 2.0,
                                'Ca': 2.0, 'Cu': 3.0},
                'valence': {'Bi': 0, 'Pb': 0, 'Sr': 0, 'Ca': 0,
                            'Cu': 0},
            },
            {
                'doi': '10.1007/s10853-008-2897-2',
                'composition': {'Ni': 1.0, 'Ti': 1.0},
                'valence': {'Ni': 0, 'Ti': 0},
            },
        ]

        self.questionable_cases = [

        ]

        self.conditions = {
            'all_oxi_states': False,
            'all_metal_oxi_states': True,
            'add_zero_valence': False,
            'add_compensator': True,
        }

    def convert_valence_to_float(self, data, inplace=False):
        if inplace:
            data_reformated = data
        else:
            data_reformated = copy.deepcopy(data)

        for case in data:
            if case['valence'] is not None:
                for k in case['valence'].keys():
                    case['valence'][k] = float(case['valence'][k])

        return data_reformated

    def test_positive_cases(self):
        test_set = self.positive_cases
        for case in test_set:
            cal_valence = get_valence_single_composition(case['composition'])
            if len(cal_valence) == 0:
                cal_valence = get_valence_single_composition(case['composition'], **self.conditions)
            print(case['composition'])
            print(case['valence'])
            print(cal_valence)
            if len(cal_valence) > 0 and 'X' in cal_valence[0]:
                del cal_valence[0]['X']
            if (len(cal_valence) != 1) or cal_valence[0] != case['valence']:
                print('False')
            print()
            self.assertEqual(1, len(cal_valence))
            self.assertEqual(case['valence'], cal_valence[0])

    def test_negative_cases(self):
        test_set = self.negative_cases
        for case in test_set:
            cal_valence = get_valence_single_composition(case['composition'])
            if len(cal_valence) == 0:
                cal_valence = get_valence_single_composition(case['composition'], **self.conditions)
            print(case['composition'])
            print(cal_valence)
            print(len(cal_valence) == 0)
            print()
            self.assertEqual(0, len(cal_valence))

    def test_questionable_cases(self):
        test_set = self.questionable_cases
        for case in test_set:
            cal_valence = get_valence_single_composition(case['composition'], **self.conditions)
            print(case['composition'])
            print(cal_valence)
            print()
