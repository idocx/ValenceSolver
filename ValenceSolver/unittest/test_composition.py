# -*- coding: utf-8 -*-
import unittest
from pprint import pprint

from ValenceSolver.core.composition_inhouse import CompositionInHouse

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

10.1039/C7TC00112F
['Cr1.9Mn0.1Al', 'Cr1.8Mn0.2Al', 'Cr1.7Mn0.3Al']
{'Cr': 1.9, 'Mn': 0.1, 'Al': 1.0}
{'Cr': 1.8, 'Mn': 0.2, 'Al': 1.0}
{'Cr': 1.7, 'Mn': 0.3, 'Al': 1.0}

10.1016/S0921-4534(00)00269-0
['Sr2Ca1-xGdxCu2Bi2.4Oy', 'Sr2Ca1-xGdxCu2Bi2.4Oy']
{'Bi': 2.4, 'Sr': 2.0, 'Gd': 1.0, 'Cu': 2.0}
{'Bi': 2.4, 'Sr': 2.0, 'Ca': 1.0, 'Cu': 2.0}

10.1016/j.matlet.2014.04.029
['Ca2.975Si2O7']
{'Ca': 2.975, 'Si': 2.0, 'O': 7.0}

10.1016/j.physc.2013.02.013
['Ba2Ca2Cu3.5O10',
 'Ba2Ca2Ti0.25Cu3.25O10',
 'Ba2Ca2Ti0.5Cu3O10',
 'Ba2Ca2Ti0.75Cu2.75O10']
{'Cu': 3.5, 'Ba': 2.0, 'Ca': 2.0, 'O': 10.0}
{'Cu': 3.25, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.25, 'O': 10.0}
{'Cu': 3.0, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.5, 'O': 10.0}
{'Cu': 2.75, 'Ba': 2.0, 'Ca': 2.0, 'Ti': 0.75, 'O': 10.0}

10.1021/ic501510u
['BaHo0.1Ce0.9O3']
{'Ba': 1.0, 'Ce': 0.9, 'O': 3.0, 'Ho': 0.1}

['Ca2.97Dy0.03Si2O7',
 'Ca2.95Dy0.05Si2O7',
 'Ca2.93Dy0.07Si2O7',
 'Ca2.9Dy0.1Si2O7']
{'Ca': 2.97, 'Si': 2.0, 'O': 7.0, 'Dy': 0.03}
{'Ca': 2.95, 'Si': 2.0, 'O': 7.0, 'Dy': 0.05}
{'Ca': 2.93, 'Si': 2.0, 'O': 7.0, 'Dy': 0.07}
{'Ca': 2.9, 'Si': 2.0, 'O': 7.0, 'Dy': 0.1}


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

10.1016/j.jpcs.2008.12.006
['Sr2Ca2-xYxCu3Pb0.3Bi1.7Oy']
{'Bi': 1.7, 'Pb': 0.3, 'Sr': 2.0, 'Ca': 2.0, 'Cu': 3.0}

10.1149/1.2988135
['In0.1Sn0.9P2O7']
{'Sn': 0.9, 'In': 0.1, 'P': 2.0, 'O': 7.0}

10.1016/j.jpowsour.2017.10.091
['Ba3Ca1.18Nb1.77Ni0.05O9']
{'Ba': 3.0, 'Ca': 1.18, 'Nb': 1.77, 'Ni': 0.05, 'O': 9.0}

10.1016/j.ijhydene.2011.10.020
['Ca0.01La0.99Ti0.01Nb0.99O4']
{'La': 0.99, 'Ca': 0.01, 'Nb': 0.99, 'Ti': 0.01, 'O': 4.0}

10.1016/j.memsci.2018.01.068
['La5.5Nb0.15W0.45Mo0.4O11.25-δ']
{'La': 5.5, 'W': 0.45, 'Nb': 0.15, 'Mo': 0.4, 'O': 11.25}

10.1016/j.matchemphys.2010.10.027
['Ba0.95La0.033Zr0.09Ti0.91O3']
{'Ba': 0.95, 'Zr': 0.09, 'Ti': 0.91, 'O': 3.0, 'La': 0.033}

negative cases
10.1016/S0167-577X(02)00433-0
['BaY2CuO45']
{'Y': 2.0, 'Ba': 1.0, 'Cu': 1.0, 'O': 45.0}

10.1007/s10853-008-2897-2
['Ti2NiC']
{'Ni': 1.0, 'Ti': 1.0}

questionable
???
10.1016/j.jmmm.2011.06.030
['BaTixFe11.8-2*xCoxO19']
{'Ba': 1.0, 'Fe': 11.8, 'O': 19.0}
f Ba(CoTi)0.9Fe11O19
x¼0.3, (b) x¼0.6, (c) x¼0.9 and (d) x¼1.2

"""

class ConceptTest(unittest.TestCase):
    def test_sinter_to_dict(self):
        sinter_action = Sinter(temperature=500.0, reaction_time=1.0)

        self.assertEqual(sinter_action.to_dict(), {'temperature': 500,
                                                   'reaction_time': 1})