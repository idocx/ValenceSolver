# -*- coding: utf-8 -*-
from __future__ import division, unicode_literals

from itertools import product

import os

from collections import defaultdict

from monty.serialization import loadfn
from six.moves import zip

from functools import total_ordering

from pymatgen.core.periodic_table import Element, Specie


from pymatgen import Composition
import  pulp

"""
The module is modified based on the original Composition class in pymatgen 
to make the function to guess oxidation states more efficiently.
This module implements a Composition class to represent compositions.
use function oxi_state_guesses or oxi_state_guesses_most_possible to solve valence faster.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 10, 2012"
__modifiedBy__ = "Tanjin He"


@total_ordering
class CompositionInHouse(Composition):
    """
    
    Represents a Composition, which is essentially a {element:amount} mapping
    type. Composition is written to be immutable and hashable,
    unlike a standard Python dict.

    Note that the key can be either an Element or a Specie. Elements and Specie
    are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Specie and
    would be put in separate keys. This differentiation is deliberate to
    support using Composition to determine the fraction of a particular Specie.

    Works almost completely like a standard python dictionary, except that
    __getitem__ is overridden to return 0 when an element is not found.
    (somewhat like a defaultdict, except it is immutable).

    Also adds more convenience methods relevant to compositions, e.g.,
    get_fraction.

    It should also be noted that many Composition related functionality takes
    in a standard string as a convenient input. For example,
    even though the internal representation of a Fe2O3 composition is
    {Element("Fe"): 2, Element("O"): 3}, you can obtain the amount of Fe
    simply by comp["Fe"] instead of the more verbose comp[Element("Fe")].

    >>> comp = Composition("LiFePO4")
    >>> comp.get_atomic_fraction(Element("Li"))
    0.14285714285714285
    >>> comp.num_atoms
    7.0
    >>> comp.reduced_formula
    'LiFePO4'
    >>> comp.formula
    'Li1 Fe1 P1 O4'
    >>> comp.get_wt_fraction(Element("Li"))
    0.04399794666951898
    >>> comp.num_atoms
    7.0
    """
    # assign X=0 with a large probability, so that X is always 0 for usual cases
    LARGE_PROBABILITY = 10000
    # might be wrong composition when valence of X is too large
    # 0.3 is an arbitrary threshold for warning
    X_valence_warning_level = 0.3

    def __init__(self, *args, **kwargs):  # allow_negative=False
        """
        Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Args:
            Any form supported by the Python built-in dict() function.

            1. A dict of either {Element/Specie: amount},

               {string symbol:amount}, or {atomic number:amount} or any mixture
               of these. E.g., {Element("Li"):2 ,Element("O"):1},
               {"Li":2, "O":1}, {3:2, 8:1} all result in a Li2O composition.
            2. Keyword arg initialization, similar to a dict, e.g.,

               Composition(Li = 2, O = 1)

            In addition, the Composition constructor also allows a single
            string as an input formula. E.g., Composition("Li2O").

            allow_negative: Whether to allow negative compositions. This
                argument must be popped from the \\*\\*kwargs due to \\*args
                ambiguity.
        """
        super().__init__(*args, **kwargs)
        # Load prior probabilities of oxidation states, used to rank solutions
        if not Composition.oxi_prob:
            module_dir = os.path.join(os.path.
                                      dirname(os.path.abspath(__file__)))
            all_data = loadfn(os.path.join(module_dir,
                                           "analysis", "icsd_bv.yaml"))
            Composition.oxi_prob = {Specie.from_string(sp): data
                                    for sp, data in
                                    all_data["occurrence"].items()}

    def get_oxid_state_guess_essentials(self,
                                        oxi_states_override=None,
                                        all_metal_oxi_states=False,
                                        all_oxi_states=False,
                                        max_sites=None,
                                        add_compensator=False,
                                        double_el_amt=False):
        """
        This is a preprocess function for oxi_state_guesses() and oxi_state_guesses()
        It will return the essential component need to guess the valence state, including
        all elements, possible oxidation state for each element, and the amount of each
        element

        :param oxi_states_override: dict. dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
        :param all_metal_oxi_states: bool. if True, besides the pymatgen
                Element.icsd_oxidation_states, the positive valence states in
                Element.oxidation_states will be appended for metals. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results. Therefore, this option is a
                trade off.
        :param all_oxi_states: bool. if True, an element defaults to
                all oxidation states in pymatgen Element.oxidation_states.
                Otherwise, default is Element.icsd_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
        :param max_sites: int. if possible, will reduce Compositions to at most
                this many many sites to speed up oxidation state guesses. Set
                to -1 to just reduce fully.
        :param add_compensator: bool. if True, add a fake element "X" to compensate
                charge. "X" represents oxygen deficiency/excess, oxidized/reduced metal,
                or other factors leading to an unsual valence environment.
        :param double_el_amt: bool. if True, double the amount of each element, because
                sometimes there is a valence skipping effect but the amount is odd.
                https://web.stanford.edu/group/fisher/research/valence_skipping_elements.html
        :return: (els, el_amt, all_oxids)
            els: list of all elements in composition
            el_amt: dict of {el: amount} as composition
            all_oxids: dict of {el: (v1, v2, ...)}. all possible valence states for each element
        """

        comp = self.copy()

        # reduce Composition if necessary
        if max_sites == -1:
            comp = self.reduced_composition

        elif max_sites and comp.num_atoms > max_sites:
            reduced_comp, reduced_factor = self. \
                get_reduced_composition_and_factor()
            if reduced_factor > 1:
                reduced_comp *= max(1, int(max_sites / reduced_comp.num_atoms))
                comp = reduced_comp  # as close to max_sites as possible
            if comp.num_atoms > max_sites:
                raise ValueError("Composition {} cannot accommodate max_sites "
                                 "setting!".format(comp))

        oxi_states_override = oxi_states_override or {}

        # assert: Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer "
                             "values in Composition!")

        el_amt = comp.get_el_amt_dict()
        els = list(el_amt.keys())
        all_oxids = {}
        for idx, el in enumerate(els):
            if oxi_states_override.get(el):
                oxids = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            elif all_metal_oxi_states:
                # icsd_oxidation_states + positive valence in all_oxi_states
                oxids = Element(el).icsd_oxidation_states
                if Element(el).is_metal:
                    all_positive = set(filter(
                        lambda x: x>0,
                        Element(el).oxidation_states
                    ))
                    oxids = tuple(set(oxids) | all_positive)
            else:
                oxids = Element(el).icsd_oxidation_states or \
                        Element(el).oxidation_states
            all_oxids[el] = oxids
        if add_compensator and 'O' in el_amt:
            el_amt['X'] = el_amt['O']
            els = list(el_amt.keys())
            all_oxids['X'] = (-1, 0, 1)
        if double_el_amt:
            for el in el_amt:
                el_amt[el] *= 2
        return els, el_amt, all_oxids

    def oxi_state_guesses(self,
                          oxi_states_override=None,
                          target_charge=0,
                          all_metal_oxi_states=False,
                          all_oxi_states=False,
                          max_sites=None,
                          add_compensator=False,
                          double_el_amt=False):
        """
        Checks if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g., if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        :param oxi_states_override: dict. dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
        :param all_metal_oxi_states: bool. if True, besides the pymatgen
                Element.icsd_oxidation_states, the positive valence states in
                Element.oxidation_states will be appended for metals. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results. Therefore, this option is a
                trade off.
        :param all_oxi_states: bool. if True, an element defaults to
                all oxidation states in pymatgen Element.oxidation_states.
                Otherwise, default is Element.icsd_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
        :param max_sites: int. if possible, will reduce Compositions to at most
                this many many sites to speed up oxidation state guesses. Set
                to -1 to just reduce fully.
        :param add_compensator: bool. if True, add a fake element "X" to compensate
                charge. "X" represents oxygen deficiency/excess, oxidized/reduced metal,
                or other factors leading to an unsual valence environment.
        :param double_el_amt: bool. if True, double the amount of each element, because
                sometimes there is a valence skipping effect but the amount is odd.
                https://web.stanford.edu/group/fisher/research/valence_skipping_elements.html
        :return: a list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
        """
        all_sols = []  # will contain all solutions
        all_scores = []  # will contain a score for each solution
        # for each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_sums = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores = defaultdict(set)  # dict of el_idx, sum -> score

        els, el_amt, all_oxids = self.get_oxid_state_guess_essentials(
            oxi_states_override=oxi_states_override,
            all_metal_oxi_states=all_metal_oxi_states,
            all_oxi_states=all_oxi_states,
            max_sites=max_sites,
            add_compensator=add_compensator,
            double_el_amt=double_el_amt
        )
        for idx, el in enumerate(els):
            el_sum_scores[idx] = {}
            el_sums.append([])
            all_sums, sum_scores = CompositionInHouse.get_possible_sums(
                el,
                all_oxids[el],
                int(el_amt[el]),
                add_compensator=add_compensator
            )
            for tmp_index, tmp_sum in enumerate(all_sums):
                el_sums[idx].append(tmp_sum)
                score = sum_scores[tmp_index]
                el_sum_scores[idx][tmp_sum] = max(el_sum_scores[idx].get(tmp_sum, 0), score)

        for x in product(*el_sums):
            # each x is a trial of one possible oxidation sum for each element
            if sum(x) == target_charge:  # charge balance condition
                el_sum_sol = dict(zip(els, x))  # element->oxid_sum
                # normalize oxid_sum by amount to get avg oxid state
                sol = {el: v / el_amt[el] for el, v in el_sum_sol.items()}
                all_sols.append(sol)  # add the solution to the list of solutions

                # determine the score for this solution
                score = 0
                for idx, v in enumerate(x):
                    score += el_sum_scores[idx][v]
                all_scores.append(score)

        # sort the solutions by highest to lowest score
        all_sols = [x for (y, x) in sorted(zip(all_scores, all_sols),
                                           key=lambda pair: pair[0],
                                           reverse=True)]

        # elementary materials are not solved but the valence should be 0
        if not all_sols and len(els) == 1:
            all_sols = [{el: 0.0 for el in els}]

        return all_sols

    def oxi_state_guesses_most_possible(self):
        """
        guess the most possible oxidation states based on the same method in pymatggen.
        linear programming is used to accelerate the computation.
        relaxation assuptions are adopted when there is no solution using default solver in pymatgen.

        :param composition: can be a plain dict or a plain string that pymatgen can interpret
        :return: (oxi_state, is_usual, comments)
            oxi_state: dict of oxidation state {el: valence}
            is_usual: bool. if True. The solution is the same as the default solution from
                pymatgen. Otherwise, there is no solution found using the default solver in
                pymatgen. Then relaxation is used by assuming the composition is unusual.
                The unusual situations include alloy, oxygen deficient/excess, more
                oxidized/reduced metal states, etc.
            comments: list of strings. Details when is_usual == False.
        """
        is_usual = True
        comments = []
        el_amt = self.get_el_amt_dict()

        # solution same as pymatgen, but much faster,
        # so that we can do relaxation if no solution found
        oxi_state = self._oxi_state_guesses_most_possible(
            all_metal_oxi_states=False,
            all_oxi_states=False,
            add_compensator=False,
            double_el_amt=False,
        )

        # deal with alloy
        if len(oxi_state) == 0 and self.is_alloy():
            oxi_state = [
                {el: 0.0 for el in el_amt}
            ]
            is_usual = False
            comments = ['is alloy']

        # solve again with relaxation
        if len(oxi_state) == 0:
            oxi_state = self._oxi_state_guesses_most_possible(
                all_metal_oxi_states=True,
                all_oxi_states=False,
                add_compensator=True,
                double_el_amt=False,
            )
            is_usual = False
            comments.append('all possible positive valence states are used for metals')

        # might be wrong composition when valence of X is too large
        # solve again for large X
        if (len(oxi_state) > 0
            and 'X' in oxi_state[0]
            and abs(oxi_state[0]['X']) > CompositionInHouse.X_valence_warning_level):
            if (oxi_state[0]['X'] == 1.0
                and el_amt['O'] == 2
                and len(el_amt) == 2
            ):
                # possibly peroxide, correct
                is_usual = True
                oxi_state[0]['O'] = -1.0
                del oxi_state[0]['X']
            else:
                # solve again by doubling the amount in case there is a valence skipping effect
                oxi_state = self._oxi_state_guesses_most_possible(
                    all_metal_oxi_states=True,
                    all_oxi_states=False,
                    add_compensator=True,
                    double_el_amt=True
                )

        if (len(oxi_state) > 0 and 'X' in oxi_state[0]):
            if oxi_state[0]['X'] > 0:
                comments.append('possibly oxygen deficient, or some elements are more oxidized than usual')
            elif oxi_state[0]['X'] < 0:
                comments.append('possibly oxygen excess, or some elements are more reduced than usual')
            if abs(oxi_state[0]['X']) > CompositionInHouse.X_valence_warning_level:
                comments.append('Warning: the input composition might be wrong')
            del oxi_state[0]['X']
        return oxi_state, is_usual, comments

    def _oxi_state_guesses_most_possible(self,
                                        oxi_states_override=None,
                                        target_charge=0,
                                        all_metal_oxi_states=False,
                                        all_oxi_states=False,
                                        max_sites=None,
                                        add_compensator=False,
                                        double_el_amt=False):
        """
        Checks if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g., if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        :param oxi_states_override: dict. dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
        :param all_metal_oxi_states: bool. if True, besides the pymatgen
                Element.icsd_oxidation_states, the positive valence states in
                Element.oxidation_states will be appended for metals. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results. Therefore, this option is a
                trade off.
        :param all_oxi_states: bool. if True, an element defaults to
                all oxidation states in pymatgen Element.oxidation_states.
                Otherwise, default is Element.icsd_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
        :param max_sites: int. if possible, will reduce Compositions to at most
                this many many sites to speed up oxidation state guesses. Set
                to -1 to just reduce fully.
        :param add_compensator: bool. if True, add a fake element "X" to compensate
                charge. "X" represents oxygen deficiency/excess, oxidized/reduced metal,
                or other factors leading to an unsual valence environment.
        :param double_el_amt: bool. if True, double the amount of each element, because
                sometimes there is a valence skipping effect but the amount is odd.
                https://web.stanford.edu/group/fisher/research/valence_skipping_elements.html
        :return: a list of dicts - the length is always 1 or 0 because only the most
                possible solution is returned. each dict reports an element symbol
                and average oxidation state across all sites in that composition.
                If the composition is not charge balanced, an empty list is returned.
        """

        all_sols = []  # will contain all solutions
        all_scores = []  # will contain a score for each solution
        els, el_amt, all_oxids = self.get_oxid_state_guess_essentials(
            oxi_states_override=oxi_states_override,
            all_metal_oxi_states=all_metal_oxi_states,
            all_oxi_states=all_oxi_states,
            max_sites=max_sites,
            add_compensator=add_compensator,
            double_el_amt=double_el_amt
        )
        solution, score = CompositionInHouse.get_most_possible_solution(
            els,
            all_oxids,
            el_amt,
            add_compensator=add_compensator,
            target_charge=target_charge,
        )
        if solution:
            all_sols = [solution]
            all_scores = [score]
        
        # elementary materials are not solved but the valence should be 0
        if not all_sols and len(els) == 1:
            all_sols = [{el:0.0 for el in els}]
            
        return all_sols

    @staticmethod
    def get_possible_sums(el, oxi_states, el_amt, add_compensator=False):
        # goal
        all_sums = []
        sum_scores = []

        oxi_names = [str(tmp_state) for tmp_state in oxi_states]
        if el in Element.__members__:
            costs = {
                oxi_names[i]: Composition.oxi_prob.get(
                    Specie(el, oxi_states[i]), -CompositionInHouse.LARGE_PROBABILITY
                ) for i in range(len(oxi_states))
            }
        else:
            costs = {}
        if add_compensator and el == 'X':
            costs = {
                '0': CompositionInHouse.LARGE_PROBABILITY,
                '1': 1,
                '-1': 1,
            }

        # definition of problem
        problem = pulp.LpProblem('score maximization problem', pulp.LpMaximize)
        # definition of variables
        oxi_vars = pulp.LpVariable.dicts('oxi_states', oxi_names, lowBound=0, cat=pulp.LpInteger)
        # objective function
        problem += pulp.lpSum([costs[i]*oxi_vars[i] for i in oxi_names]), \
                   'total probability of valence combination for one element'
        # constraints
        problem += pulp.lpSum([oxi_vars[i] for i in oxi_names]) == el_amt, 'numberOfElements'
        problem += pulp.lpSum(
            [oxi_states[i]*oxi_vars[oxi_names[i]] for i in range(len(oxi_states))]
        ) == min(oxi_states)*el_amt, 'sumOfValence'

        # change constraints and find possible sum of valence
        for tmp_sum in range(min(oxi_states)*el_amt, max(oxi_states)*el_amt+1):
            problem.constraints['sumOfValence'] = pulp.lpSum(
                [oxi_states[i]*oxi_vars[oxi_names[i]] for i in range(len(oxi_states))]
            ) == tmp_sum
            problem.solve()
            if pulp.LpStatus[problem.status] == 'Optimal':
                all_sums.append(tmp_sum)
                sum_scores.append(pulp.value(problem.objective))
                pass

        return all_sums, sum_scores

    @staticmethod
    def get_most_possible_solution(all_els,
                                   all_oxi_states,
                                   all_el_amts,
                                   add_compensator=False,
                                   target_charge=0):
        # goal
        solution = {}
        score = 0

        oxi_names = []
        costs = {}
        for el in all_els:
            oxi_names.extend([el+str(tmp_state) for tmp_state in all_oxi_states[el]])
            if el in Element.__members__:
                costs.update(
                    {
                        el+str(tmp_state): Composition.oxi_prob.get(
                            Specie(el, tmp_state), -CompositionInHouse.LARGE_PROBABILITY
                        )
                        for tmp_state in all_oxi_states[el]
                    }
                )
            if add_compensator:
                costs['X0'] = CompositionInHouse.LARGE_PROBABILITY
                costs['X1'] = 1
                costs['X-1'] = 1

        # definition of problem
        problem = pulp.LpProblem('score maximization problem', pulp.LpMaximize)
        # definition of variables
        oxi_vars = pulp.LpVariable.dicts('oxi_states', oxi_names, lowBound=0, cat=pulp.LpInteger)
        # objective function
        problem += pulp.lpSum([costs[i]*oxi_vars[i] for i in oxi_names]), \
                   'total probability of valence combination for one element'
        # constraints
        for el in all_els:
            problem += pulp.lpSum(
                [oxi_vars[el+str(tmp_state)] for tmp_state in all_oxi_states[el]]
            ) == all_el_amts[el], 'numberOfElements_'+el
        problem += pulp.lpSum(
            [tmp_state*oxi_vars[el+str(tmp_state)] for el in all_els for tmp_state in all_oxi_states[el]
             ]) == target_charge, 'sumOfValence'

        problem.solve()
        if pulp.LpStatus[problem.status] == 'Optimal':
            for el in all_els:
                solution[el] = pulp.value(
                    pulp.lpSum([tmp_state*oxi_vars[el+str(tmp_state)] for tmp_state in all_oxi_states[el]])
                )/float(all_el_amts[el])
            score = pulp.value(problem.objective)

        return solution, score

    def is_alloy(self):
        return all([Element(el).is_metal for el in self.get_el_amt_dict()])

    @staticmethod
    def get_most_possible_oxi_state_of_composition(composition):
        """
        a wrapper using the method oxi_state_guesses_most_possible to guess the most possible
        oxidation states based on the same method in pymatggen. The input can be a plain dict
        or a plain string that pymatgen can interpret

        :param composition: can be a plain dict or a plain string that pymatgen can interpret
        :return: (oxi_state, is_usual, comments)
            oxi_state: dict of oxidation state {el: valence}
            is_usual: bool. if True. The solution is the same as the default solution from
                pymatgen. Otherwise, there is no solution found using the default solver in
                pymatgen. Then relaxation is used by assuming the composition is unusual.
                The unusual situations include alloy, oxygen deficient/excess, more
                oxidized/reduced metal states, etc.
            comments: list of strings. Details when is_usual == False.
        """
        valence_comp = CompositionInHouse(composition)
        valence_comp, inte_factor = valence_comp.get_integer_formula_and_factor()
        valence_comp = CompositionInHouse(valence_comp)
        oxi_state, is_usual, comments = valence_comp.oxi_state_guesses_most_possible()
        return oxi_state, is_usual, comments


