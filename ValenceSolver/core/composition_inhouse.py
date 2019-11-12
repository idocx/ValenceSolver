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
        Composition.__init__(self, *args, **kwargs)

    def oxi_state_guesses(self, oxi_states_override=None, target_charge=0,
                          all_oxi_states=False, max_sites=None):
        """
        Checks if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g., if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many many sites to speed up oxidation state guesses. Set
                to -1 to just reduce fully.

        Returns:
            A list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
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

        # Load prior probabilities of oxidation states, used to rank solutions
        if not Composition.oxi_prob:
            module_dir = os.path.join(os.path.
                                      dirname(os.path.abspath(__file__)))
            all_data = loadfn(os.path.join(module_dir,
                                           "analysis", "icsd_bv.yaml"))
            Composition.oxi_prob = {Specie.from_string(sp): data
                                    for sp, data in
                                    all_data["occurrence"].items()}

        oxi_states_override = oxi_states_override or {}

        # assert: Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer "
                             "values in Composition!")

        # for each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_amt = comp.get_el_amt_dict()
        els = el_amt.keys()
        el_sums = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores = defaultdict(set)  # dict of el_idx, sum -> score

        for idx, el in enumerate(els):
            el_sum_scores[idx] = {}
            el_sums.append([])
            if oxi_states_override.get(el):
                oxids = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            else:
                oxids = Element(el).icsd_oxidation_states or \
                        Element(el).oxidation_states

            all_sums, all_scores = CompositionInHouse.get_possible_sums(el, oxids, int(el_amt[el]))
            for tmp_index, tmp_sum in enumerate(all_sums):
                el_sums[idx].append(tmp_sum)
                score = all_scores[tmp_index]
                el_sum_scores[idx][tmp_sum] = max(el_sum_scores[idx].get(tmp_sum, 0), score)

        all_sols = []  # will contain all solutions
        all_scores = []  # will contain a score for each solution
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
            all_sols = [{el:0.0 for ele in els}]            

        return all_sols

    def oxi_state_guesses_most_possible(self,
                                        oxi_states_override=None,
                                        target_charge=0,
                                        all_oxi_states=False,
                                        max_sites=None,
                                        add_zero_valence=False,
                                        add_env_valence=False):
        """
        Checks if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g., if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many many sites to speed up oxidation state guesses. Set
                to -1 to just reduce fully.

        Returns:
            A list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
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

        # Load prior probabilities of oxidation states, used to rank solutions
        if not Composition.oxi_prob:
            module_dir = os.path.join(os.path.
                                      dirname(os.path.abspath(__file__)))
            all_data = loadfn(os.path.join(module_dir,
                                           "analysis", "icsd_bv.yaml"))
            Composition.oxi_prob = {Specie.from_string(sp): data
                                    for sp, data in
                                    all_data["occurrence"].items()}

        oxi_states_override = oxi_states_override or {}

        # assert: Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer "
                             "values in Composition!")

        # for each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_amt = comp.get_el_amt_dict()
        els = el_amt.keys()
        el_sums = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores = defaultdict(set)  # dict of el_idx, sum -> score

        all_sols = []  # will contain all solutions
        all_scores = []  # will contain a score for each solution
        all_oxids = {}
        for idx, el in enumerate(els):
            el_sum_scores[idx] = {}
            el_sums.append([])
            if oxi_states_override.get(el):
                oxids = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            else:
                oxids = Element(el).icsd_oxidation_states or \
                        Element(el).oxidation_states

            if add_zero_valence and 0 not in oxids:
                oxids = list(oxids)
                oxids.append(0)
            all_oxids[el] = oxids

        solution, score = CompositionInHouse.get_most_possible_solution(
            els,
            all_oxids,
            el_amt,
            add_zero_valence=add_zero_valence,
        )
        if solution:
            all_sols = [solution]
            all_scores = [score]
        
        # elementary materials are not solved but the valence should be 0
        if not all_sols and len(els) == 1:
            all_sols = [{el:0.0 for ele in els}]            
            
        return all_sols

    @staticmethod
    def get_possible_sums(el, oxi_states, el_amt):
        # goal
        all_sums = []
        all_scores = []

        oxi_names = [str(tmp_state) for tmp_state in oxi_states]
        costs = {oxi_names[i]: Composition.oxi_prob.get(Specie(el, oxi_states[i]), 0) for i in range(len(oxi_states))}

        # definition of problem
        problem = pulp.LpProblem('score maximization problem', pulp.LpMaximize)
        # definition of variables
        oxi_vars = pulp.LpVariable.dicts('oxi_states', oxi_names, lowBound=0, cat=pulp.LpInteger)
        # objective function
        problem += pulp.lpSum([costs[i]*oxi_vars[i] for i in oxi_names]), 'total probability of valence combination for one element'
        # constraints
        problem += pulp.lpSum([oxi_vars[i] for i in oxi_names]) == el_amt, 'numberOfElements'
        problem += pulp.lpSum([oxi_states[i]*oxi_vars[oxi_names[i]] for i in range(len(oxi_states))]) == min(oxi_states)*el_amt, 'sumOfValence'

        # change constraints and find possible sum of valence
        for tmp_sum in range(min(oxi_states)*el_amt, max(oxi_states)*el_amt+1):
            problem.constraints['sumOfValence'] = pulp.lpSum([oxi_states[i]*oxi_vars[oxi_names[i]] for i in range(len(oxi_states))]) == tmp_sum
            problem.solve()
            if pulp.LpStatus[problem.status] == 'Optimal':
                all_sums.append(tmp_sum)
                all_scores.append(pulp.value(problem.objective))
                pass

        return all_sums, all_scores

    @staticmethod
    def get_most_possible_solution(all_els,
                                   all_oxi_states,
                                   all_el_amts,
                                   add_zero_valence=False):
        # goal
        solution = {}
        score = 0

        oxi_names = []
        costs = {}
        for el in all_els:
            oxi_names.extend([el+str(tmp_state) for tmp_state in all_oxi_states[el]])
            costs.update({el+str(tmp_state): Composition.oxi_prob.get(Specie(el, tmp_state), -10000) for tmp_state in all_oxi_states[el]})
            if add_zero_valence:
                # TODO: include Oxygen anion or not
                costs[el+"0"] = 0.1/all_el_amts[el]
        print(costs)

        # definition of problem
        problem = pulp.LpProblem('score maximization problem', pulp.LpMaximize)
        # definition of variables
        oxi_vars = pulp.LpVariable.dicts('oxi_states', oxi_names, lowBound=0, cat=pulp.LpInteger)
        # objective function
        problem += pulp.lpSum([costs[i]*oxi_vars[i] for i in oxi_names]), 'total probability of valence combination for one element'
        # constraints
        for el in all_els:
            problem += pulp.lpSum([oxi_vars[el+str(tmp_state)] for tmp_state in all_oxi_states[el]]) == all_el_amts[el], 'numberOfElements_'+el
        problem += pulp.lpSum([tmp_state*oxi_vars[el+str(tmp_state)] for el in all_els for tmp_state in all_oxi_states[el]]) == 0, 'sumOfValence'

        problem.solve()
        if pulp.LpStatus[problem.status] == 'Optimal':
            for el in all_els:
                for tmp_state in all_oxi_states[el]:
                    print(el+str(tmp_state), pulp.value(oxi_vars[el+str(tmp_state)]))
                solution[el] = pulp.value(pulp.lpSum([tmp_state*oxi_vars[el+str(tmp_state)] for tmp_state in all_oxi_states[el]]))/float(all_el_amts[el])
            score = pulp.value(problem.objective)

        return solution, score
