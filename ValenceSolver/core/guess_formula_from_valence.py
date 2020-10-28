from pymatgen.core.periodic_table import Element

import pulp

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'


def get_most_possible_solution(all_els,
                               all_oxi_states,
                               target_charge=0):
    # goal
    solution = {}
    score = 0
    valence_detail = {}

    oxi_names = []
    costs = {}
    for el in all_els:
        oxi_names.extend([el + str(tmp_state) for tmp_state in all_oxi_states[el]])
        if el in Element.__members__:
            costs.update(
                {
                    el + str(tmp_state): 1
                    for tmp_state in all_oxi_states[el]
                }
            )

    # definition of problem
    problem = pulp.LpProblem('score minimization problem', pulp.LpMinimize)
    # definition of variables
    oxi_vars = pulp.LpVariable.dicts('oxi_states', oxi_names, lowBound=1, cat=pulp.LpInteger)
    # objective function
    problem += pulp.lpSum([costs[i] * oxi_vars[i] for i in oxi_names]), \
               'total number of elements'
    # constraints
    problem += pulp.lpSum(
        [tmp_state * oxi_vars[el + str(tmp_state)] for el in all_els for tmp_state in all_oxi_states[el]
         ]) == target_charge, 'sumOfValence'

    problem.solve()
    if pulp.LpStatus[problem.status] == 'Optimal':
        for el in all_els:
            solution[el] = pulp.value(
                pulp.lpSum([oxi_vars[el + str(tmp_state)] for tmp_state in all_oxi_states[el]])
            )
            for tmp_state in all_oxi_states[el]:
                valence_detail[(el, tmp_state)] = pulp.value(oxi_vars[el + str(tmp_state)])
        score = pulp.value(problem.objective)

    return solution, score, valence_detail

def example_1():
    solution, score, valence_detail = get_most_possible_solution(
        all_els=['Fe', 'Mg', 'O'],
        all_oxi_states={
            'Fe': [+3, ],
            'Mg': [+2, ],
            'O': [-2, ],
        },
        target_charge=0
    )
    print('solution', solution)
    print('total number of elements', score)
    print('valence_detail', valence_detail)
    print()

def example_2():
    solution, score, valence_detail = get_most_possible_solution(
        all_els=['Fe', 'Al', 'O'],
        all_oxi_states={
            'Fe': [+3, ],
            'Al': [+3, ],
            'O': [-2, ],
        },
        target_charge=0
    )
    print('solution', solution)
    print('total number of elements', score)
    print('valence_detail', valence_detail)
    print()

def example_3():
    solution, score, valence_detail = get_most_possible_solution(
        all_els=['Li', 'Fe', 'P', 'O'],
        all_oxi_states={
            'Li': [+1, ],
            'Fe': [+2, ],
            'P': [+5, ],
            'O': [-2, ],
        },
        target_charge=0
    )
    print('solution', solution)
    print('total number of elements', score)
    print('valence_detail', valence_detail)
    print()

def example_4():
    solution, score, valence_detail = get_most_possible_solution(
        all_els=['Fe', 'Mg', 'O'],
        all_oxi_states={
            'Fe': Element('Fe').common_oxidation_states,
            'Mg': Element('Mg').common_oxidation_states,
            'O': Element('O').common_oxidation_states,
        },
        target_charge=0
    )
    print('solution', solution)
    print('total number of elements', score)
    print('valence_detail', valence_detail)
    print()

if __name__ == '__main__':
    example_1()
    example_2()
    example_3()
    example_4()
