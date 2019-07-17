This module is used to solve valence states of materials modified from Pymatgen default valence solver.
The principles to solve valence states is the same as Pymatgen, but it is faster.

Some packages need to be installed before use: 

    Synthepedia from CederGroupHub
    sympy
    unidecode

Install: 

    git clone git@github.com:CederGroupHub/ValenceSolver.git
    cd ValenceSolver
    pip install -e .

Example usage: 

For detailed code please refer to example/example.py. The script examples.py read in our current database and add a new field `valence` in the `composition` field for each material. There are some explanations on the returned values. 

(1) Normally, the modified data is as below, the value of `valence` is a dict similar to `elements`.

    'composition': [{'amount': '1.0',
                  'elements': {'C': '1.0',
                               'O': '3.0',
                               'Sr': '1.0'},
                  'formula': 'SrCO3',
                  'valence': {'C': 4.0, 'O': -2.0, 'Sr': 2.0}}],

(2) Sometimes, the valence state is not solved successfully. The value of `valence` would be None. 

    [{'amount': '1.0',
      'elements': {'Mn': 'y + 4', 'O': '3.0', 'Sr': '1.0', 'Ti': '-y + 1'},
      'formula': 'SrTi1-yMn4+yO3',
      'valence': None}]

(3) Sometimes, there are multiple instances corresponding to one material when there is a variable. If the valence states are the same for all the instances, the value of `valence` is still a dict as (1). Otherwise, the value of `valence` is a list including valence states for each instance. 

    [{'amount': '1.0',
      'elements': {'Bi': 'x + 2.5', 'Na': '-x + 0.5', 'Nb': '2.0', 'O': '9.0'},
      'formula': 'Na0.5-xBi2.5+xNb2O9',
      'valence': [{'amounts_vars': {'x': 0.05},
                   'elements': {'Bi': 2.55, 'Na': 0.45, 'Nb': 2.0, 'O': 9.0},
                   'valence': {'Bi': 3.0, 'Na': 1.0, 'Nb': 4.95, 'O': -2.0}},
                  {'amounts_vars': {'x': 0.1},
                   'elements': {'Bi': 2.6, 'Na': 0.4, 'Nb': 2.0, 'O': 9.0},
                   'valence': {'Bi': 3.0, 'Na': 1.0, 'Nb': 4.9, 'O': -2.0}}]}]

(4) Sometimes, valence states are not solved for every elements because the material contains a variable and the value of this variable is zero. To specify this situation, the value of `valence` is a also list including valence states and the values of variables. 

    [{'amount': '1.0',
      'elements': {'Ba': '1.0',
                   'Ge': 'y',
                   'O': '3.0',
                   'Sn': 'x',
                   'Ti': '-x - y + 1'},
      'formula': 'BaTi1-x-ySnxGeyO3',
      'valence': [{'amounts_vars': {'x': 0.0, 'y': 0.05},
                   'elements': {'Ba': 1.0, 'Ge': 0.05, 'O': 3.0, 'Ti': 0.95},
                   'valence': {'Ba': 2.0, 'Ge': 4.0, 'O': -2.0, 'Ti': 4.0}},
                  {'amounts_vars': {'x': 0.0, 'y': 0.0},
                   'elements': {'Ba': 1.0, 'O': 3.0, 'Ti': 1.0},
                   'valence': {'Ba': 2.0, 'Ge': 4.0, 'O': -2.0, 'Ti': 4.0}}]}]





