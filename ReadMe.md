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
    
Usage: 

The minimal code is as following. We first create a `GeneralComposition` object using `Synthepedia`. `GeneralComposition` is able to generate all composition instances by substituting varibles in formula and to check if the stoichiometric number is valid (should not be negative) after substitution. Then, we solve the valence by calling `get_material_valence(tmp_mat_obj, valence_cache=valence_cache)`. The valence_cache is optional, which is designed to save solved cases to avoid duplication in the future. It can be initialized as an emtpy dict and the value would be updated automatically after calling `get_material_valence`.  

    tmp_mat_obj = to_GeneralMat_obj(
        composition=composition,
        amounts_vars=amounts_vars,
        elements_vars=elements_vars
    )
    valence = get_material_valence(tmp_mat_obj, valence_cache=valence_cache)

Example: 

For detailed code please refer to example/example.py. The script examples.py read in our current database and add a new field `valence` in the `composition` field for each material. The returned value would be as following if the valence could be solved. The `valence` is a list because there could be multiple composition instances when there are variables in the formula. Each element in the `valence` filed is a dict which specify the valence state, values for variables and corresponding composition. Sometimes, the valence state is not solved successfully. The value of `valence` would be None.

    'composition': [
        { 'amount': '1.0',
          'elements': {'Bi': 'x + 2.5', 'Na': '-x + 0.5', 'Nb': '2.0', 'O': '9.0'},
          'formula': 'Na0.5-xBi2.5+xNb2O9',
          'valence': [{'amounts_vars': [{'x': 0.05}],
                       'elements': [{'Bi': 2.55, 'Na': 0.45, 'Nb': 2.0, 'O': 9.0}],
                       'valence': {'Bi': 3.0, 'Na': 1.0, 'Nb': 4.95, 'O': -2.0}d},
                      {'amounts_vars': [{'x': 0.1}],
                       'elements': [{'Bi': 2.6, 'Na': 0.4, 'Nb': 2.0, 'O': 9.0}],
                       'valence': {'Bi': 3.0, 'Na': 1.0, 'Nb': 4.9, 'O': -2.0}}]}
    ]
    