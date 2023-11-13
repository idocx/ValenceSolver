This module is used to solve valence states of materials modified from Pymatgen default valence solver.
The principles to solve valence states is the same as Pymatgen, but it is faster.

Some packages need to be installed before use: 

    sympy
    unidecode

Install: 

    git clone git@github.com:CederGroupHub/ValenceSolver.git
    cd ValenceSolver
    pip install -e .
    
Usage: 

With a composition dict or a plain formula, the minimal code is as following. (example in example_basics/example_composition.py)

    from ValenceSolver.core.composition_inhouse import CompositionInHouse

    # The composition here is a plain dict like {'ele': number} without variables. 
    oxi_state, is_usual, comments = \
        CompositionInHouse.get_most_possible_oxi_state_of_composition(composition)

