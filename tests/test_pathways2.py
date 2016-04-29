from me_scripts.hindsight.pathways import (add_heterologous_pathway,
                                           get_designs,
                                           no_route_exchanges,
                                           exchange_for_metabolite_name,
                                           add_all_heterologous_pathways)
from me_scripts.hindsight.variables import min_biomass
from theseus import setup_model, load_model
from theseus.bigg import download_model
from minime.solve.algorithms import solve_at_growth_rate, binary_search
from minime.solve.symbolic import compile_expressions
import pytest
from pytest import raises

def test_add_heterologous_pathway_me(): # TODO test all pathways
    print('loading model')
    model = setup_model(load_model('ME'), 'EX_glc__D_e', aerobic=False)
    print('adding pathway')
    model = add_heterologous_pathway(model, 'ldhL', ignore_repeats=False)
    print('solving')
    # model.change_objective('EX_lac__D_e')
    print('knockouts')
    for r in model.reactions:
        for rid in ['LDH_D', 'ACALD']:
            if r.id.startswith(rid + '_'):
                r.knock_out()
    model.expressions = compile_expressions(model)
    # sol = solve_at_growth_rate(model, 0.5, compiled_expressions=model.expressions)
    sol = binary_search(model, min_mu=0, max_mu=1.1, mu_accuracy=1e-3, compiled_expressions=model.expressions)
    assert sol.x_dict['EX_lac__L_e'] > 1

if __name__ == '__main__':
    test_add_heterologous_pathway_me()
