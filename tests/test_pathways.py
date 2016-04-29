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

DEBUG = True

models_to_test = ['iJO1366']
substrates = ['EX_glc__D_e', 'EX_xyl__D_e']
designs = get_designs()
if DEBUG:
    pathway_names = ['crotonic_acid']
    designs = {n: designs[n] for n in pathway_names}
pathways_to_test = designs.items()

@pytest.fixture(scope='session', params=pathways_to_test)
def pathway_tuple(request):
    # return the pathway key
    return request.param

@pytest.fixture(scope='function', params=models_to_test)
def model(request):
    # return the loaded model
    model_id = request.param
    return download_model(model_id)

def test_add_heterologous_pathway(model, pathway_tuple):
    additions, design = pathway_tuple
    model = add_heterologous_pathway(model, additions)
    model = setup_model(model, substrates, aerobic=False)
    # produce some biomass
    biomass = model.objective.iterkeys().next()
    biomass.lower_bound = min_biomass
    # test the exchanges
    for r in [x for x in design[3].keys()
              if 'EX_' in x and x not in no_route_exchanges]:
        print
        print 'Testing %s in %s' % (additions, model.id)
        model.optimize()
        assert model.solution.f > 1e-3
        model.change_objective(r)
        model.optimize()
        assert model.solution.f > 1e-3
        print 'Max %s: %.4f' % (r, model.solution.f)

        # look out for loops
        assert model.solution.f < 45

def test_add_heterologous_pathway_me(pathway_tuple):
    additions, design = pathway_tuple
    model = setup_model(load_model('ME'), substrates, aerobic=False,
                        sur=1000, max_our=1000)
    model = add_heterologous_pathway(model, additions)
    # test the exchanges
    for r in [x for x in design[3].keys()
              if 'EX_' in x and x not in no_route_exchanges]:
        print
        print 'Testing %s in %s' % (additions, model.id)
        model.change_objective(r)
        sol = solve_at_growth_rate(model, min_biomass,
                                   compiled_expressions=model.expressions)
        flux = sol.x_dict[r]
        assert flux > 1e-3
        print 'Max %s: %.4f' % (r, flux)
    print 'done'

def test_repeat_additions(model):
    des = get_designs().iterkeys().next()
    m = add_heterologous_pathway(model, des)
    m = add_heterologous_pathway(m, des, ignore_repeats=True)
    with raises(Exception):
        m = add_heterologous_pathway(m, des)

def test_add_all_heterologous_pathways():
    model = download_model('iJO1366')
    model = add_all_heterologous_pathways(model)
    assert 'btal_c' in model.metabolites
    assert 'HACD1' in model.reactions

def test_exchange_for_metabolite_name():
    assert exchange_for_metabolite_name('1-Butanol') == (['EX_1boh_e'], None)

if __name__ == '__main__':
    test_add_heterologous_pathway_me(pathways_to_test[0])
