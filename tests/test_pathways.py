from me_scripts.hindsight.pathways import (add_heterologous_pathway,
                                           get_designs, get_product_dictionary,
                                           get_no_route_exchanges,
                                           exchange_for_metabolite_name,
                                           add_all_heterologous_pathways)
from theseus import setup_model
from theseus.bigg import download_model

import pytest
from pytest import raises


DEBUG = False

models_to_test = ['iJO1366', 'e_coli_core', 'iJR904']
pathways_to_test = get_designs().iteritems()
if DEBUG:
    pathways_to_test = [('3-Methyl-1-Butanol-kivd-ADH2', get_designs()['3-Methyl-1-Butanol-kivd-ADH2'])]
    models_to_test = ['e_coli_core']
no_route_exchanges = get_no_route_exchanges()


@pytest.fixture(scope='session', params=pathways_to_test)
def pathway_tuple(request):
    # return the pathway key
    return request.param


@pytest.fixture(scope='function', params=models_to_test)
def model(request):
    # return the loaded model
    model_id = request.param
    return setup_model(download_model(model_id), 'EX_glc_e', aerobic=False)


def test_add_heterologous_pathway(model, pathway_tuple):
    additions, design = pathway_tuple
    model = add_heterologous_pathway(model, additions)
    # produce some biomass
    biomass = (r for r in model.reactions if r.objective_coefficient != 0).next()
    biomass.lower_bound = 0.1
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
        assert model.solution.f < 25


def test_repeat_additions(model):
    des = get_designs().iterkeys().next()
    m = add_heterologous_pathway(model, des)
    m = add_heterologous_pathway(m, des, ignore_repeats=True)
    with raises(Exception):
        m = add_heterologous_pathway(m, des)


def test_add_heterologous_pathway_core(pathway_tuple):
    model = setup_model(download_model('e_coli_core'), 'EX_glc_e', aerobic=False)
    additions, design = pathway_tuple
    if design[1] is None:
        continue
    m = add_heterologous_pathway(model.copy(), additions)


def test_add_heterologous_pathway_iJR():
    model = download_model('iJR904')
    model = setup_model(model, 'EX_glc_e', aerobic=False)
    for additions, design in designs.iteritems():
        if design[1] is None:
            continue
        m = add_heterologous_pathway(model.copy(), additions)


def test_add_all_heterologous_pathways():
    model = download_model('iJO1366')
    model = add_all_heterologous_pathways(model)
    assert 'btal_c' in model.metabolites
    assert 'HACD1' in model.reactions


def test_exchange_for_metabolite_name():
    assert exchange_for_metabolite_name('1-Butanol', get_product_dictionary()) == 'EX_1boh_e'
