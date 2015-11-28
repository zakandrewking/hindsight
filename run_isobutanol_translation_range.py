from me.problem_and_solver_classes.me_problem_w_keff_variables import EZ_Batch_Maximum_Growth
from me.problem_and_solver_classes.me_solution import LP_Solution
from me.get_db_cursor import get_db_cursor

import json
from me.parallel import run_multiprocessing
import logging
import sys
import cPickle as pickle
from os.path import join, exists

# basic setup
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='[%(process)d] %(levelname)s %(message)s')

def run_fn(simulation):
    logging.info('Solving %s' % simulation.simulation_directory_path)
    try:
        simulation.solve()
        simulation.dump_problem_and_solution()
        simulation.solution.export_to_json()
    except Exception as e:
        logging.error(e)
    
def loop():
    # get the setups
    setups = []; setup_dir = '/Users/zaking/data/hindsight'
    for j in ['14..06.30_plasmid_sims.json']:
        with open(join(setup_dir, j), 'r') as f:
            setups = setups + json.load(f)
    # set up the sims
    sims = []
    n = '14..08.08_isobutanol_translation_range'
    data_dir = join('/Users/zaking/no-backup-data/hindsight', n)
    solution_dir = join('/Users/zaking/data/hindsight', n)
    model_dir = "/Users/zaking/no-backup-data/hindsight/14..08.08_plasmids_puc19_no-maint_saved"
    
    for setup in setups:
        if setup['experiment'] != 'Atsumi2008_1_sub_exc':
            continue
        # check for the plasmid
        model_path = join(model_dir, "%s_reduced.pickle" % setup['plasmid'])
        if not exists(model_path):
            logging.warn('Skipping experiment %s with missing plasmid %s' % (setup['experiment'],
                                                                             setup['plasmid']))
            continue

        # get the relative fluxes through nonnative reactions
        sol_path = join(data_dir, 'Atsumi2008_1_translation')
        max_solution = LP_Solution.import_from_json(sol_path)
        flux = max_solution.primal_dictionary
        with open(model_path.replace('_reduced', ''), 'r') as f:
            model = pickle.load(f)
        translation = model.plasmid.get_heterologous_peptide_translation_reactions(get_db_cursor())
        elongation = [x for x in translation[0].keys() if 'elongation' in x]
        elongation_fluxes = {k, flux[k] for k in elongation}
        print elongation_fluxes
        break
        
        x_range = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.]
        for x in x_range:
            nondefault_reaction_bounds = setup['nondefault_reaction_bounds'].copy()
            # nondefault_reaction_bounds['EX_iboh_e'] = (x, x)
            exp = '%s_%s' % (setup['experiment'], x)
            sim = EZ_Batch_Maximum_Growth(data_dir,
                                          exp,
                                          solution_file_path=join(solution_dir, exp),
                                          nondefault_reaction_bounds=nondefault_reaction_bounds,
                                          use_ssadsc_keffs=True,
                                          symbolic_model_path=model_path,
                                          recover_lp=True,
                                          continue_solving_from_recover=True)
            sims.append(sim)
    run_multiprocessing(sims, run_fn, threads=1, debug_single=False)

if __name__=="__main__":
    loop()
