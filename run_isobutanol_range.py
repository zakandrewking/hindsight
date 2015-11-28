from me.problem_and_solver_classes.me_problem_w_keff_variables import EZ_Batch_Maximum_Growth
import json
from me.parallel import run_multiprocessing
import logging
import sys
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
    n = '14..08.08_isobutanol_range_no-maint'
    data_dir = join('/Users/zaking/no-backup-data/hindsight', n)
    solution_dir = join('/Users/zaking/data/hindsight', n)
    model_dir = '/Users/zaking/data/hindsight/14..07.02_plasmids_puc19_no-maint'
    for setup in setups:
        if setup['experiment'] != 'Atsumi2008_1_sub_exc':
            continue
        # check for the plasmid
        model_path = join(model_dir, '%s_reduced.pickle' % setup['plasmid'])
        if not exists(model_path):
            logging.warn('Skipping experiment %s with missing plasmid %s' % (setup['experiment'],
                                                                             setup['plasmid']))
            continue
        x_range = range(0, 60, 10) + ['max']
        for x in x_range:
            nondefault_reaction_bounds = setup['nondefault_reaction_bounds'].copy()
            if x != 'max':
                nondefault_reaction_bounds['EX_iboh_e'] = (x, x)
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
    run_multiprocessing(sims, run_fn, threads=8, debug_single=False)

if __name__=='__main__':
    loop()
