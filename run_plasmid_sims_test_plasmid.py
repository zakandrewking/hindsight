from me_scripts.me_simulation import Single_Exchange_FVA
import json
from me.parallel import run_multiprocessing
import logging
import sys
from os.path import join, exists

# basic setup
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='[%(process)d] %(levelname)s %(message)s')

def run_fn(simulation):
    logging.info('Solving %s' % simulation.experiment_name)
    try:
        simulation.solve()
    except Exception as e:
        logging.error(e)
    
def loop():
    # get the setups
    setups = []; setup_dir = '/Users/zaking/data/hindsight'
    for j in ['15..03.18_plasmid_sims_no_ko.json']:
        with open(join(setup_dir, j), 'r') as f:
            setups = setups + json.load(f)
    # set up the sims
    sims = []
    n = '15..03.18_plasmid_sims_no_ko'
    data_dir = join('/Users/zaking/no-backup-data/hindsight', n)
    solution_dir = join('/Users/zaking/data/hindsight', n)
    model_dir = "/Users/zaking/data/hindsight/15..03.18_plasmids_puc19"
    for setup in setups:
        if setup['experiment'] not in ['Atsumi2008_2_sub_exc']: continue
        # check for the solution
        batch_solution_path = join(solution_dir,
                                   setup['experiment'],
                                   '%s_batch.json' % setup['experiment'])
        if exists(batch_solution_path):
            logging.info('Skipping completed %s' % setup['experiment'])
            continue
        # check for the plasmid. Not using the reduced model
        model_path = join(model_dir, "%s.pickle" % setup['plasmid'])
        if not exists(model_path):
            logging.warn('Skipping experiment %s with missing plasmid %s' % (setup['experiment'],
                                                                             setup['plasmid']))
            continue
        
        sim = EZ_Batch_Maximum_Growth(solution_dir,
                                      setup['experiment'],
                                      solution_file_path=join(solution_dir, setup['experiment']),
                                      nondefault_reaction_bounds=setup['nondefault_reaction_bounds'],
                                      symbolic_model_path=model_path,
                                      use_ssadsc_keffs=True)
                                      
        sim = Single_Exchange_FVA(setup['experiment'],
                                  data_dir,
                                  solution_dir,
                                  {setup['objective_reaction']: 1},
                                  setup['nondefault_reaction_bounds'],
                                  symbolic_model_path=model_path)
        sims.append(sim)
    # run_multiprocessing(sims, run_fn, threads=8, debug_single=False)
    run_fn(sims[0])

if __name__=="__main__":
    loop()
