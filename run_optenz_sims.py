from me_scripts.me_simulation import Production_Envelope
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
        logging.error(type(e), e)
    
def loop():
    # get the setups
    setups = []; setup_dir = '/Users/zaking/data/me/'
    for j in ['14..08.07_optenz_sims.json']:
        with open(join(setup_dir, j), 'r') as f:
            setups = setups + json.load(f)
    # set up the sims
    sims = []
    n = '14..08.08_optenz_production_envelopes_noloops_fixed'
    data_dir = join('/Users/zaking/no-backup-data/me', n)
    solution_dir = join('/Users/zaking/data/me', n)
    for setup in setups:
        if setup['experiment'] != 'GAPD_rpiA,rpiB,mdh,pgl':
            continue
        batch_solution_path = join(solution_dir,
                                   setup['experiment'],
                                   '%s_batch.json' % setup['experiment'])
        # if exists(batch_solution_path): # TODO also check for max/min solutions
        #     logging.info('Skipping completed %s' % setup['experiment'])
        #     continue
        sim = Production_Envelope(setup['experiment'],
                                  data_dir,                                  
                                  solution_dir,
                                  objective_dictionary={setup['objective_reaction']: 1},
                                  nondefault_reaction_bounds=setup['nondefault_reaction_bounds'],
                                  gr_fraction_points=[0, 0.25, 0.5, 0.75, 0.99])
        sims.append(sim)
    run_multiprocessing(sims, run_fn, threads=len(setups), debug_single=True)

if __name__=="__main__":
    loop()
