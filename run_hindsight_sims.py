#!/usr/bin/env python

from me_scripts.me_simulation import Single_Exchange_FVA
import json
from me.parallel import run_multiprocessing
import logging
import sys
from os.path import join, exists

debug = True
max_threads = 8

# basic setup
logging.basicConfig(level=(logging.DEBUG if debug else logging.INFO),
                    stream=sys.stdout,
                    format='[%(process)d] %(levelname)s %(message)s')

def run_fn(simulation):
    logging.info('Solving %s' % simulation.experiment_name)
    try:
    # if True:
        simulation.solve()
    except Exception as e:
        logging.exception(e)
    
def loop():
    # get the setups
    setups = []; setup_dir = '/Users/zaking/data/hindsight'
    for j in ['15..03.18_native_sims.json']:
        with open(join(setup_dir, j), 'r') as f:
            setups = setups + json.load(f)
    # set up the sims
    sims = []
    n = '15..04.30_native_sims_obj-test'
    data_dir = join('/Users/zaking/no-backup-data/hindsight', n)
    solution_dir = join('/Users/zaking/data/hindsight', n)
    for setup in setups:
        if setup['experiment'] not in ['Trinh2008_sub_exc']: continue
        batch_solution_path = join(solution_dir,
                                   setup['experiment'],
                                   '%s_batch.json' % setup['experiment'])
        if exists(batch_solution_path):
            logging.info('Skipping completed %s' % setup['experiment'])
            continue
        sim = Single_Exchange_FVA(setup['experiment'],
                                  data_dir,                                  
                                  solution_dir,
                                  {setup['objective_reaction']: 1},
                                  setup['nondefault_reaction_bounds'])
        sims.append(sim)

    if debug:
        run_fn(sims[0])
    else:
        run_multiprocessing(sims, run_fn, max_threads=max_threads,
                            debug_single=False)

if __name__=="__main__":
    loop()
