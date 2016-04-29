#!/global/homes/z/zking1/me/env/bin/python

from me_scripts.sampling import generate_sampling_problems

from me.problem_and_solver_classes.me_problem_w_keff_variables import EZ_Batch_Maximum_Growth
from me.mpi.mpi_solve import mpi_solve_problems
# from me.parallel import run_multiprocessing

import json
import logging
import sys
from os.path import join, exists, expandvars

# basic setup
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='[%(process)d] %(levelname)s %(message)s')
dry_run = '--dry-run' in sys.argv

def run_fn(simulation):
    logging.info('Solving: %s' % simulation.simulation_directory_path)
    simulation.solve(dry_run=dry_run)
    simulation.dump_problem_and_solution()
    simulation.solution.export_to_json()
    
def load_problems():
    # get the setups
    setups = []; setup_dir = expandvars('$HOME')
    for j in ['14..06.30_sims.json']:
        with open(join(setup_dir, j), 'r') as f:
            setups = setups + json.load(f)
    with open(join(setup_dir,'14..08.04_keff_fits.json'), 'r') as f:
        keff_fits = json.load(f)
    # set up the sims
    n = '14..08.04_Jantama2008_sampling_lognorm_fitted'
    data_dir = join(expandvars('$SCRATCH'), 'hindsight', n)
    solution_dir = join(expandvars('$SCRATCH'), 'hindsight-solutions', n)
    for setup in setups:
        if setup['experiment'] != 'Jantama2008_sub_exc':
            continue
        batch_solution_path = join(solution_dir, setup['experiment'])
        sim = EZ_Batch_Maximum_Growth(data_dir,
                                      setup['experiment'],                                      
                                      use_ssadsc_keffs=True,
                                      use_reduced_model=False,
                                      solution_file_path=batch_solution_path,
                                      nondefault_reaction_bounds=setup['nondefault_reaction_bounds'])
        sims = generate_sampling_problems([sim], 500, distribution='lognormal',
                                          center_for_mu=True, random_seed=None,
                                          start_with=0, nondefault_fit=keff_fits)
        return sims
    raise Exception('Could not find experiment')

if __name__=="__main__":
    sys.exit(mpi_solve_problems(load_problems, run_fn))
    # problems = load_problems()
    # sys.exit(run_multiprocessing(problems, run_fn))
