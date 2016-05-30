from me.heterologous import (TranscriptionUnit,
                             Peptide,
                             Complex,
                             Plasmid,
                             CatalyzedReaction)
from me.heterologous.utilities import load_plasmid_locs, get_genome_string
from me.heterologous import parameters
from me.parallel import run_multiprocessing
from me.cobra_models import path_for_cobra_model
from me.cobra_models.reduced_model import ReducedModel
from me.get_db_cursor import get_db_cursor

from me_scripts.db.queries import m_to_me_metabolites
from me_scripts.hindsight import get_designs, get_plasmid_translation_fluxes

from multiprocessing import Lock, Manager
import cPickle as pickle
import json
import logging
import sys
import os
from os.path import join, dirname, realpath
import re

# basic setup
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='[%(process)d] %(levelname)s %(message)s')

def main():
    # make the dir
    data_directory = "/Users/zaking/data/hindsight/15..03.18_plasmids_puc19"
    logging.info('Making directory %s' % data_directory)
    try: os.makedirs(data_directory)
    except OSError: pass

    # get the sequence for the plasmid
    logging.info('Getting plasmid sequence')
    dna_sequence = get_genome_string(join(dirname(realpath(__file__)), 'puc19gbk.txt'),
                                     format='genbank')
    
    # Load the E coli genome
    logging.info('Getting genome sequence')
    genome_sequence = get_genome_string(parameters.genome_path)

    # lock for file i/o
    # https://mail.python.org/pipermail/python-list/2009-April/532771.html
    manager = Manager()
    lock = manager.Lock()
    
    # plasmid specs
    plasmid_specs = []
    # logging.info('New model copy for each design')
    for name, design in get_designs().iteritems():
        # skip finished ones
        if os.path.exists(join(data_directory, '%s_reduced.pickle' % name)):
            logging.info('Skipping finished: %s' % name)
            continue
        p = {'name': name,
             'metabolites': design[0],
             'reactions': design[1],
             'bounds': design[3],
             'directory': data_directory,
             'dna_sequence': dna_sequence,
             'genome_sequence': genome_sequence,
             'insertion_site': 397,
             'copy_number': 200,
             'has_bla': False,
             'lock': lock}
        plasmid_specs.append(p)

    logging.info('Making %d plasmids' % len(plasmid_specs))
        
    # now run the slow part
    for ps in plasmid_specs:
        run_fn(ps)
    # run_multiprocessing(plasmid_specs, run_fn, threads=8, debug_single=False)
    
def run_fn(plasmid_spec):
    logging.info("Building plasmid for %s" % plasmid_spec["name"])

    # get the cursor
    db_cursor = get_db_cursor()

    # separate out exchange reactions
    exchange_reactions = {}; reactions = {}
    for name, reaction in plasmid_spec["reactions"].iteritems():
        # conform to me naming
        reaction = {m_to_me_metabolites(k, db_cursor): v for
                    k, v in reaction.iteritems()}

        if name.startswith('EX_'):
            exchange_reactions[name] = reaction
        else:
            reactions[name] = reaction
    
    # get the insertion locs
    n_genes = len(reactions)
    plasmid_sequence, start_stops, added = load_plasmid_locs(plasmid_spec["genome_sequence"],
                                                             plasmid_spec["dna_sequence"],
                                                             n_genes,
                                                             plasmid_spec["insertion_site"])

    # make a plasmid
    plasmid = Plasmid(plasmid_spec["name"],
                      plasmid_sequence,
                      plasmid_spec["copy_number"],
                      parameters.sigma_factors,
                      maintenance_reaction=False)

    if plasmid_spec['has_bla']:
        # add bla antibiotic resistance gene
        tu = TranscriptionUnit('tu_with_bla', 1626+added, 2486+added, '-', 'RpoD_mono', True)
        plasmid.add_tu(tu)
        peptide = Peptide('bla', 1626+added, 2486+added, '-')
        plasmid.add_peptide(peptide)
        # burden of expression beta-lactamase for selection
        plasmid.couple_peptide_to_plasmid_production(peptide_id='bla', alpha=857)

    for name, reaction in exchange_reactions.iteritems():
        # add the reaction
        plasmid.add_chemical_reaction(name, reaction, (0, 1000))
    
    for pep, (name, reaction) in zip(start_stops, reactions.items()):
        # tu
        tu = TranscriptionUnit(name+'_tu', pep[0], pep[1], '+', 'RpoD_mono', True)
        plasmid.add_tu(tu)

        # peptide
        peptide = Peptide(name.replace('_', '')+'pep', pep[0], pep[1], '+')
        plasmid.add_peptide(peptide)

        # complex
        complex = Complex(name+'_complex', {peptide: 1})
        plasmid.add_complex(complex)
        
        try:
            # look for bounds
            lower_bound, upper_bound = plasmid_spec["bounds"][name]
        except (TypeError, KeyError):
            lower_bound, upper_bound = (-1000, 1000)
            
        # add forward and reverse reactions        
        def add_r(n, coeff, b):
            r = CatalyzedReaction(n, coeff)
            r.lower_bound, r.upper_bound = b
            r.set_catalyst(complex, 65)
        add_r(name+'_forward', reaction, (0, upper_bound))
        if lower_bound < 0:
            rev_reaction = {k: -1*v for k, v in reaction.iteritems()}
            add_r(name+'_reverse', rev_reaction, (0, abs(lower_bound)))

    # load model
    lock = plasmid_spec['lock']
    lock.acquire();
    logging.info('Loading model')
    with open(path_for_cobra_model(reduced_model=False)) as f:
        model = pickle.load(f)
    logging.info('Done loading model')
    lock.release();    

    # add to the model
    plasmid.add_to_model(model, save_plasmid_as_attribute=True)
    
    # write model
    outfile = join(plasmid_spec["directory"], "%s.pickle" % plasmid_spec["name"])
    with open(outfile, "wb") as f:
        pickle.dump(model, f)
                
    #reduced model
    reduced_model = ReducedModel(model)
    logging.info("%s original model has %d reactions" % (plasmid_spec["name"],
                                                         len(reduced_model.reactions)))
    reduced_model.identify_cosets()
    logging.info("%s reduced model has %d reactions" % (plasmid_spec["name"],
                                                        len(reduced_model.reactions)))
    
    outfile = join(plasmid_spec["directory"], "%s_reduced.pickle" % plasmid_spec["name"])
    with open(outfile, "wb") as f:
        pickle.dump(reduced_model, f)

if __name__=="__main__":
    main()
