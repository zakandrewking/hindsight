class NotFoundError(Exception):
    pass

class SetUpModelError(Exception):
    pass

min_biomass = 0.005
max_our = 20
default_sur = 10
carbon_uptake_rate = 60
supplement_uptake = 10 # each
growth_coupled_cutoff = 0.15
growth_coupled_h2_cutoff = 2
