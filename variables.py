class NotFoundError(Exception):
    pass

class SetUpModelError(Exception):
    pass

def get_min_biomass():
    return 0.005

def get_max_our():
    return 20
