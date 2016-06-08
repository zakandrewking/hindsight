from theseus import add_pathway
from hindsight.variables import NotFoundError
from minime.solve.algorithms import binary_search, solve_at_growth_rate

no_route_exchanges = [
    # in the byproduct secretion profile, but not in the model
    'EX_2phetoh_e',
    'EX_iamoh_e',
    'EX_2mbtoh_e',
    'EX_3hv_e',
]

EXCHANGE_NAMES = {
    'glucose':  'EX_glc__D_e',
    'xylose':   'EX_xyl__D_e',
    'D-xylose': 'EX_xyl__D_e',
    'glycerol': 'EX_glyc_e',
    'sucrose':  'EX_sucr_e',
    'LB+glucose, sorbitol, and gluconate': {
        'substrates': 'EX_glc__D_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    'glucose+betaine': 'EX_glc__D_e',
    'glucose+yeast extract': {
        'substrates': 'EX_glc__D_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    'glucose+LB': {
        'substrates': 'EX_glc__D_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    'LB': {
        'substrates': 'EX_glc__D_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    'glucose+xylose': ['EX_glc__D_e', 'EX_xyl__D_e'],
    'pyruvate+yeast extract': {
        'substrates': 'EX_pyr_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    'glucose+acetate': ['EX_glc__D_e', 'EX_ac_e'],
    'mannitol': 'EX_mnl_e',
    'D-Glucose+L-Valine+L-Isoleucine+L-Leucine': {
        'substrates': 'EX_glc__D_e',
        'supplements': ['EX_val__L_e', 'EX_leu__L_e', 'EX_ile__L_e']
    },
    '2-Ketoglutarate': 'EX_akg_e',
    '2-Phenylethanol': 'EX_2phetoh_e',
    '3-Methyl-1-butanol': 'EX_iamoh_e',
    '2-Methyl-1-butanol': 'EX_2mbtoh_e',
    '3HV': 'EX_3hv_e',
    '4HB': 'EX_4hdxbld_e',
    'Formate': 'EX_for_e',
    '1-Propanol': 'EX_1poh_e',
    'GBL': 'EX_gbl_e',
    'L-Alanine': 'EX_ala__L_e',
    '14BDO': 'EX_14btd_e',
    'Formate': 'EX_for_e',
    '1-Butanol': 'EX_1boh_e',
    '1-Hexanol': 'EX_1hex_e',
    'Fumarate': 'EX_fum_e',
    'L-Malate': 'EX_mal__L_e',
    'Isobutanol': 'EX_iboh_e',
    'Pyruvate': 'EX_pyr_e',
    '3HB': 'EX_3hb_e',
    'Acetate': 'EX_ac_e',
    'L-Lactate': 'EX_lac__L_e',
    'Ethanol': 'EX_etoh_e',
    'Lactate': 'EX_lac__D_e',
    'D-Lactate': 'EX_lac__D_e',
    'Succinate': 'EX_succ_e',
    'Butyrate': 'EX_but_e',
    'H2': 'EX_h2_e',
    'PLA': 'EX_lac__D_e',
    '3HBcoLA': 'EX_3hb_co_la_e',
    '3HBco3HV': 'EX_3hb_e',
    'PHB': 'EX_3hb_e',
    '(R,R)-2,3-BDO': 'EX_btd__RR_e',
    'meso-2,3-BDO': 'EX_btd__meso_e',
    'methyl ketones': 'EX_2ptone_e',
    'Xylitol': 'EX_xylt_e',
    'Crotonic acid': 'EX_crot_e',
}
def to_list(v):
    return v if isinstance(v, list) else [v]
def exchange_for_metabolite_name(name):
    """Returns (list of substrate exchanges, list of supplement exhanges)."""
    for k, v in EXCHANGE_NAMES.iteritems():
        if k.lower() == name.lower().strip().replace('homo', ''):
            if isinstance(v, dict):
                return to_list(v['substrates']), to_list(v['supplements'])
            else:
                return to_list(v), []
    raise NotFoundError('Could not find %s' % name)


def get_designs():
    return {
        'ldhL': # L-Lactate
        [{},
         { 'LDH_L': {'h_c': 1.0, 'lac__L_c': -1.0, 'nad_c': -1.0,
                     'nadh_c': 1.0, 'pyr_c': 1.0}},
         None, {'EX_lac__L_e': (0, 1000)}],

        # -----------
        # 1-propanol
        # -----------

        # 1-propanol SA203-pSA55-pSA62
        '1_propanol_adh2_kivd_ilva_leu': [
            {
                '1poh_c': {'formula': 'C3H8O', 'name': '1-propanol'},
                '2phetoh_c': {'formula': 'C8H10O', 'name': '2-phenylethanol'},
                'iamoh_c': {'formula': 'C5H12O', 'name': '3-Methyl-1-butanol'},
                '2mbtoh_c': {'formula': 'C5H12O', 'name': '2-Methyl-1-butanol'},
            },
            {
                '2OBUTDC': {'2obut_c': -1, 'h_c': -1, 'ppal_c': 1, 'co2_c': 1},
                '1PDH': {'ppal_c': -1, 'nadh_c': -1, 'h_c': -1, '1poh_c': 1, 'nad_c': 1},
                'EX_1poh_e': {'1poh_c': -1},
                'EX_2phetoh_e': {'2phetoh_c': -1},
                'EX_iamoh_e': {'iamoh_c': -1},
                'EX_2mbtoh_e': {'2mbtoh_c': -1},
            },
            None,
            {
                'EX_1poh_e': (0, 1000),
            },
        ],

        'ilvA, kivD, ADH2, cimA, leuBCD': # 2 pathways Shen2013
        [{'1poh_c': {'formula': 'C3H8O', 'name': '1-propanol'},
          'citmal_c': {'formula': 'C5H6O5', 'name': 'citramalate'},
          'citcon_c': {'formula': 'C5H4O4', 'name': 'citraconate'},
          'betamal_c': {'formula': 'C5H6O5', 'name': 'beta-methyl-D-malate'},},
         {'2OBUTDC': {'2obut_c': -1, 'h_c': -1, 'ppal_c': 1, 'co2_c': 1},
          '1PDH': {'ppal_c': -1, 'nadh_c': -1, 'h_c': -1, '1poh_c': 1, 'nad_c': 1},
          'CIMSY': {'pyr_c': -1, 'accoa_c': -1, 'h2o_c': -1, 'citmal_c': 1, 'coa_c': 1, 'h_c': 1},
          'CIMHYD': {'citmal_c': -1, 'citcon_c': 1, 'h2o_c': 1},
          'BMALHYD': {'citcon_c': -1, 'h2o_c': -1, 'betamal_c': 1},
          'BMALDH': {'betamal_c': -1, 'nad_c': -1, '2obut_c': 1, 'nadh_c': 1, 'co2_c': 1},
          'EX_1poh_e': {'1poh_c': -1}},
         None, {'EX_1poh_e': (0, 1000)}],

        'pdc, adhB':
        [{},
         {'PDC': {'pyr_c': -1, 'h_c': -1, 'acald_c': 1, 'co2_c': 1}},
         None, {'PDC': (0, 1000), 'EX_etoh_e': (0, 1000)}],

        'pyc':
        [{},
         {'PYC': {'pyr_c': -1, 'hco3_c': -1, 'atp_c': -1, 'adp_c': 1,
                  'oaa_c': 1, 'h_c': 1, 'pi_c': 1}},
         None, {'EX_succ_e': (0, 1000)}],

        'pepc': # Native enzyme pepc
        [{}, {}, None, {'EX_succ_e': (0, 1000)}],

        'alaD':
        [{},
         {'ALADH_L': {'ala__L_c': -1, 'nad_c': -1, 'h2o_c': -1,
                      'nh4_c': 1, 'pyr_c': 1, 'nadh_c': 1, 'h_c': 1}},
         None, {'EX_ala__L_e': (0, 1000)}],

        'Ascaris suum Malic Enzyme':  # Native enzyme
        [{}, {}, None, {'EX_succ_e': (0, 1000)}],

        'adh': # R,R butanediol
        [{'actn__R_c': {'formula': 'C4H8O2', 'name': '(R)-3-acetoin'},
          'btd__RR_c': {'formula': 'C4H10O2', 'name': '(R,R)-2,3-butanediol'}},
         {'ACLDC': {'alac__S_c': -1, 'h_c': -1, 'co2_c': 1, 'actn__R_c': 1},
          'BTDD__RR': {'actn__R_c': -1, 'nadph_c': -1, 'h_c': -1, 'nadp_c': 1, 'btd__RR_c': 1},
          'EX_btd__RR_e': {'btd__RR_c': -1}},
         None, {'EX_btd__RR_e': (0, 1000)}],

        'budC': # meso butenediol
        [{'actn__R_c': {'formula': 'C4H8O2', 'name': '(R)-3-acetoin'},
          'btd__meso_c': {'formula': 'C4H10O2', 'name': 'meso-2,3-butanediol'}},
         {'ACLDC': {'alac__S_c': -1, 'h_c': -1, 'co2_c': 1, 'actn__R_c': 1},
          'BTDD__RR': {'actn__R_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'btd__meso_c': 1},
          'EX_btd__meso_e': {'btd__meso_c': -1}},
         None, {'EX_btd__meso_e': (0, 1000)}],

        'sucD, 4hbd, sucA, cat2, 025B': # 1,4-butanediol
        [{'4hbutcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '4-hydroxybutyryl-CoA'},
          'gbl_c': {'formula': 'C4H6O2', 'name': 'gamma-Butyrolactone'},
          '4hdxbld_c': {'formula': 'C4H8O2', 'name': '4-hydroxybutaldehyde'},
          '14btd_c': {'formula': 'C4H10O2', 'name': '1,4-butanediol'}},
         {'AKGDC': {'akg_c': -1, 'h_c': -1, 'co2_c': 1, 'sucsal_c': 1},
          'SUCCALDH': {'succoa_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'coa_c': 1, 'sucsal_c': 1},
          '4HBACT': {'ghb_c': -1, 'accoa_c': -1, '4hbutcoa_c': 1, 'ac_c': 1},
          '4HBTALDDH': {'4hbutcoa_c': -1, 'nadh_c': -1, 'h_c': -1, '4hdxbld_c': 1,
                        'coa_c': 1, 'nad_c': 1},
          'EX_4hdxbld_e': {'4hdxbld_c': -1},
          'GBL_PROD': {'4hbutcoa_c': -1,
                       'gbl_c': 1, 'coa_c': 1},
          'EX_gbl_e': {'gbl_c': -1},
          'BTDP2': {'4hdxbld_c': -1, 'h_c': -1, 'nadh_c': -1, '14btd_c': 1, 'nad_c': 1},
          'EX_14btd_e': {'14btd_c': -1}},
         None, {'AKGDC': (0, 1000), 'EX_4hdxbld_e': (0, 1000),
                'EX_gbl_e': (0, 1000), 'EX_14btd_e': (0, 1000)} ],

        # -----------
        # polymers
        # -----------

        'phbCAB': # 3HBcoHV
        [{'3hbcoa__R_c': {'name': 'D(-)3-hydroxybutyryl-CoA', 'formula': 'C25H38N7O18P3S'},
          '3hb_c': {'formula': 'C4H6O2', 'name': '3HB monomer'},
          '3hv_c': {'formula': 'C5H8O2', 'name': '3HV monomer'}},
         {'PHPB': {'aacoa_c': -1, 'h_c': -1, 'nadph_c': -1, '3hbcoa__R_c': 1, 'nadp_c': 1},
          '3HB_POLYM': {'3hbcoa__R_c': -1, 'coa_c': 1, '3hb_c': 1},
          #'3HV_POLYM': {'3hv_c': 1},
          'EX_3hv_e': {'3hv_c': -1},
          'EX_3hb_e': {'3hb_c': -1}},
         None,
         {'EX_3hv_e': (0, 1000), 'EX_3hb_e': (0, 1000)}],

        'PhaC1Ps6-19, PctCp': # PLA
        [{},
         {},
         None, {'EX_lac__D_e': (0,1000)}],                     # just make D-Lactate

        'PhaC1Ps6-19, PctCp, phaABC': # P(3HB-co-LA)
        [{'3hbcoa__R_c': {'name': 'D(-)3-hydroxybutyryl-CoA',
                          'formula': 'C25H38N7O18P3S'},
          '3hb_co_la_c': {'formula': 'C7H10O4',
                          'name': '3hb_co_la_monomer'}},
         {'PHPB': {'aacoa_c': -1, 'h_c': -1, 'nadph_c': -1,
                   '3hbcoa__R_c': 1, 'nadp_c': 1},
          '3HB_LA_POLYM': {'lac__D_c': -1, '3hbcoa__R_c': -1, 'h_c': -1,
                           'h2o_c': 1, 'coa_c': 1, '3hb_co_la_c': 1},
          'EX_3hb_co_la_e': {'3hb_co_la_c': -1}},
         None, {'EX_3hb_co_la_e': (0, 1000)}],

        'phaCAB-hydrogenase-3-acsA': # TODO
        [{'3hbcoa__R_c': {'name': 'D(-)3-hydroxybutyryl-CoA',
                          'formula': 'C25H38N7O18P3S'},
          '3hb_c': {'formula': 'C4H6O2',
                    'name': '3hb_monomer'}},
         {'PHPB': {'aacoa_c': -1, 'h_c': -1, 'nadph_c': -1,
                   '3hbcoa__R_c': 1, 'nadp_c': 1},
          '3HB_POLYM': {'3hbcoa__R_c': -1, 'coa_c': 1, '3hb_c': 1},
          'EX_3hb_e': {'3hb_c': -1}},
         None, {'EX_3hb_e': (0, 1000)}],

        # -----------
        # isobutanol
        # -----------

        'isobutanol_als_ilv_kdc_adh': [
            {
                'ibald_c': {'formula': 'C4H8O', 'name': 'isobutyraldehyde'},
                'iboh_c': {'formula': 'C4H10O', 'name': 'isobutanol'},
                '2phetoh_c': {'formula': 'C8H10O', 'name': '2-phenylethanol'},
                'iamoh_c': {'formula': 'C5H12O', 'name': '3-Methyl-1-butanol'},
                '2mbtoh_c': {'formula': 'C5H12O', 'name': '2-Methyl-1-butanol'},
            },
            {
                '3MOBDC': {'3mob_c': -1, 'h_c': -1, 'ibald_c': 1, 'co2_c': 1},
                'IBDH': {'ibald_c': -1, 'nadh_c': -1, 'h_c': -1, 'iboh_c': 1, 'nad_c': 1},
                'EX_iboh_e': {'iboh_c': -1},
                'EX_2phetoh_e': {'2phetoh_c': -1},
                'EX_iamoh_e': {'iamoh_c': -1},
                'EX_2mbtoh_e': {'2mbtoh_c': -1},
            },
            None,
            {
                'EX_iboh_e': (0,1000),
            }
        ],

        'isobutanol_als_ilv_nadh_kdc_adh': [
            {
                'ibald_c': {'formula': 'C4H8O', 'name': 'isobutyraldehyde'},
                'iboh_c': {'formula': 'C4H10O', 'name': 'isobutanol'},
                '2phetoh_c': {'formula': 'C8H10O', 'name': '2-phenylethanol'},
                'iamoh_c': {'formula': 'C5H12O', 'name': '3-Methyl-1-butanol'},
                '2mbtoh_c': {'formula': 'C5H12O', 'name': '2-Methyl-1-butanol'},
            },
            {
                'KARA1x': {'h_c': 1.0, 'nad_c': -1.0, 'nadh_c': 1.0, 'alac__S_c': 1.0, '23dhmb_c': -1.0},
                '3MOBDC': {'3mob_c': -1, 'h_c': -1, 'ibald_c': 1, 'co2_c': 1},
                'IBDH': {'ibald_c': -1, 'nadh_c': -1, 'h_c': -1, 'iboh_c': 1, 'nad_c': 1},
                'EX_iboh_e': {'iboh_c': -1},
                'EX_2phetoh_e': {'2phetoh_c': -1},
                'EX_iamoh_e': {'iamoh_c': -1},
                'EX_2mbtoh_e': {'2mbtoh_c': -1},
            },
            None,
            {
                'EX_iboh_e': (0,1000),
            }
        ],

        # 'ilvIHCD, PDC6, ADH2, kivd': # Isobutanol
        # [{'ibald_c': {'formula': 'C4H8O', 'name': 'isobutyraldehyde'},
        #   'iboh_c': {'formula': 'C4H10O', 'name': 'isobutanol'},
        #   '2petoh_c': {'formula': 'C8H10O', 'name': '2-Phenylethanol'},
        #   '3m1but_c': {'formula': 'C5H12O', 'name': '3-Methyl-1-butanol'},
        #   '2m1but_c': {'formula': 'C5H12O', 'name': '2-Methyl-1-butanol'}},
        #  {'3MOBDC': {'3mob_c': -1, 'h_c': -1, 'ibald_c': 1, 'co2_c': 1},
        #   'IBDH': {'ibald_c': -1, 'nadh_c': -1, 'h_c': -1, 'iboh_c': 1, 'nad_c': 1},
        #   'EX_2petoh_e': {'2petoh_c': -1},
        #   'EX_3m1but_e': {'3m1but_c': -1},
        #   'EX_2m1but_e': {'2m1but_c': -1},
        #   'EX_iboh_e': {'iboh_c': -1}},
        #  None,
        #  {'EX_2m1but_e': (0 ,1000), 'EX_3m1but_e': (0 ,1000), 'EX_2petoh_e': (0, 1000),
        #   'EX_iboh_e': (0, 1000)}],

        # ----------
        # 3-Methyl-1-butanol
        # ----------

        '3_methyl_1_butanol_kivd_adh2': [
            {
                '3mbald_c': {'formula': 'C5H10O', 'name': '3-Methylbutanal'},
                'iamoh_c': {'formula': 'C5H12O', 'name': '3-Methyl-1-butanol'},
            },
            {
                '4MOBDC': {'4mop_c': -1, 'h_c': -1, '3mbald_c': 1, 'co2_c': 1},
                '3M1BDH': {'3mbald_c': -1, 'nadh_c': -1, 'h_c': -1, 'iamoh_c': 1, 'nad_c': 1},
                'EX_iamoh_e': {'iamoh_c': -1}
            },
            None,
            {
                'EX_iamoh_e': (0, 1000)
            }
        ],

        # -----------
        # 1-butanol, 1-hexanol
        # -----------

        '1_butanol_adh2_kivd_ilva_leu': [
            {
                '2kv_c': {'formula': 'C5H7O3', 'name': '2-ketovalerate'},
                '1boh_c': {'formula': 'C4H10O', 'name': '1-butanol'},
                '2phetoh_c': {'formula': 'C8H10O', 'name': '2-phenylethanol'},
            },
            {
                # http://ecocyc.com/ECOLI/NEW-IMAGE?type=PATHWAY&object=LEUSYN-PWY
                'ILV_PATHWAY': {'2obut_c': -1, 'h2o_c': -1,  'accoa_c': -1, 'nad_c': -1,
                                'coa_c': 1, 'h_c': 1, 'nadh_c': 1, 'co2_c': 1, '2kv_c': 1,},
                '2KVDC': {'2kv_c': -1, 'h_c': -1, 'btal_c': 1, 'co2_c': 1},
                '1BDH': {'btal_c': -1, 'nadh_c': -1, 'h_c': -1, '1boh_c': 1, 'nad_c': 1},
                'EX_1boh_e': {'1boh_c': -1},
                'EX_2phetoh_e': {'2phetoh_c': -1},
            },
            None,
            {
                'EX_1boh_e': (0, 1000),
            }
        ],

        'bcd etfAB': # 1-butanol
        [{'btcoa_c': {'formula': 'C25H38N7O17P3S', 'name': 'butyryl-CoA'}, # NOTE already in iJO1366
          '1boh_c': {'formula': 'C4H10O', 'name': '1-butanol'}},
         {'B2COAR': {'b2coa_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'btcoa_c': 1},
          'BTALDH': {'nadh_c': -1, 'h_c': -1, 'btcoa_c': -1, 'btal_c': 1,
                     'nad_c': 1, 'coa_c': 1},
          '1BDH': {'btal_c': -1, 'nadh_c': -1, 'h_c': -1, '1boh_c': 1,
                   'nad_c': 1},
          'EX_1boh_e': {'1boh_c': -1}},
         None, {'EX_1boh_e': (0,1000)}],

        'Bktb bcd etfAB': # 1-hexanol
        [{'1hex_c': {'formula': 'C6H14O', 'name': '1-hexanol'},
          'hxal_c': {'formula': 'C6H12O', 'name': '1-hexanal'},
          'hxcoa_c': {'formula': 'C27H42N7O17P3S', 'name': 'hexanoyl-CoA'}},
         {'HX2COAR': {'hx2coa_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'hxcoa_c': 1},
          'HXALDH': {'nadh_c': -1, 'h_c': -1, 'hxcoa_c': -1, 'hxal_c': 1,
                     'nad_c': 1, 'coa_c': 1},
          '1HDH': {'hxal_c': -1, 'nadh_c': -1, 'h_c': -1, '1hex_c': 1,
                   'nad_c': 1},
          'EX_1hex_e': {'1hex_c': -1}},
         None, {'EX_1hex_e': (0, 1000)}],

          'atoBEC adhE2CA crtCA bhdCA terTD fdhCB': # 1-butanol
        [{'btcoa_c': {'formula': 'C25H38N7O17P3S', 'name': 'butyryl-CoA'}, # NOTE already in iJO1366
          '1boh_c': {'formula': 'C4H10O', 'name': '1-butanol'}},
         {'B2COAR': {'b2coa_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'btcoa_c': 1},
          'BTALDH': {'nadh_c': -1, 'h_c': -1, 'btcoa_c': -1, 'btal_c': 1, 'nad_c': 1, 'coa_c': 1},
          '1BDH': {'btal_c': -1, 'nadh_c': -1, 'h_c': -1, '1boh_c': 1, 'nad_c': 1},
          'EX_1boh_e': {'1boh_c': -1}},
         None, {'EX_1boh_e': (0, 1000)}],

        # --------------------
        # butyrate, propanoate
        # --------------------

        # Lim2013-sy, Saini2014-br, Volker2014-pd
        'butyrate_crt_ter': [
             {
                 'btcoa_c': {'formula': 'C25H38N7O17P3S', 'name': 'butyryl-CoA'}, # NOTE already in iJO1366 and iAF1260(b)
             },
             {
                 'B2COAR': {'b2coa_c': -1, 'nadh_c': -1, 'h_c': -1, 'nad_c': 1, 'btcoa_c': 1},
                 'BUTTH': {'h2o_c': -1, 'btcoa_c': -1, 'but_c': 1, 'h_c': 1, 'coa_c': 1},
             },
             None,
             {
                 'EX_but_e': (0, 1000),
             }
         ],

        # Volker2014-pd 12ppd__R_c -R02376-> propionaldehyde -R09097/MMSAD2-> propionyl-coa -??-> 3-keto-pentanoate -??-> propionate
        # 'propanoate': [
        #     {},
        #     {
        #         'PPDDH': {'12ppd__R_c': -1, 'ppal_c': 1},
        #         'MMSAD2': {'coa_c': -1, 'nad_c': -1, 'ppal_c': -1,  'h_c': 1, 'nadh_c': 1, 'ppcoa_c': 1},
        #     },
        #     None,
        #     {
        #         'EX_ppa_e': (0, 1000),
        #     }
        # ],

        # -----------
        # H2, butanone
        # -----------

        # already in iJO
        'fhl': [{}, {}, None, {}],

        'cHis2A-shmks2-shmks1': [
            {
                '2ptone_c': {'formula': 'C5H10O', 'name': '2-pentatone'},
                'keto_acid_c': {'formula': 'C6H10O3', 'name': '6C keto acid'},
            },
            {
                'THE': {'3ohexACP_c': -1, 'h2o_c': -1, 'keto_acid_c': 1, 'ACP_c': 1},
                'MKS': {'keto_acid_c': -1, '2ptone_c': 1, 'co2_c': 1},
                'EX_2ptone_e': {'2ptone_c': -1},
            },
            None,
            {
                'EX_2ptone_e': (0, 1000),
            }
        ],

        # ----------------
        # Xylitol
        # ----------------

        'xylitol_xr': [
            {
                'xylt_c': {'formula': 'C5H12O5', 'name': 'Xylitol'},
            },
            {
                'XYLR': {'h_c': -1, 'nadph_c': -1, 'xyl__D_c': -1, 'nadp_c': 1, 'xylt_c': 1},
                'EX_xylt_e': {'xylt_c': -1},
            },
            None,
            {
                'EX_xylt_e': (0, 1000)
            },
        ],

        'crotonic_acid': [
            {
                'crot_c': {'formula': 'C4H6O2', 'name': 'Crotonic acid'},
            },
            {
                'CROT': {'b2coa_c': -1, 'h2o_c': -1, 'crot_c': 1, 'coa_c': 1},
                'EX_crot_e': {'crot_c': -1},
            },
            None,
            {
                'EX_crot_e': (0, 1000),
            },
        ],
    }

def get_core_designs():
    return  {
        'ldhL':
        [{'lac__L_c': {'formula': 'C3H5O3', 'name': 'L-Lactate'},
          'lac__L_e': {'formula': 'C3H5O3', 'name': 'L-Lactate'}},
         {'EX_lac__L_e': {'lac__L_e': -1.0},
          'L__LACt2rpp': {'h_c': 1.0, 'h_e': -1.0, 'lac__L_c': 1.0,
                          'lac__L_e': -1.0}},
         None, {'EX_lac__L_e': (0,1000)}],

        'pyc':
        [{'hco3_c': {'formula': 'CHO3', 'name': 'Bicarbonate'}},
         {'HCO3E': {'co2_c': -1.0, 'h2o_c': -1.0,
                    'h_c': 1.0, 'hco3_c': 1.0}},
         None, None],

        'bcd etfAB':
        [{'3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0},
          'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                    'nad_c': 1.0, 'nadh_c': -1.0}},
         None, None],

        'Bktb bcd etfAB':  # 1-hexanol
        [{'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          '3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          '3ohcoa_c': {'formula': 'C27H40N7O18P3S', 'name': '3-Oxohexanoyl-CoA'},
          '3hhcoa_c': {'formula': 'C27H42N7O18P3S', 'name': '(S)-3-Hydroxyhexanoyl-CoA'},
          'hxcoa_c': {'formula': 'C27H42N7O17P3S', 'name': 'hexanoyl-CoA'},
          'hx2coa_c': {'formula': 'C27H40N7O17P3S', 'name': 'trans-Hex-2-enoyl-CoA'},
          'fad_c': {'formula': 'C27H31N9O15P2', 'name': 'Flavin adenine dinucleotide oxidized'},
          'fadh2_c': {'formula': 'C27H33N9O15P2', 'name': 'Flavin adenine dinucleotide reduced'},
          'btcoa_c': {'formula': 'C25H38N7O17P3S', 'name': 'butyryl-CoA'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'nadh_c': -1.0},
          'ACACT2r': {'accoa_c': -1.0, 'btcoa_c': -1.0, '3ohcoa_c': 1.0, 'coa_c': 1.0},
          'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'ACOAD1f': {'btcoa_c': -1.0, 'b2coa_c': 1.0, 'fadh2_c': 1.0, 'fad_c': -1.0},
          'ACOAD2f': {'fadh2_c': 1.0, 'hxcoa_c': -1.0, 'fad_c': -1.0, 'hx2coa_c': 1.0},
          'HACD2': {'3hhcoa_c': 1.0, '3ohcoa_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'nadh_c': -1.0},
          'ECOAH2': {'3hhcoa_c': -1.0, 'hx2coa_c': 1.0, 'h2o_c': 1.0}},
         None, None],

        # Isobutanol

        'isobutanol_als_ilv_kdc_adh':
        [{'23dhmb_c': {'formula': 'C5H9O4', 'name': '(R)-2,3-Dihydroxy-3-methylbutanoate'},
          '3mob_c': {'formula': 'C5H7O3', 'name': '3-Methyl-2-oxobutanoate'},
          'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0},
          'DHAD1': {'23dhmb_c': -1.0, '3mob_c': 1.0, 'h2o_c': 1.0},
          'KARA1': {'23dhmb_c': -1.0, 'alac__S_c': 1.0, 'h_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0}},
         None, None],

        'isobutanol_als_ilv_nadh_kdc_adh':
        [{'23dhmb_c': {'formula': 'C5H9O4', 'name': '(R)-2,3-Dihydroxy-3-methylbutanoate'},
          '3mob_c': {'formula': 'C5H7O3', 'name': '3-Methyl-2-oxobutanoate'},
          'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0},
          'DHAD1': {'23dhmb_c': -1.0, '3mob_c': 1.0, 'h2o_c': 1.0},
          'KARA1': {'23dhmb_c': -1.0, 'alac__S_c': 1.0, 'h_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0}},
         None, None],

        # 3-Methyl-1-butanol

        '3_methyl_1_butanol_kivd_adh2': [
            {
                '23dhmb_c': {'formula': 'C5H9O4', 'name': '(R)-2,3-Dihydroxy-3-methylbutanoate'},
                '3mob_c': {'formula': 'C5H7O3', 'name': '3-Methyl-2-oxobutanoate'},
                'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'},
                '3c3hmp_c': {'formula': 'C7H10O5', 'name': '3-Carboxy-3-hydroxy-4-methylpentanoate'},
                '2ippm_c': {'formula': 'C7H8O4', 'name': '2-Isopropylmaleate'},
                '3c2hmp_c': {'formula': 'C7H10O5', 'name': '3-Carboxy-2-hydroxy-4-methylpentanoate'},
                '3c4mop_c': {'formula': 'C7H8O5', 'name': '3-Carboxy-4-methyl-2-oxopentanoate'},
                '4mop_c': {'formula': 'C6H9O3', 'name': '4-Methyl-2-oxopentanoate'},
            },
            {
                'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0},
                'DHAD1': {'23dhmb_c': -1.0, '3mob_c': 1.0, 'h2o_c': 1.0},
                'KARA1': {'23dhmb_c': -1.0, 'alac__S_c': 1.0, 'h_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0},
                'IPPS': {'3mob_c': -1, 'accoa_c': -1, 'h2o_c': -1, '3c3hmp_c': 1, 'coa_c': 1, 'h_c': 1},
                'IPPMIb': {'2ippm_c': -1, 'h2o_c': -1, '3c3hmp_c': 1},
                'IPPMIa': {'3c2hmp_c': -1, '2ippm_c': 1, 'h2o_c': 1},
                'IPMD': {'3c2hmp_c': -1, 'nad_c': -1, '3c4mop_c': 1, 'h_c': 1, 'nadh_c': 1},
                'OMCDC': {'3c4mop_c': -1, 'h_c': -1, '4mop_c': 1, 'co2_c': 1},
            },
            None,
            None
        ],

        'sucD, 4hbd, sucA, cat2, 025B':
        [{'sucsal_c': {'formula': 'C4H5O3', 'name': 'Succinic semialdehyde'},
          'ghb_c': {'formula': 'C4H7O3', 'name': 'gamma-hydroxybutyrate'},
          '4abut_c': {'formula': 'C4H9NO2', 'name': '4-Aminobutanoate'},
          'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'ABTA': {'akg_c': -1.0, 'glu__L_c': 1.0, 'sucsal_c': 1.0, '4abut_c': -1.0},
          'GHBDHx': {'sucsal_c': -1.0, 'nadh_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'ghb_c': 1.0}},
         None, None],

        'atoBEC adhE2CA crtCA bhdCA terTD fdhCB':
        [{'3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0},
          'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                    'nad_c': 1.0, 'nadh_c': -1.0}},
         None, None],

        'ilvA, kivD, ADH2, cimA, leuBCD':
        [{'phom_c': {'formula': 'C4H8NO6P', 'name': 'O-Phospho-L-homoserine'},
          'aspsa_c': {'formula': 'C4H7NO3', 'name': 'L-Aspartate 4-semialdehyde'},
          'asp__L_c': {'formula': 'C4H6NO4', 'name': 'L-Aspartate'},
          '2obut_c': {'formula': 'C4H5O3', 'name': '2-Oxobutanoate'},
          'hom__L_c': {'formula': 'C4H9NO3', 'name': 'L-Homoserine'},
          'thr__L_c': {'formula': 'C4H9NO3', 'name': 'L-Threonine'},
          '4pasp_c': {'formula': 'C4H6NO7P', 'name': '4-Phospho-L-aspartate'},
          'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'}},
         {'HSK': {'hom__L_c': -1.0, 'phom_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0, 'h_c': 1.0},
          'ASPK': {'asp__L_c': -1.0, '4pasp_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0},
          'HSDy': {'hom__L_c': -1.0, 'aspsa_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0},
          'ASAD': {'aspsa_c': -1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0, 'pi_c': -1.0, '4pasp_c': 1.0},
          'THRD_L': {'thr__L_c': -1.0, '2obut_c': 1.0, 'nh4_c': 1.0},
          'ASPTA': {'akg_c': -1.0, 'asp__L_c': -1.0, 'glu__L_c': 1.0, 'oaa_c': 1.0},
          'THRS': {'phom_c': -1.0, 'thr__L_c': 1.0, 'pi_c': 1.0, 'h2o_c': -1.0}},
         None, None],

        'butyrate_crt_ter':
        [{'3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'},
          'but_c': {'formula': 'C4H7O2', 'name': 'Butyrate (n-C4:0)'},
          'but_e': {'formula': 'C4H7O2', 'name': 'Butyrate (n-C4:0)'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0},
          'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                    'nad_c': 1.0, 'nadh_c': -1.0},
          'BUTtex': {'but_c': -1.0, 'but_e': 1.0},
          'EX_but_e': {'but_e': -1.0}},
         None, None],

        'alaD':
        [{'ala__L_c': {'formula': 'C3H7NO2', 'name': 'L-Alanine'},
          'ala__L_e': {'formula': 'C3H7NO2', 'name': 'L-Alanine'}},
         {'ALAt2pp': {'ala__L_c': 1.0, 'h_c': 1.0, 'h_e': -1.0, 'ala__L_e': -1.0},
          'EX_ala__L_e': {'ala__L_e': -1.0}},
         None, None],

        # 'ilvIHCD, PDC6, ADH2, kivd':
        # [{'23dhmb_c': {'formula': 'C5H9O4', 'name': '(R)-2,3-Dihydroxy-3-methylbutanoate'},
        #   '3mob_c': {'formula': 'C5H7O3', 'name': '3-Methyl-2-oxobutanoate'},
        #   'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
        #  {'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0},
        #   'DHAD1': {'23dhmb_c': -1.0, '3mob_c': 1.0, 'h2o_c': 1.0},
        #   'KARA1': {'23dhmb_c': -1.0, 'alac__S_c': 1.0, 'h_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0}},
        #  None, None],

        '1_propanol_adh2_kivd_ilva_leu':
        [{'phom_c': {'formula': 'C4H8NO6P', 'name': 'O-Phospho-L-homoserine'},
          'aspsa_c': {'formula': 'C4H7NO3', 'name': 'L-Aspartate 4-semialdehyde'},
          'asp__L_c': {'formula': 'C4H6NO4', 'name': 'L-Aspartate'},
          '2obut_c': {'formula': 'C4H5O3', 'name': '2-Oxobutanoate'},
          'hom__L_c': {'formula': 'C4H9NO3', 'name': 'L-Homoserine'},
          'thr__L_c': {'formula': 'C4H9NO3', 'name': 'L-Threonine'},
          '4pasp_c': {'formula': 'C4H6NO7P', 'name': '4-Phospho-L-aspartate'},
          'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
        {'HSK': {'hom__L_c': -1.0, 'phom_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0, 'h_c': 1.0},
         'ASPK': {'asp__L_c': -1.0, '4pasp_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0},
         'HSDy': {'hom__L_c': -1.0, 'aspsa_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0},
         'ASAD': {'aspsa_c': -1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0, 'pi_c': -1.0, '4pasp_c': 1.0},
         'THRD_L': {'thr__L_c': -1.0, '2obut_c': 1.0, 'nh4_c': 1.0},
         'ASPTA': {'akg_c': -1.0, 'asp__L_c': -1.0, 'glu__L_c': 1.0, 'oaa_c': 1.0},
         'THRS': {'phom_c': -1.0, 'thr__L_c': 1.0, 'pi_c': 1.0, 'h2o_c': -1.0}},
         None, None],

        '1_butanol_adh2_kivd_ilva_leu':
        [{'phom_c': {'formula': 'C4H8NO6P', 'name': 'O-Phospho-L-homoserine'},
          'aspsa_c': {'formula': 'C4H7NO3', 'name': 'L-Aspartate 4-semialdehyde'},
          'asp__L_c': {'formula': 'C4H6NO4', 'name': 'L-Aspartate'},
          '2obut_c': {'formula': 'C4H5O3', 'name': '2-Oxobutanoate'},
          'hom__L_c': {'formula': 'C4H9NO3', 'name': 'L-Homoserine'},
          'thr__L_c': {'formula': 'C4H9NO3', 'name': 'L-Threonine'},
          '4pasp_c': {'formula': 'C4H6NO7P', 'name': '4-Phospho-L-aspartate'},
          'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
         {'HSK': {'hom__L_c': -1.0, 'phom_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0, 'h_c': 1.0},
          'ASPK': {'asp__L_c': -1.0, '4pasp_c': 1.0, 'adp_c': 1.0, 'atp_c': -1.0},
          'HSDy': {'hom__L_c': -1.0, 'aspsa_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0},
          'ASAD': {'aspsa_c': -1.0, 'nadp_c': -1.0, 'nadph_c': 1.0, 'h_c': 1.0, 'pi_c': -1.0, '4pasp_c': 1.0},
          'THRD_L': {'thr__L_c': -1.0, '2obut_c': 1.0, 'nh4_c': 1.0},
          'ASPTA': {'akg_c': -1.0, 'asp__L_c': -1.0, 'glu__L_c': 1.0, 'oaa_c': 1.0},
          'THRS': {'phom_c': -1.0, 'thr__L_c': 1.0, 'pi_c': 1.0, 'h2o_c': -1.0}},
         None, None],

        'budC':
        [{'23dhmb_c': {'formula': 'C5H9O4', 'name': '(R)-2,3-Dihydroxy-3-methylbutanoate'},
          '3mob_c': {'formula': 'C5H7O3', 'name': '3-Methyl-2-oxobutanoate'},
          'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0},
          'DHAD1': {'23dhmb_c': -1.0, '3mob_c': 1.0, 'h2o_c': 1.0},
          'KARA1': {'23dhmb_c': -1.0, 'alac__S_c': 1.0, 'h_c': 1.0, 'nadp_c': -1.0, 'nadph_c': 1.0}},
         None, None],

        'PhaC1Ps6-19, PctCp, phaABC':
        [{'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0}},
         None, None],

        'adh':
        [{'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'ACLS': {'alac__S_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0, 'pyr_c': -2.0}},
         None, None],

        'phbCAB':
        [{'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0}},
         None, None],

        'phaCAB-hydrogenase-3-acsA':
        [{'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'}},
         {'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0}},
         None, None],

        'fhl':
        [{'h2_c': {'formula': 'H2', 'name': 'hydrogen'}},
         {'FHL': {'for_c': -1, 'h_c': -1, 'h2_c': 1, 'co2_c': 1},
          'EX_h2_e': {'h2_c': -1}},
         None, {'EX_h2_e': (0, 1000)}],

        'cHis2A-shmks2-shmks1':
        [{'accoa_c': {'formula': 'C23H34N7O17P3S', 'name': 'Acetyl-CoA'},
          'malcoa_c': {'formula': 'C24H33N7O19P3S', 'name': 'Malonyl-CoA'},
          'malACP_c': {'formula': 'C14H22N2O10PRS', 'name': 'Malonyl-[acyl-carrier protein]'},
          'actACP_c': {'formula': 'C15H25N2O9PRS', 'name': 'Acetoacetyl-ACP'},
          '3haACP_c': {'formula': 'C15H27N2O9PRS', 'name': '(3R)-3-Hydroxyacyl-[acyl-carrier protein]'},
          'but2eACP_c': {'formula': 'C15H25N2O8PRS', 'name': 'But-2-enoyl-[acyl-carrier protein]'},
          'butACP_c': {'formula': 'C15H27N2O8PRS', 'name': 'Butyryl-ACP (n-C4:0ACP)'},
          '3ohexACP_c': {'formula': 'C17H29N2O9PRS', 'name': '3-Oxohexanoyl-[acyl-carrier protein]'},
          'hco3_c': {'formula': 'CHO3', 'name': 'Bicarbonate'},
          'ACP_c': {'formula': 'C11H21N2O7PRS', 'name': 'acyl carrier protein'}},
         {'ACCOAC': {'atp_c': -1.0, 'malcoa_c': 1.0, 'adp_c': 1.0,
                     'pi_c': 1.0, 'accoa_c': -1.0, 'hco3_c': -1.0,
                     'h_c': 1.0},
          'MCOATA': {'coa_c': 1.0, 'malcoa_c': -1.0, 'malACP_c': 1.0,
                     'ACP_c': -1.0},
          'KAS15': {'accoa_c': -1.0, 'malACP_c': -1.0, 'actACP_c': 1.0,
                    'coa_c': 1.0, 'co2_c': 1.0, 'h_c': -1.0},
          '3OAR40': {'nadp_c': 1.0, '3haACP_c': 1.0, 'h_c': -1.0,
                     'nadph_c': -1.0, 'actACP_c': -1.0},
          '3HAD40': {'3haACP_c': -1.0, 'but2eACP_c': 1.0, 'h2o_c': 1.0},
          'EAR40x': {'nadh_c': -1.0, 'nad_c': 1.0, 'but2eACP_c': -1.0,
                     'butACP_c': 1.0, 'h_c': -1.0},
          '3OAS60': {'butACP_c': -1.0, 'malACP_c': -1.0,
                     '3ohexACP_c': 1.0, 'co2_c': 1.0, 'ACP_c': 1.0,
                     'h_c': -1.0},
          'HCO3E': {'h2o_c': -1.0, 'co2_c': -1.0, 'h_c': 1.0, 'hco3_c': 1.0}},
         None, {'KAS15': (0, 1000)}],

        # ----------------
        # Xylitol
        # ----------------

        'xylitol_xr': [
            {
                'xyl__D_c': {'formula': 'C5H10O5', 'name': 'D-Xylose'},
            },
            {
                'EX_xyl__D_e': {'xyl__D_c': -1},
            },
            None,
            {
                'EX_xyl__D_e': (0, 1000),
            },
        ],

        # ----------------
        # Xylitol
        # ----------------

        'crotonic_acid': [
            {
                '3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
                'aacoa_c': {'formula': 'C25H36N7O18P3S', 'name': 'Acetoacetyl-CoA'},
                'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
            },
            {
                'ACACT1r': {'aacoa_c': 1.0, 'accoa_c': -2.0, 'coa_c': 1.0},
                'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
                'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                          'nad_c': 1.0, 'nadh_c': -1.0},
            },
            None,
            None,
        ],
    }

def get_iJR904_designs():
    return {
        'bcd etfAB':
        [{'3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
         {'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                    'nad_c': 1.0, 'nadh_c': -1.0}},
         None, None],

        'Bktb bcd etfAB':  # 1-hexanol
        [{'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          '3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          '3ohcoa_c': {'formula': 'C27H40N7O18P3S', 'name': '3-Oxohexanoyl-CoA'},
          '3hhcoa_c': {'formula': 'C27H42N7O18P3S', 'name': '(S)-3-Hydroxyhexanoyl-CoA'},
          'hxcoa_c': {'formula': 'C27H42N7O17P3S', 'name': 'hexanoyl-CoA'},
          'hx2coa_c': {'formula': 'C27H40N7O17P3S', 'name': 'trans-Hex-2-enoyl-CoA'},
          'btcoa_c': {'formula': 'C25H38N7O17P3S', 'name': 'butyryl-CoA'}},
         {'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'nadh_c': -1.0},
          'ACACT2r': {'accoa_c': -1.0, 'btcoa_c': -1.0, '3ohcoa_c': 1.0, 'coa_c': 1.0},
          'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'ACOAD1f': {'btcoa_c': -1.0, 'b2coa_c': 1.0, 'fadh2_c': 1.0, 'fad_c': -1.0},
          'ACOAD2f': {'fadh2_c': 1.0, 'hxcoa_c': -1.0, 'fad_c': -1.0, 'hx2coa_c': 1.0},
          'HACD2': {'3hhcoa_c': 1.0, '3ohcoa_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'nadh_c': -1.0},
          'ECOAH2': {'3hhcoa_c': -1.0, 'hx2coa_c': 1.0, 'h2o_c': 1.0}},
         None, None],

        'sucD, 4hbd, sucA, cat2, 025B':
        [{'ghb_c': {'formula': 'C4H7O3', 'name': 'gamma-hydroxybutyrate'},
          'alac__S_c': {'formula': 'C5H7O4', 'name': '(S)-2-Acetolactate'}},
         {'GHBDHx': {'sucsal_c': -1.0, 'nadh_c': -1.0, 'h_c': -1.0,
                     'nad_c': 1.0, 'ghb_c': 1.0}},
         None, None],

        'atoBEC adhE2CA crtCA bhdCA terTD fdhCB':
        [{'3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
          'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
          'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'}},
         {'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
          'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                    'nad_c': 1.0, 'nadh_c': -1.0}},
         None, None],

        'butyrate_crt_ter': [
            {
                '3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
                'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
                'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'},
                'but_c': {'formula': 'C4H7O2', 'name': 'Butyrate (n-C4:0)'}
            },
            {
                'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
                'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                          'nad_c': 1.0, 'nadh_c': -1.0}
            },
            None,
            None
        ],

        'alaD':
        [{'ala__L_c': {'formula': 'C3H7NO2', 'name': 'L-Alanine'}},
         {},
         None, None],

        '1_propanol_adh2_kivd_ilva_leu':
        [{'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'},
          'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'}},
         {},
         None, None],

        '1_butanol_adh2_kivd_ilva_leu':
        [{'btal_c': {'formula': 'C4H8O', 'name': 'Butanal'},
          'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'}},
         {},
         None, None],

        'ilvA, kivD, ADH2, cimA, leuBCD':
        [{'ppal_c': {'formula': 'C3H6O', 'name': 'Propanal'}},
         {},
         None, None],

        'fhl':
        [{},
         {'EX_h2_e': {'h2_c': -1}},
         None, {'EX_h2_e': (0, 1000)}],

        'cHis2A-shmks2-shmks1':
        [{'3haACP_c': {'formula': 'C15H27N2O9PRS', 'name': '(3R)-3-Hydroxyacyl-[acyl-carrier protein]'},
          'but2eACP_c': {'formula': 'C15H25N2O8PRS', 'name': 'But-2-enoyl-[acyl-carrier protein]'},
          'butACP_c': {'formula': 'C15H27N2O8PRS', 'name': 'Butyryl-ACP (n-C4:0ACP)'},
          '3ohexACP_c': {'formula': 'C17H29N2O9PRS', 'name': '3-Oxohexanoyl-[acyl-carrier protein]'}},
         {'3OAR40': {'nadp_c': 1.0, '3haACP_c': 1.0, 'h_c': -1.0,
                     'nadph_c': -1.0, 'actACP_c': -1.0},
          '3HAD40': {'3haACP_c': -1.0, 'but2eACP_c': 1.0, 'h2o_c': 1.0},
          'EAR40x': {'nadh_c': -1.0, 'nad_c': 1.0, 'but2eACP_c': -1.0,
                     'butACP_c': 1.0, 'h_c': -1.0},
          '3OAS60': {'butACP_c': -1.0, 'malACP_c': -1.0,
                     '3ohexACP_c': 1.0, 'co2_c': 1.0, 'ACP_c': 1.0,
                     'h_c': -1.0}},
         None,
         None
        ],

        # ----------------
        # Xylitol
        # ----------------

        'crotonic_acid': [
            {
                '3hbcoa_c': {'formula': 'C25H38N7O18P3S', 'name': '(S)-3-Hydroxybutanoyl-CoA'},
                'b2coa_c': {'formula': 'C25H36N7O17P3S', 'name': 'Crotonoyl-CoA'},
            },
            {
                'ECOAH1': {'3hbcoa_c': -1.0, 'b2coa_c': 1.0, 'h2o_c': 1.0},
                'HACD1': {'3hbcoa_c': 1.0, 'aacoa_c': -1.0, 'h_c': -1.0,
                          'nad_c': 1.0, 'nadh_c': -1.0},
            },
            None,
            None,
        ],
    }

def get_iAF1260_designs():
    return {
        'sucD, 4hbd, sucA, cat2, 025B':
        [{'ghb_c': {'formula': 'C4H7O3', 'name': 'gamma-hydroxybutyrate'}},
         {'GHBDHx': {'sucsal_c': -1.0, 'nadh_c': -1.0, 'h_c': -1.0, 'nad_c': 1.0, 'ghb_c': 1.0}},
         None, None]
    }

def get_me_designs():
    return {
        'cHis2A-shmks2-shmks1': [
            {
                '2ptone_c': {'formula': 'C5H10O', 'name': '2-pentatone'},
                'keto_acid_c': {'formula': 'C6H10O3', 'name': '6C keto acid'},
            },
            {
                'THE': {'EG50003-MONOMER_mod_pan4p_mod_3ohex': -1, 'h2o_c': -1,
                        'keto_acid_c': 1, 'EG50003-MONOMER_mod_pan4p': 1},
                'MKS': {'keto_acid_c': -1, '2ptone_c': 1, 'co2_c': 1},
                'EX_2ptone_e': {'2ptone_c': -1},
            },
            None,
            {
                'EX_2ptone_e': (0, 1000),
            }
        ],
    }

def add_heterologous_pathway(model, additions, ignore_repeats=False, copy=False,
                             recompile_expressions=False):
    """Add the given heterologous pathway.

    recompile_expressions: If True, then recompile_expressions when new ME
    reactions are added.

    """
    if copy:
        model = model.copy()

    additions = additions.strip()
    designs = get_designs()
    core_designs = get_core_designs()
    iJR904_designs = get_iJR904_designs()
    iAF1260_designs = get_iAF1260_designs()
    me_designs = get_me_designs()
    if additions in designs:
        if model.id == 'e_coli_core' and (additions in core_designs):
            model = add_pathway(model, *core_designs[additions],
                                check_mass_balance=True,
                                check_charge_balance=False,
                                ignore_repeats=ignore_repeats,
                                recompile_expressions=recompile_expressions)
        elif model.id == 'iJR904' and (additions in iJR904_designs):
            model = add_pathway(model, *iJR904_designs[additions],
                                check_mass_balance=True,
                                check_charge_balance=False,
                                ignore_repeats=ignore_repeats,
                                recompile_expressions=recompile_expressions)
        elif 'iAF1260' in model.id and (additions in iAF1260_designs):
            model = add_pathway(model, *iAF1260_designs[additions],
                                check_mass_balance=True,
                                check_charge_balance=False,
                                ignore_repeats=ignore_repeats,
                                recompile_expressions=recompile_expressions)
        elif model.id == 'ME' and (additions in me_designs):
            return add_pathway(model, *me_designs[additions],
                               check_mass_balance=True,
                               check_charge_balance=False,
                               ignore_repeats=ignore_repeats,
                               recompile_expressions=recompile_expressions)
        return add_pathway(model, *designs[additions], check_mass_balance=True,
                           check_charge_balance=False,
                           ignore_repeats=ignore_repeats,
                           recompile_expressions=recompile_expressions)
    elif additions == 'none':
        return model
    else: # make sure we don't miss any
        raise NotFoundError()

def add_all_heterologous_pathways(model, copy=False):
    """Add all designs to the model.

    """
    if copy:
        model = model.copy()
    for additions in get_designs().iterkeys():
        model = add_heterologous_pathway(model, additions, ignore_repeats=True)
    return model
