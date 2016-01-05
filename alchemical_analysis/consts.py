r"""
Various constants for the Alchemical Analysis tool.
"""


### physical constants

# Boltzmann constant in kJ/mol/K:
# 1.3806488*10^-23 J/K (Boltzmann) * 6.02214129*10^-23 mol^-1 (Avogadro) / 1000
kB = 0.008314462145468952

# calories-to-Joules conversion factor
CAL2JOULE = 4.184


### analysis methods

# NOTE: sets do not keep order!
TI_METHODS = ('TI', 'TI-CUBIC')
BAR_METHODS = ('DEXP', 'IEXP', 'GINS', 'GDEL', 'BAR', 'UBAR', 'RBAR', 'MBAR')
ALL_METHODS = TI_METHODS + BAR_METHODS
DEFAULT_METHODS = ('TI', 'TI-CUBIC', 'DEXP', 'IEXP', 'BAR', 'MBAR')


### file names

RESULTS_FILE = 'results.txt'
RESULTS_PICKLE = 'results.pickle'
DF_T_FILE = 'dF_t.txt'

MBAR_OVERLAP_PDF = 'O_MBAR.pdf'
DF_T_PDF = 'dF_t.pdf'
DF_STATE_LONG_PDF = 'dF_state_long.pdf'
DF_STATE_PDF = 'dF_state.pdf'
DHDL_TI_PDF = 'dhdl_TI.pdf'
CFM_PDF = 'cfm.pdf'
