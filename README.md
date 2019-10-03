# Alchemical Analysis: An open tool implementing some recommended practices for analyzing alchemical free energy calculations

## Use Alchemlyb instead

We are in the process of migrating all functionality from here to instead use [`alchemlyb`](https://github.com/alchemistry/alchemlyb), which focuses on being a general Python library for performing analysis of alchemical calculations rather than a stand-alone command-line tool for analyzing these calculations.
We recommend you move to using, and contribute to the development of, `alchemlyb` instead.
It already has the following advantages, among others:
- Python 3 compatible
- Works as a library rather than a stand-alone tool, allowing easier connection with trajectory analysis
- Has automated testing/continuous integration testing
- Easily extensible

However, some of the plotting functionality and analysis features available in alchemical analysis are not yet available in `alchemlyb`.
If there is functionality which is particularly important to you which is missing there, please raise an issue on the `alchemlyb` issue tracker.

In the meantime, Alchemical Analysis is still available but **deprecated and we discourage its use. Please use alchemlyb instead.** We may occasionally perform minor/critical maintenance here, but intend to cease maintaining this package over time.

## About Alchemical Analysis

Analyze alchemical free energy calculations conducted in GROMACS, AMBER or SIRE using recommended best practices from Klimovich et al., JCAMD 29:397-411 (2015).

This tool handles analysis via a slate of free energy methods, including BAR, MBAR, TI, and the Zwanzig relationship (exponential averaging) among others, and provides a good deal of analysis of computed free energies and convergence in order to help you assess the quality of your results.

If you have problems with this tool, please use the [github issue tracker](https://github.com/mobleylab/alchemical-analysis/issues).

#### Citation [![DOI for Citing Alchemical Analysis](https://img.shields.io/badge/DOI-10.007%2Fs10822--015--9840--9-blue.svg)](http://dx.doi.org/10.1007/s10822-015-9840-9)

Alchemical Analysis is research software. If you make use of it in work which you publish, please cite it. The BibTex reference is

```
@article{Klimovich:2015er,
author = {Klimovich, Pavel V and Shirts, Michael R and Mobley, David L},
title = {Guidelines for the analysis of free energy calculations},
journal = {J Comput Aided Mol Des},
year = {2015},
volume = {29},
number = {5},
pages = {397--411},
doi = {10.1007/s10822-015-9840-9}
}
```

#### Prerequisites

Alchemical Analysis requires the `pymbar` module (we recommend
installation via conda, but other options are also available). Version
3.0.0.dev0 or above is required. The latest version of pymbar can be found
[here](https://github.com/choderalab/pymbar).

Python 2 is required.

#### Installation

    git clone https://github.com/MobleyLab/alchemical-analysis.git
    cd alchemical-analysis
    sudo python setup.py install

#### Usage

Script: `alchemical_analysis`

This implements recommended practices for analyzing alchemical free energy calculations, as described in Klimovich et al., JCAMD 29:397-411 (2015). This was motivated in part by earlier work illustrating how to apply MBAR to alchemical free energy calculations (and a comparison with other methods) in Paliwal and Shirts, J. Chem. Theory Comp, v. 7, 4115-4134 (2011).

See description in `samples/gromacs/README.md`, `samples/sire/README.md`, and `samples/amber/README.md`.


Help for `alchemical_analysis.py` (obtained with `alchemical_analysis -h`) is:

```Options:
  -h, --help            show this help message and exit
  -a SOFTWARE, --software=SOFTWARE
                        Package's name the data files come from: Gromacs,
                        Sire, Desmond, or AMBER. Default: Gromacs.
  -b FRACTION, --fraction=FRACTION
                        The fraction of the energy file will be used to
                        calculate the statistics. Default: 1.0
  -c, --cfm             The Curve-Fitting-Method-based consistency inspector.
                        Default: False.
  -d DATAFILE_DIRECTORY, --dir=DATAFILE_DIRECTORY
                        Directory in which data files are stored. Default:
                        Current directory.
  -e, --backward        Extract the energy data from the backward direction.
                        Default: False
  -f BFORWREV, --forwrev=BFORWREV
                        Plot the free energy change as a function of time in
                        both directions, with the specified number of points
                        in the time plot. The number of time points (an
                        integer) must be provided. Default: 0
  -g, --breakdown       Plot the free energy differences evaluated for each
                        pair of adjacent states for all methods, including the
                        dH/dlambda curve for TI. Default: False.
  -i UNCORR_THRESHOLD, --threshold=UNCORR_THRESHOLD
                        Proceed with correlated samples if the number of
                        uncorrelated samples is found to be less than this
                        number. If 0 is given, the time series analysis will
                        not be performed at all. Default: 50.
  -j RESULTFILENAME, --resultfilename=RESULTFILENAME
                        custom defined result filename prefix. Default:
                        results
  -k BSKIPLAMBDAINDEX, --koff=BSKIPLAMBDAINDEX
                        Give a string of lambda indices separated by '-' and
                        they will be removed from the analysis. (Another
                        approach is to have only the files of interest present
                        in the directory). Default: None.
  -m METHODS, --methods=METHODS
                        A list of the methods to esitimate the free energy
                        with. Default: [TI, TI-CUBIC, DEXP, IEXP, BAR, MBAR].
                        To add/remove methods to the above list provide a
                        string formed of the method strings preceded with +/-.
                        For example, '-ti_cubic+gdel' will turn methods into
                        [TI, DEXP, IEXP, BAR, MBAR, GDEL]. 'ti_cubic+gdel', on
                        the other hand, will call [TI-CUBIC, GDEL]. 'all'
                        calls the full list of supported methods [TI, TI-
                        CUBIC, DEXP, IEXP, GINS, GDEL, BAR, UBAR, RBAR, MBAR].
  -n UNCORR, --uncorr=UNCORR
                        The observable to be used for the autocorrelation
                        analysis; either 'dhdl_all' (obtained as a sum over
                        all energy components) or 'dhdl' (obtained as a sum
                        over those energy components that are changing;
                        default) or 'dE'. In the latter case the energy
                        differences dE_{i,i+1} (dE_{i,i-1} for the last
                        lambda) are used.
  -o OUTPUT_DIRECTORY, --out=OUTPUT_DIRECTORY
                        Directory in which the output files produced by this
                        script will be stored. Default: Same as
                        datafile_directory.
  -p PREFIX, --prefix=PREFIX
                        Prefix for datafile sets, i.e.'dhdl' (default).
  -q SUFFIX, --suffix=SUFFIX
                        Suffix for datafile sets, i.e. 'xvg' (default).
  -r DECIMAL, --decimal=DECIMAL
                        The number of decimal places the free energies are to
                        be reported with. No worries, this is for the text
                        output only; the full-precision data will be stored in
                        'results.pickle'. Default: 3.
  -s EQUILTIME, --skiptime=EQUILTIME
                        Discard data prior to this specified time as
                        'equilibration' data. Units picoseconds. Default: 0
                        ps.
  -t TEMPERATURE, --temperature=TEMPERATURE
                        Temperature in K. Default: 298 K.
  -u UNITS, --units=UNITS
                        Units to report energies: 'kJ', 'kcal', and 'kBT'.
                        Default: 'kJ'
  -v, --verbose         Verbose option. Default: False.
  -w, --overlap         Print out and plot the overlap matrix. Default: False.
  -x, --ignoreWL        Do not check whether the WL weights are equilibrated.
                        No log file needed as an accompanying input.
  -y RELATIVE_TOLERANCE, --tolerance=RELATIVE_TOLERANCE
                        Convergence criterion for the energy estimates with
                        BAR and MBAR. Default: 1e-10.
  -z INIT_WITH, --initialize=INIT_WITH
                        The initial MBAR free energy guess; either 'BAR' or
                        'zeros'. Default: 'BAR'.
```


#### License

MIT.
