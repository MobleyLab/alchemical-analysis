"""
This is the AMBER parser for the Alchemical Analysis tool.

The code can handle both sander and pmemd output.  TI gradients will be read
from the instantaneous DV/DL entries in the MDOUT file (the MDEN file is not
used) which is written every ntpr step.  MBAR values are collected but assumed
to be occuring at the same frequency as the DV/DL output (bar_intervall=ntpr).
BAR/MBAR will be switched off when the current lambda is not in the set of MBAR
lambdas.  MBAR is currently switched off also for sander output because sander
can't sample the end-points.  Component gradients are collected from the
average outputs (ntave).  GZIP/BZIP2 compression is supported.

Relevant namelist variables in MDIN: ntpr, ntave, ifmbar and bar*,
do NOT set t as the code depends on a correctly set begin time
"""

######################################################################
#
# Alchemical Analysis: An open tool implementing some recommended practices
# for analyzing alchemical free energy calculations
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: Pavel Klimovich, Michael Shirts and David Mobley
# Authors of this module: Hannes H Loeffler, Pavel Klimovich
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http://www.gnu.org/licenses/>.
#
######################################################################


import os
import re
import glob
from collections import defaultdict

import numpy as np


DVDL_COMPS = ['BOND', 'ANGLE', 'DIHED', '1-4 NB', '1-4 EEL', 'VDWAALS',
              'EELEC', 'RESTRAINT']
_FP_RE = r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
_RND_SCALE = 1e-3
_RND_SCALE_HALF = _RND_SCALE / 2.0

_MAGIC_CMPR = {
    '\x1f\x8b\x08': ('gzip', 'GzipFile'),  # last byte is compression method
    '\x42\x5a\x68': ('bz2', 'BZ2File')
}


class FEData(object):
    """A simple struct container to collect data from individual files."""

    __slots__ = ['clambda', 't0', 'dt', 'T', 'gradients',
                 'component_gradients', 'mbar_energies']

    def __init__(self):
        self.clambda = -1.0
        self.t0 = -1.0
        self.dt = -1.0
        self.T = -1.0
        self.gradients = []
        self.component_gradients = []
        self.mbar_energies = []


def _pre_gen(it, first):
    """A generator that returns first first if it exists."""

    if first:
        yield first

    while it:
        yield it.next()


class SectionParser(object):
    """
    A simple parser to extract data values from sections.
    """

    def __init__(self, filename):
        """Opens a file according to its file type."""

        self.filename = filename
        
        with open(filename, 'rb') as f:
            magic = f.read(3)   # NOTE: works because all 3-byte headers

        try:
            method = _MAGIC_CMPR[magic]
        except KeyError:
            open_it = open
        else:
            open_it = getattr(__import__(method[0]), method[1])

        try:
            self.fileh = open_it(self.filename, 'rb')
            self.filesize = os.stat(self.filename).st_size
        except IOError:
            raise SystemExit('ERROR: cannot open file %s' % filename)

        self.lineno = 0

    def skip_lines(self, nlines):
        """Skip a given number of files."""

        lineno = 0

        for line in self:
            lineno += 1

            if lineno > nlines:
                return line

        return None

    def skip_after(self, pattern):
        """Skip until after a line that matches a regex pattern."""

        for line in self:
            match = re.search(pattern, line)

            if match:
                break

        return self.fileh.tell() != self.filesize

    def extract_section(self, start, end, fields, limit=None, extra='',
                        debug=False):
        """
        Extract data values (int, float) in fields from a section
        marked with start and end regexes.  Do not read further than
        limit regex.
        """

        inside = False
        lines = []

        for line in _pre_gen(self, extra):
            if limit and re.search(limit, line):
                break

            if re.search(start, line):
                inside = True

            if inside:
                if re.search(end, line):
                    break

                lines.append(line.rstrip('\n'))

        line = ''.join(lines)
        result = []

        for field in fields:
            match = re.search(r' %s\s+=\s+(\*+|%s|\d+)'
                              % (field, _FP_RE), line)

            if match:
                value = match.group(1)

                # FIXME: assumes fields are only integers or floats
                if '*' in value:            # Fortran format overflow
                    result.append(float('Inf') )
                # NOTE: check if this is a sufficient test for int
                elif '.' not in value and re.search(r'\d+', value):
                    result.append(int(value))
                else:
                    result.append(float(value))
            else:                       # section may be incomplete
                result.append(None)

        return result

    def __iter__(self):
        return self

    def next(self):
        """Read next line of the filehandle and check for EOF."""

        self.lineno += 1
        curr_pos = self.fileh.tell()

        if curr_pos == self.filesize:
            raise StopIteration

        # NOTE: can't mix next() with seek()
        return self.fileh.readline()

    def close(self):
        """Close the filehandle."""

        self.fileh.close()

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        self.close()


def _process_mbar_lambdas(secp):
    """
    Extract the lambda points used to compute MBAR energies from an
    AMBER MDOUT file.
    """

    in_mbar = False
    mbar_lambdas = []

    for line in secp:
        if line.startswith('    MBAR - lambda values considered:'):
            in_mbar = True
            continue

        if in_mbar:
            if line.startswith('    Extra'):
                break

            if 'total' in line:
                data = line.split()
                mbar_lambdas.extend(data[2:])
            else:
                mbar_lambdas.extend(line.split())

    return mbar_lambdas


def _extrapol(x, y, scheme):
    """Simple extrapolation schemes."""

    y0 = None
    y1 = None

    if scheme == 'linear':
        if 0.0 not in x:
            y0 = (x[0] * y[1] - x[1] * y[0]) / (x[0] - x[1])

        if 1.0 not in x:
            y1 = (((x[-2] - 1.0) * y[-1] + ((1.0 - x[-1]) * y[-2])) /
                  (x[-2] - x[-1]))
    elif scheme == 'polyfit':
        nval = len(x)

        if nval < 6:
            deg = nval - 1
        else:
            deg = 6

        coeffs = np.polyfit(x, y, deg)

        if 0.0 not in x:
            y0 = coeffs[-1]

        if 1.0 not in x:
            y1 = sum(coeffs)
    else:
        raise SystemExit('ERROR: Unsupported extrapolation scheme: %s' %
                         scheme)

    return y0, y1


def _print_comps(comps, ncomp, units, conv):
    """
    Print out per-force-field-term free energy components and
    resulting dG.
    """
    
    print('\nThe correlated gradient (DV/DL) components from '
          '_every_single_ step %s:\n' % units)

    fmt = 'Lambda ' + '%10s' * ncomp
    print(fmt % tuple(DVDL_COMPS))

    fmt = '%7.5f' + ' %9.3f' * ncomp
    sep = '-' * (7 + 10 * ncomp)     # format plus ncomp spaces
    print('%s' % sep)

    x_comp = sorted(comps)
    ene_comp = []

    # comps may contain the component free energies from several time steps
    # so average them
    for clambda in x_comp:
        Enes = comps[clambda]
        ave = [0.0 for _ in range(ncomp)]

        for Ene in Enes:
            for j in range(ncomp):
                ave[j] += Ene[j] * conv

        ave = [a / len(Enes) for a in ave] # len(Enes) = no. of time steps/file
        ene_comp.append(ave)            # tansposed later to compute dG

        lstr = (clambda,) + tuple(ave)
        print(fmt % lstr)

    print('%s\n   dG =' % sep),

    for ene in zip(*ene_comp):
        x_ene = x_comp
        y_ene = ene

        if not all(y_ene):
            print(' %8.3f' % 0.0),
            continue

        ya, yb = _extrapol(x_comp, y_ene, 'polyfit')

        if ya:
            x_ene = [0.0] + x_ene
            y_ene = np.insert(y_ene, 0, ya)

        if yb:
            x_ene = x_ene + [1.0]
            y_ene = np.append(y_ene, yb)

        print(' %8.3f' % np.trapz(y_ene, x_ene)),


def any_none(sequence):
    """Check if any element of a sequence is None."""

    for element in sequence:
        if element is None:
            return True

    return False


def readDataAmber(P):
    """
    Parse free energy gradients and MBAR data from an AMBER MDOUT file. Also
    read the meta data from the header sections in the file.
    """

    P.lv_names = [r'all']               # legend for plotting
    P.bExpanded = False

    # FIXME: this should happen in the main code
    datafile_tuple = P.datafile_directory, os.sep, P.prefix, P.suffix
    filenames = glob.glob('%s%s%s*%s' % datafile_tuple)

    if not len(filenames):
        raise SystemExit("\nERROR: no files found within directory '%s' with "
                         "prefix '%s' and suffix '%s': check your inputs."
                         % datafile_tuple)

    file_data = []
    pmemd = False

    for filename in filenames:
        print('Loading in data from %s... ' % filename),

        file_datum = FEData()
        finished = False
        comps = []

        with SectionParser(filename) as secp:
            line = secp.skip_lines(5)

            if not line:
                print('  WARNING: file does not contain any useful data, '
                      'ignoring file')
                continue

            if 'PMEMD' in line:
                pmemd = True

            if not secp.skip_after('^   2.  CONTROL  DATA  FOR  THE  RUN'):
                print('  WARNING: no CONTROL DATA found, ignoring file')
                continue

            # NOTE: sections must be searched for in order
            #       some sections may not exist in MDOUT when certain flags
            #       are not set, e.g. MBAR only available if ifmbar>0
            ntpr, = secp.extract_section('^Nature and format of output:', '^$',
                                         ['ntpr'])

            nstlim, dt = secp.extract_section('Molecular dynamics:', '^$',
                                              ['nstlim', 'dt'])

            T, = secp.extract_section('temperature regulation:', '^$',
                                     ['temp0'])

            # FIXME: is it reasonable to support non-constT?
            #        probably just for BAR related methods
            if not T:
                raise SystemExit('ERROR: Non-constant temperature MD not '
                                 'currently supported')

            P.temperature = T

            clambda, = secp.extract_section('^Free energy options:', '^$',
                                            ['clambda'], '^---')

            if clambda is None:
                print('  WARNING: no free energy section found, ignoring file')
                continue

            mbar_ndata = 0
            have_mbar, mbar_ndata = secp.extract_section('^FEP MBAR options:',
                                                         '^$',
                                                         ['ifmbar',
                                                          'bar_intervall'],
                                                         '^---')

            # sander is just too cumbersome with MBAR, e.g. terminator is not
            # '^---', no end-states, etc
            if not pmemd:
                have_mbar = False

            # FIXME: what other methods depend on MBAR data?
            if 'BAR' not in P.methods and 'MBAR' not in P.methods:
                have_mbar = False

            if have_mbar:
                mbar_ndata = int(nstlim / mbar_ndata)
                mbar_lambdas = _process_mbar_lambdas(secp)
                clambda_str = '%6.4f' % clambda

                if clambda_str not in mbar_lambdas:
                    print('\n  WARNING: lambda %s not contained in set of '
                          'MBAR lambdas: %s\nNot using MBAR.' %
                          (clambda_str, ', '.join(mbar_lambdas)))

                    have_mbar = False
                else:
                    mbar_nlambda = len(mbar_lambdas)
                    mbar_lambda_idx = mbar_lambdas.index(clambda_str)

                    for _ in range(mbar_nlambda):
                        file_datum.mbar_energies.append([])

            if not secp.skip_after('^   3.  ATOMIC '):
                print('  WARNING: no ATOMIC section found, ignoring file\n')
                continue

            t0, = secp.extract_section('^ begin time', '^$', ['coords'])

            if not secp.skip_after('^   4.  RESULTS'):
                print('  WARNING: no RESULTS section found, ignoring file\n')
                continue

            file_datum.clambda = clambda
            file_datum.t0 = t0
            file_datum.dt = dt
            file_datum.T = T

            nensec = 0
            nenav = 0
            old_nstep = -1
            old_comp_nstep = -1
            high_E_cnt = 0

            in_comps = False

            for line in secp:
                if have_mbar and line.startswith('MBAR Energy analysis'):
                    mbar = secp.extract_section('^MBAR', '^ ---', mbar_lambdas,
                                                extra=line)

                    if any_none(mbar):
                        continue

                    E_ref = mbar[mbar_lambda_idx]

                    for lmbda, E in enumerate(mbar):
                        # NOTE: should be ok for pymbar because exp(-u)
                        if E > 0.0:
                            high_E_cnt += 1

                        file_datum.mbar_energies[lmbda].append(E - E_ref)

                if 'DV/DL, AVERAGES OVER' in line:
                    in_comps = True

                if line.startswith(' NSTEP'):
                    if in_comps:
                        result = secp.extract_section('^ NSTEP', '^ ---',
                                                      ['NSTEP'] + DVDL_COMPS,
                                                      extra=line)

                        if result[0] != old_comp_nstep and not any_none(result):
                            comps.append([float(E) for E in result[1:]])

                            nenav += 1
                            old_comp_nstep = result[0]

                        in_comps = False
                    else:
                        nstep, dvdl = secp.extract_section('^ NSTEP', '^ ---',
                                                           ['NSTEP', 'DV/DL'],
                                                           extra=line)

                        if nstep != old_nstep and dvdl is not None \
                               and nstep is not None:
                            file_datum.gradients.append(dvdl)
                            nensec += 1
                            old_nstep = nstep

                if line == '   5.  TIMINGS\n':
                    finished = True
                    break


            if high_E_cnt:
                print('\n  WARNING: %i MBAR energ%s > 0.0 kcal/mol' %
                      (high_E_cnt, 'ies are' if high_E_cnt > 1 else 'y is') )


        # -- end of parsing current file --

        print('%i data points, %i DV/DL averages' % (nensec, nenav))

        file_datum.component_gradients.extend(comps)
        file_data.append(file_datum)

        if not finished:
            print('  WARNING: prematurely terminated run')

        if not nensec:
            print('  WARNING: File %s does not contain any DV/DL data\n' %
                  filename)

    # -- all file parsing done --


    # NOTE: lambda sorting is not required because data is collected in dict
    #       which is later sorted on the keys
    print('\nSorting input data by starting time')

    file_data.sort(key=lambda fd: fd.t0)

    dvdl_all = defaultdict(list)
    dvdl_comps_all = defaultdict(list)
    mbar_all = {}

    t0_found = defaultdict(list)
    t0_uniq = defaultdict(set)
    dt_uniq = set()
    T_uniq = set()
    clambda_uniq = set()

    for fd in file_data:
        clambda = fd.clambda

        dvdl_all[clambda].extend(fd.gradients)
        dvdl_comps_all[clambda].extend(fd.component_gradients)
        t0_found[clambda].append(fd.t0)
        t0_uniq[clambda].add(fd.t0)

        dt_uniq.add(fd.dt)
        T_uniq.add(fd.T)
        clambda_uniq.add(clambda)

        if have_mbar:
            if clambda not in mbar_all:
                mbar_all[clambda] = [[] for _ in range(mbar_nlambda)]

            for mbar, data in zip(mbar_all[clambda], fd.mbar_energies):
                mbar.extend(data)


    lvals = sorted(clambda_uniq)

    if not dvdl_all:
        raise SystemExit('ERROR: No DV/DL data found')

    if not have_mbar:
        # FIXME: needs to be handled by main code
        if 'BAR' in P.methods:
            P.methods.remove('BAR')

        if 'MBAR' in P.methods:
            P.methods.remove('MBAR')

        print('\nNote: BAR/MBAR results are not computed.')
    elif len(dvdl_all) != len(mbar_lambdas):
        ndvdl = len(dvdl_all)
        raise SystemExit('ERROR: gradient samples have been found for %i '
                         'lambda%s:\n%s\n       but MBAR data has %i:\n%s\n' %
                         (ndvdl, 's' if ndvdl > 1 else '',
                          ', '.join([str(l) for l in lvals]),
                          len(mbar_lambdas),
                          ', '.join([str(float(l)) for l in mbar_lambdas]) ) )

    for found, uniq in zip(t0_found.values(), t0_uniq.values() ):
        if len(found) != len(uniq):
            raise SystemExit('ERROR: Same starting time occurs multiple '
                             'times.')

    if len(dt_uniq) != 1:
        raise SystemExit('ERROR: Not all files have the same time step (dt).')

    if len(T_uniq) != 1:
        raise SystemExit('ERROR: Not all files have the same temperature (T).')

    # P.beta has been computed with P.temperature from the command line
    T = T_uniq.pop()
    P.beta *= P.temperature / T
    P.beta_report *= P.temperature / T
    P.temperature = T

    # FIXME: check if ntpr is consistent across all input files
    P.snap_size = [dt_uniq.pop() * ntpr]

    start_from = int(round(P.equiltime / (ntpr * float(dt))))
    print('\nSkipping first %d steps (= %f ps)\n' % (start_from, P.equiltime))

    # FIXME: compute maximum number of MBAR energy sections
    K = len(lvals)
    nsnapshots = [len(dvdl_all[clambda]) - start_from for clambda in lvals]
    nsnapshots = [int(i*P.fraction) for i in nsnapshots]
    maxn = max(nsnapshots)

    # AMBER has currently only one global lambda value, hence 2nd dim = 1
    dhdlt = np.zeros([K, 1, maxn], np.float64)

    if have_mbar:
        assert K == len(mbar_all)
        u_klt = np.zeros([K, K, maxn], np.float64)
    else:
        u_klt = None

    # FIXME: needs centralised solution
    if P.units == '(kcal/mol)':
        Econv = 1.0
    elif P.units == '(kJ/mol)':
        Econv = 4.184
    elif P.units == '(k_BT)':
        Econv = P.beta
    else:
        raise SystemExit('ERROR: unknown units %s' % P.units)

    ave = []
    std = []

    for i, clambda in enumerate(lvals):
        vals = dvdl_all[clambda][start_from:]
        if P.backward:
            end_position = len(vals)-int(P.fraction*len(vals))
            vals = vals[end_position:]
        else:
            end_position = int(len(vals)*P.fraction)
            vals = vals[:end_position]
        ave.append(np.average(vals) * Econv)
        std.append(np.std(vals) * Econv)

        dhdlt[i][0][:len(vals)] = np.array(vals)

        if have_mbar:
            for j, ene in enumerate(mbar_all[clambda]):
                enes = ene[start_from:]
                l_enes = len(enes)
                u_klt[i][j][:l_enes] = enes

    if have_mbar:
        u_klt = P.beta * u_klt


    # sander does not sample end-points...
    y0, y1 = _extrapol(lvals, ave, 'polyfit')

    if y0:
        print('Note: extrapolating missing gradient for lambda = 0.0: %f' % y0)

        K += 1
        lvals.insert(0, 0.0)
        nsnapshots.insert(0, maxn)
        ave.insert(0, y0)

        # we need a little bit of noise to get the statistics, otherwise
        # covariance will be zero and program will terminate with error
        frand = y0 + _RND_SCALE * np.random.rand(maxn) - _RND_SCALE_HALF
        dhdlt = np.insert(dhdlt, 0, frand, 0)

    if y1:
        print('Note: extrapolating missing gradient for lambda = 1.0: %f' % y1)

        K += 1
        lvals.append(1.0)
        nsnapshots.append(maxn)
        ave.append(y1)

        frand = y1 + _RND_SCALE * np.random.rand(maxn) - _RND_SCALE_HALF
        dhdlt = np.append(dhdlt, [[frand]], 0)

    print('\nThe gradients (DV/DL) from the correlated samples %s:\n\n'
          'Lambda   gradient' % P.units)

    sep = '-' * (7 + 9 + 1)             # format plus 1 space
    print('%s' % sep)

    for clambda, dvdl in zip(lvals, ave):
        print('%7.5f %9.3f' % (clambda, dvdl))

    print('%s\n   dG = %9.3f' % (sep, np.trapz(ave, lvals)))

    _print_comps(dvdl_comps_all, len(DVDL_COMPS), P.units, Econv)

    print('\n')

    # FIXME: make this a parser dependent option
    with open(os.path.join(P.output_directory, 'grads.dat'), 'w') as gfile:
        gfile.write('# gradients (DV/DL) for lambdas: %s\n' %
                    ' '.join('%s' % l for l in lvals) )

        for i in range(maxn):
            gfile.write('%f ' % (i * dt) )

            for clambda in lvals:
                try:
                    val = dvdl_all[clambda][i]
                except IndexError:
                    val = 0.0

                gfile.write('%f ' % val)

            gfile.write('\n')


    return (np.array(nsnapshots), np.array(lvals).reshape(K, 1), dhdlt, u_klt)
