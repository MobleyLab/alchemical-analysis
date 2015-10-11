"""
This is the AMBER parser for Alchemical Analysis.  The code can handle both
sander and pmemd output.  TI gradients will be read from the instantaneous
DV/DL entries at every ntpr.  MBAR values are collected too but assumed to be
occuring at the same frequence as the energy output.  BAR/MBAR will be switched
off when the current lambda is not in the set of MBAR lambdas and when the code
detects an overflow (value is all '*') in the energy output.  Component
gradients are collected from the average outputs.
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
        self.filename = filename

        if (filename[-2:] == 'gz'):
          import gzip
          self.fileh = gzip.GzipFile(self.filename, 'rb')
        else:
          self.fileh = open(self.filename, 'rb')

        self.filesize = os.fstat(self.fileh.fileno()).st_size
        self.lineno = 0
        self.last_pos = 0

    def skip_after(self, pattern):
        """Skip until after a line that matches pattern."""

        for line in self:
            match = re.search(pattern, line)

            if match:
                break

        return self.fileh.tell() != self.filesize

    def extract_section(self, start, end, fields, limit=None, extra=''):
        """
        Extract data values (int, float) in fields from a section
        marked with start and end.  Do not read further than limit.
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

                # FIXME: assumes fields are only integers and floats
                if '*' in value:            # Fortran format overflow
                    result.append(None)
                # NOTE: check if this is a sufficient test for int
                elif '.' not in value and re.search(r'\d+', value):
                    result.append(int(value))
                else:
                    result.append(float(value))
            else:                       # section may be incomplete
                result.append(None)

        return result

    def pushback(self):
        """Reposition to one line back.  Works only once, see next()."""
        self.lineno -= 1
        self.fileh.seek(self.last_pos)

    def __iter__(self):
        return self

    def next(self):
        """Read next line of filehandle and remember current position."""
        self.lineno += 1
        curr_pos = self.fileh.tell()

        if curr_pos == self.filesize:
            raise StopIteration

        self.last_pos = curr_pos

        # NOTE: can't mix next() with seek()
        return self.fileh.readline()

    def close(self):
        """Close the filehandle."""
        self.fileh.close()

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        self.close()


class OnlineAvVar(object):
    '''Online algorithm to compute mean (and variance).'''

    def __init__(self):
        self.step = 0
        self.mean = 0.0
        # self.M2 = 0.0

    def accumulate(self, val):
        '''Accumulate data points to compute mean and variance on-the-fly.'''

        self.step += 1

        delta = val - self.mean

        self.mean += delta / self.step
        # self.M2 += delta * (val - self.mean)


def process_mbar_lambdas(secp):
    """
    Extract the lambda points used to compute MBAR energies from
    AMBER MDOUT.
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
    """Simple extrapolation scheme."""

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
        raise SystemExit('Unsupported extrapolation scheme: %s' % scheme)

    return y0, y1


def _print_comps(dvdl_comps_all, ncomp):
    """Print out per force field term free energy components."""
    
    print("\nThe correlated DV/DL components from "
          "_every_single_ step (kcal/mol):")

    ene_comp = []
    x_comp = sorted(dvdl_comps_all.keys())

    for comp in sorted(dvdl_comps_all.items()):
        ene_comp.append(comp[1:])

    fmt = 'Lambda ' + '%10s' * ncomp
    print(fmt % tuple(DVDL_COMPS))

    fmt = '%7.5f' + ' %9.3f' * ncomp

    outstr = {}

    for clambda in x_comp:
        lstr = (clambda,) + tuple(dvdl_comps_all[clambda])
        print(fmt % lstr)

    print('   TI ='),

    for ene in np.transpose(ene_comp):
        x_ene = x_comp
        y_ene = ene[0]

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


def readDataAmber(P):
    """
    Parse free energy gradients and MBAR data from AMBER MDOUT file based on
    sections in the file.
    """

    # To suppress unwanted calls in __main__.
    P.lv_names = ['']

    datafile_tuple = P.datafile_directory, P.prefix, P.suffix
    filenames = glob.glob('%s/%s*%s' % datafile_tuple)

    if not len(filenames):
        raise SystemExit("\nERROR!\nNo files found within directory '%s' with "
                         "prefix '%s' and suffix '%s': check your inputs."
                         % datafile_tuple)

    dvdl_all = defaultdict(list)
    dvdl_comps_all = defaultdict(list)
    mbar_all = {}

    nsnapshots = []

    ncomp = len(DVDL_COMPS)
    global_have_mbar = True
    pmemd = False

    for filename in filenames:
        print('Loading in data from %s... ' % filename),

        in_comps = False
        finished = False

        dvdl_data = []
        dvdl_comp_data = []

        for _ in DVDL_COMPS:
            dvdl_comp_data.append(OnlineAvVar())

        with SectionParser(filename) as secp:
            lineno = 0
            line = ''

            for line in secp:
                lineno += 1

                if lineno > 5:
                    break

            if 'PMEMD' in line:
                pmemd = True

            if not secp.skip_after('^   2.  CONTROL  DATA  FOR  THE  RUN'):
                print('WARNING: no control data found, ignoring file')
                continue

            # NOTE: sections must be searched for in order!
            ntpr, = secp.extract_section('^Nature and format of output:', '^$',
                                         ['ntpr'])

            nstlim, dt = secp.extract_section('Molecular dynamics:', '^$',
                                              ['nstlim', 'dt'])

            # FIXME: check if temp0, dt, etc. are consistent between files
            P.temperature, = secp.extract_section('temperature regulation:',
                                                  '^$', ['temp0'])

            # FIXME: file may end just after "2. CONTROL DATA" so vars will
            #        be all None

            # NOTE: some sections may not have been created
            clambda, = secp.extract_section('^Free energy options:', '^$',
                                            ['clambda'], '^---')

            if clambda is None:
                print('WARNING: no free energy section found, ignoring file')
                continue

            mbar_ndata = 0
            have_mbar, mbar_ndata = secp.extract_section('^FEP MBAR options:',
                                                         '^$',
                                                         ['ifmbar',
                                                          'bar_intervall'],
                                                         '^---')

            # sander is just too cumbersome with MBAR, e.g. terminator is not
            # '^---', no end-points, etc
            if not pmemd:
                have_mbar = False

            if have_mbar:
                mbar_ndata = int(nstlim / mbar_ndata)
                mbar_lambdas = process_mbar_lambdas(secp)
                clambda_str = '%6.4f' % clambda

                # FIXME: case when lambda is contained in mbar_lambdas but
                #        mbar_lambdas has additional entries
                if clambda_str not in mbar_lambdas:
                    if global_have_mbar:
                        print('\nWARNING: lambda %s not contained in set of '
                              'MBAR lambdas: %s\nNot using MBAR.' %
                              (clambda_str, ', '.join(mbar_lambdas)))

                    global_have_mbar = have_mbar = False
                else:
                    mbar_nlambda = len(mbar_lambdas)
                    mbar_lambda_idx = mbar_lambdas.index(clambda_str)
                    mbar_data = []

                    for _ in range(mbar_nlambda):
                        mbar_data.append([])
            else:
                global_have_mbar = False

            if not secp.skip_after('^   4.  RESULTS'):
                print('WARNING: no results found, ignoring file\n')
                continue

            nensec = 0
            nenav = 0
            old_nstep = -1
            old_comp_nstep = -1
            incomplete = False

            for line in secp:
                if have_mbar and line.startswith('MBAR Energy analysis'):
                    mbar = secp.extract_section('^MBAR', '^ ---', mbar_lambdas,
                                                extra=line)

                    if not all(mbar):
                        if global_have_mbar:
                            print('\nWARNING: some MBAR energies cannot be '
                                  'read. Not using MBAR.')

                        global_have_mbar = have_mbar = False
                        continue

                    E_ref = mbar[mbar_lambda_idx]

                    for lmbda, E in enumerate(mbar):
                        mbar_data[lmbda].append(E - E_ref)

                if 'DV/DL, AVERAGES OVER' in line:
                    in_comps = True

                if line.startswith(' NSTEP'):
                    if in_comps:
                        result = secp.extract_section('^ NSTEP', '^ ---',
                                                      ['NSTEP'] + DVDL_COMPS,
                                                      extra=line)

                        for res in result:
                            if res is None:
                                incomplete = True

                        if result[0] != old_comp_nstep and not incomplete:
                            for i, E in enumerate(DVDL_COMPS):
                                dvdl_comp_data[i].accumulate(
                                    float(result[i+1]))

                            nenav += 1
                            old_comp_nstep = result[0]
                            incomplete = False

                        in_comps = False
                    else:
                        nstep, dvdl = secp.extract_section('^ NSTEP', '^ ---',
                                                           ['NSTEP', 'DV/DL'])

                        for res in nstep, dvdl:
                            if res is None:
                                incomplete = True

                        if nstep != old_nstep and not incomplete:
                            dvdl_data.append(dvdl)
                            nensec += 1
                            old_nstep = nstep
                            incomplete = False

                if line == '   5.  TIMINGS\n':
                    finished = True
                    break

        # -- end of parsing current file --

        print('%i data points, %i DV/DL averages' % (nensec, nenav))

        if not finished:
            print('WARNING: prematurely terminated run\n')
            continue

        if not nensec:
            print('WARNING: File %s does not contain any DV/DL data\n' %
                  filename)
            continue

        if have_mbar:
            if clambda not in mbar_all:
                mbar_all[clambda] = []

                for _ in range(mbar_nlambda):
                    mbar_all[clambda].append([])

            for mbar, data in zip(mbar_all[clambda], mbar_data):
                mbar.extend(data)

        dvdl_all[clambda].extend(dvdl_data)
        dvdl_comps_all[clambda] = [Enes.mean for Enes in dvdl_comp_data]

    # -- all file parsing finished --

    if not dvdl_all:
        raise SystemExit('No DV/DL data found')

    if not global_have_mbar:
        if 'BAR' in P.methods:
            P.methods.remove('BAR')

        if 'MBAR' in P.methods:
            P.methods.remove('MBAR')

        print('\nWARNING: BAR/MBAR have been switched off.')

    ave = []
    start_from = int(round(P.equiltime / (ntpr * float(dt))))

    # FIXME: compute maximum number of MBAR energy sections
    lvals = sorted(dvdl_all)
    K = len(dvdl_all)
    nsnapshots = [len(e) - start_from for e in dvdl_all.values()]
    maxn = max(nsnapshots)

    # AMBER has currently only one global lambda value, hence 2nd dim = 1
    dhdlt = np.zeros([K, 1, maxn], np.float64)

    if have_mbar and global_have_mbar:
        assert K == len(mbar_all)
        u_klt = np.zeros([K, K, maxn], np.float64)
    else:
        u_klt = None

    for i, clambda in enumerate(lvals):
        vals = dvdl_all[clambda][start_from:]
        ave.append(np.average(vals))

        dhdlt[i][0][:len(vals)] = np.array(vals)

        if u_klt is not None:
            for j, ene in enumerate(mbar_all[clambda]):
                enes = ene[start_from:]
                l_enes = len(enes)
                u_klt[i][j][:l_enes] = enes

    if u_klt is not None:
        u_klt = P.beta * u_klt

    # sander does not sample end-points...
    y0, y1 = _extrapol(lvals, ave, 'polyfit')

    if y0:
        print('Note: adding missing lambda = 0.0: %f' % y0)

        K += 1
        lvals.insert(0, 0.0)
        nsnapshots.insert(0, maxn)
        ave.insert(0, y0)

        # we need a little bit of noise to get the statistics, otherwise
        # covariance will be zero and error termination
        frand = y0 + _RND_SCALE * np.random.rand(maxn) - _RND_SCALE_HALF
        dhdlt = np.insert(dhdlt, 0, frand, 0)

    if y1:
        print('Note: adding missing lambda = 1.0: %f' % y1)

        K += 1
        lvals.append(1.0)
        nsnapshots.append(maxn)#
        ave.append(y1)

        frand = y1 + _RND_SCALE * np.random.rand(maxn) - _RND_SCALE_HALF
        dhdlt = np.append(dhdlt, [[frand]], 0)

    print("\nThe gradients from the correlated samples (kcal/mol):")

    for clambda, dvdl in zip(lvals, ave):
        print('%7.5f %9.3f' % (clambda, dvdl) )

    print('   TI = %9.3f' % np.trapz(ave, lvals))

    _print_comps(dvdl_comps_all, ncomp)

    print('\n\n')

    return (np.array(nsnapshots), np.array(lvals).reshape(K, 1), dhdlt, u_klt)
