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



import os, re, glob
from collections import defaultdict

import numpy



DVDL_COMPS = ['BOND', 'ANGLE', 'DIHED', '1-4 NB', '1-4 EEL', 'VDWAALS', 'EELEC',
              'RESTRAINT']
_FP_RE = r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'


class SectionParser(object):
    """
    A simple parser to extract data values from sections.
    """

    def __init__(self, filename):
        self.filename = filename
        self.fh = open(self.filename, 'rb')
        self.filesize = os.fstat(self.fh.fileno()).st_size
        self.lineno = 0
        self.last_pos = 0


    def skip_after(self, pattern):
        """Skip until after a line that matches pattern."""

        for line in self:
            match = re.search(pattern, line)

            if match:
                break


    def extract_section(self, start, end, fields, limit = None):
        """
        Extract data values (int, float) in fields from a section
        marked with start and end.  Do not read farther than limit.
        """

        inside = False
        result = []

        for line in self:
            if limit and re.search(limit, line):
                break

            if re.search(start, line):
                inside = True

            if inside:
                if re.search(end, line):
                    break

                for field in fields:
                    # FIXME: a) assumes fields are only integers and floats
                    #        b) return fields in order of requested fields and not
                    #           in order of occurance in input (accumulate and
                    #           then parse?)
                    match = re.search(' %s\s+=\s+(\*+|%s|\d+)'
                                     % (field, _FP_RE), line)
                    if match:
                        if re.search('\*+', match.group(1) ):
                            raise SystemExit('Cannot parse value: %s' %
                                             match.group(1))
                        elif re.search('%s' % _FP_RE, match.group(1) ):
                            result.append(float(match.group(1) ) )
                        else:
                            result.append(int(match.group(1) ) )

        return result


    def pushback(self):
       """Reposition one line back."""
       self.lineno -= 1
       self.fh.seek(self.last_pos)


    def __iter__(self):
        return self


    def next(self):
        self.lineno += 1
        curr_pos = self.fh.tell()

        # FIXME: just a bad hack?
        if curr_pos == self.filesize:
            raise StopIteration

        self.last_pos = curr_pos

        # NOTE: can't mix next() with seek()
        return self.fh.readline()


    def __enter__(self):
        return self


    def __exit__(self, typ, value, traceback):
        self.fh.close()


class OnlineAvVar(object):
   '''Online algorithm to compute mean and variance.'''

   def __init__(self):
      self.step = 0
      self.mean = 0.0
      #self.M2 = 0.0

   def accumulate(self, x):
      '''Accumulate data points to compute mean and variance on-the-fly.'''

      self.step += 1

      delta = x - self.mean

      self.mean += delta / self.step
      #self.M2 += delta * (x - self.mean)


def process_mbar_lambdas(sp):
   """Extract the lambda points used to compute MBAR energies."""

   mbar_nlambda = 0
   in_mbar = False
   mbar_lambdas = []
            
   for line in sp:
      if line.startswith('    MBAR - lambda values considered:'):
         in_mbar = True
         continue

      if in_mbar:
         if line.startswith('    Extra'):
            break

         if 'total' in line:
            data = line.split()
            mbar_nlambda = data[0]
            mbar_lambdas.extend(data[2:])
         else:
            mbar_lambdas.extend(line.split())

   # FIXME: report error more appropriately
   if len(mbar_lambdas) != int(mbar_nlambda):
      raise ValueError


   return mbar_lambdas


def getG(A_n, mintime = 3):
   """Computes 'g', the statistical inefficiency, as defined in eq (20)
      in J D Chodera, W C Swope, J W Pitera, C Seok, and K A Dill.
      Use of the weighted histogram analysis method for the analysis
      of simulated and parallel tempering simulations.
      JCTC 3(1):26-41, 2007."""

   N = A_n.size
   dA_n = A_n - A_n.mean()

   g = denominator = N * dA_n.var() / 2

   if g == 0:
      raise SystemExit('The variance is 0; cannot proceed.')

   for t in range(1, mintime + 1):
      N_t = N - t
      g += N_t * (dA_n[0:N_t] * dA_n[t:N]).mean()

   t += 1

   while t < N  - 1:
      N_t = N - t
      C = (dA_n[0:N_t] * dA_n[t:N]).mean()

      if C <= 0:
         break

      g += C * N_t
      t += 1

   g /= denominator

   return g if g > 1 else 1


def uncorrelateAmber(dhdl_k, uncorr_threshold):
   """Retain every 'g'th sample of the original array 'dhdl_k'."""

   g = getG(dhdl_k)

   if g == 1:
      return dhdl_k

   N = dhdl_k.size  # Number of correlated samples.
   N_k = 1 + int(N / g) # Number of uncorrelated samples.

   if int(round(N_k * g - g) ) > N - 1:
      N_k -= 1

   if N_k < uncorr_threshold:
      print('WARNING:\nOnly %s uncorrelated samples found;\n'
            'proceeding with analysis using correlated samples...' % N_k)
      return dhdl_k

   indices = numpy.rint(g*numpy.arange(N_k)).astype(int)

   return dhdl_k[indices]


def _extrapol(x, y, scheme):
   """Simple extrapolation scheme."""

   y0 = None
   y1 = None

   if scheme == 'linear':
      if 0.0 not in x:
         y0 = (x[0] * y[1] - x[1] * y[0]) / (x[0] - x[1])

      if 1.0 not in x:
         y1 = ( (x[-2] - 1.0) * y[-1] + ((1.0 - x[-1]) * y[-2]) ) / (x[-2] - x[-1])
   elif scheme == 'polyfit':
      nf = len(x)
      
      if nf < 6:
         deg = nf - 1
      else:
         deg = 6

      coeffs = numpy.polyfit(x, y, deg)

      if 0.0 not in x:
         y0 = coeffs[-1]

      if 1.0 not in x:
         y1 = sum(coeffs)
   else:
      raise SystemExit('Unsupported extrapolation scheme: %s' % scheme)

   return y0, y1


def readDataAmber(P):

   # To suppress unwanted calls in __main__.
   P.lv_names = ['']

   datafile_tuple = P.datafile_directory, P.prefix, P.suffix
   filenames = glob.glob('%s/%s*%s' % datafile_tuple) # will be sorted later

   if not len(filenames):
      raise SystemExit("\nERROR!\nNo files found within directory '%s' with "
                       "prefix '%s' and suffix '%s': check your inputs."
                       % datafile_tuple)

   mbar_all = defaultdict(list)
   dvdl_all = defaultdict(list)
   dvdl_comps_all = {}
   ncomp = len(DVDL_COMPS)
   
   for filename in filenames:
      print('Loading in data from %s... ' % filename),

      in_comps = False
      finished = False

      mbar_data = []
      dvdl_data = []
      dvdl_comp_data = []

      for cmp in DVDL_COMPS:
         dvdl_comp_data.append(OnlineAvVar() )


      with SectionParser(filename) as sp:
         # NOTE: some sections may not exist
         sp.skip_after('^   2.  CONTROL  DATA  FOR  THE  RUN')
         ntpr, = sp.extract_section('^Nature and format of output:', '^$',
                                    ['ntpr'])
         nstlim, dt = sp.extract_section('Molecular dynamics:', '^$',
                                         ['nstlim', 'dt'])

         try:
            ifsc, clambda = sp.extract_section('^Free energy options:', '^$',
                                               ['ifsc', 'clambda'])
         except ValueError:
            raise SystemExit('File %s does not contain a free energy section' %
                             filename)

         try:
            have_mbar, mbar_ndata = sp.extract_section('^FEP MBAR options:', '^$',
                                                       ['ifmbar',
                                                        'bar_intervall'],
                                                       '^---')
            mbar_ndata = int(nstlim / mbar_ndata)
         except ValueError:
            have_mbar = False
            mbar_ndata = 0


         if have_mbar:
            mbar_lambdas = process_mbar_lambdas(sp)
            clambda_str = '%6.4f' % clambda

            if clambda_str not in mbar_lambdas:
               raise SystemExit('Lambda %s not contained in MBAR lambdas: %s' %
                                (clambda_str, ', '.join(mbar_lambdas) ) )

            mbar_nlambda = len(mbar_lambdas)
            mbar_lambda_idx = mbar_lambdas.index(clambda_str)

         sp.skip_after('^   4.  RESULTS')

         mbar = []

         nenergy = int(nstlim / ntpr)
         nensec = 0
         nenav = 0
         old_nstep = -1
         old_comp_nstep = -1

         for line in sp:
            if have_mbar and line.startswith('MBAR Energy analysis'):
               sp.pushback()
               mbar = sp.extract_section('^MBAR', '^ ---', mbar_lambdas)
               Eref = mbar[mbar_lambda_idx]
               mbar_data.append([E - Eref for E in mbar])

            if 'DV/DL, AVERAGES OVER' in line:
                in_comps = True

            if line.startswith(' NSTEP'):
               sp.pushback()

               if in_comps:
                   result = sp.extract_section('^ NSTEP', '^ ---',
                                               ['NSTEP'] + DVDL_COMPS)

                   if result[0] != old_comp_nstep:
                       for i, E in enumerate(DVDL_COMPS):
                           dvdl_comp_data[i].accumulate(float(result[i+1]) )

                       nenav += 1
                       old_comp_nstep = result[0]

                   in_comps = False
               else:
                   nstep, dvdl = sp.extract_section('^ NSTEP', '^ ---',
                                                    ['NSTEP', 'DV/DL'])

                   if nstep != old_nstep:
                       dvdl_data.append(dvdl)
                       nensec += 1
                       old_nstep = nstep

            if line == '   5.  TIMINGS\n':
               finished = True
               break

      print('%i data points, %i DV/DL averages' % (nensec, nenav) )

      if not finished:
         print('WARNING: prematurely terminated run\n')
         next
           
      if not nensec:
          raise SystemExit('File %s does not contain DV/DL data' %
                           filename)

      mbar_all[clambda] = numpy.append(mbar_all[clambda], mbar_data)
      dvdl_all[clambda] = numpy.append(dvdl_all[clambda], dvdl_data)
      dvdl_comps_all[clambda] = [Es.mean for Es in dvdl_comp_data]


   lv = []
   ave = []
   std = []

   start_from = int(round(P.equiltime / (ntpr * float(dt) ) ) )

   print('\nThe average and standard error of the mean (kcal/mol):')

   if start_from:
       print('(first %s ignored)' %
             ('%i data points' % start_from if start_from > 1 else 'data point') )

   print ('%6s %12s %12s %12s %12s %12s' %
          ('State', 'Lambda', 'N', '(Total N)', '<dv/dl>', 'SEM') )

   # FIXME: do not store data again?
   for i, clambda in enumerate(sorted(dvdl_all.keys() ) ):
      dhdl_k = dvdl_all[clambda][start_from:]
      N = dhdl_k.size

      if P.uncorr_threshold:
         dhdl_k = uncorrelateAmber(dhdl_k, P.uncorr_threshold)

      N_k = dhdl_k.size
      ave_dhdl = numpy.average(dhdl_k)
      std_dhdl = numpy.std(dhdl_k) / numpy.sqrt(N_k - 1)

      print('%6s %12.5f %12s %12s %12.6f %12.6f' %
            (i, clambda, N_k, '(' + str(N) + ')', ave_dhdl, std_dhdl) )

      lv.append(clambda)
      ave.append(ave_dhdl)
      std.append(std_dhdl)

   print


   # sander does not sample end-points...
   y0, y1 = _extrapol(lv, ave, 'polyfit')

   if y0:
      print('Note: adding missing lambda = 0.0: %f' % y0)
      lv.insert(0, 0.0)
      ave.insert(0, y0)
      std.insert(0, 0.0)

   if y1:
      print('Note: adding missing lambda = 1.0: %f' % y1)
      lv.append(1.0)
      ave.append(y1)
      std.append(0.0)


   print("\nThe DV/DL components from gradients of "
         "_every_single_ step (kcal/mol):")

   ene_comp = []
   x_comp = sorted(dvdl_comps_all.keys() )

   for en in sorted(dvdl_comps_all.items() ):
      ene_comp.append(en[1:])

   fmt = 'Lambda ' + '%10s' * ncomp
   print(fmt % tuple(DVDL_COMPS) )

   fmt = '%7.5f' + ' %9.3f' * ncomp


   for clambda in x_comp:
       l = (clambda,) + tuple(dvdl_comps_all[clambda])
       print(fmt % l)

   print('   TI ='),

   for ene in numpy.transpose(ene_comp):
      y_comp = ene[0]

      if not all(y_comp):
         print(' %8.3f' % 0.0),
         continue

      y0, y1 = _extrapol(x_comp, y_comp, 'polyfit')

      if y0:
         x_comp.insert(0, 0.0)
         y_comp.insert(0, y0)

      if y1:
         x_comp.append(1.0)
         y_comp.append(y1)

      print(' %8.3f' % numpy.trapz(y_comp, x_comp) ),
      
   print('\n')


   K = len(lv)

   return (numpy.array(lv).reshape(K, 1),
           P.beta * numpy.array(ave).reshape(K, 1),
           P.beta * numpy.array(std).reshape(K, 1) )
