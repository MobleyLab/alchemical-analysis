######################################################################
#
# Alchemical Analysis: An open tool implementing some recommended practices
# for analyzing alchemical free energy calculations
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: Pavel Klimovich, Hannes H Loeffler, Michael Shirts and David Mobley
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



DVDL_COMPS = ('BOND', 'ANGLE', 'DIHED', '1-4 NB', '1-4 EEL', 'VDWAALS', 'EELEC',
              'RESTRAINT')
RE_FP_NUM= '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
MAX_N_AVER = 100000                     # FIXME: arbitrary
N_AVER = 10                             # FIXME: arbitrary


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

   dvdl = defaultdict(list)
   comps = {}
   comp_lines = []
   ncomp = len(DVDL_COMPS)
   nwarn = 0

   for filename in filenames:
      print 'Loading in data from %s... ' % filename,

      data = []
      Ecomps = []

      for cmp in DVDL_COMPS:
         Ecomps.append(OnlineAvVar() )

      ctrl_data = False
      in_dvdl = False
      in_comps = False
      finished = False
      reduced = False
      old_nsteps = ''
      mden_file = ''
      tot_nsteps = 0
      ndvdl = 0
      nstlim = 0
      dt = 0 

      # FIXME: change strategy:
      #        MDEN: file is not in sync with  MDCRD and MDVEL but one step
      #          behind, manual says that's why it is rarely written,
      #          instantenous values(?)
      #        MDOUT: writes DV/DL every ntpr steps, instantenous values
      #          setting ntave=N computes averages over the last N steps
      #          (ene_avg_sampling should not be set)
      with open(filename, 'rb') as out:
         for line in out:
            if ' MDEN:' in line:        # pmemd adds another space...
               mden_file = line.split()[2]

            if line == '   2.  CONTROL  DATA  FOR  THE  RUN\n':
               ctrl_data = True

            if ctrl_data:
               if line.startswith('     clambda ='):
                  clambda = float(line.split()[2][:-1])

               if 'nstlim  =' in line:
                  nstlim = float(line.split()[2][:-1])

               if 'dt      =' in line:
                  dt = float(line.split()[5][:-1])

            if line.startswith('Summary of dvdl values over'):
               in_dvdl = True
               continue

            if in_dvdl:
               if line.startswith('End of dvdl summary'):
                  in_dvdl = False
                  continue

               data.append(float(line) )
               ndvdl += 1

            if 'DV/DL, AVERAGES OVER' in line:
               nsteps = re.search('OVER\s+(.*?)\s+STEPS', line).group(1)

               if nsteps == old_nsteps or not old_nsteps:
                  in_comps = True

                  # Fortran format may be exceeded
                  try:
                     tot_nsteps += int(nsteps)
                  except ValueError:
                     pass
               else:
                  in_comps = False

               old_nsteps = nsteps
               continue

            if in_comps and 'DV/DL  =' in line:
               in_comps = False
               lines = ''.join(comp_lines)

               # NOTE: assumes all terms are present
               for i, ene_comp in enumerate(DVDL_COMPS):
                  E = re.search('%s\s*=\s*(%s)' % (ene_comp, RE_FP_NUM), lines)
                  Ecomps[i].accumulate(float(E.group(1) ) )

               comp_lines = []

            if in_comps:
               comp_lines.append(line)

            if line == '   5.  TIMINGS\n':
               finished = True
               break

      if not ndvdl:
         print ('no DV/DL summary found!')
         mden_file = os.path.join(os.path.split(filename)[0], mden_file)

         with open(mden_file, 'rb') as en:
            print 'Loading in data from %s... ' % mden_file,

            for line in en:
               if line.startswith('L9') and not 'dV/dlambda' in line:
                  try:
                     data.append(float(line.split()[5]) )
                  except IndexError:
                     print ('WARNING: %s does not contain gradients' %
                            mden_file)
                     nwarn += 1
                     next

                  ndvdl += 1

      print '%i data points' % ndvdl

      # FIXME: check for exceeded Fortran format
      if tot_nsteps < nstlim and tot_nsteps > 0:
         print ('WARNING: DV/DL components for %i steps only' % tot_nsteps)
         nwarn += 1

      if not finished:
         print ('WARNING: prematurely terminated run')
         nwarn += 1
         next

      # reduce data and store in numpy array to save on memory
      ndata = len(data)

      # FIXME: arbitrary sizes
      if ndata > MAX_N_AVER:
         reduced = True
         res = []
         overhang = ndata % N_AVER

         for i in range(0, ndata - overhang, N_AVER):
            res.append(reduce(lambda x, y: (x+y), data[i:i+N_AVER]) / N_AVER)

         # FIXME: not really a clean way
         if overhang:
            res.append(reduce(lambda x, y: (x+y), data[i+N_AVER:]) / overhang)

         del(data)
      else:
         res = data   
      
      dvdl[clambda] = numpy.append(dvdl[clambda], res)

      comps[clambda] = [Es.mean for Es in Ecomps]

   if reduced:
      print '\nSome or all data have been reduced\n'

   lv = []
   ave = []
   std = []

   start_from = int(P.equiltime / float(dt) )

   if reduced:
      start_from = int(start_from / N_AVER)


   print('\nThe average and standard error of the mean in raw data units:\n'
         '(first %i data points ignored)' % start_from)
   print ('%6s %12s %12s %12s %12s %12s' %
          ('State', 'Lambda', 'N', '(Total N)', '<dv/dl>', 'SEM') )

   # FIXME: do not store data again!
   for i, clambda in enumerate(sorted(dvdl.keys() ) ):
      dhdl_k = dvdl[clambda][start_from:]
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

   y0, y1 = _extrapol(lv, ave, 'polyfit')

   if y0:
      print 'Note: adding missing lambda = 0.0: %f' % y0
      lv.insert(0, 0.0)
      ave.insert(0, y0)
      std.insert(0, 0.0)

   if y1:
      print 'Note: adding missing lambda = 1.0: %f' % y1
      lv.append(1.0)
      ave.append(y1)
      std.append(0.0)


   if old_nsteps:
      print "\nThe DV/DL components:"

      fmt = 'Lambda ' + '%10s' * ncomp
      print (fmt % DVDL_COMPS)

      fmt = '%7.5f' + ' %9.3f' * ncomp

      for clambda in sorted(comps.keys() ):
         l = (clambda,) + tuple(comps[clambda])
         print fmt % l

      print

   K = len(lv)

   if nwarn:
      print('\nWARNING: %i warning%s been issued, check integrity of '
            'data and files\n' % (nwarn, ' has' if nwarn == 1 else 's have') )


   return (numpy.array(lv).reshape(K, 1),
           P.beta * numpy.array(ave).reshape(K, 1),
           P.beta * numpy.array(std).reshape(K, 1) )
