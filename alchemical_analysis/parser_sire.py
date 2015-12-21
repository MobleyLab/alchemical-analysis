######################################################################
# Alchemical Analysis: An open tool implementing some recommended practices for analyzing alchemical free energy calculations
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: Pavel Klimovich, Michael Shirts and David Mobley
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2.1 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, see <http://www.gnu.org/licenses/>.
######################################################################

import numpy
from glob import glob # for pathname matching

import unixlike       # some implemented unixlike commands

#===================================================================================================
# FUNCTIONS: This is the Sire lambda gradient file parser.
#===================================================================================================

def readDataSire(P):
   """Read in Sire's output files; return nsnapshots, lv, dhdlt, and u_klt."""

   class F:
      """This is the object to be built on the filename."""

      def __init__(self, filename):
         self.filename   = filename
         self.skip_lines = 0  # Number of lines from the top that are to be skipped.
         snap_size       = [] # Time from first two snapshots to determine snapshot's size.

         print "Reading metadata from %-*s" % (len_fstring+1, self.filename+';'),
         with open(self.filename,'r') as infile:
            for line in infile:

               if line.startswith('#'):
                  self.skip_lines += 1
                  elements = line.split()
                  if 'lambba_val.val' in elements:
                     self.lv = elements[-1]
                     lv.append(elements[-1:])
               else:
                  snap_size.append(float(line.split()[0]))
                  if len(snap_size) > 1:
                     self.snap_size = numpy.diff(snap_size)[0]
                     break
            equilsnapshots  = int(P.equiltime/self.snap_size)
            self.skip_lines += equilsnapshots
            nsnapshots.append(unixlike.wcPy(infile) + 2 - equilsnapshots)
            print "first %s ps (%s snapshots) will be discarded due to equilibration..." % (P.equiltime, equilsnapshots)

      def loadtxtSire(self, state):
         print "Loading in data from %-*s (state %d) ..." % (len_fstring, self.filename, state)
         try: #Numpy 1.1 removes skiprows
            dhdlt[state, :, :nsnapshots[state]] = numpy.genfromtxt(self.filename, dtype=float, skiprows=self.skip_lines, usecols=1)
         except TypeError:
            dhdlt[state, :, :nsnapshots[state]] = numpy.genfromtxt(self.filename, dtype=float, skip_header=self.skip_lines, usecols=1)
         return

   # Preliminaries I-III: Sort the dhdl.xvg files; read in the header; count up the equilibrated snapshots.
   datafile_tuple = P.datafile_directory, P.prefix, P.suffix
   fs = glob('%s/%s*%s' % datafile_tuple)
   K = len(fs)
   if not K:
      raise SystemExit("\nERROR!\nNo files found within directory '%s' with prefix '%s' and suffix '%s': check your inputs." % datafile_tuple)
   len_fstring = max([len(i) for i in fs])

   n_components = 1
   lv           = []
   nsnapshots   = []
   P.lv_names   = ['']
   fs = [ F(filename) for filename in fs ]
   print ''
   fs = sorted(fs, key=lambda f: f.lv)
   lv, nsnapshots = zip(*sorted(zip(lv, nsnapshots)))
   lv = numpy.array(lv, float)

   # Preliminaries IV: Load in equilibrated data.
   maxn  = max(nsnapshots)                                   # maximum number of the equilibrated snapshots from any state
   dhdlt = numpy.zeros([K,n_components,int(maxn)], float)    # dhdlt[k,n,t] is the derivative of energy component n with respect to state k of snapshot t
   u_klt = None

   for nf, f in enumerate(fs):
      f.loadtxtSire(nf)
   print ''

   return nsnapshots, lv, dhdlt, u_klt
