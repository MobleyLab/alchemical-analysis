###==========================###
#   Desmond parser module
#   Module adapted from Gromacs parser
#   Adapted by Nathan M. Lim
###==========================###
import numpy
import os                       # for os interface
import re                       # for regular expressions
from glob import glob           # for pathname matching
from collections import Counter # for counting elements in an array

import unixlike                 # some implemented unixlike commands

#===================================================================================================
# FUNCTIONS: This is the Desmond gibbs.N.dE file parser.
#===================================================================================================

def readDataDesmond(P):

   class F:

      def __init__(self, filename):
         self.filename = filename

      def sortedHelper(self):
         meat = os.path.basename(self.filename).replace(P.prefix, '').replace(P.suffix, '')
         l = [i for i in re.split('\.|-|_', meat) if i]
         try:
            self.state = l[0] = int(l[0]) # Will be of use for selective MBAR analysis.
         except:
            print("\nERROR!\nFile's prefix should be followed by a numerical character. Cannot sort the files.\n")
            raise
         return tuple(l)

      def get_snapsize(self):
         self.skip_lines = 0
         self.lv_names   = ()
         snap_size       = [] # Time from first two snapshots to determine snapshot's size.
         self.lv         = [] # Lambda vectors, e.g. (0, 0), (0.2, 0), (0.5, 0).


         with open(self.filename,'r') as infile:
            for line in infile:
               snap_size.append(float(line.split()[0]))
               if len(snap_size) > 1:
                  self.snap_size = numpy.diff(snap_size)[0]
                  P.snap_size.append(self.snap_size)
                  break

      def iter_loadtxt(self, state):
         def iter_func():
            with open(self.filename, 'r') as infile:
               for _ in range(self.skip_lines):
                  next(infile)
               for line in infile:
                  line = line.split()
                  for item in line:
                     yield item

         def slice_data(data, state=state):
            #Energies stored in:
            #   Reverse: data[1,:]
            #   Forward: data[2,:]
            #Desmond unit input: kcal/mol, conversion factor 4.184kJ/kcal
            #P.beta from alchemical_analysis.py in kJ/mol/K
            #Return: u_klt contains energies of adjacent lambdas only

            data = data.T
            if state == 0:
               u_klt[state, state+1 , :nsnapshots[state]] = data[ 2 , : ]*4.184*P.beta
            elif state == K:
               u_klt[state, state-1 , :nsnapshots[state]] = data[ 2 , : ]*4.184*P.beta
            else:
               u_klt[state, state-1, :nsnapshots[state]] = data[ 1 , :]*4.184*P.beta
               u_klt[state, state+1, :nsnapshots[state]] = data[ 2 , :]*4.184*P.beta
            return

         print("Loading in data from %s (%s) ...") % (self.filename, 'state %d' % state)
         data = numpy.fromiter(iter_func(), dtype=float)
         if not self.len_first == self.len_last:
            data = data[: -self.len_last]
         data = data.reshape((-1, self.len_first))

         slice_data(data)
   #===================================================================================================
   # Preliminaries I: Get LV,Snapsize,consistency check, and skip frames
   #===================================================================================================

   datafile_tuple = P.datafile_directory, P.prefix, P.suffix
   fs = [ F(filename) for filename in glob( '%s/%s*%s' % datafile_tuple ) ]
   n_files = len(fs)

   if not n_files:
      raise SystemExit("\nERROR!\nNo files found within directory '%s' with prefix '%s' and suffix '%s': check your inputs." % datafile_tuple)
   if n_files > 1:
      fs = sorted(fs, key=F.sortedHelper)

   ###Set lambda vector and get snapsize
   lv = []
   P.snap_size = []
   for nf, f in enumerate(fs):
      lv.append( [nf,0] )
      f.get_snapsize()

      P.lv_names = lv_names = f.lv_names
      n_components = len(lv_names)

   lv = numpy.array(lv, float) # *** Lambda vectors.
   K  = len(lv)                # *** Number of lambda states.
   equiltime = P.equiltime
   nsnapshots = numpy.zeros(K, int)

   for nf, f in enumerate(fs):
      ###Check for consistent timestep???
      f.len_first, f.len_last = (len(line.split()) for line in unixlike.tailPy(f.filename, 2))
      bLenConsistency = (f.len_first != f.len_last)

      ###Skip N snapshots
      equilsnapshots  = int(equiltime/f.snap_size)
      f.skip_lines   += equilsnapshots
      nsnapshots[nf] += unixlike.wcPy(f.filename) - f.skip_lines - 1*bLenConsistency
      print("First %s ps (%s snapshots) will be discarded due to equilibration from file %s...") % (equiltime, equilsnapshots, f.filename)

   #===================================================================================================
   # Preliminaries: Load in equilibrated data.
   #===================================================================================================

   maxn  = max(nsnapshots)                                   # maximum number of the equilibrated snapshots from any state
   u_klt = numpy.zeros([K,K+1,int(maxn)], numpy.float64)       # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m
   for nf, f in enumerate(fs):
      f.iter_loadtxt(nf)
   return nsnapshots, lv, u_klt
