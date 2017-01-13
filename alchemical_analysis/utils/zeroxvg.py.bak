import numpy,os,sys

def zero_output(K,P):
   """Adapted from totalEnergies() function to write out all zeros in results.txt file"""
   print "Generating results.txt with all 0"
   df_allk = list(); ddf_allk = list()
   for k in range(K-1):
      df=dict(); ddf=dict()
      for name in P.methods:
         df[name] = 0.0
         ddf[name] = 0.0
      df_allk = numpy.append(df_allk,df)
      ddf_allk = numpy.append(ddf_allk,ddf)
   segments      = ['Coulomb'  , 'vdWaals'  , 'TOTAL']
   def printLine(str1, str2, d1=None, d2=None):
      text = str1
      for name in P.methods:
         if d1 == 'plain':
            text += ' ' + str2
         if d1 == 'name':
            text += ' ' + str2 % (name, P.units)
         if d1 and d2:
            text += ' ' + str2 % (d1[name]/P.beta_report, d2[name]/P.beta_report)
      outtext.append(text + '\n')
      return

   d = P.decimal
   str_dash  = (d+7 + 6 + d+2)*'-'
   str_dat   = ('X%d.%df  +-  X%d.%df' % (d+7, d, d+2, d)).replace('X', '%')
   str_names = ('X%ds X-%ds' % (d+6, d+8)).replace('X', '%')
   outtext = []
   printLine(12*'-', str_dash, 'plain')
   printLine('%-12s' % '   States', str_names, 'name')
   printLine(12*'-', str_dash, 'plain')
   for k in range(K-1):
      printLine('%4d -- %-4d' % (k, k+1), str_dat, df_allk[k], ddf_allk[k])
   printLine(12*'-', str_dash, 'plain')
   w = 12 + (1+len(str_dash))*len(P.methods)
   str_align = '{:I^%d}' % w
   if len(P.lv_names)>1:
      for i in range(len(segments)):
         printLine('%9s:  ' % segments[i], str_dat, df_allk[i], ddf_allk[i])
   else:
      printLine('%9s:  ' % segments[-1], str_dat, 0.000, 0.000)
   # Store results.
   outfile = open(os.path.join(P.output_directory, 'results.txt'), 'w')
   outfile.write('# Command line was: %s\n\n' % ' '.join(sys.argv) )
   outfile.writelines(outtext)
   outfile.close()

def zero_dFt(K,P,nsnapshots):

   print "Generating dF_t.txt with all 0"

   # Define a list of bForwrev equidistant time frames at which the free energy is to be estimated; count up the snapshots embounded between the time frames.
   n_tf = P.bForwrev + 1
   nss_tf = numpy.zeros([n_tf, K], int)
   increment = 1./(n_tf-1)
   if P.bExpanded:
      from collections import Counter # for counting elements in an array
      tf = numpy.linspace(0.0, 1.0, n_tf)*(numpy.sum(nsnapshots)-1)+1
      tf[0] = 0
      for i in range(n_tf-1):
         nss = Counter(extract_states[tf[i]:tf[i+1]])
         nss_tf[i+1] = numpy.array([nss[j] for j in range(K)])
   else:
      tf = numpy.linspace(0.0, 1.0, n_tf)*(max(nsnapshots)-1)+1
      tf[0] = 0
      for i in range(n_tf-1):
         nss_tf[i+1] = numpy.array([min(j, tf[i+1]) for j in nsnapshots]) - numpy.sum(nss_tf[:i+1],axis=0)
   # Define the real time scale (in ps) rather than a snapshot sequence.
   ts = ["%.1f" % ((i-(i!=tf[0]))*P.snap_size[0] + P.equiltime) for i in tf]
   # Initialize arrays to store data points to be plotted.
   F_df  = numpy.zeros(n_tf-1, float)
   F_ddf = numpy.zeros(n_tf-1, float)
   R_df  = numpy.zeros(n_tf-1, float)
   R_ddf = numpy.zeros(n_tf-1, float)

   outtext = ["%12s %10s %-10s %17s %10s %s\n" % ('Time (ps)', 'Forward', P.units, 'Time (ps)', 'Reverse', P.units)]
   outtext+= ["%10s %11.3f +- %5.3f %18s %11.3f +- %5.3f\n" % (ts[1:][i], F_df[i], F_ddf[i], ts[:-1][i], R_df[i], R_ddf[i]) for i in range(len(F_df))]
   outfile = open(os.path.join(P.output_directory, 'dF_t.txt'), 'w'); outfile.writelines(outtext); outfile.close()

