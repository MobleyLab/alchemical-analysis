import numpy,os,sys

def zero_output(K,P):
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
      """Fills out the results table linewise."""
      #print str1,
      text = str1
      for name in P.methods:
         if d1 == 'plain':
            #print str2,
            text += ' ' + str2
         if d1 == 'name':
            #print str2 % (name, P.units),
            text += ' ' + str2 % (name, P.units)
         if d1 and d2:
            #print str2 % (d1[name]/P.beta_report, d2[name]/P.beta_report),
            text += ' ' + str2 % (d1[name]/P.beta_report, d2[name]/P.beta_report)
      #print ''
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
