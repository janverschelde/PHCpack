# This scripts runs the modified Gram-Schmidt method 
# for the largest dimensions: n = 256 for double,
# n = 128 for double double, and n = 80 for quad double.

from os import system

def run_maxdim():
   """
   Runs the modified Gram-Schmidt method for the largest
   dimensions with double (256), double double (128),
   and quad double (85). 
   """
   rd = '../run_mgs_d'    # executable for double precision
   rdd = '../run_mgs_dd'  # executable for double doubles
   rqd = '../run_mgs_qd'  # executable for quad doubles
   nd = 256               # maximal dimension for doubles
   ndd = 128              # maximal dimension for double doubles
   nqd = 85               # maximal dimension for quad doubles
   freq = 10000           # frequency
   print 'running %s for dimension %d ...' % (rd,nd)
   cmd = 'time ' + rd + ' ' + str(nd) + ' ' + str(nd) \
       + ' ' + str(freq) + ' 0'
   print 'launching', cmd
   system(cmd)
   print 'running %s for dimension %d ...' % (rdd,ndd)
   cmd = 'time ' + rdd + ' ' + str(ndd) + ' ' + str(ndd) \
       + ' ' + str(freq) + ' 0'
   print 'launching', cmd
   system(cmd)
   print 'running %s for dimension %d ...' % (rqd,nqd)
   cmd = 'time ' + rqd + ' ' + str(nqd) + ' ' + str(nqd) \
       + ' ' + str(freq) + ' 0'
   print 'launching', cmd
   system(cmd)

run_maxdim()
