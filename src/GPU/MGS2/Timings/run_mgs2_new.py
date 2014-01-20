"""
This script contains the timings in dictionary formats for runs
-- see the files run_mgs2_dd_new_k20c.txt and run_mgs2_qd_new_k20c.txt,
for the raw outputs of time and the command line arguments of run_mgs2_d ---
with run_mgs2_dd and run_mgs2_qd for running Newton's method on 
the chandra system for dimension 1024.
Running the script computes the speedups and
displays the results in a format that is easy to cut and paste
into a LaTeX table.
"""

def seconds(strtime):
   """
   Splits of the 'm' and the 's' in the time string
   to return the number of seconds as a float.
   """
   m = strtime.split('m')
   s = int(m[0])*60 + float(m[1][0:-1])
   return s

Cdd = {}
Cdd[1024] = {'real':'42m30.192s', 'user':'42m26.101s', 'sys':'0m0.039s'}
Gdd = {}
Gdd[(1024,64)] = {'real':'1m1.815s', 'user':'0m52.883s', 'sys':'0m8.676s'}
Cqd = {}
Cqd[1024] = {'real':'256m41.262s', 'user':'256m18.870s', 'sys':'0m0.074s'}
Gqd = {}
Gqd[(1024,64)] = {'real':'17m34.150s', 'user':'14m4.030s', 'sys':'3m28.215s'}

print (1024,'dd cpu'), Cdd[1024]['real'], Cdd[1024]['user'], Cdd[1024]['sys']
Sdd = list(Gdd.keys())
Sdd.sort()
for L in Sdd:
   spdup = seconds(Cdd[L[0]]['real'])/seconds(Gdd[L]['real'])
   print L, Gdd[L]['real'], Gdd[L]['user'], Gdd[L]['sys'], ':' , \
         seconds(Cdd[L[0]]['real']) , '/' , seconds(Gdd[L]['real']), \
         '=', spdup
print (1024,'qd cpu'), Cqd[1024]['real'], Cqd[1024]['user'], Cqd[1024]['sys']
Sqd = list(Gqd.keys())
Sqd.sort()
for L in Sqd:
   spdup = seconds(Cqd[L[0]]['real'])/seconds(Gqd[L]['real'])
   print L, Gqd[L]['real'], Gqd[L]['user'], Gqd[L]['sys'], ':' , \
         seconds(Cqd[L[0]]['real']) , '/' , seconds(Gqd[L]['real']), \
         '=', spdup
