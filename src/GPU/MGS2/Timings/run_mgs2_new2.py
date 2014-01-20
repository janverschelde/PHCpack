"""
This script contains the timings in dictionary formats for runs
-- see the files run_mgs2_new2_k20c.txt for the raw outputs of time
and the command line arguments of run_mgs2_dd and run_mgs2_qd ---
with run_mgs2_dd and run_mgs2_qd for running Newton's method on 
the chandra system for dimension 1024.
Running the script computes the speedups and
displays the results in a format that is easy to cut and paste
into a LaTeX table.
$p$ &     & $n$  &     real~~  &    user ~~  &      sys~~ & speedup \\ \hline
DD  & CPU & 1024 &  42m41.480s &  42m37.692s &     0.038s & \\
    & GPU & 1024 &     35.220s &     23.000s &    11.961s & \\ \cline{2-7}
    & CPU & 2048 & 341m47.998s & 341m18.009s &     0.362s & \\
    & GPU & 2048 &   4m18.193s &   2m18.193s &  1m33.281s & \\ \hline
QD  & CPU & 1024 & 253m51.126s & 253m24.180s &     4.802s & \\
    & GPU & 1024 &  15m11.362s &   9m28.399s &  5m41.532s & \\ \cline{2-7}
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
Cdd[1024] = {'real':'42m41.480s', 'user':'42m37.692s', 'sys':'0m0.038s'}
Cdd[2048] = {'real':'341m47.998', 'user':'341m18.009s', 'sys':'0m0.362s'}
Gdd = {}
Gdd[(1024,64)] = {'real':'0m35.220s', 'user':'0m23.000s', 'sys':'0m11.961s'}
Gdd[(2048,64)] = {'real':'4m18.193s', 'user':'2m18.193s', 'sys':'1m33.281s'}
Cqd = {}
Cqd[1024] = {'real':'253m51.126s', 'user':'253m24.180s', 'sys':'0m4.802s'}
Gqd = {}
Gqd[(1024,64)] = {'real':'15m11.362s', 'user':'9m28.399s', 'sys':'5m41.532s'}

print (1024,'dd cpu'), Cdd[1024]['real'], Cdd[1024]['user'], Cdd[1024]['sys']
print (2048,'dd cpu'), Cdd[2048]['real'], Cdd[2048]['user'], Cdd[2048]['sys']
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
