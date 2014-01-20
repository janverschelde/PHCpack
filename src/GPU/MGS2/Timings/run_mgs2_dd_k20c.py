"""
This script contains the timings in dictionary formats for runs
-- see the files run_mgs2_dd_1024_k20c.txt, run_mgs2_dd_2048_k20c.txt,
   run_mgs2_dd_3072_k20c.txt, and run_mgs2_dd_4096_k20c.txt
   for the raw outputs of time and the command line arguments of run_mgs2_dd ---
with run_mgs2_dd for dimensions 1024, 2048, 3072, and 4096,
on the host and accelerated for various block sizes.
The dimension 4096 is an outlier since we only have accelerated times.
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

C = {}
C[1024] = {'real':'7m1.064s', 'user':'7m0.199s', 'sys':'0m0.055s'}
G = {}
G[(1024,32)] = {'real':'0m11.846s', 'user':'0m8.563s', 'sys':'0m3.050s'}
G[(1024,64)] = {'real':'0m6.608s', 'user':'0m4.437s', 'sys':'0m1.967s'}
G[(1024,96)] = {'real':'0m5.359s', 'user':'0m3.440s', 'sys':'0m1.665s'}
G[(1024,128)] = {'real':'0m4.270s', 'user':'0m2.749s', 'sys':'0m1.320s'}
C[2048] = {'real':'56m2.627s', 'user':'55m56.677s', 'sys':'0m0.586s'}
G[(2048,64)] = {'real':'0m44.798s', 'user':'0m27.813s', 'sys':'0m16.720s'}
G[(2048,96)] = {'real':'0m34.929s', 'user':'0m21.765s', 'sys':'0m12.860s'}
G[(2048,128)] = {'real':'0m27.039s', 'user':'0m16.838s', 'sys':'0m9.922s'}
C[3072] = {'real':'189m22.950s', 'user':'189m4.729s', 'sys':'0m0.575s'}
G[(3072,96)] = {'real':'1m49.947s', 'user':'1m6.248s', 'sys':'0m43.375s'}
G[(3072,128)] = {'real':'1m26.581s', 'user':'0m51.433s', 'sys':'0m34.724s'}
C[4096] = {'real':'452m53.113s', 'user':'452m12.340s', 'sys':'0m2.046s'}
G[(4096,128)] = {'real':'3m21.074s', 'user':'1m57.756s', 'sys':'1m22.789s'}

print (1024,'cpu'), C[1024]['real'], C[1024]['user'], C[1024]['sys']
print (2048,'cpu'), C[2048]['real'], C[2048]['user'], C[2048]['sys']
print (3072,'cpu'), C[3072]['real'], C[3072]['user'], C[3072]['sys']
print (4096,'cpu'), C[4096]['real'], C[4096]['user'], C[4096]['sys']
S = list(G.keys())
S.sort()
for L in S:
   spdup = seconds(C[L[0]]['real'])/seconds(G[L]['real'])
   print L, G[L]['real'], G[L]['user'], G[L]['sys'], ':' , \
         seconds(C[L[0]]['real']) , '/' , seconds(G[L]['real']), \
         '=', spdup
