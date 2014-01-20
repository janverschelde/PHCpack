"""
This script contains the timings in dictionary formats for runs
-- see the files run_mgs2_d_1024_k20c.txt, run_mgs2_d_2048_k20c.txt,
   run_mgs2_d_3072_k20c.txt, run_mgs2_d_4096_k20c.txt, and
   run_mgs2_d_5120_k20c.txt for the raw outputs of time
   and the command line arguments of run_mgs2_d ---
with run_mgs2_d for dimensions 1024, 2048, 3072, 4096, 5120,
on the host and accelerated for various block sizes.
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
C[1024] = {'real':'0m39.872s', 'user':'0m39.649s', 'sys':'0m0.028s'}
G = {}
G[(1024,32)] = {'real':'0m5.067s', 'user':'0m3.138s', 'sys':'0m1.501s'}
G[(1024,64)] = {'real':'0m3.045s', 'user':'0m1.747s', 'sys':'0m1.099s'}
G[(1024,96)] = {'real':'0m2.544s', 'user':'0m1.434s', 'sys':'0m0.966s'}
G[(1024,128)] = {'real':'0m2.159s', 'user':'0m1.116s', 'sys':'0m0.801s'}
G[(1024,160)] = {'real':'0m2.140s', 'user':'0m1.083s', 'sys':'0m0.859s'}
G[(1024,192)] = {'real':'0m1.971s', 'user':'0m1.039s', 'sys':'0m0.782s'}
G[(1024,224)] = {'real':'0m1.862s', 'user':'0m0.874s', 'sys':'0m0.794s'}
G[(1024,256)] = {'real':'0m1.712s', 'user':'0m0.807s', 'sys':'0m0.744s'}

C[2048] = {'real':'5m23.402s', 'user':'5m22.850s', 'sys':'0m0.075s'}
G[(2048,64)] = {'real':'0m18.112s', 'user':'0m10.892s', 'sys':'0m6.950s'}
G[(2048,96)] = {'real':'0m14.568s', 'user':'0m8.866s', 'sys':'0m5.516s'}
G[(2048,128)] = {'real':'0m11.414s', 'user':'0m6.841s', 'sys':'0m4.342s'}
G[(2048,160)] = {'real':'0m10.610s', 'user':'0m6.335s', 'sys':'0m4.042s'}
G[(2048,192)] = {'real':'0m9.595s', 'user':'0m5.717s', 'sys':'0m3.634s'}
G[(2048,224)] = {'real':'0m9.244s', 'user':'0m5.530s', 'sys':'0m3.461s'}
G[(2048,256)] = {'real':'0m8.098s', 'user':'0m4.768s', 'sys':'0m3.147s'}

C[3072] = {'real':'18m19.975s', 'user':'18m18.073s', 'sys':'0m0.289s'}
G[(3072,96)] = {'real':'0m44.692s', 'user':'0m26.928s', 'sys':'0m17.518s'}
G[(3072,128)] = {'real':'0m35.435s', 'user':'0m21.003s', 'sys':'0m14.198s'}
G[(3072,160)] = {'real':'0m33.395s', 'user':'0m19.848s', 'sys':'0m13.251s'}
G[(3072,192)] = {'real':'0m28.310s', 'user':'0m17.109s', 'sys':'0m11.002s'}
G[(3072,224)] = {'real':'0m26.342s', 'user':'0m15.780s', 'sys':'0m10.390s'}
G[(3072,256)] = {'real':'0m24.411s', 'user':'0m14.616s', 'sys':'0m9.542s'}

C[4096] = {'real':'43m18.139s', 'user':'43m13.184s', 'sys':'0m0.328s'}
G[(4096,128)] = {'real':'1m20.761s', 'user':'0m47.926s', 'sys':'0m32.490s'}
G[(4096,160)] = {'real':'1m14.422s', 'user':'0m44.555s', 'sys':'0m29.596s'}
G[(4096,192)] = {'real':'1m6.437s', 'user':'0m43.248s', 'sys':'0m22.877s'}
G[(4096,224)] = {'real':'1m1.094s', 'user':'0m39.699s', 'sys':'0m21.001s'}
G[(4096,256)] = {'real':'0m55.140s', 'user':'0m32.937s', 'sys':'0m21.898s'}

C[5120] = {'real':'83m59.885s', 'user':'83m50.417s', 'sys':'0m0.464s'}
G[5120,160] = {'real':'2m20.774s', 'user':'1m20.936s', 'sys':'0m59.490s'}
G[5120,192] = {'real':'2m4.366s', 'user':'1m13.106s', 'sys':'0m50.782s'}
G[5120,224] = {'real':'1m53.530s', 'user':'1m7.168s', 'sys':'0m45.975s'}
G[5120,256] = {'real':'1m45.579s', 'user':'1m3.446s', 'sys':'0m41.842s'}

print (1024,'cpu'), C[1024]['real'], C[1024]['user'], C[1024]['sys']
print (2048,'cpu'), C[2048]['real'], C[2048]['user'], C[2048]['sys']
print (3072,'cpu'), C[3072]['real'], C[3072]['user'], C[3072]['sys']
print (4096,'cpu'), C[4096]['real'], C[4096]['user'], C[4096]['sys']
print (5120,'cpu'), C[5120]['real'], C[5120]['user'], C[5120]['sys']
S = list(G.keys())
S.sort()
for L in S:
   spdup = seconds(C[L[0]]['real'])/seconds(G[L]['real'])
   print L, G[L]['real'], G[L]['user'], G[L]['sys'], ':' , \
         seconds(C[L[0]]['real']) , '/' , seconds(G[L]['real']), \
         '=', spdup
