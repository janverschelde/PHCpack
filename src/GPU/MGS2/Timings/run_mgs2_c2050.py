"""
This script processes timings with the C2050 GPU,
with the CPU frequency scaled at 2.66GHz to match
the 2.60GHz clock speed of the host of the K20C.
See the file run_mgs2_c2050.txt for the raw data.
The C2050 has fewer cores per streaming multiprocessor
than the K20C, but they run at higher clock speed.
"""

def seconds(strtime):
   """
   Splits of the 'm' and the 's' in the time string
   to return the number of seconds as a float.
   """
   m = strtime.split('m')
   s = int(m[0])*60 + float(m[1][0:-1])
   return s

# for dimension 1024 only
C2050d = {}
C2050d[32] = {'real':'0m5.012s', 'user':'0m2.654s', 'sys':'0m1.861s'}
C2050d[64] = {'real':'0m3.494s', 'user':'0m1.533s', 'sys':'0m1.451s'}
C2050d[96] = {'real':'0m3.154s', 'user':'0m1.501s', 'sys':'0m1.596s'}
C2050d[128] = {'real':'0m2.902s', 'user':'0m1.072s', 'sys':'0m1.247s'}
C2050d[160] = {'real':'0m3.125s', 'user':'0m1.229s', 'sys':'0m1.236s'}
C2050d[192] = {'real':'0m3.126s', 'user':'0m1.450s', 'sys':'0m1.427s'}
C2050d[224] = {'real':'0m3.184s', 'user':'0m1.534s', 'sys':'0m1.594s'}
C2050d[256] = {'real':'0m3.029s', 'user':'0m1.430s', 'sys':'0m1.542s'}
# for double double arithmetic
C2050dd = {}
C2050dd[32] = {'real':'0m10.866s', 'user':'0m7.260s', 'sys':'0m3.468s'}
C2050dd[64] = {'real':'0m6.834s', 'user':'0m4.042s', 'sys':'0m2.339s'}
C2050dd[96] = {'real':'0m6.259s', 'user':'0m3.421s', 'sys':'0m2.128s'}
C2050dd[128] = {'real':'0m5.880s', 'user':'0m3.543s', 'sys':'0m2.218s'}
# for quad double arithmetic
C2050qd = {}
C2050qd[32] = {'real':'3m47.729s', 'user':'2m58.814s', 'sys':'0m48.703s'}
C2050qd[64] = {'real':'2m34.034s', 'user':'1m50.693s', 'sys':'0m41.865s'}
# entries copied from run_mgs2_d_k20c.py for speedups
K20d = {}
K20d[(1024,32)] = {'real':'0m5.067s', 'user':'0m3.138s', 'sys':'0m1.501s'}
K20d[(1024,64)] = {'real':'0m3.045s', 'user':'0m1.747s', 'sys':'0m1.099s'}
K20d[(1024,96)] = {'real':'0m2.544s', 'user':'0m1.434s', 'sys':'0m0.966s'}
K20d[(1024,128)] = {'real':'0m2.159s', 'user':'0m1.116s', 'sys':'0m0.801s'}
K20d[(1024,160)] = {'real':'0m2.140s', 'user':'0m1.083s', 'sys':'0m0.859s'}
K20d[(1024,192)] = {'real':'0m1.971s', 'user':'0m1.039s', 'sys':'0m0.782s'}
K20d[(1024,224)] = {'real':'0m1.862s', 'user':'0m0.874s', 'sys':'0m0.794s'}
K20d[(1024,256)] = {'real':'0m1.712s', 'user':'0m0.807s', 'sys':'0m0.744s'}
# entries copied from run_mgs2_dd_k20c.py for speedups
K20dd = {}
K20dd[(1024,32)] = {'real':'0m11.846s', 'user':'0m8.563s', 'sys':'0m3.050s'}
K20dd[(1024,64)] = {'real':'0m6.608s', 'user':'0m4.437s', 'sys':'0m1.967s'}
K20dd[(1024,96)] = {'real':'0m5.359s', 'user':'0m3.440s', 'sys':'0m1.665s'}
K20dd[(1024,128)] = {'real':'0m4.270s', 'user':'0m2.749s', 'sys':'0m1.320s'}
# entries copied from run_mgs2_qd_k20c.py for speedups
K20qd = {}
K20qd[(1024,32)] = {'real':'4m12.664s', 'user':'3m22.920s', 'sys':'0m49.078s'}
K20qd[(1024,64)] = {'real':'2m30.463s', 'user':'1m53.848s', 'sys':'0m36.221s'}

S = list(C2050d.keys())
S.sort()
for L in S:
   LL = (1024,L)
   sys_spdup = seconds(C2050d[L]['sys'])/seconds(K20d[LL]['sys'])
   real_spdup = seconds(C2050d[L]['real'])/seconds(K20d[LL]['real'])
   print L, C2050d[L]['real'], C2050d[L]['user'], C2050d[L]['sys'], ':' , \
         seconds(C2050d[L]['sys']), '/', seconds(K20d[LL]['sys']), \
         '=', sys_spdup, ':', \
         seconds(C2050d[L]['real']), '/', seconds(K20d[LL]['real']), \
         '=', real_spdup

S = list(C2050dd.keys())
S.sort()
for L in S:
   LL = (1024,L)
   sys_spdup = seconds(C2050d[L]['sys'])/seconds(K20dd[LL]['sys'])
   real_spdup = seconds(C2050d[L]['real'])/seconds(K20dd[LL]['real'])
   print L, C2050dd[L]['real'], C2050dd[L]['user'], C2050dd[L]['sys'], ':', \
         seconds(C2050d[L]['sys']), '/', seconds(K20dd[LL]['sys']), \
         '=', sys_spdup, \
         seconds(C2050d[L]['real']), '/', seconds(K20dd[LL]['real']), \
         '=', real_spdup

S = list(C2050qd.keys())
S.sort()
for L in S:
   LL = (1024,L)
   sys_spdup = seconds(C2050qd[L]['sys'])/seconds(K20qd[LL]['sys'])
   real_spdup = seconds(C2050qd[L]['real'])/seconds(K20qd[LL]['real'])
   print L, C2050qd[L]['real'], C2050qd[L]['user'], C2050qd[L]['sys'], ':', \
         seconds(C2050qd[L]['sys']) , '/' , seconds(K20qd[LL]['sys']), \
         '=', sys_spdup, \
         seconds(C2050qd[L]['real']) , '/' , seconds(K20qd[LL]['real']), \
         '=', real_spdup
