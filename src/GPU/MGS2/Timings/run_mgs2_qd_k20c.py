def seconds(strtime):
   """
   Splits of the 'm' and the 's' in the time string
   to return the number of seconds as a float.
   """
   m = strtime.split('m')
   s = int(m[0])*60 + float(m[1][0:-1])
   return s


C = {}
C[1024] = {'real':'41m9.521s', 'user':'41m5.699s', 'sys':'0m0.092s'}
C[2048] = {'real':'329m4.188s', 'user':'328m35.558s', 'sys':'0m0.346s'}
G = {}
G[(1024,32)] = {'real':'4m12.664s', 'user':'3m22.920s', 'sys':'0m49.078s'}
G[(1024,64)] = {'real':'2m30.463s', 'user':'1m53.848s', 'sys':'0m36.221s'}
G[(2048,64)] = {'real':'19m20.893s', 'user':'11m48.509s', 'sys':'7m30.749s'}

print (1024,'cpu'), C[1024]['real'], C[1024]['user'], C[1024]['sys']
print (2048,'cpu'), C[2048]['real'], C[2048]['user'], C[2048]['sys']
S = list(G.keys())
S.sort()
for L in S:
   spdup = seconds(C[L[0]]['real'])/seconds(G[L]['real'])
   print L, G[L]['real'], G[L]['user'], G[L]['sys'], ':' , \
         seconds(C[L[0]]['real']) , '/' , seconds(G[L]['real']), \
         '=', spdup
