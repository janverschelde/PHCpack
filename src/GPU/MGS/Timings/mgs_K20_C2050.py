# times on K20C for double, double double, and quad double 
k20c_real = [ ['2m44.308s', '2m23.335s', '9m23.748s'], \
              ['2m44.368s', '2m23.286s', '9m23.687s'], \
              ['2m44.353s', '2m23.308s', '9m23.668s'], \
              ['2m44.432s', '2m23.360s', '9m24.233s'], \
              ['2m44.286s', '2m23.279s', '9m23.741s'], \
              ['2m44.304s', '2m23.307s', '9m23.812s'], \
              ['2m44.353s', '2m23.284s', '9m23.704s'], \
              ['2m44.326s', '2m23.349s', '9m23.805s'], \
              ['2m44.325s', '2m23.361s', '9m23.670s'], \
              ['2m44.348s', '2m23.288s', '9m23.765s'] ]
k20c_user = [ ['1m36.782s', '1m24.326s', '5m33.092s'], \
              ['1m38.062s', '1m23.755s', '5m29.009s'], \
              ['1m37.999s', '1m23.742s', '5m23.443s'], \
              ['1m38.193s', '1m24.422s', '5m28.469s'], \
              ['1m38.510s', '1m24.647s', '5m23.170s'], \
              ['1m38.315s', '1m24.661s', '5m23.802s'], \
              ['1m37.976s', '1m24.632s', '5m23.448s'], \
              ['1m38.197s', '1m23.963s', '5m28.907s'], \
              ['1m39.346s', '1m24.103s', '5m29.003s'], \
              ['1m39.116s', '1m24.181s', '5m28.569s'] ]
k20c_sys = [  ['1m7.016s', '0m58.493s', '3m49.601s'], \
              ['1m5.663s', '0m59.125s', '3m53.849s'], \
              ['1m5.957s', '0m59.148s', '3m59.116s'], \
              ['1m5.768s', '0m58.631s', '3m54.889s'], \
              ['1m5.386s', '0m58.213s', '3m59.720s'], \
              ['1m5.546s', '0m58.283s', '3m58.971s'], \
              ['1m5.770s', '0m58.343s', '3m59.010s'], \
              ['1m5.514s', '0m58.839s', '3m53.797s'], \
              ['1m4.389s', '0m58.718s', '3m53.721s'], \
              ['1m4.746s', '0m58.754s', '3m54.248s'] ]
C2050_real = [['3m56.248s', '2m41.758s', '11m6.346s'], \
              ['3m56.420s', '2m41.546s', '11m7.573s'], \
              ['3m58.219s', '2m41.130s', '11m5.448s'], \
              ['3m56.545s', '2m41.467s', '11m6.106s'], \
              ['3m57.021s', '2m41.110s', '11m6.459s'], \
              ['3m56.862s', '2m40.925s', '11m6.819s'], \
              ['3m56.518s', '2m41.218s', '11m6.898s'], \
              ['3m56.720s', '2m41.925s', '11m8.060s'], \
              ['3m56.425s', '2m41.213s', '11m7.226s'], \
              ['3m56.367s', '2m41.549s', '11m6.196s'] ]
C2050_user = [['1m33.960s', '1m0.690s', '4m3.236s'], \
              ['1m31.900s', '1m1.897s', '4m7.095s'], \
              ['1m32.641s', '1m0.739s', '4m5.208s'], \
              ['1m32.583s', '1m0.682s', '4m5.856s'], \
              ['1m35.224s', '1m0.955s', '4m8.324s'], \
              ['1m30.009s', '1m0.425s', '4m0.998s'], \
              ['1m31.834s', '1m1.346s', '4m2.226s'], \
              ['1m30.745s', '1m1.570s', '4m4.052s'], \
              ['1m30.961s', '1m0.754s', '4m3.007s'], \
              ['1m31.640s', '1m0.613s', '4m4.975s'] ]
C2050_sys = [ ['2m21.025s', '1m39.732s', '7m1.296s'], \
              ['2m23.286s', '1m38.332s', '6m58.066s'], \
              ['2m23.955s', '1m39.214s', '6m58.573s'], \
              ['2m23.042s', '1m39.268s', '6m58.444s'], \
              ['2m19.801s', '1m38.840s', '6m56.435s'], \
              ['2m25.312s', '1m39.340s', '7m3.937s'], \
              ['2m23.463s', '1m38.600s', '7m3.397s'], \
              ['2m24.639s', '1m38.829s', '7m2.126s'], \
              ['2m24.244s', '1m39.322s', '7m2.619s'], \
              ['2m23.395s', '1m39.305s', '6m59.156s'] ]
print ''
print 'times on the K20C :'
print ''
sum = 0
for L in k20c_real: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_real_d = sum/len(k20c_real)
print 'average real time for double precision :', k20c_avg_real_d, 'seconds'
sum = 0
for L in k20c_real: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_real_dd = sum/len(k20c_real)
print '  average real time for double doubles :', k20c_avg_real_dd, 'seconds'
sum = 0
for L in k20c_real: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_real_qd = sum/len(k20c_real)
print '    average real time for quad doubles :', k20c_avg_real_qd, 'seconds'
sum = 0
for L in k20c_user: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_user_d = sum/len(k20c_user)
print 'average user time for double precision :', k20c_avg_user_d, 'seconds'
sum = 0
for L in k20c_user: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_user_dd = sum/len(k20c_user)
print '  average user time for double doubles :', k20c_avg_user_dd, 'seconds'
sum = 0
for L in k20c_user: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_user_qd = sum/len(k20c_user)
print '    average user time for quad doubles :', k20c_avg_user_qd, 'seconds'
sum = 0
for L in k20c_sys: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_sys_d = sum/len(k20c_sys)
print ' average sys time for double precision :', k20c_avg_sys_d, 'seconds'
sum = 0
for L in k20c_sys: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_sys_dd = sum/len(k20c_sys)
print '   average sys time for double doubles :', k20c_avg_sys_dd, 'seconds'
sum = 0
for L in k20c_sys: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
k20c_avg_sys_qd = sum/len(k20c_sys)
print '     average sys time for quad doubles :', k20c_avg_sys_qd, 'seconds'
print ''
print 'times on the C2050 :'
print ''
sum = 0
for L in C2050_real: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_real_d = sum/len(C2050_real)
print 'average real time for double precision :', C2050_avg_real_d, 'seconds'
sum = 0
for L in C2050_real: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_real_dd = sum/len(C2050_real)
print '  average real time for double doubles :', C2050_avg_real_dd, 'seconds'
sum = 0
for L in C2050_real: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_real_qd = sum/len(C2050_real)
print '    average real time for quad doubles :', C2050_avg_real_qd, 'seconds'
sum = 0
for L in C2050_user: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_user_d = sum/len(C2050_user)
print 'average user time for double precision :', C2050_avg_user_d, 'seconds'
sum = 0
for L in C2050_user: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_user_dd = sum/len(C2050_user)
print '  average user time for double doubles :', C2050_avg_user_dd, 'seconds'
sum = 0
for L in C2050_user: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_user_qd = sum/len(C2050_user)
print '    average user time for quad doubles :', C2050_avg_user_qd, 'seconds'
sum = 0
for L in C2050_sys: 
   t = L[0]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_sys_d = sum/len(C2050_sys)
print ' average sys time for double precision :', C2050_avg_sys_d, 'seconds'
sum = 0
for L in C2050_sys: 
   t = L[1]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_sys_dd = sum/len(C2050_sys)
print '   average sys time for double doubles :', C2050_avg_sys_dd, 'seconds'
sum = 0
for L in C2050_sys: 
   t = L[2]
   s = t.split('m')
   sum = sum + int(s[0])*60 + float(s[1][0:-1])
C2050_avg_sys_qd = sum/len(C2050_sys)
print '     average sys time for quad doubles :', C2050_avg_sys_qd, 'seconds'
print ''
print 'speedups, system times :'
print ''
print '       double precision :', C2050_avg_sys_d/k20c_avg_sys_d
print 'double double precision :', C2050_avg_sys_dd/k20c_avg_sys_dd
print '  quad double precision :', C2050_avg_sys_qd/k20c_avg_sys_qd
