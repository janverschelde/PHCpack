n = range(16,257,16)
c = [2.26,16.48,53.03,123.80,238.83,409.30,649.30,962.17,1363.57,1867.36,\
     2480.18, 3221.12, 4080.77, 5099.98, 6253.65, 7575.85]
g = [3.36, 5.58,10.26, 15.16, 22.45, 28.64,36.59,45.29,59.48,70.54,\
     84.26, 97.62, 112.19, 128.11, 146.31, 164.33]
for k in range(0,len(c)):
   s = '%4d' % n[k]
   s = s + ' ' + ('%5.2f' % c[k]).rjust(8)
   s = s + ' ' + ('%5.2f' % g[k]).rjust(7)
   s = s + ' ' + ('%5.2f' % (c[k]/g[k])).rjust(6)
   print s

# Below is the screen output of the Python script run_mgs_timings.py
# comparing the running times of the modified Gram-Schmidt method 
# for the NVIDIA Tesla K20C GPU with a 2.60GHz CPU,
# done with complex double precision arithmetic.

# $ python run_mgs_timings.py
# running the modified Gram-Schmidt method
# type precision level (d, dd, or qd) : d
# give start dimension (>= 2) : 16
# give stop dimension (<= 256) : 256
# give the step : 16
# give frequency per dimension (>= 1) : 10000
# give the name of the output file : run_mgs_d_times.txt
# running for dimension 16 and frequency 10000 on CPU
# ['real 2.25', 'user 2.23', 'sys 0.02']
# running for dimension 16 and frequency 10000 on GPU
# ['real 1.23', 'user 0.45', 'sys 0.63']
# running for dimension 32 and frequency 10000 on CPU
# ['real 16.31', 'user 16.20', 'sys 0.08']
# running for dimension 32 and frequency 10000 on GPU
# ['real 6.51', 'user 3.88', 'sys 2.50']
# running for dimension 48 and frequency 10000 on CPU
# ['real 53.12', 'user 52.88', 'sys 0.15']
# running for dimension 48 and frequency 10000 on GPU
# ['real 10.24', 'user 6.40', 'sys 3.67']
# running for dimension 64 and frequency 10000 on CPU
# ['real 123.70', 'user 123.19', 'sys 0.31']
# running for dimension 64 and frequency 10000 on GPU
# ['real 2.55', 'user 1.66', 'sys 0.71']
# running for dimension 80 and frequency 10000 on CPU
# ['real 239.62', 'user 238.34', 'sys 0.88']
# running for dimension 80 and frequency 10000 on GPU
# ['real 22.41', 'user 13.94', 'sys 8.20']
# running for dimension 96 and frequency 10000 on CPU
# ['real 409.30', 'user 407.17', 'sys 1.47']
# running for dimension 96 and frequency 10000 on GPU
# ['real 28.64', 'user 17.94', 'sys 10.50']
# running for dimension 112 and frequency 10000 on CPU
# ['real 649.00', 'user 645.68', 'sys 2.29']
# running for dimension 112 and frequency 10000 on GPU
# ['real 36.59', 'user 22.38', 'sys 13.87']
# running for dimension 128 and frequency 10000 on CPU
# ['real 962.17', 'user 957.79', 'sys 2.71']
# running for dimension 128 and frequency 10000 on GPU
# ['real 45.29', 'user 28.02', 'sys 16.99']
# running for dimension 144 and frequency 10000 on CPU
# ['real 1363.57', 'user 1358.10', 'sys 3.30']
# running for dimension 144 and frequency 10000 on GPU
# ['real 59.48', 'user 36.73', 'sys 22.45']
# running for dimension 160 and frequency 10000 on CPU
# ['real 1867.36', 'user 1858.74', 'sys 5.56']
# running for dimension 160 and frequency 10000 on GPU
# ['real 70.54', 'user 42.69', 'sys 27.56']
# running for dimension 176 and frequency 10000 on CPU
# ['real 2480.18', 'user 2469.74', 'sys 6.24']
# running for dimension 176 and frequency 10000 on GPU
# ['real 84.26', 'user 51.40', 'sys 32.52']
# running for dimension 192 and frequency 10000 on CPU
# ['real 3221.12', 'user 3206.63', 'sys 8.81']
# running for dimension 192 and frequency 10000 on GPU
# ['real 97.62', 'user 58.97', 'sys 38.17']
# running for dimension 208 and frequency 10000 on CPU
# ['real 4080.77', 'user 4064.61', 'sys 9.46']
# running for dimension 208 and frequency 10000 on GPU
# ['real 112.19', 'user 67.24', 'sys 44.56']
# running for dimension 224 and frequency 10000 on CPU
# ['real 5099.98', 'user 5077.51', 'sys 13.65']
# running for dimension 224 and frequency 10000 on GPU
# ['real 128.11', 'user 76.88', 'sys 50.85']
# running for dimension 240 and frequency 10000 on CPU
# ['real 6253.65', 'user 6230.79', 'sys 12.75']
# running for dimension 240 and frequency 10000 on GPU
# ['real 146.31', 'user 88.93', 'sys 56.93']
# running for dimension 256 and frequency 10000 on CPU
# ['real 7575.85', 'user 7549.78', 'sys 14.52']
# running for dimension 256 and frequency 10000 on GPU
# ['real 164.33', 'user 99.73', 'sys 64.22']
# $ 
