n = [  16,    32,    48,     64,     80,     96,    112,    128,     144,     160,     176,     192,     208,     224,     240]
c = [2.01, 14.61, 47.80, 112.60, 217.52, 373.06, 589.35, 876.11, 1243.26, 1701.57, 2260.07, 2932.15, 3722.77, 4641.71, 5703.77]
g = [4.11,  6.52, 11.11,  15.38,  22.89,  30.43,  40.82,  49.10,   67.41,   80.42,   99.94,  116.90,  149.45,  172.30,  211.30]
for k in range(0,len(n)):
   print n[k], c[k], g[k], c[k]/g[k]

# the original times are below:
# [jan@dezon MGS]$ python run_mgs_timings.py
# running the modified Gram-Schmidt method
# type precision level (d, dd, or qd) : d
# give start dimension (>= 2) : 16
# give stop dimension (<= 254) : 255
# give the step : 16
# give frequency per dimension (>= 1) : 10000
# give the name of the output file : run_mgs_d_times
# running for dimension 16 and frequency 10000 on CPU
# ['real 2.01', 'user 1.99', 'sys 0.01']
# running for dimension 16 and frequency 10000 on GPU
# ['real 4.11', 'user 1.04', 'sys 1.89']
# running for dimension 32 and frequency 10000 on CPU
# ['real 14.61', 'user 14.55', 'sys 0.05']
# running for dimension 32 and frequency 10000 on GPU
# ['real 6.52', 'user 2.34', 'sys 2.94']
# running for dimension 48 and frequency 10000 on CPU
# ['real 47.80', 'user 47.61', 'sys 0.15']
# running for dimension 48 and frequency 10000 on GPU
# ['real 11.11', 'user 4.10', 'sys 5.65']
# running for dimension 64 and frequency 10000 on CPU
# ['real 112.60', 'user 112.29', 'sys 0.21']
# running for dimension 64 and frequency 10000 on GPU
# ['real 15.38', 'user 6.48', 'sys 7.17']
# running for dimension 80 and frequency 10000 on CPU
# ['real 217.52', 'user 216.96', 'sys 0.33']
# running for dimension 80 and frequency 10000 on GPU
# ['real 22.89', 'user 10.38', 'sys 11.38']
# running for dimension 96 and frequency 10000 on CPU
# ['real 373.06', 'user 372.20', 'sys 0.50']
# running for dimension 96 and frequency 10000 on GPU
# ['real 30.43', 'user 13.31', 'sys 15.21']
# running for dimension 112 and frequency 10000 on CPU
# ['real 589.35', 'user 588.07', 'sys 0.75']
# running for dimension 112 and frequency 10000 on GPU
# ['real 40.82', 'user 18.59', 'sys 21.22']
# running for dimension 128 and frequency 10000 on CPU
# ['real 876.11', 'user 874.32', 'sys 0.99']
# running for dimension 128 and frequency 10000 on GPU
# ['real 49.10', 'user 22.46', 'sys 25.50']
# running for dimension 144 and frequency 10000 on CPU
# ['real 1243.26', 'user 1240.84', 'sys 1.22']
# running for dimension 144 and frequency 10000 on GPU
# ['real 67.41', 'user 30.40', 'sys 35.86']
# running for dimension 160 and frequency 10000 on CPU
# ['real 1701.57', 'user 1698.37', 'sys 1.58']
# running for dimension 160 and frequency 10000 on GPU
# ['real 80.42', 'user 35.74', 'sys 43.11']
# running for dimension 176 and frequency 10000 on CPU
# ['real 2260.07', 'user 2256.19', 'sys 1.86']
# running for dimension 176 and frequency 10000 on GPU
# ['real 99.94', 'user 44.00', 'sys 54.77']
# running for dimension 192 and frequency 10000 on CPU
# ['real 2932.15', 'user 2926.93', 'sys 2.50']
# running for dimension 192 and frequency 10000 on GPU
# ['real 116.90', 'user 51.93', 'sys 63.43']
# running for dimension 208 and frequency 10000 on CPU
# ['real 3722.77', 'user 3716.35', 'sys 2.94']
# running for dimension 208 and frequency 10000 on GPU
# ['real 149.45', 'user 66.72', 'sys 81.44']
# running for dimension 224 and frequency 10000 on CPU
# ['real 4641.71', 'user 4634.09', 'sys 3.33']
# running for dimension 224 and frequency 10000 on GPU
# ['real 172.30', 'user 75.59', 'sys 94.71']
# running for dimension 240 and frequency 10000 on CPU
# ['real 5703.77', 'user 5694.41', 'sys 3.96']
# running for dimension 240 and frequency 10000 on GPU
# ['real 211.30', 'user 92.85', 'sys 117.15']
