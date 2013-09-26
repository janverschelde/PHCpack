n = [   16,     32,     48,     64,      80,      96,     112]
c = [17.17, 125.06, 408.20, 952.35, 1841.07, 3159.93, 4994.48]
g = [11.85,  22.44,  35.88,  55.18,   79.11,  105.67,  143.60]
for k in range(0,len(n)):
   print n[k], c[k], g[k], c[k]/g[k]

#running for dimension 16 and frequency 10000 on CPU
#['real 17.17', 'user 17.13', 'sys 0.02']
#running for dimension 16 and frequency 10000 on GPU
#['real 11.85', 'user 4.84', 'sys 5.94']
#running for dimension 32 and frequency 10000 on CPU
#['real 125.06', 'user 124.85', 'sys 0.10']
#running for dimension 32 and frequency 10000 on GPU
#['real 22.44', 'user 8.29', 'sys 12.57']
#running for dimension 48 and frequency 10000 on CPU
#['real 408.20', 'user 407.61', 'sys 0.26']
#running for dimension 48 and frequency 10000 on GPU
#['real 35.88', 'user 14.74', 'sys 19.88']
#running for dimension 64 and frequency 10000 on CPU
#['real 952.35', 'user 951.06', 'sys 0.48']
#running for dimension 64 and frequency 10000 on GPU
#['real 55.18', 'user 23.03', 'sys 30.26']
#running for dimension 80 and frequency 10000 on CPU
#['real 1841.07', 'user 1838.74', 'sys 0.81']
#running for dimension 80 and frequency 10000 on GPU
#['real 79.11', 'user 33.64', 'sys 44.20']
#running for dimension 96 and frequency 10000 on CPU
#['real 3159.93', 'user 3156.02', 'sys 1.27']
#running for dimension 96 and frequency 10000 on GPU
#['real 105.67', 'user 45.65', 'sys 58.37']
#running for dimension 112 and frequency 10000 on CPU
#['real 4994.48', 'user 4988.46', 'sys 1.83']
#running for dimension 112 and frequency 10000 on GPU
#['real 143.60', 'user 60.37', 'sys 81.75']
