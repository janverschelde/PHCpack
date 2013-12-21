# This scripts allows to take timings for various dimensions.

from commands import getoutput

def get_times(s):
   """
   Processes the timings in the string s.
   """
   R = []
   L = s.split('\n')
   for e in L:
      # K = e.split('\t')
      R.append(e)
   return R

def time_d_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 256) : ')
   step = input('give the step : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   name = raw_input('give the name of the output file : ')
   file = open(name,'w')
   p = 'time -p ../run_mgs_d '
   for i in range(stard,stopd+1,step):
      s =  'running for dimension %d and frequency %d on CPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 1'
      K = getoutput(cmd)
      t = get_times(K)
      print t
      file.write(str(t) + '\n')
      s = 'running for dimension %d and frequency %d on GPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 0'
      L = getoutput(cmd)
      t = get_times(L)
      print t
      file.write(str(t) + '\n')
   file.close()

def time_dd_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for double double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 127) : ')
   step = input('give the step : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   name = raw_input('give the name of the output file : ')
   file = open(name,'w')
   p = 'time -p ../run_mgs_dd '
   for i in range(stard,stopd+1,step):
      s =  'running for dimension %d and frequency %d on CPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 1'
      K = getoutput(cmd)
      t = get_times(K)
      print t
      file.write(str(t) + '\n')
      s = 'running for dimension %d and frequency %d on GPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 0'
      L = getoutput(cmd)
      t = get_times(L)
      print t
      file.write(str(t) + '\n')
   file.close()

def time_qd_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for quad double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 63) : ')
   step = input('give the step : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   name = raw_input('give the name of the output file : ')
   file = open(name,'w')
   p = 'time -p ../run_mgs_qd '
   for i in range(stard,stopd+1,step):
      s =  'running for dimension %d and frequency %d on CPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 1'
      K = getoutput(cmd)
      t = get_times(K)
      print t
      file.write(str(t) + '\n')
      s = 'running for dimension %d and frequency %d on GPU' % (i,freq)
      print s
      file.write(s + '\n')
      cmd = p + str(i) + ' ' + str(i) + ' ' + str(freq) + ' 0'
      L = getoutput(cmd)
      t = get_times(L)
      print t
      file.write(str(t) + '\n')
   file.close()

def main():
   """
   Prompts the user for a precision level.
   """
   print 'running the modified Gram-Schmidt method'
   p = raw_input("type precision level (d, dd, or qd) : ")
   if p == "d": time_d_cmplx()
   if p == "dd": time_dd_cmplx()
   if p == "qd": time_qd_cmplx()

main()
