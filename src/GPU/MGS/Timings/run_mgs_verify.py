# This scripts allows to verify the correctness for various dimensions.

from os import system

def test_d_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 127) : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   p = './run_mgs_d '
   for i in range(stard,stopd):
      print 'running for dimension %d ...' % i
      cmd = p + str(i) + ' ' + str(i) + ' 1 2'
      for j in range(freq): system(cmd)

def test_dd_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for double double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 127) : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   p = './run_mgs_dd '
   for i in range(stard,stopd):
      print 'running for dimension %d ...' % i
      cmd = p + str(i) + ' ' + str(i) + ' 1 2'
      for j in range(freq): system(cmd)

def test_qd_cmplx():
   """
   Prompts the user for a range of dimensions, a frequency,
   and runs the modified Gram-Schmidt method on both the CPU
   and the GPU for quad double complex arithmetic.
   """
   stard = input('give start dimension (>= 2) : ')
   stopd = input('give stop dimension (<= 63) : ')
   freq = input('give frequency per dimension (>= 1) : ') 
   p = './run_mgs_qd '
   for i in range(stard,stopd):
      print 'running for dimension %d ...' % i
      cmd = p + str(i) + ' ' + str(i) + ' 1 2'
      for j in range(freq): system(cmd)

def main():
   """
   Prompts the user for a precision level.
   """
   print 'verification of modified Gram-Schmidt method'
   p = raw_input("type precision level (d, dd, or qd) : ")
   if p == "d": test_d_cmplx()
   if p == "dd": test_dd_cmplx()
   if p == "qd": test_qd_cmplx()

main()
