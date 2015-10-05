#############################   Phcpy   ###############################
### To test Phcpy
# polysys = Print_Python_Eq(name+".poly")
# print polysys
# sols = Phcpy_Solve(polysys)
# Print_Python_Sols(sols)

# This part is not used.
def Phcpy_Solve(polysys):
   import sys
   sys.path.append("../../PHCv2_3p/Python")
   import phcpy
   sols = phcpy.solve(polysys)
   return sols

def Print_Python_Eq(f):
   file = open(f,'r')
   line = file.readline()
   polysys = []
   symbols = ['\r','\n']
   neweq = "";
   for line in file:
      for i in range(0,len(line)):
         if line[i] not in symbols:
            neweq += line[i]
         if line[i] == ';':
            polysys.append(neweq+'')
            neweq = ""
   return polysys

def Print_Python_Sols(sols):
   for line in sols:
      print line + "</p>"
