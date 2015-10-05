def check_file_exist(file_name):
   from os import path
   if path.isfile(file_name):
      return 1
   else:
      from phc_config import phc_debug
      if phc_debug:
         print "<p><font color='red'>The file %s doesn't exist. Try it again later.</font></p>"%file_name
      return 0

def get_file_id(Uid, poly_id):
   return "../users/%s/%s"%(Uid, poly_id)

def get_file_pid(file_id):
   pid = 0
   with open('%s.pid'%file_id,'r') as f:
      pid = int(f.read())
   return pid
      

def scp_upload(filename, destination):
   from os import waitpid
   from subprocess import Popen
#   myfile = "phc_html.py"
#   destination = "yxc@131.193.178.177:Desktop/"
   p = Popen(["scp",  filename, destination])
   sts = waitpid(p.pid, 0)
   print "scp uploading done"

def rm_exts(file_id, exts):
   """
   Delete all files_id.ext files
   """
   from os import remove, path
   
   for ext in exts:
      file_name = "%s.%s"%(file_id, ext)
      if path.isfile(file_name):
         remove(file_name)
 
def rm_all_exts(file_id):
   from os import system
   """
   Delete all files_id.*
   """
   r = 'rm -f %s.*'%file_id
   s = system(r)

def rm_without_exts(file_id, exts):
   """
   Delete all files_id.* except .ext files
   """
   from os import remove, path, listdir
   
   l = len(file_id) - 1
   
   for i in xrange(l+1):
      if i == l:
         path ='.'
         file_name = file_id
         break;
      elif file_id[l-i] == '/':
         file_name = file_id[-i:]
         path = file_id[:l-i]
         break;
   
   l = len(file_name)
   for f in listdir(path):
      if f[:l] == file_name:
         if f[l+1:] not in exts:
            remove(path+'/'+f)

def poly_file_id(Uid, poly_id):
   return "../users/%s/%s"%(Uid, poly_id)


# This part need to read every file, needs to be improve
def Deg(poly_id):
   file = open(poly_id+".poly",'r')
   line = file.next()
   file.close()
   try:
      deg = int(line)
   except ValueError:
      deg = -1
   return deg

def check_int(s):
   text_sym = ['\n','\r', '\t']
   l = len(s)-1
   for i in xrange(l):
      if s[i] > '9' or s[i] < '0':
         return 0

   if s[l] > '9' or s[l] < '0':
      if s[l] not in text_sym:
         return 0

   return 1

def Deg_Sols(poly_id):
   file = open(poly_id+".sol",'r')
   line = file.readline().split(' ')
   if len(line) == 2 and check_int(line[0]) and check_int(line[1]):
      sols = int(line[0])
      deg = int(line[1])
   else:
      for line in file:
         if line == "THE SOLUTIONS :\n":
            line = file.next()
            break

      for i in range(0, len(line)):
         if line[i] == ' ': 
            sols = int(line[:i])
            deg = int(line[i+1:])
            break
   file.close()
   return sols, deg

def Sols(Uid, poly_id):
   file_id = poly_file_id(Uid,poly_id)
   file = open("%.phc"%file_id,'r')
   for line in file:
      if line[0] == 'T':
         line = file.next()
         break
   for i in range(0, len(line)):
      if line[i] == ' ': 
         sols = int(line[:i])
         deg = int(line[i+1:])
         break
   file.close()
   return sols

#This part need to read every file, needs to be improve
def Slvtime(Uid, poly_id):
   file_id = poly_file_id(Uid,poly_id)
   file = open("%.phc"%file_id,'r')
   file.seek(-1024, 2)
   last = file.readlines()[-1].decode()
   file.close()
   slvtime = last[25:]
   for i in range(0, len(slvtime)):
      if slvtime[i] == 's':
         slvtime = slvtime[:i+1]
         break
   return slvtime

#This part need to read every file, needs to be improve
def Solve_Time(poly_id):
   file = open(poly_id+".phc",'r')
   file.seek(-1024, 2)
   last = file.readlines()[-1].decode()
   file.close()
   solvetime = last[25:]
   for i in range(0, len(solvetime)):
      if solvetime[i] == 's':
         solvetime = solvetime[:i+1]
         break
   return solvetime

def solve_time(poly_id):
   """
   read time information for phc file
   """
   file = open(poly_id+".phc",'r')
   n_steps = 3
   solvetime = [0.0 for i in xrange(n_steps + 1)]

   time_line = ["TIMING INFORMATION for Root Counting\n",
                "TIMING INFORMATION for Construction of Start System\n",
                "TIMING INFORMATION for continuation\n",
                "TIMING INFORMATION for Solving the polynomial system\n"]
   time_mark = "The elapsed time in seconds was";


   time_record = -1
   for line in file:
      if time_record != -1:
         line = line[len(time_mark):].strip()
         solvetime[time_record] = float(line.split(" =")[0])
         time_record = -1
      elif line in time_line:
         time_record = time_line.index(line)

   if solvetime[n_steps] == 0.0:
      for i in xrange(n_steps):
         solvetime[n_steps] += solvetime[i]

   return solvetime[n_steps]


def phc_start(startfile):
   """
   Get the start system from the solution file
   """
   from os import remove, path
   
   
   if path.isfile('%s.start'% startfile):
      remove('%s.start'% startfile)

   Start_System_Line = ["START SYSTEM : \n", " START SYSTEM :\n"]
   
   Start_System_Line_Len = len(Start_System_Line[0])
   
   f = open('%s.phc'% startfile,'r')
   
   g = open('%s.start'% startfile,'w')
   
   print startfile

   f.readline()
   
   breaknextloop = 0
   for line in f:
      if breaknextloop == 1:
         for x in line:
            if x != " " and x != "\n":
               breaknextloop = 2
               break
         if breaknextloop == 2:
            break
      if len(line)>Start_System_Line_Len and line[-Start_System_Line_Len:] in Start_System_Line:
         breaknextloop = 1
   
   if breaknextloop != 2:
      g.close()
      remove('%s.start'% startfile)
      return

   g.write(line)
   
   Solution_Start_Line = ["THE START SOLUTIONS : \n", "START SOLUTIONS : \n"]
   Solution_End_Line = ["HOMOTOPY PARAMETERS :\n", "TIMING INFORMATION for continuation\n"]
   
   for line in f:
      if line in Solution_Start_Line :
         g.write("THE SOLUTIONS : \n")
      elif line in Solution_End_Line:
         break
      else:
         g.write(line)
    
   f.close()
   g.close()
