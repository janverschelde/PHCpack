def phccheck(polyname):
   from subprocess import Popen
   from phc_config import phc, phc_check
   from os import system, path, waitpid
   
   filename = "%s.poly"% polyname
   result = "%s.err"% polyname
   if path.isfile(result):
      r = 'rm -f %s'% result  # XXX remove outfile?
      s = system(r)

   p = Popen([phc, "-g", filename, result])
   s = waitpid(p.pid, 0)

   f = open(result, 'r')

   success_mark = 'PHC parsed input successfully on'
   success_mark_l = len(success_mark)
   
   poly_err = 1
   for line in f:
      if line[:success_mark_l] == success_mark:
         poly_err = 0
         break;

   if not poly_err:
      f.seek(0)
      degree = f.readline()
      f.close()

      r = 'rm -f %s'% result  # XXX remove outfile?
      s = system(r)
      
      for i in xrange(len(degree)):
         if degree[i] < '0' or degree[i] > '9':
            if i != len(degree)-1:
               return 0
            else:
               return int(degree[:i])
   else:
      f.close()
      return 0
   

def PHC_Solve(Uid,data):
   from subprocess import PIPE, Popen
   from os import setsid
   from phc_config import phc_cmd
   
   with open(Uid+'.poly', 'w') as f:
      f.write(data)
   f.close()
   with open(Uid+'.sol', 'w') as f:
      f.write(data)
   f.close()
   cmd = phc_cmd + " %s.sol %s.phc > %s &" %(Uid,Uid,Uid)
   proc = Popen(cmd, stdout=PIPE, shell=True, preexec_fn=setsid)
   return proc.pid

def PHC_Kill(Pid, file_id):
    import signal, os
    if Pid != '':
       #if os.path.getsize(file_id + ".poly") == os.path.getsize(file_id + ".sol"):
        print "Kill Pid = %s"%Pid
        try:
           os.kill(int(Pid), signal.SIGKILL)
        except OSError as ex:
           print "Process is finished"
        r = 'rm -f %s.sol %s.phc %s.pid %s' %(file_id, file_id, file_id, file_id)
        s = os.system(r)
          
## check whether a polynomial has been solved
#def PHC_Status(name):
#   from os import path, system
#   # solved, solution is writing
#   if path.isfile(name+".sol") and path.getsize(name+".poly") != path.getsize(name+".sol"): 
#      # XXX solution file writting might not be complete
#      if path.isfile(name+".pid"): # remove pid file
#         r = 'rm -f %s.pid'% name
#         s = system(r)
#      return 1
#   # not solved
#   else:
#      return 0

def processes():
   import subprocess
   import shlex
   proc1 = subprocess.Popen(['ps', 'aux'],stdout=subprocess.PIPE)
   proc2 = subprocess.Popen(['grep', '[p]hc64'],stdin=proc1.stdout,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
   return proc2.communicate()[0].split('\n')
