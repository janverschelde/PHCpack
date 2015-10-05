import threading
import Queue
from subprocess import PIPE, Popen
from phc_config import client_id_len, phc_cmd, phc_python, phc_maple
from phc_phc import PHC_Kill, processes
from phc_file import get_file_id, get_file_pid
import os

class solver_local(threading.Thread):
   def __init__(self, queue, Usersolver, name):
      threading.Thread.__init__(self)
      self.queue = queue
      self.Usersolver = Usersolver
      self.name = name
      # self.queue_file = queue_file

   #----------------------------------------------------------------------
   def run(self):
      while True:
         # gets the polynomial system from the queue
         [opt, Uid, msg] = self.queue.get()
         if opt == '1':
            # solves the polynomial system
            self.polysolver(Uid, msg)
            self.queue.task_done()
         elif opt == '2':
            # kill polynomial
            self.phc_kill(Uid,msg)
            self.queue.task_done()
         elif opt == '6':
            # generate print different formats
            self.phcprint(msg)
            self.queue.task_done()
         elif opt == '7':
            # generate homotopy
            self.phc_homotopy(Uid,msg)
            self.queue.task_done()


   #----------------------------------------------------------------------
   def polyfile(self, Uid, msg, poly_id):
      from phc_file import rm_exts
      file_id = "../users/%s/%s" %(Uid, poly_id)
      
      exts = ['cmd','sol','phc','err']
      rm_exts(file_id,exts)
      
      with open( file_id+'.sol', 'w') as f:
         f.write(msg)
         f.close()
      
      return file_id

   #----------------------------------------------------------------------
   def polysolver(self, Uid, msg, poly_id = "current"):
      file_id = self.polyfile(Uid, msg, poly_id)
      #cmd = phc_cmd + " %s.sol %s.phc > %s.outfile &" %(file_id, file_id, file_id)
      #print cmd
      #proc = Popen(cmd, stdout=PIPE, shell=True, preexec_fn=os.setsid)

      cmd = phc_cmd + ['%s.sol'%file_id, '%s.phc'%file_id]
      print 'cmd = %s'%cmd
      outf = open('%s.outfile'%file_id,'w')
      proc = Popen(cmd, stdout=outf)

      Pid = str(proc.pid)
      print "Pid = %s"%Pid
      with open('%s.pid'%file_id,'w') as f:
         f.write(Pid)
         f.close()
      return Pid

   #----------------------------------------------------------------------
   def polyfile_homotopy(self, Uid, msg):
      from phc_file import rm_exts, phc_start
      from os import path
      
      [poly_id, poly_id_start, poly_data] = msg.split('|')
      
      # Generate poly file
      file_id = "../users/%s/%s" %(Uid, poly_id)

      exts = ['cmd', 'sol','phc']
      rm_exts(file_id,exts)

      print "msg = %s"% msg
   
      # Generate start system file from phc file         
      file_id_start= "../users/%s/%s" %(Uid, poly_id_start)
      
      rm_exts(file_id_start, ['cmd'])
      
      if not path.isfile(file_id_start+".start"):
         phc_start(file_id_start)
      
      with open('%s.cmd'%file_id_start, 'w') as f:
         f.write("y\n%s.sol\n%s.start\n0\n0\n0\n0\n"%(file_id, file_id_start))
      f.close()
      
      return file_id, file_id_start

   #----------------------------------------------------------------------
   def phc_homotopy(self, Uid, msg):
      """
      Path Tracking by homotypy in 1 parameter
      """ 
      from subprocess import PIPE, Popen
      from os import setsid
      from phc_config import phc_hom

      file_id, file_id_start = self.polyfile_homotopy(Uid, msg)
      cmd = phc_hom + ['%s.poly'%file_id, '%s.phc'%file_id]
      print 'cmd = %s'%cmd
      inf = open('%s.cmd'%file_id_start)
      outf = open('output','w')
      proc = Popen(cmd,stdin=inf,stdout=outf)
      Pid = str(proc.pid)
      print "Pid = %s"%Pid
      with open('%s.pid'%file_id,'w') as f:
         f.write(Pid)
         f.close()

      return Pid

   def phc_kill(self, Uid, msg, poly_id = ''):
      if poly_id == '':
         poly_id = 'current'
      print 'Before kill:\n'
      print processes()
      file_id = get_file_id(Uid, poly_id)
      pid = self.get_pid(Uid, poly_id, file_id)
      if pid !='':
         PHC_Kill(pid, file_id)
         print 'After kill:\n'
         print processes()
      else:
         print "No process to kill. The process has not started."

   def get_pid(self, Uid, poly_id, file_id):
      from os import path
      pid = ''
      job_key = Uid+'/'+poly_id
      if job_key in self.Usersolver:
         pid = self.Usersolver[job_key][1]
         print 'aaa, pid = %d'%pid
      elif path.isfile('%s.pid'%file_id):
         pid = get_file_pid(file_id)
         print 'bbb, pid = %d'%pid
      return pid
      

   #----------------------------------------------------------------------
   def phcprint(self, file_id):
      # send a signal to the queue that the job is done
      cmd = phc_python + " %s.phc > %s.python"% (file_id,file_id)
      s = os.system(cmd)
      cmd = phc_maple + " %s.phc > %s.maple"% (file_id,file_id)
      s = os.system(cmd)
