# L-28 MCS 507 Wed 31 Oct 2012 : clockserver.py

# Illustration of using SocketServer module.
# The server will send the time to clients.

from SocketServer import StreamRequestHandler, ThreadingTCPServer
from SocketServer import TCPServer
from collections import deque
import Queue
import threading
from os import path, system, rename, mkdir
from time import ctime, time
from shutil import move, copyfile

# PHC config file
from phc_config import Sid_len, client_id_len, portnumber, file_ext_all
from phc_config import buffer, poly_small_threads, poly_small

# local sover
from solver_local import solver_local
from phc_phc import phccheck
from phc_sql import PHCSQL

from phc_recover import recover_queue, file_recover, save_queue

from phc_file import phc_start, get_file_id, rm_all_exts, rm_exts

"""
Client Type and Id

client_id:  User Id
      client to submit new polynomial systems

Jid:  Job Id

client_id: Solver Id
     Two kind of solvers: local solver and remote solver
     


Message

Polyclient -> Server type 0
1. Submit New poly:  0 + 1 + self.client_id + New poly
2. kill computing poly: 0 + 2 + self.client_id (+ Jid)
3. Register: 0 + 3 + self.client_id

Solverclient -> Server

1. request New poly:  2 + 1 + self.client_id + Message
   when the solver is not full loaded, it can 
   request new poly, but it still first check
   whether there is any poly to kill

2. kill poly: 2 + 2 + self.client_id
   when the solver is full loaded, it only 
   check whether there is any poly to kill

3. status poly: 2 + 3 + self.client_id + update information
   update the status of polynomial system
   if polynomial is solved, send the solutions
   back to the server

Server -> Solverclient
0. Nothing: 0
1. New poly:  1 + self.client_id + Message
2. kill poly: 2 + self.client_id

"""

class PHC_Server(StreamRequestHandler):
   """
   PHC_Server handles jobs form polyclient and solverclient.
   """
   def handle(self):
      """
      Handler handles jobs form polyclient and solverclient.
      """
      start_time = time()
      #return empty message by default
      self.msg = ''

      # identify client_type, client_opt, client_id, client_msg 
      self.read_client_type()
 
      # handle different client 
      if self.client_type == "0":   # Client type 0 -- polyclient to submit poly or kill process
         self.polyclient()
      elif self.client_type == "2": # Client type 2 -- large poly remote solver
         self.solverclient()
      
      # return message
      self.wfile.write(self.msg)
      
      # write log
      self.phc_server_log()
      
      end_time = time()
      time_cost = end_time - start_time
      print time_cost
   
   def phc_server_log(self):
      self.print_queue()
   
   def print_queue(self):
      print "Polynomials in the large queue: "
      print queue
      print "Polynomials in the small queue: "
      print queue_local.queue
      print "Polynomials to kill: "
      print queue_kill
      print list_kill
      print "return message:"
      print self.msg
      print "******************************"
      
   def polyclient(self):
   
      self.read_poly_client_msg()
      
      # polyclient opt 1 -- submit new poly
      if self.client_opt =='1':
         self.polyclient_newpoly()
         
      # polyclient opt 7 -- submit new poly by homotopy
      elif self.client_opt =='7':
         self.polyclient_newhom()

      # polyclient opt 2 -- kill solving poly
      elif self.client_opt == "2":
         self.polyclient_kill()

      # polyclient opt 3 -- create folder for the user
      elif self.client_opt == "3":
         self.polyclient_mkdir()

      # polyclient opt 4 -- rename files
      elif self.client_opt == "4":
         self.polyclient_rename()

      # polyclient opt 5 -- delete files
      elif self.client_opt == "5":
         self.polyclient_delete()
         
      # polyclient opt 6 -- print maple python files from phc files
      elif self.client_opt == '6':
         self.polyclient_print()
 
      # polyclient opt others: no such option
      else:
         print "Message Error: no such option " + self.client_opt

   def polyclient_newhom(self):
      # check whether polynomial system is efficient
      rm_all_exts(self.file_id)
      
      newpoly = [self.client_opt,self.client_id,self.client_msg]
      self.client_msg = self.client_msg.split('|')

      if len(self.client_msg) == 1:        
         file_id_start= get_file_id(self.client_id, self.client_msg[0]) 
         if not path.isfile(file_id_start+".start"):
            phc_start(file_id_start)
         self.msg = 'homotopy start system generated'

      else:
         with open('%s.poly'% self.file_id, 'w') as f:
            f.write(self.client_msg[2])
            f.close()


         # distribute job to small or large queue
         worker_ind, deg = self.job_distributor()
         
         if worker_ind == 0:
            phcdb.delpoly(self.client_id, self.poly_id)
            with open('%s.err'% self.file_id, 'r') as f:
               self.msg = f.read() + "received"
         else:
            print "**********%s"%deg
            rm_exts(self.file_id, 'err')
            
            with open('%s.sta'% self.file_id,'w') as f:
               f.write('submitted')
               f.close()

         queue_local.put(newpoly)

         self.msg = 'homotopy received'

         # update userstatus to computing
         phcdb.delpoly(self.client_id, self.poly_id)
         phcdb.newpoly(self.client_id, self.poly_id, deg, path.getctime('%s.poly'% self.file_id))

   def polyclient_newpoly(self):
      """
      polyclient to submit new poly
      """
      if self.client_id in Usersolver:
         self.msg = "one polynomial is being solved"
      else:
         # check whether polynomial system is efficient
         rm_all_exts(self.file_id)
         phcdb.delpoly(self.client_id, self.poly_id)
         
         with open('%s.poly'% self.file_id, 'w') as f:
            f.write(self.client_msg)
            f.close()

         start_time = time()
         # distribute job to small or large queue
         worker_ind, deg = self.job_distributor()
         end_time = time()
         time_cost = end_time - start_time
         print "check_time = %f"%time_cost
         poly_status = '1'
         if worker_ind == 0:
            with open('%s.err'% self.file_id, 'r') as f:
               self.msg = f.read() + "received"
            poly_status = '9'
         else:
            rm_exts(self.file_id, 'err')
            with open('%s.sta'% self.file_id,'w') as f:
               f.write('submitted')
               f.close()
            
            newpoly = [self.client_opt,self.client_id,self.client_msg]
            if worker_ind == 1:
               queue_local.put(newpoly)               
               self.msg = 'small polynomial received'
            else:
               queue.append(newpoly)
               self.msg = 'large polynomial received'
            
         # update userstatus to computing
         phcdb.newpoly(self.client_id, self.poly_id, deg, path.getctime('%s.poly'% self.file_id), poly_status)
   
   def job_distributor(self):
      # check polynomial
      deg = phccheck(self.file_id)
      
      print 'degree = %d' %deg

      if deg == 0:
         worker_ind = 0
      elif deg < poly_small:
         worker_ind = 1
      else:
         worker_ind = 2

      return worker_ind, deg

   def polyclient_kill(self):
      """
      polyclient to kill solving poly
      """
      # request has been sent, ask solverclient to stop
      if self.client_id in Usersolver:
         self.solver_id = Usersolver[self.client_id]
         if Usersolver[self.client_id]+'/':
            queue_local.put(['2',self.client_id,''])
         else:
            queue_kill[self.solver_id].append(self.client_id)

      # request hasn't been sent, block it before sending
      else:
         print 'abccccc'
         queue_local.put(['2',self.client_id,''])
#         list_kill.append(self.client_id)
         #del Usersolver[self.client_id]
      
      phcdb.killpoly(self.client_id, self.poly_id)

      # update message to return
      self.msg = self.client_id + ' Kill request received'
      
   def polyclient_mkdir(self): 
      """
      polyclient to create folder for the user
      """
      mkdir("../users/"+ self.client_id,0755)
      mkdir("../users/"+ self.client_id + "/.Deleted",0755)
      
      sample_files = ['sample.poly', 'sample.phc', 'sample.sta', 'sample.sol']
      
      for sample_file in sample_files:
         copyfile("../Sample/%s"%(sample_file), "../users/%s/%s"%(self.client_id,sample_file))
      
      poly_create_time = path.getctime("../users/%s/sample.poly"%(self.client_id))
      phcdb.newpoly(self.client_id, 'sample', 6, poly_create_time)
      
      # update message to return
      self.msg = self.client_id + ' Mkdir request received'
      
   def polyclient_rename(self):
      """
      polyclient to rename file for the user
      """
      name = self.client_msg.split('|')
      file_path = "../users/" + self.client_id + "/"
      file_src = file_path + name[0]
      file_dst = file_path + name[1]
      
      for ext in file_ext_all:
         if path.isfile(file_src + ext):
            rename(file_src + ext, file_dst + ext)
      
      phcdb.delpoly(self.client_id, name[1])
      phcdb.renamepoly(self.client_id, name[0],name[1])
      
      # update message to return
      self.msg = self.client_id + ' Rename request received'
      
   def polyclient_delete(self):
      """
      polyclient to delete files
      
      Uid       = self.client_id
      file name = self.client_msg
      """
      filesrc = "../users/" + self.client_id + "/"
      filedst = filesrc + "/.Deleted/"
      # add deleting time to avoid repeated file name
      deltime = "("+ctime()[4:]+")"
      filesrc = filesrc + self.client_msg
      filedst = filedst + self.client_msg + deltime
      
      for ext in file_ext_all:
         if path.isfile(filesrc + ext):
            move(filesrc + ext, filedst + ext)
      
      # update sql datebase
      phcdb.delpoly(self.client_id, self.client_msg)

      # update message to return
      self.msg = self.client_id + ' Delete request received'
      
   def polyclient_print(self):
      """
      polyclient to maple python files from phc files
      
      put job to the queue of local solver
      """
      newprint = [self.client_opt,self.client_id,self.client_msg]
      queue_local.put(newprint)

      # update message to return
      self.msg = self.client_id + ' View request received'
         
   def solverclient(self):
      """
      solverclient check new jobs on the server after certain amount of time
      """   
      self.read_solver_client_msg()
      
      # Generate a kill queue for each solver, when the solver connects the first time
      if self.client_id not in queue_kill:
         queue_kill[self.client_id]=deque()
         print "Initialize Solver " + self.client_id
      
      # solverclient opt 1 -- solverclient to get more polys to solve
      if self.client_opt == '1':
         self.solverclient_solve()

      # solverclient opt 2 -- solverclient to check solving polys to kill
      elif self.client_opt == '2':
         self.solverclient_kill()
   
      # solverclient opt 3 -- solverclient to return results of phc
      elif self.client_opt == '3':
         self.solverclient_results()
         
   def solverclient_solve(self):
      # no large polynomial to solve.
      if not queue: # solver do nothing
         self.msg = '0'
         print 'writing \"' + self.msg + '\" to client'
      else:
         # remove killed process in the queue
         while queue:
            newpoly = 1
            [opt, Uid, poly] = queue.popleft()
            if Uid in list_kill:
               newpoly = 0
               list_kill.remove(Uid)
               print "request %s Killed" % Uid
            else:
               break
         if newpoly:
               print "new poly request from " + self.client_id
               self.msg = opt + Uid + poly
               Usersolver[Uid] = self.client_id
         else:
            self.msg = '0'
         print Usersolver
            
   def solverclient_kill(self):
      if not queue_kill[self.client_id]: # process to kill
         self.msg = '0'
         print 'writing \"' + self.msg + '\" to client'
      else:
         self.client_id = queue_kill[self.client_id].popleft()
         self.msg = '2'+ self.client_id
   
   def solverclient_results(self):
      u_id = self.client_msg[0:client_id_len]
      print "Receiving the solution for %s"%self.client_id
      del Usersolver[u_id]
      queue_local.put(['6', u_id, self.poly_id])
      
      # update sql datebase -- solve imformation
      phcdb.slvpoly(u_id, self.poly_id)
      
      self.msg = 'r'
      
   def read_client_type(self):
      # Use 1st digit to identify client type
      self.client_type = self.rfile.read(1)
   
   def read_poly_client_msg(self):
      # Use 2st digit to identify client self.client_opt
      self.client_opt = self.rfile.read(1)
      # identify client id
      self.client_id = self.rfile.read(client_id_len)  # buffer size to be changed     
      # Save the message
      messg = self.rfile.read(buffer-client_id_len-2).strip()
      if '||' in messg:
         [self.poly_id, self.client_msg] = messg.split('||')
      else:
         self.poly_id = 'current'
         self.client_msg = messg
      self.job_key = self.client_id + '/' + self.poly_id
      
      self.file_id = get_file_id(self.client_id, self.poly_id)
      
      # print client information and message
      self.print_poly_client_msg()
      
   def print_poly_client_msg(self):
      print "client_addr = ", self.client_address
      print "client_type = %s" % self.client_type
      print "client_opt  = %s" % self.client_opt
      print "client_id   = %s" % self.client_id
      print "  poly_id   = %s" % self.poly_id
      print "client_msg  = %s" % self.client_msg
      print "   job_key  = %s" % self.job_key
      print "  file_id   = %s" % self.job_key
   
   def read_solver_client_msg(self):
      # Use 2st digit to identify client self.client_opt
      self.client_opt = self.rfile.read(1)
      # identify client id
      self.client_id = self.rfile.read(client_id_len)  # buffer size to be changed     
      # Save the message
      messg = self.rfile.read(buffer-client_id_len-2).strip()
      if '||' in messg:
         [self.poly_id, self.client_msg] = messg.split('||')
      else:
         self.poly_id = 'current'
         self.client_msg = messg
      self.job_key = self.client_id + '/' + self.poly_id
      
      # print client information and message
      self.print_solver_client_msg()
         
   def print_solver_client_msg(self):
      print "client_addr = ", self.client_address
      print "client_type = %s" % self.client_type
      print "client_opt  = %s" % self.client_opt
      print "client_id   = %s" % self.client_id
      print "  poly_id   = %s" % self.poly_id
      print "client_msg  = %s" % self.client_msg
      print "   job_key  = %s" % self.job_key
      

def small_solvers(queue_local, Usersolver, poly_small_threads):
   # create a thread pool and give them a queue
   smsolver = []
   for i in range(poly_small_threads):
      smsolver.append(solver_local(queue_local, Usersolver, "sslv%d"%i))
      smsolver[i].setDaemon(False)
      smsolver[i].start()
      
   return smsolver

def main():

   smsolver = small_solvers(queue_local, Usersolver, poly_small_threads)

   # Reuse the port to debug
   TCPServer.allow_reuse_address = True
   ss = TCPServer(('0.0.0.0',portnumber),PHC_Server)

   print 'server is listening to', portnumber

   try:
      print 'press ctrl c to stop server'
      ss.serve_forever()
   except KeyboardInterrupt:
      print ' ctrl c pressed, closing server'

      if queue:
         print 'writing recover file'
         save_queue(queue)
      
      print 'closing solver threads'
      for i in range(poly_small_threads):
         smsolver[i]._Thread__stop()

      print 'server closed'
      ss.socket.close()

if __name__ == "__main__":
    queue_local = Queue.Queue()         # small polynomial system job queue
    queue = recover_queue(file_recover) # large polynomial system job queue
    print queue
    queue_kill = {}                     # large polynomial system job queue

    list_kill = ['Kill']
    Usersolver={}
    phcdb =PHCSQL()
    main()
