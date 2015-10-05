# L-28 MCS 507 Wed 31 Oct 2012 : clockclient.py

# Code for client corresponding to clockserver.py.
# It prints the data received from the server.

from socket import *
import os,signal
from subprocess import PIPE, Popen
from time import sleep

from phc_config import sleeptime, Solver_Remote_Folder, usersfolder
from phc_phc import PHC_Solve, PHC_Kill, PHC_Status, processes
from phc_file import scp_upload
from phc_client import retry_exchange_msg

cores = 2
Userproc = {}

computer_id = 'T430s'
solver_id = computer_id + (client_id_len-len(computer_id))*' '

def solver(full,Uid):
   if Uid:
      print "upload result"
      solver_sleep = solver_result(Uid)
   elif not full:
      print "request any running job to kill"
      solver_sleep = solver_kill()
      if solver_sleep == 1:
         print "request any new job"
         solver_sleep = solver_solve()
   else:
      print "request any running job to kill"
      solver_sleep = solver_kill()

   return solver_sleep # 0: Empty; 1: kill the process
   
def solver_solve():
   data = '21%s' % solver_id
   
   data = retry_exchange_msg(data)
   
   solver_sleep = 1
   Uid = data[1:1+client_id_len]
   print "Uid = " + Uid
   if data[0] =='1':
      # new poly system to solve
      if Uid!= 'Empty':
         pid = PHC_Solve(Uid,data[1+client_id_len:])
         Userproc[Uid]= pid
         with open(Uid+'.pid', 'w') as f:
            f.write(data)
         f.close()
      solver_sleep = 0
   return solver_sleep

def solver_kill():
   # job to kill
   data = '22%s' % solver_id

   data = retry_exchange_msg(data)
   
   solver_sleep = 1
   if data[0] =='2':
      Uid = data[1:1+client_id_len]
      print "Uid = " + Uid
      print Uid + " to kill"
      if Uid in Userproc:
         Pid = Userproc[Uid]
         PHC_Kill(Pid,Uid)
         del Userproc[Uid]
      else:
         print "no such a process"
      solver_sleep = 0
   return solver_sleep
   
def scp_sol(Uid):
   print "sending sol and phc files"
   scp_upload(Uid+'.sol',usersfolder+Uid+'/current.sol')
   scp_upload(Uid+'.phc',usersfolder+Uid+'/current.phc')
   print "sending files done"

def solver_result(Uid):
"""
   copy solution by scp and report to server
"""
   # copy solution files
   scp_sol(Uid)

   # report the poly has been solved
   data = "23%s%s"%(solver_id,Uid)
   
   data = retry_exchange_msg(data)
   
   solver_sleep = 0
   return solver_sleep
   
def open_solver_folder():   
   Open_Folder = os.path.isdir(Solver_Remote_Folder)
   if not Open_Folder:
      print "Create Remote Solver Temprary Solution Folder"
      os.mkdir(Solver_Remote_Folder)
   os.chdir(Solver_Remote_Folder)
   print "Get to Remote Solver Temprary Solution Folder"

def main():
   # change directory

   open_solver_folder()
   
   while True:
      full = 0
      Procs = processes()
      print "Procs"
      comps = len(Procs)-1
      for i in xrange(comps):
         print Procs[i]
      if comps >= cores:
         full = 1
         Uid = 0
         print "No more Cores to use"
      else: # XXX to check the status, check by proccess Uid but not 
         if comps < len(Userproc): # some poly has been solved
            for Jid in Userproc:
               if PHC_Status(Jid):
                  del Userproc[Jid]
                  Uid = Jid
                  print "$$$$ Job done: " + Jid 
                  break
         else:
            Uid = 0
      solver_sleep = solver(full, Uid)
      if solver_sleep:
         # No polynomial to solve and no request
         print "No request"
         sleep(sleeptime)
      print "###########################################"
 
if __name__ == "__main__":
   main()
