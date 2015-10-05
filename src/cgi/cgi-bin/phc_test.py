def rm_without_exts(file_id, exts):
   """
   Delete all files_id.ext files
   """
   from os import remove, path, listdir
   
   l = len(file_id) - 1
   
   for i in xrange(l):
      if file_id[l-i] == '/':
         file_name = file_id[-i:]
         if i < l-1:
            path = file_id[:l-i]
         else:
            path = '.'
         break;
   
   l = len(file_name)
   for f in listdir(path):
      if f[:l] == file_name:
         if f[l+1:] not in exts:
            remove(path+'/'+f)
            
file_id = 'current'
exts = 'poly'
rm_without_exts(file_id,exts)

#def arg_test(a, b= '', c = ''):
#   print a,b,c

#a = '123'

#def arg_test2():
#   c = '4'
#   arg_test(a, c )

#arg_test2()

# ----------------------------------------------------------
#import socket
#import sys
#from time import sleep

#from phc_config import client_id_len, hostname, portnumber, sleeptime

#server_address = (hostname, portnumber)
#client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#while True:
#   try:
#      client.connect(server_address)
#      break
#   except socket.error as msg:
#      print "sleeping"
#      sleep(sleeptime)

#-------------------------------------------------
#from phc_sql import PHCSQL

#phcdb = PHCSQL()


#Uid = "0126964600478804a826e8d0f76d2db92dadfffe"
#Status = 1

#phcdb.sql_userstatus(Uid, Status)

#-------------------------------------------------

#from phc_client import PHC_Server

#Uid = "1d8a4fe98efd9fb16353a611c8b79bce3b3198cf"
#Opt = "7"

#poly_id = "current"
#poly_id_start = "bb"

#file_id =  "../users/%s/%s.poly" %(Uid, poly_id)

#Msg = "%s/%s/"%(poly_id,poly_id_start)

#with open(file_id,'r') as f:
#   Msg += f.read()
#   
#f.close()

#PHC_Server(Uid,Opt,Msg)

#-------------------------------------------------

#from collections import deque

#def read_next_symbol(data_str, pos_start, symbol):
#   pos_start += 1
#   while pos_start < len(data_str):
#      if data_str[pos_start] == symbol:
#         break;
#      pos_start += 1
#   return pos_start

#def read_element(data_str, symbol = ','):
#   element = []
#   data = data_str.split(', ')
#   for i in xrange(len(data)):
#      element.append(data[i][1:-1])
#   return element

#def read_queue(file_id):
#   queue = deque()
#   with open(file_id,'r') as f:
#      data_str = f.read()
#      pos_start = 0
#      pos_start = read_next_symbol(data_str,pos_start,'[')

#      while True:
#         if data_str[pos_start+1] == ']':
#            break;
#         pos_start = read_next_symbol(data_str,pos_start,'[')
#         pos_end = read_next_symbol(data_str,pos_start,']')
#         element = read_element(data_str[pos_start+1:pos_end])
#         queue.append(element)
#         pos_start = pos_end

#   return queue

#queue = deque()
#queue.append(["1","2","3"])
#print queue
#queue.append(['4','5','6'])
#print queue

#with open("phc_test.test",'w') as f:
#   f.write(str(queue))
#   f.close()

#print read_queue("phc_test.test")




