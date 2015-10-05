from socket import *
from phc_config import server_address, buffer,  client_id_len, phc_debug


def PHC_Server(Uid,Opt,Msg):
   client = try_connect(server_address)
   connect_status = 1
   if client:
      data = '0' + Opt + Uid + Msg
      data = exchange_msg(client, data)
      client.close()
   else:
      print '<font color="red"> Sorry. The server is down. Please try again later. </font>'
      connect_status = 0
   
   return connect_status
    
def try_connect(server_address):  
   from socket import error
   client = socket(AF_INET, SOCK_STREAM)
   try:
      client.connect(server_address)
      return client;
   except error as msg:
      return 0

def retry_connect(server_address): 
   from socket import error 
   client = socket(AF_INET, SOCK_STREAM)
   while True:
      try:
         client.connect(server_address)
         print "client is connected"
         return client;
      except error as msg:
         print "connect refused, sleeping, try later"
         sleep(sleeptime)


def exchange_msg(client, data):
   if phc_debug: print "<p>send msg = %s</p>"% data
   client.send(data + (buffer-len(data))*' ')
   
   data = client.recv(buffer).strip()
   if phc_debug: print "<p>receive msg = %s</p>"% data
   
   return data
   
def retry_exchange_msg(data):
   client = retry_connect(client, server_address)
   data = send_rec(client, data)
   client.close()

   return data
