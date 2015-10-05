from collections import deque
file_recover = 'phcserver.recover'

def read_next_symbol(data_str, pos_start, symbol):
   pos_start += 1
   while pos_start < len(data_str):
      if data_str[pos_start] == symbol:
         break;
      pos_start += 1
   return pos_start

def read_element(data_str, symbol = ', '):
   element = []
   data = data_str.split(symbol)
   for i in xrange(len(data)):
      element.append(data[i])
   return element

def read_queue(file_id):
   queue = deque()
   with open(file_id,'r') as f:
      data_str = f.read()
      pos_start = 0
      pos_start = read_next_symbol(data_str,pos_start,'[')
      print pos_start

      while True:
         if data_str[pos_start+1] == ']':
            break;
         pos_start = read_next_symbol(data_str,pos_start,'[')
         pos_end = read_next_symbol(data_str,pos_start,']')
         element = read_element(data_str[pos_start+1:pos_end])
         queue.append(element)
         pos_start = pos_end

   return queue

#read_queue(file_recover)

def recover_queue(file_id):
   """
   Recover a queue from a file
   """
   from os import path, remove
   if path.isfile(file_id):
      queue = read_queue(file_id)
      remove(file_id)
      return queue
   else:
      queue = deque()
      return queue

def element2str(element, symbol = ', '):
   """
   Change a list to a string, elements are splitted by symbol
   """
   element_str = '['
   for i in xrange(len(element)-1):
      element_str += element[i]
      element_str += symbol
   element_str +=  element[-1] + ']'
   return element_str

def save_queue(queue, file_id):
   """
   Save a queue into a file
   """
   with open(file_id, 'w') as f:
      f.write('queue = [')
      while queue:
         element = queue.popleft()
         f.write(element2str(element))
      f.write(']')
      f.close()
