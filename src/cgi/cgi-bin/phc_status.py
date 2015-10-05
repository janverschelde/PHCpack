def PHC_File_Status(file_id):
   from os import path, system, remove
   status = 0
   if path.isfile(file_id+'.sta'):
      f = open(file_id+'.sta','r')
      if f.readline() == 'submitted':
         status = 1
         if path.isfile(file_id+".sol") and  path.isfile(file_id+".poly"):
            sol_file_size  = path.getsize(file_id+".sol")
            poly_file_size = path.getsize(file_id+".poly")
            if sol_file_size > 3 and sol_file_size != poly_file_size: 
                status = 2
      f.close()
   return status
            
def PHC_Status(Uid, poly_id, file_id, phcdb):
   status = phcdb.check_poly_status(Uid, poly_id)
   if status == 1:
      return PHC_File_Status(file_id)
   else:
      return status
