from phc_client import PHC_Server

def PHC_Kill(Uid,form,phcdb):
   from os import path
   phcdb.userstatus(Uid, "NULL")
   if path.isfile("current.sol"):
      if path.getsize("current.poly") != path.getsize("current.sol"):
         data = PHC_Server(Uid,'2','')
   else:
      data = PHC_Server(Uid,'2','')

   if data != 0:
      Print_Poly_File_New(Uid, "current")
