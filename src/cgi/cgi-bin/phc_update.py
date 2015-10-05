from phc_html import *
from phc_file import get_file_id, check_file_exist
from phc_status import PHC_Status

def PHC_Update(Uid,Pid,poly_id,phcdb, sol_page=0):
   from os import path
   if poly_id == '':
      Print_Poly_Empty(Uid)
   else:
      file_id = get_file_id(Uid, poly_id)
      if check_file_exist(file_id+".poly"):
         Print_Notice_Top()
         status = PHC_Status(Uid, poly_id, file_id, phcdb)
         Print_Area_Eq(file_id+".poly")
         if path.isfile(file_id+".err"):
            phcdb.userstatus(Uid, "NULL")
            Print_Button_Status(status, Uid, Pid, poly_id)
            Print_Notice_Mid(status, Uid, poly_id)

         elif status == 1: # being computed
            Print_Button_Status(status, Uid, Pid, poly_id)
            Print_Notice_Mid(status, Uid, poly_id)
            
         elif status == 0: # killed
            Print_Button_Status(status, Uid, Pid, poly_id)
            Print_Notice_Mid(status, Uid, poly_id)
            Print_Report(status, file_id, Uid, poly_id)
            
         else: # Just Solved
            phcdb.userstatus(Uid, "NULL")
            phcdb.slvpoly2(Uid, poly_id)
            Print_Button_Status(status, Uid, Pid, poly_id)
            Print_Notice_Mid(status, Uid, poly_id)
            Print_Report(status, file_id, Uid, poly_id, sol_page)

      else:
         Print_Poly_Empty(Uid)
         if poly_id == "current":
            phcdb.userstatus(Uid, "NULL")

def Print_Report(status, file_id, Uid, poly_id, sol_page):
   if status == 4 or status == 0:
      Print_Area_Report(file_id+".phc","report")
   elif status != 9:
      Print_Area_Sol(file_id+".sol",Uid, poly_id,"sols", sol_page)
   else:
      Print_Area_Report(file_id+".phc","report")
