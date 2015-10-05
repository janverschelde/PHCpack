from phc_edit import *
from phc_update import *

from time import sleep
from phc_config import solve_sleep_time

def hom_check(Uid, file_id_start, poly):
    """
    homotopy check, not finished
    """
    return 1


def PHC_Solve(Uid, Status, form, phcdb):
   if Status != None:
      PHC_Update(Uid,'0',"current",phcdb)
   else:
      if not form.has_key('poly'):
         PHC_Edit(Uid, form)
      else:
         Poly = form['poly'].value.strip()
         if form.has_key('homotopy'):
#            error_hom_check = hom_check(Uid, form['homotopy'].value, Poly)
#            
#            if not error_hom_check:
            data = PHC_Server(Uid,'7', 'current|' + form['homotopy'].value + '|'+ Poly)
#            else:
#               data = 0
#               PHC_Edit_Hom(Uid,form)
#               print "<p><font color='red'><b>Error:</b> These two polynomials doesn't match.</font></p>"
               
            if data != 0:
               print "<p>Start polynomial systems = %s</p>"%form['homotopy'].value
         else:
            data = PHC_Server(Uid,'1',Poly)

         if data != "error" and data != 0:
            phcdb.userstatus(Uid, "1")
            sleep(solve_sleep_time)
            PHC_Update(Uid,'0',"current",phcdb)
