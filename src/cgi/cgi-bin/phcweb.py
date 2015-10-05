#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug
if phc_debug:
    import cgitb
    cgitb.enable()

import cgi, Cookie, os, sha, datetime

from sql_login import PHC_login
from phc_client import PHC_Server
from phc_html import *

from phc_config import style_web

from phc_login_err import *
from phc_solve import *
from phc_edit import *
from phc_update import *
from phc_preview import *
from phc_kill import *
from phc_mypoly import *

#####################     File Management Section    ######################

def PHC_Rename(Uid, s,t):
   return PHC_Server(Uid,'4',s+'|'+t)

def PHC_Delete(Uid,name):  
   return PHC_Server(Uid,'5',name)

def polypid(Uid):
   f = "poly.pid"
   file = open(f,'r')
   Pid = file.readline()
   file.close()
   return '1'

############################  PHC Function Section ############################

def PHCWEB(form,Uid, Status, phcdb): ### combine all these key to phcaction
   from os import path

   # Job 0 and 1: Solve polynomial by homotopy
   if form.has_key('Solve'):
      if phc_debug: print "<p>Action: Solve</p>"
      PHC_Solve(Uid, Status, form, phcdb)

   # Job 2: Kill
   elif form.has_key('Kill'):
      if phc_debug: print "<p>Action: Kill</p>"
      PHC_Kill(Uid, form, phcdb)

   # Manage 0: New Polynomial
   elif form.has_key('New'):
      if phc_debug: print "<p>Action: New</p>"
      Print_Poly_Empty(Uid)

   # Manage 1: View Polynomials
   elif form.has_key('print'):
      if phc_debug: print "<p>Action: Print</p>"
      PHC_Print(Uid,form)

   # Manage 2: Save Polynomials XXX rename need to be done
   elif form.has_key('Save'):
      if phc_debug: print "<p>Action: Save</p>"
      poly_id = form['poly_id'].value
      new_poly_id = form['newname'].value
      data = PHC_Rename(Uid, poly_id, new_poly_id)
      if data != 0:
         form['poly_id'].value = new_poly_id
         PHC_Update(Uid, '0', new_poly_id, phcdb)
      else:
         PHC_Update(Uid, '0', poly_id, phcdb)

   # Manage 3: Delete Polynomials
   elif form.has_key('Delete'):
      if phc_debug: print "<p>Action: Delete</p>"
      if form.has_key('poly_id'):
         poly_id = form['poly_id'].value
      else:
         poly_id = ''
      poly_id_del = form['Delete'].value
      data = PHC_Delete(Uid,poly_id_del)
      PHC_Update(Uid, '0', poly_id, phcdb)

   # Manage 6: Preview
   elif form.has_key('Preview'):
      if phc_debug: print "<p>Action: Preview</p>"
      PHC_Preview(Uid, form)

   # Manage 4: Edit
   elif form.has_key('Edit') or form.has_key('homotopy'):
      if phc_debug: print "<p>Action: Edit\n"
      PHC_Edit(Uid, form)

   # Manage 5: Update
   elif form.has_key('Update') :
      if phc_debug: print "<p>Action: Update</p>"
      Pid = form['Pid'].value
      PHC_Update(Uid,Pid,"current",phcdb)

   # Manage 7: view poly
   elif form.has_key('poly_id'):
      if phc_debug: print "<p>Action: Update2</p>"
      poly_id = form['poly_id'].value
      sol_page = 0
      if form.has_key('sol_page'):
         sol_page = int(form['sol_page'].value)
      PHC_Update(Uid,'0',poly_id,phcdb, sol_page)

   # Manage 8: User Come back
   elif Status != None:
      if phc_debug: print "<p>Action: Update</p>"
      PHC_Update(Uid,'0',"current",phcdb)

   # Manage 0: New Polynomial
   else:
      if phc_debug: print "<p>Action: New</p>"
      Print_Poly_Empty(Uid)
   
   PHC_Mypoly(form, Uid, phcdb)


########################   Main Function    #############################

def main():
   """
   Form to process Email.
   """
   form = cgi.FieldStorage()
   # Processing cookie and user pass ### put this part into cookie_login.py
   from phc_sql import PHCSQL
   phcdb = PHCSQL()
   error, Uid, Name_First, Status = PHC_login(form,phcdb)
   # PHCpack Web Interface Header
   PrintHeader(style_web)
   #print "Uid = %s</p>" % Uid
   if not error or error == 7: # error 7: folder not exist, recovery
      Print_User_Bar(Name_First)
      PHCWEB(form,Uid, Status, phcdb)
      # Need to close database
   else:
      handle_error(error)
   phcdb.db.close()
   Footer_Bar()
   print "</div></body></html>\n"

main()

