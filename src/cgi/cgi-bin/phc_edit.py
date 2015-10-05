from phc_html import *
from phc_file import get_file_id, check_file_exist
from os import path
from phc_client import PHC_Server

def PHC_Edit_Hom(Uid,form):  
   poly_id = form['poly_id'].value
   file_id = get_file_id(Uid, poly_id)
   if not path.isfile(file_id+".start"):
      data = PHC_Server(Uid, '7', poly_id)
      if data != 0 and path.isfile(file_id+".start"):
         Print_Poly_File_New(Uid, poly_id, 1)
      else:
         Print_Poly_File_New(Uid, poly_id, 0)
         print """<p><font color='red'><b>Error:</b> the original system doesn't have a start system. Please just solve it directly.</font></p>"""
         
   else:
      Print_Poly_File_New(Uid, poly_id, 1)

def PHC_Edit(Uid, form):
   if form.has_key('homotopy'):
      PHC_Edit_Hom(Uid,form)

   elif form.has_key('poly'):
      poly = form['poly'].value
      Print_Poly_New(Uid, poly)

   elif form.has_key('poly_id'):
      poly_id = form['poly_id'].value
      Print_Poly_File_New(Uid, poly_id)

   else:
      Print_Poly_Empty(Uid)
