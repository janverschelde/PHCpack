from phc_client import PHC_Server
from phc_file import get_file_id
from phc_html import PHC_View_Bar,Print_Text_File
import os

####################### Polynomial Management ################################   
def PHC_Mypoly(form, Uid,phcdb):
   import time
   if form.has_key('poly_id'):
      poly_id = form['poly_id'].value
   else:
      poly_id = ''
   
   s = phcdb.allpoly(Uid)
   print """<p><h2><a name="mypoly" id="mypoly">All Polynomial Systems</a></h2>"""
   print """<table class="mypolys" rules ="rows" border =1><tr><th>Name</th><th>Sols</th><th  class="hidden_column">Dim</th><th  class="hidden_column">Create</th><th>Report</th><th>Actions</th></tr>"""
   for poly in s:
      # Name
      print """<tr align="center"><td><a href="phcweb.py?poly_id=%s">%s </a></td>"""% (poly[2],poly[2])
      
      # Solutions
      if poly[6] != None:
         print """<td><a href="phcweb.py?print=%s&r=1&m=0">%d </a></td>"""% (poly[2], poly[6])
      elif poly[5] == 0:     # in the pipeline, it doesn't work
         print "<td>Killed</td>"
      elif poly[5] == 1:
         print "<td>Solving</td>"
      else:
         print "<td>Error</td>"
      
      # Dim
      print """<td class="hidden_column">%d</td>"""% poly[3]
         
      # Create Time
      print """<td class="hidden_column">%s</td>"""% poly[4]

      # PHC Report
      if poly[6] != None:
         print """<td align="right"><a href="phcweb.py?print=%s&r=2&m=0">%.2f s</a></td>"""% (poly[2], float(poly[7]))
      else:
         print "<td></td>"
         
      # Actions
      if poly[5] != 1:
         if poly_id == '' or poly_id == poly[2]:
            print """<td><a href="phcweb.py?Delete=%s#mypoly">Del </a></td></tr>"""% (poly[2])
         else:
            print """<td><a href="phcweb.py?Delete=%s&poly_id=%s#mypoly">Del </a></td></tr>"""% (poly[2],poly_id)
            
      else:
         print "<td></td></tr>"
   
   print "</table>"
   print """<div><a href="#top">Back to Top</a> </div>"""

def PHC_Print(Uid, form):
   if form.has_key('r') and form.has_key('m'):
      f = form['print'].value
      report = form['r'].value
      formt = form['m'].value
      PHC_View(Uid, f, report, formt)
   else:
      print """ Welcome to PHC Web Interface. </p>
                You request can't be served. Please sign in with your Email and password. </p>
               Redirecting to
               <a href='login.py'>
               Sign in PHC Web Interface</a>
					<meta http-equiv="REFRESH" content="0;url=login.py">"""


def PHC_View(Uid, poly_id, report, formt):
   file_id = get_file_id(Uid, poly_id)
   if not os.path.isfile(file_id+".maple") or not os.path.isfile(file_id+".python") :
      data = PHC_Server(Uid,'6',file_id)
   if report == '0':
      Print_Text_File(file_id+".poly")
   elif report == '1':
      PHC_View_Bar(poly_id, formt)
      if formt == '0':
         Print_Text_File(file_id+".sol")
      elif formt == '1':
         if not os.path.isfile(file_id+".python"):
            print "Not Ready, Please Refresh."
         else:
            Print_Text_File(file_id+".python")
      elif formt == '2':
         if not os.path.isfile(file_id+".maple"):
            print "Not Ready, Please Refresh."
         else:
            Print_Text_File(file_id+".maple")
   elif report == '2':
      Print_Text_File(file_id+".phc")
   else:
      print report + ": not such option"
