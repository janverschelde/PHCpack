#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug, style_forgetpwd
if phc_debug:
   import cgitb; cgitb.enable()

import cgi, Cookie, os, sha, datetime, time

from phc_html import PrintHeader, Footer_Bar
from phc_sql import sql_forgetpwd

"""
infos[0]: Email
infos[1]: password
infos[2]: repeated password
infos[3]: Name_First
infos[4]: Name_Last
infos[5]: Org

error[0]: 1. general error
          2. password don't match
error[1]: infos[0]: Email
error[2]: infos[1]: password
error[3]: infos[2]: repeated password
error[4]: infos[3]: Name_First
error[5]: infos[4]: Name_Last
error[6]: infos[5]: Org
"""

def ProcessName(form):
   """
   Processes name of Email form.
   Returns True if error, else False.
   """
   error = [0,0,0,0,0]
   infos = ["","","",""]
   form_list = ['login','Name_First','Name_Last','Org']

   for i in xrange(len(form_list)):
      try:
         infos[i] = form[form_list[i]].value
      except KeyError:
         #print "Please enter your email"
         error[i+1] = 1
         error[0] = 1

   return infos, error

def handle_error(error):
   if error[0] == 1:
      print "Please enter these information"
      info = ["Email", "First Name", "Last Name", "Organization" ]
      for i in xrange(0, len(info)):
         if error[i+1]==1:
            print ", " + info[i]
   elif error[0] == 2:
      print "Sorry. Your password don't match."
   elif error[0] == 3:
      print """Sorry. Your email has been registered. 
         <a href='Find_Password.html'>Find my password</a>"""
   elif error[0] == 4:
      print """Sorry. We couldn't send the email to your mailbox. Your information doesn't match our database. Please check your email address and other information."""

def main():
   """
   Form to process Email.
   """
   PrintHeader(style_forgetpwd)
   form = cgi.FieldStorage()
   infos, error = ProcessName(form)
   
   Email_Address = infos[0]
   Name_First    = infos[1]
   Name_Last     = infos[2]
   Org           = infos[3]
   
   if not error[0]:
      error[0] = sql_forgetpwd(Email_Address, Name_First, Name_Last, Org)
      if not error[0]:
         print """
               The reset link has been sent to your mailbox.</p>
               Please check your mailbox to reset your password."""
         
   if error[0]:
      print """<html><body>
               <div>
               <h1> Forget password </h1>
               <form method = "post" 
                     action = "forgetpwd.py">"""
      print """<table border="0">"""

      # Email Address
      print """<tr><td class="reg_form">Email</td>"""

      if Email_Address == '':
         print """<td><input type="text" name="login" spellcheck="false">
         </td><td></td>"""
      else:
         if error[0] == 4:
            print """<td><input type = "text" name = "login"  
                                 spellcheck="false"  value =%s ></td><td><font color="red">  X</font></td>"""% infos[0]
         else:
            print """<td><input type = "text" name = "login"  
                                 spellcheck="false"  value =%s></td><td>V</td>"""% Email_Address

      print "</tr>"

      form_list = [[Name_First, 'First Name', 'Name_First'],\
                   [Name_Last , 'Last Name','Name_Last'],\
                   [Org, 'Organization','Org']]
      
      # First Name
      for i in xrange(len(form_list)):
         print """<tr><td class="reg_form">%s</td>"""%form_list[i][1]
         if form_list[i][0] == '':
            print """<td><input type="text" name="%s" spellcheck="false" >
            </td></tr>"""%form_list[i][2]
         else:
            print """<td><input type = "text" name = "%s"  
                   spellcheck="false"  value =%s></td><td>V</td>"""\
                  %(form_list[i][2],form_list[i][0])
         print "</tr>"

      print "</table>"

      # Submit bottum
      print """<input type="submit" value="Submit">
         </form>
         </div>"""
   handle_error(error)
   Footer_Bar()
   print "</body></html>\n"

main()

