#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug
if phc_debug:
   import cgitb; cgitb.enable()

import cgi, Cookie, os, sha, datetime, time

from phc_html import PrintHeader, Footer_Bar
from phc_sql import sql_checkticket, sql_resetpwd

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
   error = [0,0,0]
   infos = ["","","",""]
   
   # get email from form
   try:
      infos[0] = form['login'].value
   except KeyError:
      #print "No email"
      error[0] = 3
   # get ticket from form
   try:
      infos[1] = form['ticket'].value
   except KeyError:
      error[0] = 3
   
   # password
   try:
      infos[2] = form['passw'].value
   except KeyError:
      #print "Please enter your password"
      error[1] = 1
      error[0] = 1
   try:
      infos[3] = form['passw_check'].value
   except KeyError:
      #print "Please enter your password"
      error[2] = 1
      error[0] = 1
      
   # check if password match
   if error[0] == 0:
      if infos[2] != infos[3]:
         error[0] = 2
   
   return infos, error

def handle_error(error):
   if error[0] == 1:
      print "Please enter these information"
      error[1]=1
      error[2]=1
      info = ["Password","Confirm Password"]
      for i in xrange(0, len(info)):
         if error[i+1]==1:
            print ", " + info[i]
   elif error[0] == 2:
      print "Sorry. Your password don't match."
   elif error[0] == 3:
      print """Sorry. This reset link is invalid."""
   elif error[0] == 4:
      print """Sorry. We couldn't send the activation email to your mailbox. Please check your email address"""

def main():
   """
   Form to process Email.
   """
   PrintHeader('login')
   form = cgi.FieldStorage()
   infos, error = ProcessName(form)
   if not error[0]:
      error[0] = sql_checkticket(infos[0],infos[1])
      if not error[0]:
         sql_resetpwd(infos[1],infos[2])
         print """<p>Welcome to PHC Web Interface. Your password has been reset.</p>  <a href='login.py'>
               Sign in PHC Web Interface</a>"""
         
   if error[0]:
      print """<html><body>
               <div>
               <h1> Reset password </h1>
               <form method = "post" 
                     action = "resetpwd.py">"""

      # Password
      print """<table border="0"><tr>
      <td>Password</td>
      <td><input type = "password" name = "passw"
             size = "30" maxlength = "15">
      </td>
      </tr>
      <tr>
      <td>Confirm Password</td>
      <td>
      <input type = "password" name = "passw_check"
             size = "30" maxlength = "15">
      </td>
      </tr>
      </table>"""
      
      # Submit bottum
      print """
         <input type = "hidden" name = "login" value ="%s">
         <input type = "hidden" name = "ticket" value ="%s">
         <input type="submit" value="Reset">
         </form>
         </div>""" %(infos[0],infos[1])
   
   handle_error(error)
   Footer_Bar()
   print "</body></html>\n"

main()
