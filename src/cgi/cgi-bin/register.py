#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug, style_reg
if phc_debug:
    import cgitb
    cgitb.enable()
   
import cgi, Cookie, os, sha, datetime, time

from phc_html import PrintHeader, Footer_Bar
from phc_sql import sql_reg

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
   error = [0,0,0,0,0,0,0]
   infos = ["","","","","",""]
   try:
      infos[0] = form['login'].value
   except KeyError:
      #print "Please enter your email"
      error[1] = 1
      error[0] = 1

   # password
   try:
      infos[1] = form['passw'].value
   except KeyError:
      #print "Please enter your password"
      error[2] = 1
      error[0] = 1
   try:
      infos[2] = form['passw_check'].value
   except KeyError:
      #print "Please enter your password"
      error[3] = 1
      error[0] = 1
   try:
      infos[3] = form['Name_First'].value
   except KeyError:
      error[4] = 1
      error[0] = 1
   try:
      infos[4] = form['Name_Last'].value
   except KeyError:
      error[5] = 1
      error[0] = 1
   try:
      infos[5] = form['Org'].value
   except KeyError:
      error[6] = 1
      error[0] = 1
   # check if password match
   if error[0] == 0:
      if infos[1] != infos[2]:
         error[0] = 2
   return infos, error

def handle_error(error):
   if error[0] == 1:
      print "Please enter these information"
      info = ["Email", "Password","Password Confirm",  "First Name", "Last Name", "Organization" ]
      for i in xrange(0, len(info)):
         if error[i+1]==1:
            print ", " + info[i]
   elif error[0] == 2:
      print "Sorry. Your password don't match."
   elif error[0] == 3:
      print """Sorry. Your email has been registered. 
         <a href='Find_Password.html'>Find my password</a>"""
   elif error[0] == 4:
      print """Sorry. We couldn't send the activation email to your mailbox. Please check your email address"""

def main():
   """
   Form to process Email.
   """
   PrintHeader(style_reg)
   form = cgi.FieldStorage()
   infos, error = ProcessName(form)
   if not error[0]:
      error[0], ticket = sql_reg(infos)
      if not error[0]:
         print """
               Congratulations! Register successful.</p>
               Please check your mailbox to activate your account."""
         
   if error[0]:
      print """<html><body>
               <div>
               <h1> Register </h1>
               <form method = "post" autocomplete="off"
                     action = "register.py">"""
      print """<table border="0">"""


      form_list = [[infos[0], 'Email'           , 'login'      ],\
                   [infos[1], 'Password'        , 'passw'      ],\
                   [infos[2], 'Confirm Password', 'passw_check'],\
                   [infos[3], 'First Name'      , 'Name_First' ],\
                   [infos[4], 'Last Name'       , 'Name_Last'  ],\
                   [infos[5], 'Organization'    , 'Org'        ]]
      
      # First Name
      for i in xrange(len(form_list)):
         print """<tr><td class="reg_form">%s</td>"""%form_list[i][1]
         if i==1 or i==2:
            print """<td><input type="password" name="%s" maxlength="15">
                     </td></tr>"""%form_list[i][2]
         else:
            if form_list[i][0] == '':
               print """<td><input type="text" name="%s" spellcheck="false" >
                        </td></tr>"""%form_list[i][2]
            else:
               print """<td><input type = "text" name = "%s"  
                      spellcheck="false"  value =%s></td><td>V</td>"""\
                     %(form_list[i][2],form_list[i][0])
            print "</tr>"

      print """</table>
         <input type="submit" value="Submit">
         </form>
         </div>"""
   handle_error(error)
   Footer_Bar()
   print "</body></html>\n"

main()

