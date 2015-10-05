#!/usr/bin/python

"""
This script activates a new account.
"""

# Enter debug mode
from phc_config import phc_debug, style_activate
if phc_debug:
    import cgitb
    cgitb.enable()

import cgi, os

from phc_html import PrintHeader, Footer_Bar

from phc_sql import sql_activate

def ProcessName(form):
    """
    Processes name of Email form.
    Returns True if error, else False.
    """
    error = 0
    email = ""
    ticket = ""
    # get email from form
    try:
        email = form['login'].value
    except KeyError:
        #print "No email"
        error = 1
    # get ticket from form
    try:
        ticket = form['ticket'].value
    except KeyError:
        error = 1
    return email, ticket, error

def Activate_Link():
    """
    Prints the form that will be handled by phcweb.py.
    """
    print """<form method = "post" action = "phcweb.py">
    <p><b>Email</b></p>
    <input type="text" name="login" spellcheck="false" size="23" ></p>
    <p><b>Password</b></p>
    <input type = "password" name = "passw"
           size = "23" maxlength = "10"></p>
    <input type="submit" value="Send">
    </form>"""

def main():
   """
   Form to process Email.
   """
   PrintHeader(style_activate)
   form = cgi.FieldStorage()
   email, ticket, error = ProcessName(form)
   if not error:
      error = sql_activate(email, ticket)
   if error == 0:
      print """<p>Welcome to PHC Web Interface. Your account has been activated.</p>  <a href='login.py'>
               Sign in PHC Web Interface</a>"""
   elif error == 1:
      print " Welcome to PHC Web Interface. You haven't registered. Please register here."
   elif error == 2:
      print """<p> Welcome to PHC Web Interface. There is something wrong with the activation link. </p> Send me the link again?"""
      Activate_Link()
      print """ <p> Haven't registered? Click here. <a href='register.py'>Create a new account.</a>"""
   elif error == 3:
      print """<p>Welcome to PHC Web Interface. You have activated your account before.</p>  <a href='login.py'>
               Sign in PHC Web Interface</a>"""
   else:
      print """ Sorry. There is some mistake in your account. Error = %d. </p> Report here. <a href='contact.py'>Contact Us.</a>""" % error
   Footer_Bar()
   print "</body></html>\n"

main()
