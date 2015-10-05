#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug, style_contact
if phc_debug:
   import cgitb; cgitb.enable()
   
import cgi, Cookie, os, sha, datetime

from phc_html import PrintHeader, Footer_Bar
from phc_email import phc_email

from sql_login import phc_login_simple

def GetCookie():
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   if c.has_key('login'):
      v = c['login'].value
   else:
      v = ""
   return v

def Contact_Process(form):
   """
   Processes name of Email form.
   Returns True if error, else False.
   """
   error = [0,0,0,0,0,0,0]
   infos = ["","","","","",""]
   try:
      infos[0] = form['login'].value
   except KeyError:
      infos[0] = GetCookie()
      if infos[0] == "":
         #print "Please enter your email"
         error[1] = 1
         error[0] = 1
   try:
      infos[1] = form['Name_First'].value
   except KeyError:
      error[2] = 1
      error[0] = 1
   try:
      infos[2] = form['Name_Last'].value
   except KeyError:
      error[3] = 1
      error[0] = 1
   try:
      infos[3] = form['Org'].value
   except KeyError:
      error[4] = 1
      error[0] = 1
   try:
      infos[4] = form['Subject'].value
   except KeyError:
      error[5] = 1
      error[0] = 1
   try:
      infos[5] = form['Message'].value
   except KeyError:
      error[6] = 1
      error[0] = 1
   return infos, error

def handle_error(error):
   if error[0] == 1:
      print "Please enter these information"
      info = ["Email","First Name", "Last Name", "Organization", "Subject","Message" ]
      for i in xrange(0, len(info)):
         if error[i+1]==1:
            print ", " + info[i]
   elif error[0] == 0:
      print """ <p>Thank you. Your message has been sent to us. A copy has been forwarded to your mailbox.</p> Redirect to
               <a href='cookie_login.py'>PHC Web Interface</a>"""
   elif error[0] == 4:
      print """ Sorry. We couldn't verify your email address. Please double check your email address."""

def Contact_Email(infos):
   addrs = infos[0]
   msg_cont = " Welcome to PHC Web Interfact. You just left a message to us. We will reply to you as soon as possible. This is a copy of your message.\n\n From: %s\n First Name: %s \n Last Name: %s \n Org: %s \n\n Subject: %s\n\n Message: %s" %  (infos[0], infos[1],infos[2], infos[3],infos[4], infos[5])
   subject = "%s : %s" % (infos[1], infos[4]);
   error = phc_email(addrs, subject, msg_cont);

def Print_JS():
    # empty notification is not finished
    print"""<script type="text/javascript">
       function clearTextArea() {
          document.contact.Message.value='';
       }
       function emptyTextArea() {
          document.contact.Message.value='please enter the polynomial system';
       }
       </script>"""

def main():
   """
   Form to process login.
   """
   PrintHeader(style_contact)
   form = cgi.FieldStorage()
   infos, error = Contact_Process(form)
   send_gmail = 0
   if not error[0]:
      error[0] = Contact_Email(infos)
      if not error[0]:
         send_gmail = 1

   phc_login_simple()

   print """<html><body>
            <p><h2>Contact Us</h2>
            <form name = "contact" method = "post" 
                  action = "contact.py">
            <table border="0" width = "100%"><tr>"""

   form_list = [['Email','login'],\
                ['First Name', 'Name_First'],\
                ['Last Name','Name_Last'],\
                ['Organization','Org'],\
                ['Subject','Subject']]

   for i in xrange(len(form_list)):
      print """<td class="reg_form">%s</td>"""%form_list[i][0]
      if infos[i] == '':
         print """<td><input type="text" name="%s" spellcheck="false"  >
         </td></tr>"""%form_list[i][1]
      else:
         print """<td><input type = "text" name = "%s"  
                              spellcheck="false"  value ="%s" ></td><td>V</td>
         </tr>"""%(form_list[i][1],infos[i])

   # Message
   print """<tr><td class="reg_form" valign="top">Message</td>"""
   Print_JS()
   if infos[5] == '':
      print """
      <td><textarea name="Message" rows = 8  onfocus="clearTextArea();">Please enter your information</textarea></td>
      </tr>
      </table>"""
   else:
      print """<td><textarea name="Message" rows = 8 >%s</textarea></td><td valign="top">V</td>
      </tr>
      </table>"""% infos[5]

   if send_gmail == 1:
      print """<p>Thank you. Your message has been sent to us. A copy of message has been forwarded to your own email address.</p>
    <p>We will reply to you as soon as possible.</p>"""
   else:
      # Submit bottum
      print """
         <input type="submit" value="Submit">
         </form>"""
   handle_error(error)
   Footer_Bar()
   print "</body></html>\n"

main()
