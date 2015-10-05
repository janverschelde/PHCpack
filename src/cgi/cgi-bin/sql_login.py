import cgi

from os import chdir, path 

from phc_client import PHC_Server
from phc_cookie import UpdateCookie, CookieUid, RemoveUid

def ProcessName(form):
   """
   Processes name of Email form.
   Returns True if error, else False.
   """
   error = 0
   usermail = ''
   userpwd = ''
   try:
      usermail = form['login'].value
   except KeyError:
      #print "Please enter your email"
      error = 1
   try:
      userpwd = form['passw'].value
   except KeyError:
      #print "Please enter your password"
      error = error + 2
   return error, usermail, userpwd

def AccessFolder(Folder):
   error = 0

   ## The account hasn't be activated
   if Folder == "None":
      error = 5
      return error
   elif len(Folder) != 40:
      error = 6
      return error
   Upath="../users/"+Folder
   Open_Folder = path.isdir(Upath)
   if not Open_Folder:
      PHC_Server(Folder,'3',Folder)
      error=7
   return error

def PHC_login(form,phcdb):
   import Cookie
   error = 0
   Folder = ''
   Name_First=''
   Status = None
   if form.has_key('Signout'):
      error = 10  # Error 10 is Sign out
      UpdateCookie(form, '', '', error)
   else:
      if form.has_key('phcaction'):
         error, usermail, userpwd = ProcessName(form)
         if not error:
            error, Folder, Name_First, Status= phcdb.check(usermail, userpwd)
            UpdateCookie(form, Folder, Name_First, error)
      else:
         Folder, Name_First = CookieUid()
         if Folder == '':
            error = 5 # haven't signin
         else:
            Status = phcdb.check_userstatus(Folder)
   if not error:
      error = AccessFolder(Folder)
   return error, Folder, Name_First, Status
   
def phc_login_simple():
   from sql_login import PHC_login
   from phc_sql import PHCSQL
   from phc_html import Print_User_Bar
   form = cgi.FieldStorage()
   # Processing cookie and user pass ### put this part into cookie_login.py
   phcdb = PHCSQL()
   error, Uid, Name_First, Status = PHC_login(form,phcdb)
   #print "Uid = %s</p>" % Uid
   if not error or error == 7: # error 7: folder not exist, recovery
      Print_User_Bar(Name_First, 1)
   else:      
      Print_User_Bar(Name_First, 0)
   phcdb.db.close()

