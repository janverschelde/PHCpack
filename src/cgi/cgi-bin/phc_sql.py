import MySQLdb

class PHCSQL():
   def __init__(self):
      """
      Start SQL Datebase and Return 
      """
      from phc_config import sql_host, sql_user, sql_passwd, sql_datebase
      
      self.db = MySQLdb.connect(sql_host, sql_user, sql_passwd, sql_datebase)
      self.d = self.db.cursor()
      
   def try_excute(self, commd):
      try:
         self.d.execute(commd)
      except MySQLdb.Error,e:
         print "# SQL Error = %s, %s"%(e[0],e[1])
         self.d.close()
         self.db.close()
         self.__init__()
         self.d.execute(commd)
   
   def run_commd(self, commd):
      self.try_excute(commd)
      self.db.commit()
      return 0

   def query(self,commd):
      self.try_excute(commd)
      s = self.d.fetchall()
#      self.db.commit()
      return s
   
   def close(self):
      self.db.close()

 ##################### account management #######################
   def userstatus(self, Uid, Status):
      commd = "Update users SET Status = %s WHERE Folder = '%s';" %(Status, Uid)
      self.run_commd(commd)
      return 0
      
   def check_userstatus(self, Uid):
      commd = "SELECT Status FROM users WHERE Folder='%s'"%(Uid)
      self.d.execute(commd)
      s = self.d.fetchall()
      Status = s[0][0]
      return Status
      
   def check_poly_status(self, Uid, poly_id):
      commd = "SELECT Status FROM polys WHERE Uid='%s' and Name='%s'"%(Uid, poly_id)
      self.d.execute(commd)
      s = self.d.fetchall()
      if len(s)==0:
         Status = 9
      else:
         Status = s[0][0]
      
      return Status

   def check(self, usermail, userpwd):
      """
      Check usermail and userpwd.
      If match, return the user's Folder and Name_First
      """
      error = 0
      Folder = ''
      Name_First = ''
      Status = None
      
      if len(userpwd) != 40:
         import sha
         userpwd = sha.new(userpwd).hexdigest()

      commd = "SELECT '" + userpwd + \
         "'= passwd FROM users WHERE Email='" + usermail +"'"
      self.d.execute(commd)
      s = self.d.fetchall()

      if len(s)==0 or s[0][0]==0:
         error = 4
      else:
         commd = "SELECT Folder, Name_First, Status FROM users WHERE Email='" + usermail +"'"
         self.d.execute(commd)
         s = self.d.fetchall()
         Folder = str(s[0][0])
         Name_First = str(s[0][1])
         Status = s[0][2]
         
      return error, Folder, Name_First, Status 
      
##################### polynomials management #######################
   def newpoly(self, Uid, Name, Deg, Ctime, Status='1'):
      # pop a windows to let user check same name polynomial
#      commd = "DELETE FROM polys WHERE Uid = '%s' and Name = '%s';"%(Uid,Name)
#      self.run_commd(commd)
      from datetime import datetime
      file_ctime = datetime.fromtimestamp(Ctime)
      commd = "INSERT INTO polys (Uid,Name, Deg, Ctime, Status) VALUES ('%s','%s','%d','%s','%s');" %(Uid, Name,Deg,file_ctime,Status)
      self.run_commd(commd)
      return 0

   def killpoly(self,Uid, poly_id):
      commd = "Update polys SET Status = 0 WHERE Uid = '%s' and NAME = '%s';" %(Uid, poly_id)
      self.run_commd(commd)
      return 0

   def slvpoly(self, Uid, Name):
      from phc_file import Sols, Slvtime
      sols = Sols(Uid, Name)
      slvtime = Slvtime(Uid, Name)
      commd = "Update polys SET Status = 2, Sols = %d, Time ='%s' WHERE Uid = '%s' and Name = '%s';" %(sols, slvtime, Uid, Name)
      self.run_commd(commd)
      return 0

   def slvpoly2(self, Uid, poly_id):
      from phc_file import Deg_Sols, solve_time, get_file_id
      file_id = get_file_id(Uid, poly_id)
      sols, deg = Deg_Sols(file_id)
      slvtime = solve_time(file_id)
      commd = "Update polys SET Status = 2, Sols = %d, Time ='%s' WHERE Uid = '%s' and Name = '%s';" %(sols, slvtime, Uid, poly_id)
      self.run_commd(commd)
      return 0

   def renamepoly(self, Uid, Name, NewName):
      commd = "Update polys SET Name = '%s' WHERE Uid = '%s' and Name = '%s';" %(NewName, Uid, Name)
      self.run_commd(commd)
      return 0
      
   def delpoly(self,Uid,Name):
      commd = "DELETE FROM polys WHERE Uid = '%s' and Name = '%s';" %(Uid, Name)
      self.run_commd(commd)
      return 0
      
   def allpoly(self,Uid):
      query = "SELECT * FROM polys WHERE Uid = '%s' ORDER BY Polyid DESC;" %(Uid)
      s = self.query(query)
      return s
 
def StartSQL():
   """
   Start SQL Datebase and Return 
   """
   import MySQLdb
   from phc_config import sql_host, sql_user, sql_passwd, sql_datebase

   db = MySQLdb.connect(sql_host, sql_user, sql_passwd, sql_datebase)
   d = db.cursor()
   return d,db 
 
def sql_reg(infos):
   """
   Processes password of Email form.
   """
   import sha, datetime
   from time import time
   from phc_email import Activate_Gmail
   d,db = StartSQL()
   np = sha.new(infos[1]).hexdigest()
   commd = "SELECT Uid FROM users WHERE Email='" + infos[0] +"'"
   d.execute(commd)
   s = d.fetchall()
   if len(s)==0:
      ticket = sha.new(str(time())).hexdigest()
      error = Activate_Gmail(infos[0], infos[3] , ticket)
      if error:
         db.close()
         return 4,"" # Email address is not right
      else:
         commd = "INSERT INTO users VALUES(NULL,'%s','%s','%s','%s','%s','%s','%s','%s',NULL)" % (infos[3],infos[4],infos[0],infos[5],str(np),str(datetime.date.today()),ticket,"")
         d.execute(commd)
         s = d.fetchall()
         db.commit()
         db.close()
         return 0, ticket
   else:
      db.close()
      return 3, ""

def sql_activate(email, ticket):
   """
   Processes email and ticket and generate the folder for the user.
   """
   import os, sha, time
   d,db = StartSQL()
   commd = "SELECT '%s'= Ticket FROM users WHERE Email='%s'" %(ticket, email)
   d.execute(commd)
   s = d.fetchall()
   if len(s)==0 or s[0][0]==0:
      return 2
   commd = "SELECT Folder FROM users WHERE Email='%s'" % email
   d.execute(commd)
   s = d.fetchall()
   if s[0][0] =="":
      Folder = sha.new(str(time.time())).hexdigest()
      from phc_client import PHC_Server
      PHC_Server(Folder,'3',Folder)
      commd = "UPDATE users SET Folder = '%s' WHERE Email = '%s'" % (Folder, email)
      d.execute(commd)
      s = d.fetchall()
      db.commit()
      db.close()
      return 0

   elif len(s[0][0]) == 40:
      db.close()
      Upath="../users/" + s[0][0]   
      Open_Folder = os.path.isdir(Upath)
      if Open_Folder:    
         return 3 # The account has been activated, and the folder exists.
      else:
         return 4 # The account has been activated, but the folder disappears.
   else:
      db.close()
      return 5 # The Folder data in SQL is wrong.
      
def sql_forgetpwd(Email_Address, Name_First, Name_Last, Org):
   import sha, time
   from phc_email import Resetpwd_Gmail
   
   d,db = StartSQL()
   error = 0
   
   commd = "SELECT Name_First='%s' AND Name_Last='%s' AND Org='%s' FROM users WHERE Email='%s'"%(Name_First, Name_Last, Org, Email_Address)
   d.execute(commd)
   
   s = d.fetchall()
   ticket = ''
   
   if len(s)==0 or s[0][0]==0:
      error = 4
   else:
      ticket = sha.new(str(time.time())).hexdigest()
      error = Resetpwd_Gmail(Email_Address, Name_First, ticket)
      if not error:
         commd = "UPDATE users SET Ticket = '%s', passwd = NULL WHERE Email = '%s';" % (ticket, Email_Address)
         d.execute(commd)
         s = d.fetchall()
         db.commit()
   
   db.close()
      
   return error

def sql_checkticket(email, ticket):
   """
   Processes email and ticket and generate the folder for the user.
   """
   import os, sha, time
   d,db = StartSQL()
   commd = "SELECT '%s'= Ticket FROM users WHERE Email='%s'" %(ticket, email)
   d.execute(commd)
   s = d.fetchall()
   if len(s)==0 or s[0][0]==0:
      return 1
   else:
      return 0


def sql_resetpwd(ticket, pwd):
   import sha, time
   pwd = sha.new(pwd).hexdigest()
   
   d,db = StartSQL()
   new_ticket = sha.new(str(time.time())).hexdigest()
   commd = "UPDATE users SET passwd = '%s', Ticket = '%s' WHERE Ticket = '%s';" % (pwd, new_ticket, ticket)
   d.execute(commd)
   s = d.fetchall()
   db.commit()
   
   return 0

