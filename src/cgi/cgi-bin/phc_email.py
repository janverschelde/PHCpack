def  phc_email(addrs, subject, msg_cont):
   import smtplib
   from email.mime.text import MIMEText

   from phc_config import phcstmp, phcmail, phcmailps
   
   msg = MIMEText(msg_cont)
   msg['Subject'] = subject
   msg['To'] = addrs

   # The actual mail send  
   server = smtplib.SMTP(phcstmp)  
   server.starttls()  
   server.login(phcmail, phcmailps)
   try:
      server.sendmail(phcmail, [addrs,phcmail], msg.as_string())
   except smtplib.SMTPRecipientsRefused:
      return 4
   server.quit()
   return 0 # There might be some errors when send email, I haven't checked that.

def Activate_Gmail(addrs, first_name, Ticket):
   import smtplib
   from email.mime.text import MIMEText
   from phc_config import hostadd
   from phc_email import phc_email

   subject = "Welcome %s to PHC Web Interface" % first_name
   msg_cont = """Hello %s,\n\n   Welcome to PHC Web Interface. Please click the following link to activate your account.\n\n
%s/cgi-bin/activate.py?login=%s&ticket=%s""" % (first_name,hostadd,addrs,Ticket)
   
   phc_email(addrs, subject, msg_cont);

def Resetpwd_Gmail(addrs, first_name, Ticket):
   import smtplib
   from email.mime.text import MIMEText
   from phc_config import hostadd
   from phc_email import phc_email
   
   subject = "Your PHCWEB password has been reset."
   msg_cont = """Hello %s,\n\n   Welcome to PHC Web Interface. Please click the following link to reset your password.\n\n
%s/cgi-bin/resetpwd.py?login=%s&ticket=%s""" % (first_name,hostadd,addrs,Ticket)
   
   phc_email(addrs, subject, msg_cont);
