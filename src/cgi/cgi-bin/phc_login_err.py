########################   Handle Login Error    #############################
def handle_error(error):
   if error == 1:
      print "Please enter your email"
   elif error == 2:
      print "Please enter your password"
   elif error == 3:
      print "Please enter your email and password"
   elif error == 4:
      print """<p>Sorry. Email and password don't match. Please choose from the following options:</p>
      <p><a href="login.py" id="login_retry">Try it again?</a></p>
      <p>or if you can't remember the password,</p>
      <p><a href="forgetpwd.py" id="login-forget-password">Forgot my password?</a></p>
      <p>or if you haven't registered,</p>
      <p><a href="register.py" id="register">Create my account?</a></p>"""
   elif error == 5:
      print """ Welcome to PHCpack Web Interface.  </p>
               Please sign in with your Email and password. </p>
               Redirecting to
               <a href='../index.html'>
               Sign in PHCpack Web Interface</a>
					<meta http-equiv="REFRESH" content="0;url=login.py">"""
   elif error == 6:
      print "Your account haven't been activated. Please check your email and use the link to activate your account."
   elif error == 10:
      print """ Thank you for using PHCpack Web Interface. You have successfully signed out.</p>
               Redirecting to
               <a href='../index.html'>
               Sign in PHCpack Web Interface</a>
               <meta http-equiv="REFRESH" content="2;url=login.py">"""
   else:
      print """ Sorry. There is some mistake in your account. Error = %d. </p> Report here. <a href='contact.py'>Contact Us.</a>""" % error
