#!/usr/bin/python

"""
This is the first script launched via index.html.
"""

from phc_config import phc_debug, style_login
if phc_debug:
    import cgitb
    cgitb.enable()

import cgi, Cookie, os, datetime
from phc_cookie import GetCookie
from phc_html import Footer_Bar, PrintHeader

def AskName():
    """
    Form to enter user name, using cookie c
    to show user name.
    """
    print """<div class="signin-box">
    <div  style = "max-width:220px; margin-left:auto; margin-right:auto" >
    <h2> Sign in </h2>
    <form method = "post"
    action = "phcweb.py">
    """
    v, w = GetCookie()
    if v == '':
        print """<p>
        <b>Email</b></p>
        <p><input type="email" name="login" spellcheck="false" >"""
    else:
        print """<p>
        <b>Email</b></p>
        <p><input type="email" name = "login"
                  spellcheck="false"  value =%s >"""% v
    if w == '':
        print """<p>
        <b>Password</b></p>
        <p><input type = "password" name = "passw">"""
    else:
        print """<p>
        <b>Password</b></p>
        <p><input type = "password" name = "passw" value = %s>""" % w
    if v != '' and w != '':
        print """
           </p>
           <input type="checkbox" name="rememberme" 
                checked="checked" value="yes">Remember me"""
    else:
        print """
        </p>
        <p><input type="checkbox" name="rememberme" value="yes">Remember me"""
    print """<input type = "hidden" name = "phcaction" value ="signin">
    <input style="float:right" type="submit" class="g-button g-button-submit" value="Sign in"></p>
    </form>
    <center>
    <a href="register.py" id="login-create-an-account">Create an account</a>
    </center>
    <center>
    <a href="forgetpwd.py" id="login-forget-password">Forgot your password?</a>
    </center></div>
    </div>"""

def main():
    """
    Form to process login.
    """
    PrintHeader(style_login)
    print """<div class="center1">"""
    AskName()
    print"""<div class="photo"><a href="demo.py" id="login-create-an-account"><img src="../imag/Path.jpg" alt="PHC Web Interface"></a></div>"""
    print "</div>"
    Footer_Bar()
    print "</div></body></html>\n"

main()
