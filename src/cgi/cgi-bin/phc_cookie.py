#!/usr/bin/python

import cgi, Cookie, os, sha, datetime

def CookieUid():
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   if c.has_key('Uid'):
      Uid = c['Uid'].value
      if c.has_key('Name_First'):
         Name_First = c['Name_First'].value
   else:
      Uid = ''
      Name_First = ''
   return Uid, Name_First

def GetCookie():
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   if c.has_key('login'):
      v = c['login'].value
   else:
      v = ''
   if c.has_key('passw'):
      w = c['passw'].value
   else:
      w = ''
   return v,w

def WriteCookie(form,Uid):
   """
   Retrieves cookie and uses form to update.
   """
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   expires = datetime.datetime.now() + datetime.timedelta(days=7)
   if form.has_key('login'):
      c['login'] = form['login'].value
      c['login']['expires'] = expires.strftime("%a, %d %b %Y %H:%M:%S GMT")
   if form.has_key('passw'):
      if len(form['passw'].value) != 40:
         p = form['passw'].value
         d = sha.new(p).hexdigest()
         c['passw'] = d
      else:
         c['passw'] = form['passw'].value
      c['passw']['expires'] = expires.strftime("%a, %d %b %Y %H:%M:%S GMT")
   return c

def RemoveCookie():
   """
   Retrieves cookie and uses form to update.
   """
   expires = datetime.datetime.now() - datetime.timedelta(days=100)
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   if c.has_key('login'):
      c['login']['expires'] = expires.strftime("%a, %d %b %Y 23:59:59 GMT")
   if c.has_key('passw'):
      c['passw']['expires'] = expires.strftime("%a, %d %b %Y 23:59:59 GMT")
   return c

def UpdateCookie(form,Uid,Name_First,error):
   if form.has_key('rememberme') and not error:
      c = WriteCookie(form,Uid)
   else:
      c = RemoveCookie()
   c['Uid']=Uid
   c['Name_First']=Name_First # don't store name in cookie
   print c

def RemoveUid():
   """
   Retrieves cookie and uses form to update.
   """
   update = 0
   expires = datetime.datetime.now() - datetime.timedelta(days=1)
   if os.environ.has_key('HTTP_COOKIE'):
      c = Cookie.Cookie(os.environ['HTTP_COOKIE'])
   else:
      c = Cookie.Cookie()
   if c.has_key('Uid'):
      c['Uid'] = ''
   if c.has_key('Name_First'):
      c['Name_First'] = ''
   print c
