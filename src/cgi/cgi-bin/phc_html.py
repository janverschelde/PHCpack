def PrintHeader(style):
   """
   writes title and header of page
   """
   from phc_config import title, header, small_screen

   print """Content-type: text/html;charset=utf-8

<!doctype html>
<head>
<title>%s</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" media="screen and (max-width:%dpx)" href="/style/%s.css" >
<link rel="stylesheet" media="screen and (min-width:%dpx)" href="/style/%s.css" >
<link rel="stylesheet" href="/style/%s.css" >
<link rel="stylesheet" href="/style/%s.css" >
</head> 
<body>
<div class="content">
<div class="header"><h1><a href="phcweb.py" name="top" id="top" style="text-decoration:none; color:black">%s</a></h1></div>
""" % (title, small_screen, style[0], small_screen+1, style[1], style[2], style[3], header)

def Print_User_Bar(Name_First, user_login = 1):
   if user_login:
      print """ <div class='userbar'><table style="width:100%%"><tr>
      <td>Hello, <a href="phcweb.py#mypoly">%s</a></td>
      <td style="text-align:right" ><a href="phcweb.py?Signout=1">Sign Out</a></td></tr></table></div>""" % Name_First
   else:
      print """ <div class='userbar'><table style="width:100%%"><tr>
      <td>Wecome PHCpack user, </td>
      <td><a href="register.py">Register</a></td>
      <td style="text-align:right" ><a href="login.py">Login</a></td>
      </tr></table></div>"""

def Footer_Bar():
   print """<div class="footer"><table class="footer_link"><tr>
   <td style="width:25%%"><a href="demo.py">Demo</a></td>
   <td style="width:25%%"><a href="about.py">About</a></td>
   <td style="text-align:center; width:25%%"><a href="http://homepages.math.uic.edu/~jan/download.html">PHCpack</a></td>
   <td style="text-align:right; width:25%%"><a href="contact.py">Contact</a></td></tr></table>
   <div class="grant"> <p> This material is based upon work supported by the National Science Foundation under Grant No. 1115777 and 1440534.</p>
   <p>The computer hosting the web server was purchased through a UIC LAS Science award.</p></div></div>
   """

def Print_JS():
    # empty notification is not finished
    print"""<script type="text/javascript">
       function clearFocus() {
         if (document.textform.poly.value=='please enter the polynomial system') {
            document.textform.poly.value='';
         }
       }
       function showstartpoly(){
         if (document.getElementById("start_poly").style.display !="none"){
            document.getElementById("start_poly").style.display="none";
         }
         else{
            document.getElementById("start_poly").style.display="inline";
         }
       }
       </script>"""
       
def Print_Text_File(f):
   from phc_file import check_file_exist
   if check_file_exist(f):
      print """<div class="report_text">"""
      file = open(f,'r')
      cont =""
      for line in file:
         cont += line + '<br>'
      print cont
      print '</div>'
      file.close()

def Print_Area(poly,name):
   print """<div><textarea name=%s rows = 10 autofocus="autofocus">%s</textarea></div>"""% (name, poly)
   print "</td></tr></table>"

def Print_Area_File(f,name):
    print """<div><textarea name=%s rows = 10 autofocus="autofocus">"""% name
    if f=="empty":
         print "please enter the polynomial system"
    else:
       file = open(f,'r')
       cont =""
       for line in file:
          cont += line
       print cont
       file.close()
    print """</textarea></div>"""

def Print_Area_Eq(file_id, style = ''):
   from os import path
   error_poly_file = 0
   
   if path.isfile(file_id):
      file = open(file_id,'r')
      poly = ""
      end_mark = "THE SOLUTIONS :"
      for line in file:
         if line[:len(end_mark)] == end_mark:
            break
         else:
            poly += line
      file.close()
      error_poly_file = Print_Area_Preview(poly, "poly_preview", style)
   else:
      print '<p>Polynomial file does not exist.</p>'
      
   return error_poly_file
    

def Print_Area_Preview(poly, name, style = ''):
   if poly == "":
      print """<textarea name=%s rows = 10 autofocus="autofocus"  style="width:100%%">"""% name
      print "please enter the polynomial system"
      print """</textarea>"""
   else:
      for i in xrange(len(poly)):
         if poly[i] == '\n':
            dim = poly[:i]
            break;
      poly = poly[i:].split('\n')
      try:
         dim = int(dim)
      except:
         print "<font color='red'>Miss the Dimension of the polynomial system in the first line.</font>"
         return 1;
      print """<div %s><fieldset><legend>Dim= %d</legend><table>"""% (style,dim)
      cont = ""
      lastfuns =0      
      funs = 1
      lastchar = '\0'
      symbols = [' ',';','+','-','*','^','.','(',')', chr(13)]
      sup = 0 
      sub = 0
      find_all_eq = 0

      for line in poly:
         if find_all_eq == 1:
            break;
         if lastfuns < funs:
            cont += """<tr></tr><tr></tr><tr><th valign="top">f<sub>"""+ str(funs) +"""</sub></th><th valign="top"> = </th><td>"""
            lastfuns = funs
         len_line = len(line)
         for i in range(0,len_line):
            if line[i] == ';':
               sup = 0 
               sub = 0
               funs += 1
               if funs == dim+1:
                  find_all_eq = 1;
                  break;
               cont += "</td></tr>\n"
               lastchar = '\0'
            elif line[i] != ' ':
               if sub == 0: # subscript for index
                  if lastchar < '0' or lastchar > '9': ## not good
                     if line[i] >= '0' and line[i] <= '9':
                        if lastchar != '\0' and lastchar not in symbols:
                           cont += "<sub>"
                           sub = 1
               else:
                  if line[i] < '0' or line[i] > '9':
                     cont += "</sub>" 
                     sub = 0    
               if sup == 0:  # superscript for power
                  if line[i] =='*' and lastchar == '*':
                        cont += "<sup>"
                        sup = 1
                  elif line[i] == '^':
                     cont += "<sup>"
                     sup = 1
                  elif line[i] == '*':
                     if i < len_line-1 and line[i+1] !='*':
                        cont += "  "
                  elif line[i] == '+':
                     if lastchar != '\0' and lastchar != 'E' and  lastchar != 'e':
                        cont += "  +  "
                     else:
                        cont += '+'
                  elif line[i] == '-':
                     if lastchar != '\0' and lastchar != 'E' and  lastchar != 'e':
                        cont += "  -  "
                     else:
                        cont += "-"
                  elif line[i] == 'i':
                     cont += 'I'
                  else:
                     cont += line[i]
               else:
                  if line[i] < '0' or line[i] > '9':
                     cont += "</sup>"
                     sup = 0
                  if line[i] == '*':
                     cont += "  "
                  elif line[i] == '+':
                     if lastchar != '\0' and lastchar != 'E' and  lastchar != 'e':
                        cont += "  +  "
                     else:
                        cont += '+'
                  elif line[i] == '-':
                     if lastchar != '\0' and lastchar != 'E' and  lastchar != 'e':
                        cont += "  -  "
                     else:
                        cont += '-'
                  elif line[i] == 'i':
                     cont += 'I'
                  else:
                     cont += line[i]                        
               lastchar = line[i]
      if sup == 1:
         cont += "</sup>"
      elif sub == 1:
         cont += "</sub>"
      cont += "</th></tr>\n</table></fieldset></div>"
      print cont
      
   return 0

def Print_Notice_Top():
   print """<div class="notice"><table><tr>
      <td>Solve the polynomial system:</td>
      <td style="text-align:right" ><form action="phcweb.py" style="display: inline"><input type="submit" value = "New"></form></td>
      </tr></table></div>"""

def Print_Notice_Mid(status, Uid, name):
   print """<div class="notice">"""
   Print_Status(status, Uid, name)
   print """</div>"""

def Start_Poly_Button():
   print """<button type="button" onclick="showstartpoly()">On/Off</button>"""

def Print_Area_Report(f,name):
    from os import path
    if path.isfile(f):
       print """<p><textarea name=%s rows = 30 readonly="readonly">"""% name
       file = open(f,'r')
       cont =""
       for line in file:
          cont += line
       print cont
       file.close()
       print """</textarea></p>"""
    else:
       print"<p>Still computing mix volume. No report yet.</p>"

def Print_Sol_Page(filename, poly_id, sol_page):
    from phc_sol import read_sol_number_file
    from phc_config import n_sol_page
    n_sol = read_sol_number_file(filename)
    print """<div>
    <table>
    """
    if sol_page > 0:
        last_page_link = "phcweb.py?poly_id=%s&sol_page=%s"%(poly_id, sol_page-1)
        print """
        <td style="text-align:left;"><a href="%s"> Last Page </a></td>
        """%(last_page_link)
    if (sol_page+1)*n_sol_page < n_sol:
        next_page_link = "phcweb.py?poly_id=%s&sol_page=%s"%(poly_id, sol_page+1)
        print """
        <td  style="text-align:right;"><a href="%s"> Next Page </a></td>
        """%(next_page_link)
    print """
    </table>
    </div>"""

def Print_Area_Sol(filename, Uid, poly_id,name, sol_page):
    from phc_sol import read_sol_file
    Print_Sol_Page(filename, poly_id, sol_page)
    read_sol_file(filename, sol_page)
    Print_Sol_Page(filename, poly_id, sol_page)

def Print_Poly_New(Uid, poly, homotopy = 0):
   Print_Notice_Top()
   print """<form name="textform" method="post" action="phcweb.py">"""

   Print_Area(poly,'poly')

   if homotopy == 0:
      Print_Solve_Button(Uid);
   else:
      Print_Solve_Button(Uid, hom = poly_id) #XXX
      Print_Start_Eq(file_id)

   print "</form>"

# Print form.poly in the box with solve bottum

def Print_Poly_Preview(Uid, poly, homotopy = 0, poly_id=''):
   Print_Notice_Top()
   print """<form name="textform" method="post" action="phcweb.py">"""

   Print_Area_Preview(poly,"poly_preview")

   if homotopy == 0:
      Print_Solve_Button(Uid,poly);
   else:
      from phc_file import get_file_id
      file_id = get_file_id(Uid, poly_id)
      Print_Solve_Button(Uid, poly= poly, hom = poly_id)
      Print_Start_Eq(file_id)

   print "</form>"

# Print current.poly in the box with solve bottum

def Print_Poly_File_New(Uid, poly_id, homotopy = 0):
   from phc_file import get_file_id
   file_id = get_file_id(Uid, poly_id)
   Print_Notice_Top()
   print """<form name="textform" method="post" action="phcweb.py">"""

   Print_Area_File("%s.poly"%file_id,"poly")

   if homotopy == 0:
      Print_Solve_Button(Uid);
   else:
      Print_Solve_Button(Uid, hom = poly_id)
      Print_Start_Eq(file_id)

   print "</form>"

def Print_Start_Eq(file_id):
   poly_id = file_id.split('/')[-1]
   print "<p>Starting system of %s used in homotopy:"%poly_id
   Start_Poly_Button()
   print "</p>"
   Print_JS()
   error_start_poly = Print_Area_Eq(file_id+".start", 'id ="start_poly" style="display:none; width:100%%"')
   if error_start_poly:
      print "<p><font color = 'red'>The original polynomial doesn't have a starting system.</font></p>"
   else:
      print "<p>Important Notice: in order to use homotopy, you can only change the coefficients of the new system, but don't change any variable or term</p>"

# Print empty box with solve bottum

def Print_Poly_Empty(Uid):
   Print_JS()
   Print_Notice_Top()
   print """<form name="textform" method="post" action="phcweb.py">
   <div style="width:100%%">
   <textarea name="poly" rows = 10 style="width:100%%" onfocus="clearFocus();">please enter the polynomial system</textarea>
   </div>"""

   Print_Solve_Button(Uid)

   print """</form>"""

def Print_Solve_Button(Uid, poly = '', hom =''):
   print "<div>"
   if hom != '':
       print """<input type = "hidden" name = "homotopy" value =%s >"""% hom
       print """<input type = "hidden" name = "poly_id" value =%s >"""% hom

   print """<input type = "hidden" name = "Uid" value ="""+Uid+ ">"

   
   if poly != '':
      print """<input type = "hidden" name = "poly" value ="%s">"""% poly
      button_list = [("submit", "Edit"     ,"Edit"     ,1),\
                     ("submit", "Solve"    ,"Solve"       ,1)]
   else:
      if hom != '':
         button_list = [("submit", "Preview"  ,"Preview"     ,1),\
                        ("submit", "Solve"    ,"Solve by homotopy" ,1)]
      else:
         button_list = [("submit", "Preview"  ,"Preview"     ,1),\
                        ("submit", "Solve"    ,"Solve"       ,1)]
         
      
   Print_Button_Table(button_list)
   print "</div>"
   
def PHC_View_Bar(f,formt):
   option=["","",""]
   option[int(formt)]="selected"
   print """<form action="phcweb.py">
              <select name="m">
                <option value="0" %s >PHC Format</option>
                <option value="1" %s >Python Format</option>
                <option value="2" %s >Maple Format</option>
              </select>
              <input type = "hidden" name = "print" value = %s>
              <input type = "hidden" name = "r" value = "1">
              <input type= "submit" value = "View">
            </form>""" % (option[0], option[1], option[2],f)

def Print_Status(status, Uid, name):
   from os import path
   from phc_file import get_file_id
   file_id = get_file_id(Uid, name)
   print """<table><tr><td>Status:</td>"""
   if status == 0:
       print """<td><b>Killed</b></td></tr><td></td>
               <td>Check the following phc report before killed."""
       file_name = "../users/%s/%s" %(Uid, name+".phc")
       if path.isfile(file_name):
          print """Download PHC detailed report <a href="%s" target="_blank"> %s </a>."""%(file_name, name+".phc")
       print "</td></tr>"
   elif status == 1:
      print """<td><b>Computing</b></td></tr><tr><td></td><td>
               <p>Your polynomial has been submitted and is being solved.</p>
            Notice: You can only solve one polynomial each time. If you want to solve a new one, you have to kill current computation first."""
       
   elif status == 2:
      if name == "current":
         print """<td><b>Solved</b>"""
      else:
         print """<td><b>Saved as %s</b>"""% name

      print """</td></tr><tr><td></td>
               <td><div>Download solutions <a href="../users/%s/%s" target="_blank"> %s </a>. </div>
               <div>Download PHC detailed report <a href="../users/%s/%s" target="_blank"> %s </a>.</div>
               </td></tr>"""%(Uid,name+".sol",name+".sol",Uid,name+".phc",name+".phc")
   elif status == 9:
      print "<td><b>Error:</b></td></tr> <tr><td></td><td>The input polynomial has something wrong. Please check the following report:</td></tr><tr><td></td><td>"
      if path.isfile(file_id+".err"):
         f = open(file_id+".err", 'r')
         for line in f:
            if line == '\n':
               print "</br>"
            else:
               print line
         print "</td></tr>"
   else:
      print """<td>Status report does not exist.</td>"""
   print "</table>"

def Print_Button_Status(status, Uid, Pid, poly_id):
   # start form
   print """<form method="post" action="phcweb.py">
   <input type = "hidden" name = "Uid" value ="""+Uid+ """>
   <input type = "hidden" name = "poly_id" value ="""+ poly_id + """>"""

   button_list = []

   if status == 1:
      print """<input type = "hidden" name = "Pid" value ="""+ Pid + """>"""
      button_list = [("submit", "Kill"      ,"Kill"     ,1),\
#                     ("submit", "Update"    ,"Email Me" ,1),\
                     ("submit", "Update"    ,"Update"   ,0)]

   elif status == 0 or (status == 2 and poly_id == 'current'):
      button_list = [("submit", "Edit"    ,"Edit"         ,1),\
                     ("submit", "Save"    ,"Save as"      ,0),\
                     ("text"  , "newname" ,"New name"     ,1)]

   elif status == 2:
      button_list = [("submit", "Edit"    ,"Edit"         ,1),\
                     ("submit", "Save"    ,"Rename"       ,0),\
                     ("text"  , "newname" ,"New name"     ,1),\
                     ("submit", "homotopy","Similar",1)]

   elif status == 9:
      button_list = [("submit", "Edit"    ,"Edit"         ,1),\
                     ("submit", "Save"    ,"Rename"       ,0),\
                     ("text"  , "newname" ,"New name"     ,1)]

   if button_list != []:
      Print_Button_Table(button_list)
   else:
      print "<p>Error: Status = %s Status is not right</p>"%status
    
   # end form   
   print """</form>"""

def Print_Button_Table(button_list):

   n_buttons = len(button_list)-1
   button_width = 1.0/(n_buttons+1)*100

   cont = """<div style ="margin-top:10px"><table><tr>\n"""

   for i in xrange(n_buttons):

      #if button_list[i][3] == 0:
      #   cont += """<td>"""
      #else:
      if button_list[i-1][3] != 0:
         cont += """<td>"""
      cont += """<input type="%s" name="%s" value = "%s" size = "9">
      """%(button_list[i][0], button_list[i][1],button_list[i][2])

      if button_list[i][3]:
         cont += "</td>\n"

   if button_list[n_buttons-1][3]:
      cont += """<td style="text-align:right">"""
   cont += """<input type="%s" name="%s" value = "%s" size = "9">\n</td>\n"""%(button_list[n_buttons][0], button_list[n_buttons][1],button_list[n_buttons][2])
   
   cont += """</tr></table></div>"""
   print cont
