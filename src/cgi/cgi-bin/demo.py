#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug, style_demo
if phc_debug:
   import cgitb; cgitb.enable()
   
import cgi
from phc_html import PrintHeader, Footer_Bar

from sql_login import phc_login_simple

def main():
   """
   About
   """
   
   PrintHeader(style_demo)

   phc_login_simple()
   
   print """<p><h2>A brief introduction to PHCWeb</h2>"""
   
   print """<div class='notice'>
<p>PHCpack Web Interface intends to create an intuitive user interface to solve polynomial systems.</p>"""
   
   print """<div style="text-align:center;"><img src="../imag/Path.jpg" alt="PHC Web Interface" width="400"></div>"""
      
   print"""<p><b>Benefits: </b></p><ul>
<li>No installation required</li>
<li>We provide a computational server, hosted by fast computers</li>
<li>Manage the solved polynomial systems</li>
<li>Solve polynomial systems anywhere from cell phone to tablet</li>
</ul>"""
   
   print"""<p><b>Main functions of PHCpack Web Interface:</b></p>
   <ul>
<li><a href="#solve">How to solve a polynomial system by PHCWeb?</a></li>
<li><a href="#mypoly">How to manage my polynomial systems in PHCWeb?</a></li>
</ul>"""

   print """</div>"""
   
   
   print """<h2><a name="solve" id="solve">How to solve a polynomial system by PHCWeb?</a></h2>"""
   
   print"""<ol>
<li style="font-size:120%">Enter your polynomial system in the textarea.</li>"""
   
   print """<div style="text-align:center; margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/poly.png" alt="PHC Web Interface"></div>"""
   
   print """<li style="font-size:120%">Click "Solve" button to get solutions. If the solution is large, it will take some time and you might need to click update to check the newest result.</li>"""
   
   print """<div style="text-align:center;margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/Solution.png" alt="PHC Web Interface"></div>"""
   
   print """<li style="font-size:120%">Give your system a name. Click on "Save as" to save it.</li>"""
   
   print """<div style="text-align:center;margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/saveas.png" alt="PHC Web Interface"></div>"""
   
   print """<li style="font-size:120%">If you have a similar solution with different efficients. Click on "Solve Similar", change the coefficents and use homotopy to solve your new system.</li>"""
   
   print """<div style="text-align:center;margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/startsystem.png" alt="PHC Web Interface"></div></ol>"""
   
   print """<h2><a name="mypoly" id="mypoly">How to manage my polynomial systems in PHCWeb?</a></h2>"""
   
   print"""<ol>
<li style="font-size:120%">Click "My Polynomial Systems" on the top.</li>"""
   print """<div style="text-align:center;margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/topbar.png" alt="PHC Web Interface"></div>"""
   print """<li style="font-size:120%">You have entire table of all polynomials you have solved. You can check solutions and phc report. You can also transfer solution to maple or python format.</li>"""
   print """<div style="text-align:center;margin-top:15px; margin-bottom: 25px"><img src="../imag/Intro/mypoly.png" alt="PHC Web Interface"></div></ol>"""
   
   print """</ol>"""
   
   
   print """</ol>"""
   
   Footer_Bar()
   
   print "</body></html>\n"

main()
