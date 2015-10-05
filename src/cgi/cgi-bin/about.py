#!/usr/bin/python

# Enter debug mode
from phc_config import phc_debug, style_about
if phc_debug:
   import cgitb; cgitb.enable()
   
import cgi
from phc_html import PrintHeader, Footer_Bar

from sql_login import phc_login_simple

def print_about():
   
   print """<div><p>This is a beta version of a web interface to PHCpack,
developed by Xiangcheng Yu.</p>

<p>Its current version exports the blackbox solver (phc -b) and
tracking paths defined by a homotopy in one parameter (phc -p).</p>

<p>The blackbox solver computes the mixed volume of the Newton polytopes
of the system, via an integrated version of Algorithm 846: MixedVol:
A software package for mixed volume computation, by T. Gao, T.Y. Li,
and M. Wu, published in ACM Trans. Math. Softw. 31(4):555-560, 2005.</p>
   
<div style="text-align:center;"><img src="../imag/Path.jpg" alt="PHC Web Interface" width="400"></div>
<p>The web interface logo results from running HomLab in Matlab:

The Numerical Solution of Systems of Polynomials Arising in Engineering and Science, by Andrew J. Sommese and Charles W. Wampler, II, World Scientific, 2005.</p>

<div>"""

def main():
   """
   About
   """
   # PHCpack Web Interface Header
   
   PrintHeader(style_about)
   
   phc_login_simple()
   
   print_about()
   
   Footer_Bar()
   
   print "</body></html>\n"

main()
