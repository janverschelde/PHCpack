from phc_html import *

def PHC_Preview(Uid, form):
   from phc_html import Print_Area_Preview
   poly = form['poly'].value
   if form.has_key('homotopy') and form.has_key('poly_id'):
      Print_Poly_Preview(Uid, poly, 1, form['poly_id'].value)
   else:
      Print_Poly_Preview(Uid, poly)

