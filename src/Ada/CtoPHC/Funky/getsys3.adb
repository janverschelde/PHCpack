with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with C_Integer_io,C_Double_io;          use C_Integer_io,C_Double_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;

procedure getsys3 ( n,m : in integer32; mc : in C_intarrs.Pointer;
                    ns : in integer32; s : in C_intarrs.Pointer;
                    nc : in integer32; c : in C_dblarrs.Pointer ) is

-- DESCRIPTION :
--   The number of monomials in the i-th equation is mc[i], where
--   i runs from 0 to m-1.  The support vector s has range 0..ns-1
--   and the coefficient vector c has range 0..nc-1.
--   The number of variables in every equation is n.

-- ON ENTRY :
--   n         number of variables;
--   m         number of polynomials in the system;
--   mc        number of monomials in every polynomial, range 0..m-1;
--   ns        number of exponents;
--   s         support vector for the system, range 0..ns-1;
--   nc        number of coefficients;
--   c         coefficient vector for the system, range 0..nc-1.

  mva : C_Integer_Array(0..Interfaces.C.size_T(m-1))
      := C_intarrs.Value(mc,Interfaces.C.ptrdiff_T(m));
  sva : C_Integer_Array(0..Interfaces.C.size_T(ns-1))
      := C_intarrs.Value(s,Interfaces.C.ptrdiff_T(ns));
  cva : C_Double_Array(0..Interfaces.C.size_T(nc-1)) 
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_T(nc));
  p : Poly_Sys(1..m);

begin
  put("Number of variables : "); put(n,1); new_line;
  put("Number of equations : "); put(m,1); new_line;
  put("Number of monomials in every equation : ");
  for i in 0..m-1 loop
    put(" "); put(mva(Interfaces.C.size_T(i)),1);
  end loop;
  new_line;
  put("Number of exponents in support : "); put(ns,1); new_line;
  put("The support :");
  for i in 0..ns-1 loop
    put(" "); put(sva(Interfaces.C.size_T(i)),1);
  end loop;
  new_line;
  put("Number of doubles in coefficient vector : "); put(nc,1); new_line;
  put_line("The coefficients : ");
  for i in 0..nc-1 loop
    put(" "); put(cva(Interfaces.C.size_T(i)));
    new_line;
  end loop;
  p := Create(natural32(n),mva,cva,sva);
  put_line("The polynomial system defined by coefficients and support : ");
  put_line(p);
end getsys3;
