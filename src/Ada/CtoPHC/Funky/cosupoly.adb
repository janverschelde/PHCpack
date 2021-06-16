with text_io;                          use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with C_Integer_io,C_Double_io;         use C_Integer_io,C_Double_io;
with C_Integer_Arrays;                 use C_Integer_Arrays;
with C_Double_Arrays;                  use C_Double_Arrays;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;  use Standard_Complex_Polynomials_io;
with Coefficient_Support_Polynomials;  use Coefficient_Support_Polynomials;

procedure cosupoly ( n : in integer32; s : in C_intarrs.Pointer;
                     m : in integer32; c : in C_dblarrs.Pointer ) is

-- DESCRIPTION :
--   The support vector s has range 0..n-1 and the coefficient
--   vector c has range 0..m-1.

  sva : C_Integer_Array(0..Interfaces.C.size_T(n-1))
      := C_intarrs.Value(s,Interfaces.C.ptrdiff_T(n));
  cva : C_Double_Array(0..Interfaces.C.size_T(m-1)) 
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_T(m));
  numvars : constant natural32 := natural32(n)/(natural32(m)/2);
  p : Poly := Create(numvars,cva,sva);

begin
  put("The n is "); put(n,1);
  put(" and m is "); put(m,1);
  put(" number of variables : "); put(numvars,1); new_line;
 -- put("sva'first : "); put(integer(sva'first),1); new_line;
 -- put("sva'last : "); put(integer(sva'last),1); new_line;
 -- put("cva'first : "); put(integer(cva'first),1); new_line;
 -- put("cva'last : "); put(integer(cva'last),1); new_line;
  put("The support :");
  for i in 0..n-1 loop
    put(" "); put(sva(Interfaces.C.size_T(i)),1);
  end loop;
  new_line;
  put_line("The coefficients : ");
  for i in 0..m-1 loop
    put(" "); put(cva(Interfaces.C.size_T(i)));
    new_line;
  end loop;
  put_line("The polynomial defined by coefficients and support : ");
  put_line(p);
end cosupoly;
