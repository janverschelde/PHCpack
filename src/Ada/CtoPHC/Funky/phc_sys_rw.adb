with system;
with text_io,integer_io;                use text_io,integer_io;
with Interfaces.C;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Integer_Arrays_io;               use C_Integer_Arrays_io;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;

function phc_sys_rw ( rw,size_p : integer; p : C_DblArrs.Pointer )
                    return C_DblArrs.Pointer is

  res : C_DblArrs.Pointer;

  procedure Coefficient_Support_Representation ( s : in Poly_Sys ) is

    n : constant natural := Number_of_Unknowns(s(s'first));
    mon : constant C_Integer_Array := Monomial_Count(s);
    moncnt : constant natural := Sum(mon);
    sup : constant C_Integer_Array := Support(n,moncnt,mon,s);
    cff : constant C_Double_Array := Coefficients(moncnt,mon,s);
    cct : C_Double_Array(0..Interfaces.C.size_T(size_p-1))
        := C_DblArrs.Value(p,Interfaces.C.ptrdiff_T(size_p));

  begin
    cct := Concat(n,mon,cff,sup);
    res := cct(0)'unchecked_access;
  end Coefficient_Support_Representation;

  procedure Write_Polynomial_System is

    cct : constant C_Double_Array
        := C_DblArrs.Value(p,Interfaces.C.ptrdiff_T(size_p));
    size : constant natural := integer(cct(cct'first));
    x : C_Double_Array(0..Interfaces.C.size_t(size))
      := cct(0..Interfaces.C.size_t(size));

  begin
    if size > size_p
     then put_line("DANGER: MORE ELEMENTS THAN SIZE OF BUFFER !!!");
          put(size,1); put(" > "); put(size_p,1); new_line;
    end if;
    declare
      n : constant natural := Dimension(x);
      m : constant C_Integer_Array := Monomial_Count(cct);
      s : constant C_Integer_Array := Support(x);
      c : constant C_Double_Array := Coefficients(x);
      ps : constant Poly_Sys := Create(n,m,c,s);
      k : natural := integer(m'last)-1;
    begin
      put(ps);
    end;
  end Write_Polynomial_System;

begin
  if rw = 0
   then new_line;
        declare
          lp : Link_to_Poly_Sys;
        begin
          get(lp);
          Coefficient_Support_Representation(lp.all);
        end;
   else Write_Polynomial_System;
  end if;
  return res;
end phc_sys_rw;
