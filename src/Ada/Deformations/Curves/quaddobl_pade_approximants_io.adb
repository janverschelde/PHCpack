with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Poly_Strings;

package body QuadDobl_Pade_Approximants_io is

  function to_Poly ( c : Vector ) return Poly is

    res : Poly := Null_Poly;
    trm : Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..1);
    for k in c'range loop
      trm.dg(1) := natural32(k);
      trm.cf := c(k);
      Add(res,trm);
    end loop;
    Clear(trm);
    return res;
  end to_Poly;

  function Write ( c : Vector; s : Symbol ) return string is

    sa : Array_of_Symbols(1..1);
    p : Poly := to_Poly(c);

  begin
    sa(1) := s;
    declare
      res : constant string
          := QuadDobl_Complex_Poly_Strings.Write(p,sa);
    begin
      QuadDobl_Complex_Polynomials.Clear(p);
      return res;
    end;
  end Write;

  function Write ( p : Pade; s : Symbol ) return string is

    numcff : constant QuadDobl_Complex_Vectors.Vector
           := Numerator_Coefficients(p);
    numstr : constant string := Write(numcff,s);
    dencff : constant QuadDobl_Complex_Vectors.Vector
           := Denominator_Coefficients(p);
    denstr : constant string := Write(dencff,s);
    res : constant string
        := "(" & numstr(1..numstr'last-1) & ")/("
               & denstr(1..denstr'last-1) & ")";

  begin
    return res;
  end Write;

end QuadDobl_Pade_Approximants_io;
