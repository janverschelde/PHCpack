with Interfaces.C;                      use Interfaces.C;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;

package body Coefficient_Support_Polynomials is

  function Support ( p : Poly ) return C_Integer_Array is

    m : constant natural32 := Number_of_Terms(p);
    n : constant natural32 := Number_of_Unknowns(p);
    numexp : constant size_T := size_T(n*m-1);
    res : C_Integer_Array(0..numexp);
    ind : Size_T := 0;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      for i in t.dg'range loop
        res(ind) := Interfaces.C.int(t.dg(i));
        ind := ind + 1;
      end loop;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Support;

  function Coefficients ( p : Poly ) return C_Double_Array is

    m : constant natural32 := Number_of_Terms(p);
    numcff : constant size_T := size_T(2*m-1);
    res : C_Double_Array(0..numcff);
    ind : Size_T := 0;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      res(ind) := Interfaces.C.double(REAL_PART(t.cf));
      ind := ind + 1;
      res(ind) := Interfaces.C.double(IMAG_PART(t.cf));
      ind := ind + 1;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Coefficients;

  function Create ( n : natural32; c : C_Double_Array; s : C_Integer_Array ) 
                  return Poly is

    res : Poly := Null_Poly;
    indcff : size_T := 0;
    indsup : size_T := 0;
    t : Term;

  begin
   -- put("c'first : "); put(integer(c'first),1); new_line;
   -- put("c'last : "); put(integer(c'last),1); new_line;
   -- put("s'first : "); put(integer(s'first),1); new_line;
   -- put("s'last : "); put(integer(s'last),1); new_line;
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while (indcff < c'last) and (indsup <= s'last) loop
      t.cf := Create(double_float(c(indcff)),double_float(c(indcff+1)));
      indcff := indcff + 2;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(s(indsup));
        indsup := indsup + 1;
        exit when (indsup > s'last);
      end loop;
      Add(res,t);
    end loop;
    Clear(t);
    return res;
  end Create;

end Coefficient_Support_Polynomials;
