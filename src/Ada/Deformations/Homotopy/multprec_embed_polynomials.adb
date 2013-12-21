with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;         use Multprec_Complex_Numbers;
with Standard_Natural_Vectors;         use Standard_Natural_Vectors;

package body Multprec_Embed_Polynomials is

  function Add_Variables ( p : Poly; k : natural32 ) return Poly is

    res : Poly := Null_Poly;
    ek : constant integer32 := integer32(k);

    procedure Add_Variable_in_Term ( t : in Term; continue : out boolean ) is

      ext : Term;

    begin
      Copy(t.cf,ext.cf);
      ext.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+ek);
      ext.dg(t.dg'range) := t.dg.all;
      ext.dg(t.dg'last+1..t.dg'last+ek) := (t.dg'last+1..t.dg'last+ek => 0);
      Add(res,ext);
      Clear(ext);
      continue := true;
    end Add_Variable_in_Term;
    procedure Add_Variables_in_Terms is
      new Visiting_Iterator(Add_Variable_in_Term);

  begin
    Add_Variables_in_Terms(p);
    return res;
  end Add_Variables;

  function Add_Variables ( p : Poly_Sys; k : natural32 ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
    end loop;
    return res;
  end Add_Variables;

  function Add_Variables ( p : Jaco_Mat; k : natural32 ) return Jaco_Mat is

    res : Jaco_Mat(p'range(1),p'range(2));

  begin
    for i in p'range(1) loop
      for j in p'range(2) loop
        res(i,j) := Add_Variables(p(i,j),k);
      end loop;
    end loop;
    return res;
  end Add_Variables;

end Multprec_Embed_Polynomials;
