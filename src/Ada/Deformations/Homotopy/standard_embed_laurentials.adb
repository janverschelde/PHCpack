with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Vectors;         use Standard_Integer_Vectors;

package body Standard_Embed_Laurentials is

  function Add_Variables ( p : Poly; k : natural32 ) return Poly is

    res : Poly := Null_Poly;
    ek : constant integer32 := integer32(k);

    procedure Add_Variable_in_Term ( t : in Term; continue : out boolean ) is

      ext : Term;

    begin
      ext.cf := t.cf;
      ext.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last+ek);
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

  function Add_Variables ( p : Laur_Sys; k : natural32 ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
    end loop;
    return res;
  end Add_Variables;

end Standard_Embed_Laurentials;
