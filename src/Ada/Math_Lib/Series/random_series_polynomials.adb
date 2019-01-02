with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Random_Series;
with DoblDobl_Complex_Random_Series;
with QuadDobl_Complex_Random_Series;

package body Random_Series_Polynomials is

  function Standard_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return Standard_CSeries_Polynomials.Term is

    res : Standard_CSeries_Polynomials.Term;
    deg : constant integer32 := integer32(degcff);
    pos,rnd : integer32;

  begin
    res.cf := Standard_Complex_Random_Series.Random_Series(deg);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(nbvar) => 0);
    for i in 1..integer32(degtrm) loop
      pos := Standard_Random_Numbers.Random(1,integer32(nbvar));
      rnd := Standard_Random_Numbers.Random(0,1);
      res.dg(pos) := res.dg(pos) + natural32(rnd);
    end loop;
    return res;
  end Standard_Random_Term;

  function DoblDobl_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return DoblDobl_CSeries_Polynomials.Term is

    res : DoblDobl_CSeries_Polynomials.Term;
    deg : constant integer32 := integer32(degcff);
    pos,rnd : integer32;

  begin
    res.cf := DoblDobl_Complex_Random_Series.Random_Series(deg);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(nbvar) => 0);
    for i in 1..integer32(degtrm) loop
      pos := Standard_Random_Numbers.Random(1,integer32(nbvar));
      rnd := Standard_Random_Numbers.Random(0,1);
      res.dg(pos) := res.dg(pos) + natural32(rnd);
    end loop;
    return res;
  end DoblDobl_Random_Term;

  function QuadDobl_Random_Term
             ( nbvar,degtrm,degcff : natural32 )
             return QuadDobl_CSeries_Polynomials.Term is

    res : QuadDobl_CSeries_Polynomials.Term;
    deg : constant integer32 := integer32(degcff);
    pos,rnd : integer32;

  begin
    res.cf := QuadDobl_Complex_Random_Series.Random_Series(deg);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(nbvar) => 0);
    for i in 1..integer32(degtrm) loop
      pos := Standard_Random_Numbers.Random(1,integer32(nbvar));
      rnd := Standard_Random_Numbers.Random(0,1);
      res.dg(pos) := res.dg(pos) + natural32(rnd);
    end loop;
    return res;
  end QuadDobl_Random_Term;

  function Standard_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return Standard_CSeries_Polynomials.Poly is

    use Standard_CSeries_Polynomials;

    res : Poly := Null_Poly;

  begin
    for i in 1..nbterms loop
      declare
        ranterm : Term := Standard_Random_Term(nbvar,degpol,degcff);
      begin
        Add(res,ranterm);
        Clear(ranterm);
      end;
    end loop;
    return res;
  end Standard_Random_Polynomial;

  function DoblDobl_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return DoblDobl_CSeries_Polynomials.Poly is

    use DoblDobl_CSeries_Polynomials;

    res : Poly := Null_Poly;

  begin
    for i in 1..nbterms loop
      declare
        ranterm : Term := DoblDobl_Random_Term(nbvar,degpol,degcff);
      begin
        Add(res,ranterm);
        Clear(ranterm);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Polynomial;

  function QuadDobl_Random_Polynomial
             ( nbvar,nbterms,degpol,degcff : natural32 )
             return QuadDobl_CSeries_Polynomials.Poly is

    use QuadDobl_CSeries_Polynomials;

    res : Poly := Null_Poly;

  begin
    for i in 1..nbterms loop
      declare
        ranterm : Term := QuadDobl_Random_Term(nbvar,degpol,degcff);
      begin
        Add(res,ranterm);
        Clear(ranterm);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Polynomial;

  function Standard_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(1..integer32(nbequ));

  begin
    for i in res'range loop
      res(i) := Standard_Random_Polynomial(nbvar,nbterms,degpol,degcff);
    end loop;
    return res;
  end Standard_Random_System;

  function DoblDobl_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..integer32(nbequ));

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Polynomial(nbvar,nbterms,degpol,degcff);
    end loop;
    return res;
  end DoblDobl_Random_System;

  function QuadDobl_Random_System
             ( nbequ,nbvar,nbterms,degpol,degcff : natural32 )
             return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..integer32(nbequ));

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Polynomial(nbvar,nbterms,degpol,degcff);
    end loop;
    return res;
  end QuadDobl_Random_System;

end Random_Series_Polynomials;
