with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Dense_Series;             use Standard_Dense_Series;
with Standard_Dense_Series_io;          use Standard_Dense_Series_io;

package body Series_and_Polynomials is

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  The first variable in t
    --   is considered as the variable in the truncated power series.
    --   The other variables are copied to variable of the same index
    --   minus one in the term of the series polynomial.

      rtm : Standard_Series_Polynomials.Term;
      ord : constant integer32 := integer32(t.dg(t.dg'first));
      dim : constant integer32 := t.dg'last-1;
      rcf : Series := Create(0.0,ord);
      cnt : natural32 := 0;

    begin
      rcf.cff(ord) := t.cf;
      rtm.cf := rcf;
      rtm.dg := new Standard_Natural_Vectors.Vector(t.dg'first..dim);
      for i in rtm.dg'range loop
        rtm.dg(i) := t.dg(i+1);
      end loop;
      cnt := cnt + 1;
      if verbose then
        put("Adding term "); put(cnt,1); put_line(" with coefficient :");
        put(rtm.cf);
        put("order : "); put(ord,1);
        put(" and degrees : "); put(rtm.dg.all); new_line;
      end if;
      Standard_Series_Polynomials.Add(res,rtm);
      Standard_Series_Polynomials.Clear(rtm);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : Standard_Series_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Series_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : Series := t.cf;
      zero : constant Complex_Number := Create(0.0);
      rtpc : Complex_Number;
      dim1 : constant integer32 := t.dg'last+1;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of order "); put(cffs.order,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      for k in 0..cffs.order loop
        rtpc := cffs.cff(k);
        if not Equal(rtpc,zero) then
          declare
            rt : Standard_Complex_Polynomials.Term;
          begin
            rt.cf := rtpc;
            rt.dg := new Standard_Natural_Vectors.Vector(1..dim1);
            rt.dg(1) := natural32(k);
            for i in t.dg'range loop
              rt.dg(i+1) := t.dg(i);
            end loop;
            Standard_Complex_Polynomials.Add(res,rt);
            Standard_Complex_Polynomials.Clear(rt);
          end;
        end if;
      end loop;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;
 
end Series_and_Polynomials;
