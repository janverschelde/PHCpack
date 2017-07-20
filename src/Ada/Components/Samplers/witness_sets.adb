with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with DoblDobl_Random_Matrices;           use DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;           use QuadDobl_Random_Matrices;
with Standard_Complex_Substitutors;      use Standard_Complex_Substitutors;
with Standard_Embed_Polynomials;         use Standard_Embed_Polynomials;
with Standard_Embed_Laurentials;         use Standard_Embed_Laurentials;
with DoblDobl_Embed_Polynomials;         use DoblDobl_Embed_Polynomials;
with DoblDobl_Embed_Laurentials;         use DoblDobl_Embed_Laurentials;
with QuadDobl_Embed_Polynomials;         use QuadDobl_Embed_Polynomials;
with QuadDobl_Embed_Laurentials;         use QuadDobl_Embed_Laurentials;
with Planes_and_Polynomials;             use Planes_and_Polynomials;

package body Witness_Sets is

  function Random_Hyperplanes ( k,n : natural32 )
                              return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Vectors.Vector'
                      (Random_Vector(0,integer32(n)));
    end loop;
    return res;
  end Random_Hyperplanes;

  function Random_Hyperplanes ( k,n : natural32 )
                              return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'
                      (Random_Vector(0,integer32(n)));
    end loop;
    return res;
  end Random_Hyperplanes;

  function Random_Hyperplanes ( k,n : natural32 )
                              return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'
                      (Random_Vector(0,integer32(n)));
    end loop;
    return res;
  end Random_Hyperplanes;

  function Random_Hyperplanes ( k,n,size : natural32 )
                              return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (Random_Vector(0,integer32(n),size));
    end loop;
    return res;
  end Random_Hyperplanes;

-- MANIPULATION OF POLYNOMIALS :

  function Add_Dummy ( n,k,i : natural32 )
                     return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(1.0);
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummy ( n,k,i : natural32 )
                     return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Numbers,Standard_Complex_Laurentials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(1.0);
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummy ( n,k,i : natural32 )
                     return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(integer(1));
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummy ( n,k,i : natural32 )
                     return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(integer(1));
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummy ( n,k,i : natural32 )
                     return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(integer(1));
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummy ( n,k,i : natural32 )
                     return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n+k) => 0);
    t.cf := Create(integer(1));
    t.dg(integer32(n+i)) := 1;
    Add(res,t);
    Clear(t);
    return res;
  end Add_Dummy;

  function Add_Dummies
             ( n,k : natural32 ) return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n+k) => 0);
    for i in 1..integer32(k) loop
      t.cf := Standard_Random_Numbers.Random1;
      t.dg(integer32(n)+i) := 1;
      Add(res,t);
      t.dg(integer32(n)+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Dummies;

  function Add_Embedding ( p : Standard_Complex_Polynomials.Poly;
                           k : natural32 )
                         return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := Standard_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : Standard_Complex_Laurentials.Poly;
                           k : natural32 )
                         return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := Standard_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : DoblDobl_Complex_Polynomials.Poly;
                           k : natural32 )
                         return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := DoblDobl_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : DoblDobl_Complex_Laurentials.Poly;
                           k : natural32 )
                         return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := DoblDobl_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : QuadDobl_Complex_Polynomials.Poly;
                           k : natural32 )
                         return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := QuadDobl_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : QuadDobl_Complex_Laurentials.Poly;
                           k : natural32 )
                         return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Add_Variables(p,k);
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in 1..integer32(k) loop
      t.cf := QuadDobl_Random_Numbers.Random1;
      t.dg(n+i) := 1;
      Add(res,t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                           k : natural32 )
                         return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Add_Embedding(p(i),k);
    end loop;
    return res;
  end Add_Embedding;

  function Add_Embedding ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                           k : natural32 )
                         return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Add_Embedding(p(i),k);
    end loop;
    return res;
  end Add_Embedding;

-- OPERATIONS ON POLYNOMIAL SYSTEMS :

  function Make_Square ( f : Standard_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(k));
    r : Standard_Complex_Numbers.Complex_Number;

  begin
    if f'last <= integer32(k) then
      for i in f'range loop
        Copy(f(i),res(i));
      end loop;
    else
      for i in 1..integer32(k) loop
        Copy(f(i),res(i));
      end loop;
      for i in integer32(k)+1..f'last loop
        for j in res'range loop
          r := Standard_Random_Numbers.Random1;
          declare
            rp : Poly := r*f(i);
          begin
            Add(res(j),rp);
            Clear(rp);
          end;
        end loop;
      end loop;
    end if;
    return res;
  end Make_Square;

  function Make_Square ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(k));
    r : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if f'last <= integer32(k) then
      for i in f'range loop
        Copy(f(i),res(i));
      end loop;
    else
      for i in 1..integer32(k) loop
        Copy(f(i),res(i));
      end loop;
      for i in integer32(k)+1..f'last loop
        for j in res'range loop
          r := DoblDobl_Random_Numbers.Random1;
          declare
            rp : Poly := r*f(i);
          begin
            Add(res(j),rp);
            Clear(rp);
          end;
        end loop;
      end loop;
    end if;
    return res;
  end Make_Square;

  function Make_Square ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(k));
    r : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if f'last <= integer32(k) then
      for i in f'range loop
        Copy(f(i),res(i));
      end loop;
    else
      for i in 1..integer32(k) loop
        Copy(f(i),res(i));
      end loop;
      for i in integer32(k)+1..f'last loop
        for j in res'range loop
          r := QuadDobl_Random_Numbers.Random1;
          declare
            rp : Poly := r*f(i);
          begin
            Add(res(j),rp);
            Clear(rp);
          end;
        end loop;
      end loop;
    end if;
    return res;
  end Make_Square;

  function Add_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                       hyp : Standard_Complex_Vectors.Vector )
                     return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    for i in p'range loop
      Standard_Complex_Polynomials.Copy(p(i),res(i));
    end loop;
    res(res'last) := Hyperplane(hyp);
    return res;
  end Add_Slice;

  function Add_Slice ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                       hyp : DoblDobl_Complex_Vectors.Vector )
                     return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    for i in p'range loop
      DoblDobl_Complex_Polynomials.Copy(p(i),res(i));
    end loop;
    res(res'last) := Hyperplane(hyp);
    return res;
  end Add_Slice;

  function Add_Slice ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                       hyp : QuadDobl_Complex_Vectors.Vector )
                     return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    for i in p'range loop
      QuadDobl_Complex_Polynomials.Copy(p(i),res(i));
    end loop;
    res(res'last) := Hyperplane(hyp);
    return res;
  end Add_Slice;

  function Remove_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                        return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Natural_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := Standard_Complex_Numbers.Create(1.0);
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Remove_Slice ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                        return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laurentials;
    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Integer_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := Standard_Complex_Numbers.Create(1.0);
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Remove_Slice ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Natural_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := DoblDobl_Complex_Numbers.Create(integer(1));
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Remove_Slice ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laurentials;
    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Integer_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := DoblDobl_Complex_Numbers.Create(integer(1));
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Remove_Slice ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Natural_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := QuadDobl_Complex_Numbers.Create(integer(1));
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Remove_Slice ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laurentials;
    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    t : Term;

  begin
    for i in p'first..p'last-1 loop
      Copy(p(i),res(i));
    end loop;
    t.dg := new Standard_Integer_Vectors.Vector'(p'range => 0);
    t.dg(t.dg'last) := 1;
    t.cf := QuadDobl_Complex_Numbers.Create(integer(1));
    res(res'last) := Create(t);
    Clear(t.dg);
    return res;
  end Remove_Slice;

  function Eliminate_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                             k,i : natural32 )
                           return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
    n : constant integer32
      := integer32(Number_of_Unknowns(p(integer32(k))));
    s,srv : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.dg(integer32(i)) := 1;
    t.cf := Coeff(p(integer32(k)),t.dg);
    s := p(integer32(k))-t;
    Min(s);
    srv := Remove_Variable(s,integer32(i));
    for j in res'range loop
      if j < integer32(k)
       then res(j) := Substitute(integer32(i),srv,p(j)); 
       else res(j) := Substitute(integer32(i),srv,p(j+1));
      end if;
    end loop;
    Clear(t); Clear(s); Clear(srv);
    return res;
  end Eliminate_Slice;

  function Embed ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                 return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                 return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                 return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                 return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                 return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                 return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
    end loop;
    return res;
  end Embed;

  function Embed ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                   gamma : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.dg(n+1) := 1;
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
      t.cf := gamma(i);
      Add(res(i),t);
    end loop;
    Clear(t);
    return res;
  end Embed;

  function Embed ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                   gamma : DoblDobl_Complex_Vectors.Vector )
                 return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.dg(n+1) := 1;
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
      t.cf := gamma(i);
      Add(res(i),t);
    end loop;
    Clear(t);
    return res;
  end Embed;

  function Embed ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                   gamma : QuadDobl_Complex_Vectors.Vector )
                 return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.dg(n+1) := 1;
    for i in res'range loop
      res(i) := Add_Variables(p(i),1);
      t.cf := gamma(i);
      Add(res(i),t);
    end loop;
    Clear(t);
    return res;
  end Embed;

  function Slice_and_Embed1
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               hyp : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Adds one random hyperplane to the system and does an embedding
  --   with one new variable.

    sli : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Add_Slice(p,hyp); 
    gam : constant Standard_Complex_Vectors.Vector(p'first..p'last+1)
        := Random_Vector(1,p'last+1);
    res : constant Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Embed(sli,gam);

  begin
    Standard_Complex_Polynomials.Clear(sli(sli'last));
    return res;
  end Slice_and_Embed1;

  function Slice_and_Embed1
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               hyp : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Adds one random hyperplane to the system and does an embedding
  --   with one new variable.

    sli : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Add_Slice(p,hyp); 
    gam : constant DoblDobl_Complex_Vectors.Vector(p'first..p'last+1)
        := Random_Vector(1,p'last+1);
    res : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Embed(sli,gam);

  begin
    DoblDobl_Complex_Polynomials.Clear(sli(sli'last));
    return res;
  end Slice_and_Embed1;

  function Slice_and_Embed1
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               hyp : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Adds one random hyperplane to the system and does an embedding
  --   with one new variable.

    sli : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Add_Slice(p,hyp); 
    gam : constant QuadDobl_Complex_Vectors.Vector(p'first..p'last+1)
        := Random_Vector(1,p'last+1);
    res : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1)
        := Embed(sli,gam);

  begin
    QuadDobl_Complex_Polynomials.Clear(sli(sli'last));
    return res;
  end Slice_and_Embed1;

  procedure Store_Random_Hyperplanes
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 ) is

    use Standard_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : Standard_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  procedure Store_Random_Hyperplanes
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 ) is

    use Standard_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : Standard_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  procedure Store_Random_Hyperplanes
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 ) is

    use DoblDobl_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : DoblDobl_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  procedure Store_Random_Hyperplanes
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 ) is

    use DoblDobl_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : DoblDobl_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  procedure Store_Random_Hyperplanes
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 ) is

    use QuadDobl_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : QuadDobl_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  procedure Store_Random_Hyperplanes
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 ) is

    use QuadDobl_Complex_Matrices;

    dim : constant integer32 := integer32(n);
    nbr : constant integer32 := integer32(k);
    lastidx : constant integer32 := dim+nbr;
    cff : Matrix(1..lastidx+1,1..nbr);

  begin
    if k = 1
     then cff := Random_Matrix(n+k+1,k);
     else cff := Random_Orthogonal_Matrix(n+k+1,k);
    end if;
    for i in 1..nbr loop
      declare
        hyp : QuadDobl_Complex_Vectors.Vector(0..lastidx);
      begin
        for j in hyp'range loop
          hyp(j) := cff(j+1,i);
        end loop;
        p(dim+i) := Hyperplane(hyp);
      end;
    end loop;
  end Store_Random_Hyperplanes;

  function Slice_and_Embed
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               k : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;

    dim : constant integer32 := p'last + integer32(k);
    res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := Standard_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Slice_and_Embed
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               k : natural32 )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laurentials;

    dim : constant integer32 := p'last + integer32(k);
    res : Standard_Complex_Laur_Systems.Laur_Sys(p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := Standard_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Slice_and_Embed
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               k : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;

    dim : constant integer32 := p'last + integer32(k);
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys (p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := DoblDobl_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Slice_and_Embed
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               k : natural32 )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laurentials;

    dim : constant integer32 := p'last + integer32(k);
    res : DoblDobl_Complex_Laur_Systems.Laur_Sys (p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := DoblDobl_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Slice_and_Embed
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               k : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;

    dim : constant integer32 := p'last + integer32(k);
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := QuadDobl_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Slice_and_Embed
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               k : natural32 )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laurentials;

    dim : constant integer32 := p'last + integer32(k);
    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'first..dim);
    n : constant integer32 := p'length;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(k) => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),k);
      for j in n+1..n+integer32(k) loop
        t.cf := QuadDobl_Random_Numbers.Random1;
        t.dg(j) := 1;
        Add(res(i),t);
        t.dg(j) := 0;
      end loop;
    end loop;
    Clear(t);
    Store_Random_Hyperplanes(res,natural32(n),k);
    return res;
  end Slice_and_Embed;

  function Embed_with_Dummies
              ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    nvars : constant integer32 := p'length+integer32(k);
    res : Standard_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nvars => 0);
    t.cf := Standard_Complex_Numbers.Create(1.0);
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Embed_with_Dummies
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laurentials;
    nvars : constant integer32 := p'length+integer32(k);
    res : Standard_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..nvars => 0);
    t.cf := Standard_Complex_Numbers.Create(1.0);
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Embed_with_Dummies
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    nvars : constant integer32 := p'length+integer32(k);
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nvars => 0);
    t.cf := DoblDobl_Complex_Numbers.Create(integer(1));
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Embed_with_Dummies
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laurentials;
    nvars : constant integer32 := p'length+integer32(k);
    res : DoblDobl_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..nvars => 0);
    t.cf := DoblDobl_Complex_Numbers.Create(integer(1));
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Embed_with_Dummies
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    nvars : constant integer32 := p'length+integer32(k);
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nvars => 0);
    t.cf := QuadDobl_Complex_Numbers.Create(integer(1));
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Embed_with_Dummies
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laurentials;
    nvars : constant integer32 := p'length+integer32(k);
    res : QuadDobl_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last+integer32(k))
        := Slice_and_Embed(p,k);
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..nvars => 0);
    t.cf := QuadDobl_Complex_Numbers.Create(integer(1));
    for i in 0..integer32(k)-1 loop
      t.dg(nvars-i) := 1;
      Clear(res(p'last-i));
      res(p'last-i) := Create(t);
      t.dg(nvars-i) := 0;
    end loop;
    return res;
  end Embed_with_Dummies;

  function Slice_and_Embed
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : natural32 ) 
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys is

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;

    res : Array_of_Poly_Sys(0..integer32(k));
    slices : Matrix(1..p'last+integer32(k)+1,1..integer32(k));

  begin
    if k = 1
     then slices := Random_Matrix(natural32(p'last)+k+1,1);
     else slices := Random_Orthogonal_Matrix(natural32(p'last)+k+1,k);
    end if;
    res(0) := new Poly_Sys'(p);
    for i in 1..integer32(k) loop
      declare
        hyp : Standard_Complex_Vectors.Vector(0..res(i-1)'last);
      begin
        for j in hyp'range loop
          hyp(j) := slices(j+1,i);
        end loop;
        res(i) := new Poly_Sys'(Slice_and_Embed1(res(i-1).all,hyp));
      end;
    end loop;
    return res;
  end Slice_and_Embed;

  function Copy_Last_and_Embed
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Copies the last equation and embeds everything.

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);
    gam : constant Standard_Complex_Vectors.Vector(p'first..p'last+1)
        := Random_Vector(1,p'last+1);

  begin
    res(p'range) := p;
    res(res'last) := p(p'last);
    return Embed(res,gam);
  end Copy_Last_and_Embed;

  function Slice_and_Embed
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               k,l : natural32 )
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..integer32(k));
    use Standard_Complex_Poly_Systems;

  begin
    res(0) := new Poly_Sys'(p);
    for i in 1..integer32(l) loop
      declare
        sys : constant Poly_Sys(p'first..p'last+i)
            := Copy_Last_and_Embed(res(i-1).all);
      begin
        res(i) := new Poly_Sys'(sys);
      end;
    end loop;
    for i in integer32(l)+1..integer32(k) loop
      declare
        hyp : constant Standard_Complex_Vectors.Vector(0..res(i-1)'last)
            := Random_Vector(0,res(i-1)'last);
      begin
        res(i) := new Poly_Sys'(Slice_and_Embed1(res(i-1).all,hyp));
      end;
    end loop;
    return res;
  end Slice_and_Embed;

  function Slices ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Slices ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new Standard_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Slices ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Slices ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Slices ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Slices ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..integer32(k));

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'
                      (Polynomial(p(p'last-integer32(k)+i)));
    end loop;
    return res;
  end Slices;

  function Square ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials,Standard_Complex_Poly_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Poly_Sys(1..m);
        hyp : Standard_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Poly_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := Standard_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Square ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laurentials,Standard_Complex_Laur_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Laur_Sys(1..m);
        hyp : Standard_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Laur_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := Standard_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Square ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials,DoblDobl_Complex_Poly_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Poly_Sys(1..m);
        hyp : DoblDobl_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Poly_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := DoblDobl_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Square ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laurentials,DoblDobl_Complex_Laur_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Laur_Sys(1..m);
        hyp : DoblDobl_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Laur_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := DoblDobl_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Square ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials,QuadDobl_Complex_Poly_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Poly_Sys(1..m);
        hyp : QuadDobl_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Poly_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := QuadDobl_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Square ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laurentials,QuadDobl_Complex_Laur_Systems;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if n < m then
      declare
        res : Laur_Sys(1..m);
        hyp : QuadDobl_Complex_Vectors.Vector(0..m);
      begin
        res(1..n) := p;
        for i in 1..m-n loop
          hyp := Random_Vector(0,m);
          res(n+i) := Hyperplane(hyp);
        end loop;
        return res;
      end;
    elsif n > m then
      declare
        res : Laur_Sys(1..n);
        k : constant natural32 := natural32(n-m);
        t : Term;
      begin
        t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
        for i in 1..n loop
          res(i) := Add_Variables(p(i),k);
          for j in 1..integer32(k) loop
            t.dg(m+j) := 1; 
            t.cf := QuadDobl_Random_Numbers.Random1;
            Add(res(i),t);
            t.dg(m+j) := 0;
          end loop;
        end loop;
        Clear(t);
        return res;
      end;
    else
      return p;
    end if;
  end Square;

  function Is_Added_Term
             ( t : Standard_Complex_Polynomials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Is_Added_Term
             ( t : Standard_Complex_Laurentials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Is_Added_Term
             ( t : DoblDobl_Complex_Polynomials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Is_Added_Term
             ( t : DoblDobl_Complex_Laurentials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Is_Added_Term
             ( t : QuadDobl_Complex_Polynomials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Is_Added_Term
             ( t : QuadDobl_Complex_Laurentials.Term; dim : natural32 )
             return boolean is

  -- DESCRIPTION :
  --   This routine tells whether the term was added in the embedding,
  --   i.e.: whether it has only nonzero degrees in the last dim positions.

  begin
    for i in t.dg'first..t.dg'last-integer32(dim) loop
      if t.dg(i) /= 0
       then return false;        -- term has nonzero original variable
      end if;
    end loop;                    -- term has no nonzero original variable
    for i in t.dg'last-integer32(dim)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;         -- moreover, term has added variable
      end if;
    end loop;
    return false;                -- term is constant term
  end Is_Added_Term;

  function Remove_Embedding ( p : Standard_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Natural_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Remove_Embedding ( p : Standard_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Integer_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Remove_Embedding ( p : DoblDobl_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Natural_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Remove_Embedding ( p : DoblDobl_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Integer_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Remove_Embedding ( p : QuadDobl_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Natural_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Remove_Embedding ( p : QuadDobl_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,dim) then
        rt.cf := t.cf;
        rt.dg := new Standard_Integer_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(dim)));
        Add(res,rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Embedding;

  function Number_of_Zero_Equations
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return natural32 is

    res : natural32 := 0;
    use Standard_Complex_Polynomials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Number_of_Zero_Equations
             ( p : Standard_Complex_Laur_Systems.Laur_Sys ) return natural32 is

    res : natural32 := 0;
    use Standard_Complex_Laurentials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Number_of_Zero_Equations
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return natural32 is

    res : natural32 := 0;
    use DoblDobl_Complex_Polynomials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Number_of_Zero_Equations
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys ) return natural32 is

    res : natural32 := 0;
    use DoblDobl_Complex_Laurentials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Number_of_Zero_Equations
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural32 is

    res : natural32 := 0;
    use QuadDobl_Complex_Polynomials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Number_of_Zero_Equations
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys ) return natural32 is

    res : natural32 := 0;
    use QuadDobl_Complex_Laurentials;

  begin
    while p(p'last-integer32(res)) = Null_Poly loop
      res := res + 1;
    end loop;
    return res;
  end Number_of_Zero_Equations;

  function Remove_Embedding1 ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Remove_Embedding1 ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Remove_Embedding1 ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Remove_Embedding1 ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Remove_Embedding1 ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Remove_Embedding1 ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys
            (p'first..p'last-integer32(dim));

  begin
    for i in res'range loop
      res(i) := Remove_Embedding(p(i),dim);
    end loop;
    return res;
  end Remove_Embedding1;

  function Complete ( n,k : natural32;
                      p : Standard_Complex_Poly_Systems.Poly_Sys )
                    return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n-k));
    use Standard_Complex_Polynomials;

  begin
    if p'last = integer32(n-k) then
      Standard_Complex_Poly_Systems.Copy(p,res);
    else
      declare
        rnd : Standard_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := Standard_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

  function Complete ( n,k : natural32;
                      p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                    return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n-k));
    use DoblDobl_Complex_Polynomials;

  begin
    if p'last = integer32(n-k) then
      DoblDobl_Complex_Poly_Systems.Copy(p,res);
    else
      declare
        rnd : DoblDobl_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := DoblDobl_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

  function Complete ( n,k : natural32;
                      p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                    return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n-k));
    use QuadDobl_Complex_Polynomials;

  begin
    if p'last = integer32(n-k) then
      QuadDobl_Complex_Poly_Systems.Copy(p,res);
    else
      declare
        rnd : QuadDobl_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := QuadDobl_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

  function Complete ( n,k : natural32;
                      p : Standard_Complex_Laur_Systems.Laur_Sys )
                    return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n-k));
    use Standard_Complex_Laurentials;

  begin
    if p'last = integer32(n-k) then
      Standard_Complex_Laur_Systems.Copy(p,res);
    else
      declare
        rnd : Standard_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := Standard_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

  function Complete ( n,k : natural32;
                      p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                    return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n-k));
    use DoblDobl_Complex_Laurentials;

  begin
    if p'last = integer32(n-k) then
      DoblDobl_Complex_Laur_Systems.Copy(p,res);
    else
      declare
        rnd : DoblDobl_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := DoblDobl_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

  function Complete ( n,k : natural32;
                      p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                    return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n-k));
    use QuadDobl_Complex_Laurentials;

  begin
    if p'last = integer32(n-k) then
      QuadDobl_Complex_Laur_Systems.Copy(p,res);
    else
      declare
        rnd : QuadDobl_Complex_Numbers.Complex_Number;
        extra : Poly;
      begin
        for i in 1..integer32(n-k) loop
          Copy(p(i),res(i));
        end loop;
        for i in integer32(n-k)+1..p'last loop
          for j in res'range loop
            rnd := QuadDobl_Random_Numbers.Random1;
            extra := rnd*p(i);
            Add(res(j),extra);
            Clear(extra);
          end loop;
        end loop;
      end;
    end if;
    return res;
  end Complete;

-- OPERATIONS ON SOLUTION LISTS :

  function Add_Component
             ( s : Standard_Complex_Solutions.Solution;
               c : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    res : Solution(s.n+1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'range) := s.v;
    res.v(res.n) := c;
    return res;
  end Add_Component;
 
  function Add_Component
             ( sols : Standard_Complex_Solutions.Solution_List;
               c : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Add_Component(Head_Of(tmp).all,c));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Add_Component;

  function Add_Embedding
             ( s : Standard_Complex_Solutions.Solution; k : natural32 )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    res : Solution(s.n+integer32(k));

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'range) := s.v;
    for i in s.n+1..res.n loop
      res.v(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    return res;
  end Add_Embedding;

  function Add_Embedding
             ( s : DoblDobl_Complex_Solutions.Solution; k : natural32 )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    res : Solution(s.n+integer32(k));

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'range) := s.v;
    for i in s.n+1..res.n loop
      res.v(i) := DoblDobl_Complex_Numbers.Create(integer(0));
    end loop;
    return res;
  end Add_Embedding;

  function Add_Embedding
             ( s : QuadDobl_Complex_Solutions.Solution; k : natural32 )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    res : Solution(s.n+integer32(k));

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'range) := s.v;
    for i in s.n+1..res.n loop
      res.v(i) := QuadDobl_Complex_Numbers.Create(integer(0));
    end loop;
    return res;
  end Add_Embedding;
 
  function Add_Embedding
             ( sols : Standard_Complex_Solutions.Solution_List;
               k : natural32 )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Add_Embedding(Head_Of(tmp).all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Add_Embedding;
 
  function Add_Embedding
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Add_Embedding(Head_Of(tmp).all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Add_Embedding;
 
  function Add_Embedding
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Add_Embedding(Head_Of(tmp).all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Add_Embedding;
 
  function Remove_Component
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    res : Solution(s.n-1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'first..s.v'last-1) := s.v(s.v'first..s.v'last-1);
    return res;
  end Remove_Component;
 
  function Remove_Component
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    res : Solution(s.n-1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'first..s.v'last-1) := s.v(s.v'first..s.v'last-1);
    return res;
  end Remove_Component;
 
  function Remove_Component
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    res : Solution(s.n-1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(s.v'first..s.v'last-1) := s.v(s.v'first..s.v'last-1);
    return res;
  end Remove_Component;

  procedure Remove_Component
              ( sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        remsol : constant Solution(ls.n-1) := Remove_Component(ls.all);
      begin
        Clear(ls);
        ls := new Solution'(remsol);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Remove_Component;
 
  function Remove_Component
             ( sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Remove_Component(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Component;
 
  function Remove_Component
             ( sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Remove_Component(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Component;
 
  function Remove_Component
             ( sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      Append(res,res_last,Remove_Component(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Component;

  function Remove_Embedding
             ( s : Standard_Complex_Solutions.Solution;
               k : natural32 )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    res : Solution(s.n-integer32(k));

  begin
    res.m := s.m;
    res.t := s.t;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1..s.n-integer32(k)) := s.v(1..s.n-integer32(k));
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( s : DoblDobl_Complex_Solutions.Solution;
               k : natural32 )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    res : Solution(s.n-integer32(k));

  begin
    res.m := s.m;
    res.t := s.t;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1..s.n-integer32(k)) := s.v(1..s.n-integer32(k));
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( s : QuadDobl_Complex_Solutions.Solution;
               k : natural32 )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    res : Solution(s.n-integer32(k));

  begin
    res.m := s.m;
    res.t := s.t;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1..s.n-integer32(k)) := s.v(1..s.n-integer32(k));
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( sols : Standard_Complex_Solutions.Solution_List;
               k : natural32 )
             return Standard_Complex_Solutions.Solution_List is
   
    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Remove_Embedding(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List is
   
    use DoblDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Remove_Embedding(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List is
   
    use QuadDobl_Complex_Solutions;
    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Remove_Embedding(ls.all,k));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Embedding;
 
end Witness_Sets;
