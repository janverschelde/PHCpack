with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with Witness_Sets;                      use Witness_Sets;
with Planes_and_Polynomials;            use Planes_and_Polynomials;

package body QuadDobl_Diagonal_Polynomials is

  function Create
             ( n,i : integer32 ) return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    res : Term;
 
  begin
    res.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    res.dg(i) := 1;
    res.cf := Create(integer(1));
    return res;
  end Create;

  function Create
             ( n,i : integer32 ) return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    res : Term;
 
  begin
    res.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    res.dg(i) := 1;
    res.cf := Create(integer(1));
    return res;
  end Create;

  function Create
             ( n,i : integer32 ) return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res_term : constant Term := Create(n,i);
    res : constant Poly := Create(res_term);
 
  begin
    return res;
  end Create;

  function Create
             ( n,i : integer32 ) return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res_term : constant Term := Create(n,i);
    res : constant Poly := Create(res_term);
 
  begin
    return res;
  end Create;

  function Insert_Variables
             ( n : integer32; t : QuadDobl_Complex_Polynomials.Term )
             return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..t.dg'last+n);
    for i in 1..n loop
      res.dg(i) := 0;
    end loop;
    for i in t.dg'range loop
      res.dg(i+n) := t.dg(i);
    end loop;
    return res;
  end Insert_Variables;

  function Insert_Variables
             ( n : integer32; t : QuadDobl_Complex_Laurentials.Term )
             return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(1..t.dg'last+n);
    for i in 1..n loop
      res.dg(i) := 0;
    end loop;
    for i in t.dg'range loop
      res.dg(i+n) := t.dg(i);
    end loop;
    return res;
  end Insert_Variables;

  function Insert_Variables
             ( n : integer32; p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Add_Inserted_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Insert_Variables(n,t); 

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Add_Inserted_Term;
    procedure Add_Inserted_Terms is new Visiting_Iterator(Add_Inserted_Term);

  begin
    Add_Inserted_Terms(p);
    return res;
  end Insert_Variables;

  function Insert_Variables
             ( n : integer32; p : QuadDobl_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Add_Inserted_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Insert_Variables(n,t); 

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Add_Inserted_Term;
    procedure Add_Inserted_Terms is new Visiting_Iterator(Add_Inserted_Term);

  begin
    Add_Inserted_Terms(p);
    return res;
  end Insert_Variables;

  function Insert_Variables ( n : integer32; p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Insert_Variables(n,p(i));
    end loop;
    return res;
  end Insert_Variables;

  function Insert_Variables ( n : integer32; p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Insert_Variables(n,p(i));
    end loop;
    return res;
  end Insert_Variables;

  function Append_Variables
             ( n : integer32; t : QuadDobl_Complex_Polynomials.Term )
             return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..t.dg'last+n);
    res.dg(t.dg'range) := t.dg.all;
    for i in 1..n loop
      res.dg(t.dg'last+i) := 0;
    end loop;
    return res;
  end Append_Variables;

  function Append_Variables
             ( n : integer32; t : QuadDobl_Complex_Laurentials.Term )
             return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(1..t.dg'last+n);
    res.dg(t.dg'range) := t.dg.all;
    for i in 1..n loop
      res.dg(t.dg'last+i) := 0;
    end loop;
    return res;
  end Append_Variables;

  function Append_Variables
             ( n : integer32; p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Add_Appended_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Append_Variables(n,t);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Add_Appended_Term;
    procedure Add_Appended_Terms is new Visiting_Iterator(Add_Appended_Term);

  begin
    Add_Appended_Terms(p);
    return res;
  end Append_Variables;

  function Append_Variables
             ( n : integer32; p : QuadDobl_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Add_Appended_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Append_Variables(n,t);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Add_Appended_Term;
    procedure Add_Appended_Terms is new Visiting_Iterator(Add_Appended_Term);

  begin
    Add_Appended_Terms(p);
    return res;
  end Append_Variables;

  function Append_Variables ( n : integer32; p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Append_Variables(n,p(i));
    end loop;
    return res;
  end Append_Variables;

  function Append_Variables ( n : integer32; p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Append_Variables(n,p(i));
    end loop;
    return res;
  end Append_Variables;

  function Is_Slack
             ( t : QuadDobl_Complex_Polynomials.Term; n : integer32 )
             return boolean is

  -- DESCRIPTION :
  --   Returns true if there are any variables with indices > n.

  begin
    for i in n+1..t.dg'last loop
      if t.dg(i) > 0
       then return true;
      end if;
    end loop;
    return false;
  end Is_Slack;

  function Is_Slack
             ( t : QuadDobl_Complex_Laurentials.Term; n : integer32 )
             return boolean is

  -- DESCRIPTION :
  --   Returns true if there are any variables with indices > n.

  begin
    for i in n+1..t.dg'last loop
      if t.dg(i) > 0
       then return true;
      end if;
    end loop;
    return false;
  end Is_Slack;

  function Truncate
             ( t : QuadDobl_Complex_Polynomials.Term; n : integer32 )
             return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..n);
    res.dg(1..n) := t.dg(1..n);
    return res;
  end Truncate;

  function Truncate
             ( t : QuadDobl_Complex_Laurentials.Term; n : integer32 )
             return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(1..n);
    res.dg(1..n) := t.dg(1..n);
    return res;
  end Truncate;

  function Truncate
             ( p : QuadDobl_Complex_Polynomials.Poly; n : integer32 )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Truncate_Term ( t : in Term; continue : out boolean ) is

      tt : Term;

    begin
      if not Is_Slack(t,n) then
        tt := Truncate(t,n);
        Add(res,tt);
        Clear(tt);
      end if;
      continue := true;
    end Truncate_Term;
    procedure Truncate_Terms is new Visiting_Iterator(Truncate_Term);

  begin
    Truncate_Terms(p);
    return res;
  end Truncate;

  function Truncate
             ( p : QuadDobl_Complex_Laurentials.Poly; n : integer32 )
             return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Truncate_Term ( t : in Term; continue : out boolean ) is

      tt : Term;

    begin
      if not Is_Slack(t,n) then
        tt := Truncate(t,n);
        Add(res,tt);
        Clear(tt);
      end if;
      continue := true;
    end Truncate_Term;
    procedure Truncate_Terms is new Visiting_Iterator(Truncate_Term);

  begin
    Truncate_Terms(p);
    return res;
  end Truncate;

  function Collapse
             ( t : QuadDobl_Complex_Polynomials.Term; n : integer32 )
             return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    rt : Term;
    done : boolean := false;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      if t.dg(i) /= 0 then
        done := true;
        rt.dg(i) := t.dg(i);
      end if;
    end loop;
    if not done
     then rt.dg(1..n) := t.dg(n+1..2*n);
    end if;
    return rt;
  end Collapse;

  function Collapse
             ( t : QuadDobl_Complex_Laurentials.Term; n : integer32 )
             return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    rt : Term;
    done : boolean := false;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      if t.dg(i) /= 0 then
        done := true;
        rt.dg(i) := t.dg(i);
      end if;
    end loop;
    if not done
     then rt.dg(1..n) := t.dg(n+1..2*n);
    end if;
    return rt;
  end Collapse;

  function Collapse
             ( p : QuadDobl_Complex_Polynomials.Poly; n : integer32 )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Collapse_Term ( t : in Term; continue : out boolean ) is

      ct : Term;

    begin
      if not Is_Slack(t,2*n) then
        ct := Collapse(t,n);
        Add(res,ct);
        Clear(ct);
      end if;
      continue := true;
    end Collapse_Term;
    procedure Collapse_Terms is new Visiting_Iterator(Collapse_Term);

  begin
    Collapse_Terms(p);
    return res;
  end Collapse;

  function Collapse
             ( p : QuadDobl_Complex_Laurentials.Poly; n : integer32 )
             return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Collapse_Term ( t : in Term; continue : out boolean ) is

      ct : Term;

    begin
      if not Is_Slack(t,2*n) then
        ct := Collapse(t,n);
        Add(res,ct);
        Clear(ct);
      end if;
      continue := true;
    end Collapse_Term;
    procedure Collapse_Terms is new Visiting_Iterator(Collapse_Term);

  begin
    Collapse_Terms(p);
    return res;
  end Collapse;

  function Collapse ( p : Poly_Sys; n : integer32 ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Collapse(p(i),n);
    end loop;
    return res;
  end Collapse;

  function Collapse ( p : Laur_Sys; n : integer32 ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Collapse(p(i),n);
    end loop;
    return res;
  end Collapse;

  function Collapse
             ( t : QuadDobl_Complex_Polynomials.Term;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    rt : Term;
    first : boolean := false;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      if t.dg(i) /= 0 then
        first := true;   -- monomial belongs to first group of variables
        rt.dg(i) := t.dg(i);
      end if;
    end loop;
    if not first then
      for i in 1..n loop -- place variables of second group
        rt.dg(q(i)) := t.dg(i+n);
      end loop;
     end if;
   -- for i in 2*n+1..t.dg'last loop  -- copy embedding variables
   --   rt.dg(i-n) := t.dg(i);
   -- end loop;
    return rt;
  end Collapse;

  function Collapse
             ( t : QuadDobl_Complex_Laurentials.Term;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    rt : Term;
    first : boolean := false;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      if t.dg(i) /= 0 then
        first := true;   -- monomial belongs to first group of variables
        rt.dg(i) := t.dg(i);
      end if;
    end loop;
    if not first then
      for i in 1..n loop -- place variables of second group
        rt.dg(q(i)) := t.dg(i+n);
      end loop;
     end if;
   -- for i in 2*n+1..t.dg'last loop  -- copy embedding variables
   --   rt.dg(i-n) := t.dg(i);
   -- end loop;
    return rt;
  end Collapse;

  function Collapse
             ( p : QuadDobl_Complex_Polynomials.Poly;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Collapse_Term ( t : in Term; continue : out boolean ) is

      ct : Term;

    begin
      if not Is_Slack(t,2*n) then
        ct := Collapse(t,n,q);
        Add(res,ct);
        Clear(ct);
      end if;
      continue := true;
    end Collapse_Term;
    procedure Collapse_Terms is new Visiting_Iterator(Collapse_Term);

  begin
    Collapse_Terms(p);
    return res;
  end Collapse;

  function Collapse
             ( p : QuadDobl_Complex_Laurentials.Poly;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Collapse_Term ( t : in Term; continue : out boolean ) is

      ct : Term;

    begin
      if not Is_Slack(t,2*n) then
        ct := Collapse(t,n,q);
        Add(res,ct);
        Clear(ct);
      end if;
      continue := true;
    end Collapse_Term;
    procedure Collapse_Terms is new Visiting_Iterator(Collapse_Term);

  begin
    Collapse_Terms(p);
    return res;
  end Collapse;

  function Collapse ( p : Poly_Sys; n : integer32; q : Permutation )
                    return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Collapse(p(i),n,q);
    end loop;
    return res;
  end Collapse;

  function Collapse ( p : Laur_Sys; n : integer32; q : Permutation )
                    return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Collapse(p(i),n,q);
    end loop;
    return res;
  end Collapse;

  function Diagonal ( n : integer32 ) return Poly_Sys is

    use QuadDobl_Complex_Polynomials;

    res : Poly_Sys(1..n);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2*n => 0);
    t.cf := Create(integer(1));
    for i in 1..n loop
      t.dg(i) := 1;
      res(i) := Create(t);
      t.dg(i) := 0;
      t.dg(n+i) := 1;
      Sub(res(i),t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Diagonal;

  function Diagonal ( n : integer32 ) return Laur_Sys is

    use QuadDobl_Complex_Laurentials;

    res : Laur_Sys(1..n);
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..2*n => 0);
    t.cf := Create(integer(1));
    for i in 1..n loop
      t.dg(i) := 1;
      res(i) := Create(t);
      t.dg(i) := 0;
      t.dg(n+i) := 1;
      Sub(res(i),t);
      t.dg(n+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Diagonal;

  function Product ( n1,n2 : integer32; p1,p2 : Poly_Sys ) return Poly_Sys is

    n : constant integer32 := p1'length + p2'length;
    res : Poly_Sys(1..n);
    ind : integer32 := 0;

  begin
    for i in p1'range loop
      ind := ind + 1;
      res(ind) := Append_Variables(n2,p1(i));
    end loop;
    for i in p2'range loop
      ind := ind + 1;
      res(ind) := Insert_Variables(n1,p2(i));
    end loop;
    return res;
  end Product;

  function Product ( n1,n2 : integer32; p1,p2 : Laur_Sys ) return Laur_Sys is

    n : constant integer32 := p1'length + p2'length;
    res : Laur_Sys(1..n);
    ind : integer32 := 0;

  begin
    for i in p1'range loop
      ind := ind + 1;
      res(ind) := Append_Variables(n2,p1(i));
    end loop;
    for i in p2'range loop
      ind := ind + 1;
      res(ind) := Insert_Variables(n1,p2(i));
    end loop;
    return res;
  end Product;

  function Product ( n,k : integer32; hyp1,hyp2 : VecVec ) return VecVec is

    res : VecVec(1..k);
    n2 : constant integer32 := 2*n;

  begin
    for i in 1..k loop
      declare
        v : QuadDobl_Complex_Vectors.Vector(0..n2);
      begin
        v(0) := hyp1(k)(0) + hyp2(k)(0);
        for j in 1..n loop
          v(j) := hyp1(k)(j);
          v(j+n) := hyp2(k)(j);
        end loop;
        res(i) := new QuadDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Product;

end QuadDobl_Diagonal_Polynomials;
