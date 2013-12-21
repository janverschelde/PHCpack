-- io added for debugging
--with text_io,integer_io; use text_io,integer_io;

with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Witness_Sets;                      use Witness_Sets;
with Planes_and_Polynomials;            use Planes_and_Polynomials;

package body Extrinsic_Diagonal_Homotopies is

-- OPERATIONS ON POLYNOMIALS :

  function Create ( n,i : integer32 ) return Term is

    res : Term;
 
  begin
    res.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    res.dg(i) := 1;
    res.cf := Create(1.0);
    return res;
  end Create;

  function Create ( n,i : integer32 ) return Poly is

    res_term : constant Term := Create(n,i);
    res : constant Poly := Create(res_term);
 
  begin
    return res;
  end Create;

  function Insert_Variables ( n : integer32; t : Term ) return Term is

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

  function Insert_Variables ( n : integer32; p : Poly ) return Poly is

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

  function Append_Variables ( n : integer32; t : Term ) return Term is

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

  function Append_Variables ( n : integer32; p : Poly ) return Poly is

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

  function Is_Slack ( t : Term; n : integer32 ) return boolean is

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

  function Truncate ( t : Term; n : integer32 ) return Term is

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..n);
    res.dg(1..n) := t.dg(1..n);
    return res;
  end Truncate;

  function Truncate ( p : Poly; n : integer32 ) return Poly is

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

  function Collapse ( t : Term; n : integer32 ) return Term is

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

  function Collapse ( p : Poly; n : integer32 ) return Poly is

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

  function Collapse ( t : Term; n : integer32; q : Permutation ) return Term is

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

  function Collapse ( p : Poly; n : integer32; q : Permutation ) return Poly is

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

  function Diagonal ( n : integer32 ) return Poly_Sys is

    res : Poly_Sys(1..n);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2*n => 0);
    t.cf := Create(1.0);
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

  function Product ( n,k : integer32; hyp1,hyp2 : VecVec ) return VecVec is

    res : VecVec(1..k);
    n2 : constant integer32 := 2*n;

  begin
    for i in 1..k loop
      declare
        v : Standard_Complex_Vectors.Vector(0..n2);
      begin
        v(0) := hyp1(k)(0) + hyp2(k)(0);
        for j in 1..n loop
          v(j) := hyp1(k)(j);
          v(j+n) := hyp2(k)(j);
        end loop;
        res(i) := new Standard_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Product;

-- OPERATIONS ON SOLUTIONS :

  function Product ( s1,s2 : Solution ) return Solution is

    res : Solution(s1.n+s2.n);

  begin
    res.m := s1.m;
    res.t := s1.t;
    res.err := s1.err;
    res.rco := s1.rco;
    res.res := s1.res;
    res.v(s1.v'range) := s1.v;
    res.v(s1.v'last+1..res.v'last) := s2.v;
    return res;
  end Product;

  function Product ( s1,s2 : Solution_List ) return Solution_List is

    res,res_last,ptr1,ptr2 : Solution_List;

  begin
    ptr1 := s1;
    while not Is_Null(ptr1) loop
      ptr2 := s2;
      while not Is_Null(ptr2) loop
        Append(res,res_last,Product(Head_Of(ptr1).all,Head_Of(ptr2).all));
        ptr2 := Tail_Of(ptr2);
      end loop;
      ptr1 := Tail_Of(ptr1);
    end loop;
    return res;
  end Product;

  function Truncate ( s : Solution; n : integer32 ) return Solution is

    res : Solution(n);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := s.v(1..n);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Truncate;

  function Truncate ( s : Solution_List; n : integer32 ) return Solution_List is

    res,res_last,tmp : Solution_List;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      Append(res,res_last,Truncate(Head_Of(tmp).all,n));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Truncate;

-- CREATION OF THE CASCADE OF HOMOTOPIES :

  function Cascade_Dimension ( n1,n2,a,b : natural32 ) return natural32 is

    k : constant natural32 := n1-a;
    res : natural32;

  begin
    if a + b >= k
     then res := 3*k - a;
     else res := 2*k + b;
    end if;
    return res;
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : Poly_Sys; a,b : natural32 ) return natural32 is

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  procedure Cascade1 ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                       start,target : out Poly_Sys ) is

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    cdi : constant Poly_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade2 ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                       start,target : out Poly_Sys ) is

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : Standard_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(0.0);
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                start,target : out Poly_Sys ) is

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  function Extrinsic_Product
               ( a,b : natural32; s1,s2 : Solution ) return Solution is

    k : constant natural32 := natural32(s1.n);
    s : constant Solution(s1.n+s2.n) := Product(s1,s2);

  begin
    if a+b < k
     then return Add_Embedding(s,b);
     else return Add_Embedding(s,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b,k : natural32; sols1,sols2 : Solution_List )
               return Solution_List is

    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k
     then return Add_Embedding(sols,b);
     else return Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Product;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                sols1,sols2 : in Solution_List;
                start,target : out Poly_Sys; esols : out Solution_List ) is

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

end Extrinsic_Diagonal_Homotopies;
