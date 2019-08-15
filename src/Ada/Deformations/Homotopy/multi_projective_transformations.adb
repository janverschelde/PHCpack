with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Degrees_in_Sets_of_Unknowns;

package body Multi_Projective_Transformations is

  function Make_Partition
             ( n,m : natural32; p : Standard_Natural_Vectors.Vector )
             return Partition is

    res : Partition(1..m);

  begin
    for i in res'range loop
      res(i) := Sets_of_Unknowns.Create(n); -- initialize each set
    end loop;
    for i in p'range loop
      Sets_of_Unknowns.Add(res(p(i)),natural32(i));
    end loop;
    return res;
  end Make_Partition;

  function Multiset_Degrees
             ( p : in Standard_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in DoblDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in QuadDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Make_Homogeneous
             ( t : Standard_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : DoblDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : QuadDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in Standard_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : Standard_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DoblDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return DoblDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in QuadDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return QuadDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      QuadDobl_Complex_Polynomials.Add(res,rt);
      QuadDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Standard_Random_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Random_Linear_Term;

  function DoblDobl_Random_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;

  begin
    res.cf := DoblDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DoblDobl_Random_Linear_Term;

  function QuadDobl_Random_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;

  begin
    res.cf := QuadDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end QuadDobl_Random_Linear_Term;

  function Standard_Start_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Complex_Numbers.Create(1.0);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Start_Linear_Term;

  function DoblDobl_Start_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;
    one : constant double_double := create(1.0);

  begin
    res.cf := DoblDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DoblDobl_Start_Linear_Term;

  function QuadDobl_Start_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    one : constant quad_double := create(1.0);

  begin
    res.cf := QuadDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end QuadDobl_Start_Linear_Term;

  function Standard_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : Standard_Complex_Polynomials.Term
            := Standard_Random_Linear_Term(n,i);
        begin
          Standard_Complex_Polynomials.Add(res,t);
          Standard_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end Standard_Random_Linear_Polynomial;

  function DoblDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : DoblDobl_Complex_Polynomials.Term
            := DoblDobl_Random_Linear_Term(n,i);
        begin
          DoblDobl_Complex_Polynomials.Add(res,t);
          DoblDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end DoblDobl_Random_Linear_Polynomial;

  function QuadDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : QuadDobl_Complex_Polynomials.Term
            := QuadDobl_Random_Linear_Term(n,i);
        begin
          QuadDobl_Complex_Polynomials.Add(res,t);
          QuadDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end QuadDobl_Random_Linear_Polynomial;

  function Standard_Start_Linear_Polynomial
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Poly is

    trm : Standard_Complex_Polynomials.Term
        := Standard_Start_Linear_Term(n,i);
    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    Standard_Complex_Polynomials.Sub(res,trm);
    Standard_Complex_Polynomials.Clear(trm);
    return res;
  end Standard_Start_Linear_Polynomial;

  function DoblDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Poly is

    trm : DoblDobl_Complex_Polynomials.Term
        := DoblDobl_Start_Linear_Term(n,i);
    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    DoblDobl_Complex_Polynomials.Sub(res,trm);
    DoblDobl_Complex_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Start_Linear_Polynomial;

  function QuadDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Poly is

    trm : QuadDobl_Complex_Polynomials.Term
        := QuadDobl_Start_Linear_Term(n,i);
    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    QuadDobl_Complex_Polynomials.Sub(res,trm);
    QuadDobl_Complex_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Start_Linear_Polynomial;

  function Standard_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : Standard_Complex_Polynomials.Term;
    ztm : Standard_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := Standard_Random_Linear_Polynomial(dim,z(i));
      cst.cf := Standard_Random_Numbers.Random1;
      ztm.cf := Standard_Random_Numbers.Random1;
      Standard_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      Standard_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    Standard_Complex_Polynomials.Clear(cst);
    Standard_Complex_Polynomials.Clear(ztm);
    return res;
  end Standard_Random_Linear_Polynomials;

  function DoblDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : DoblDobl_Complex_Polynomials.Term;
    ztm : DoblDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := DoblDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := DoblDobl_Random_Numbers.Random1;
      ztm.cf := DoblDobl_Random_Numbers.Random1;
      DoblDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      DoblDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    DoblDobl_Complex_Polynomials.Clear(cst);
    DoblDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end DoblDobl_Random_Linear_Polynomials;

  function QuadDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : QuadDobl_Complex_Polynomials.Term;
    ztm : QuadDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := QuadDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := QuadDobl_Random_Numbers.Random1;
      ztm.cf := QuadDobl_Random_Numbers.Random1;
      QuadDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      QuadDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    QuadDobl_Complex_Polynomials.Clear(cst);
    QuadDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end QuadDobl_Random_Linear_Polynomials;

  function Standard_Start_Linear_Polynomials
             ( n,m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := Standard_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end Standard_Start_Linear_Polynomials;

  function DoblDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := DoblDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end DoblDobl_Start_Linear_Polynomials;

  function QuadDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := QuadDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end QuadDobl_Start_Linear_Polynomials;

  function Add_Ones ( s : Standard_Complex_Solutions.Solution;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : Standard_Complex_Solutions.Solution(dim+integer32(m));

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := Standard_Complex_Numbers.Create(1.0);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : DoblDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : DoblDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant double_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := DoblDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : QuadDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : QuadDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant quad_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := QuadDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : Standard_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution_List is

    res,res_last : Standard_Complex_Solutions.Solution_List;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : DoblDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution_List is

    res,res_last : DoblDobl_Complex_Solutions.Solution_List;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : QuadDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution_List is

    res,res_last : QuadDobl_Complex_Solutions.Solution_List;
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  procedure Add_Ones ( sols : in out Standard_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant Standard_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        Standard_Complex_Solutions.Clear(ls);
        ls := new Standard_Complex_Solutions.Solution'(sol);
        Standard_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant DoblDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        DoblDobl_Complex_Solutions.Clear(ls);
        ls := new DoblDobl_Complex_Solutions.Solution'(sol);
        DoblDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant QuadDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        QuadDobl_Complex_Solutions.Clear(ls);
        ls := new QuadDobl_Complex_Solutions.Solution'(sol);
        QuadDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

end Multi_Projective_Transformations;
