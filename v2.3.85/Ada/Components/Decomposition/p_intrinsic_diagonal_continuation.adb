with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Increment_and_Fix_Continuation;    use Increment_and_Fix_Continuation;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Affine_Sampling_Machine;           use Affine_Sampling_Machine;

package body P_Intrinsic_Diagonal_Continuation is

  function Create ( n : natural; p,q : Poly ) return Eval_Jaco_Mat is

    res : Eval_Jaco_Mat(1..2,1..n);
    pq : Poly_Sys(1..2);
    jm : Jaco_Mat(1..2,1..n);

  begin
    pq(1) := p;
    pq(2) := q;
    jm := Create(pq);
    res := Create(jm);
    Clear(jm);
    return res;
  end Create;

  function Eval ( n1,n2 : natural; p,q : Eval_Poly;
                  b : Vector; v : VecVec; c : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the evaluation of b + c*v in the polynomials p and q.

    res : Vector(1..2);
    x : constant Vector := Affine_Expand(c,b,v); -- Combine(c,b,v);
    x1,x2 : Vector(1..n1);

  begin
    for i in 1..n1 loop
      x1(i) := x(i);
      x2(i) := x(i+n1);
    end loop;
    res(1) := Eval(p,x1);
    res(2) := Eval(q,x2);
    return res;
  end Eval;

  function Diff ( n1,n2 : natural; jm : Eval_Jaco_Mat;
                  b : Vector; v : VecVec; c : Vector ) return Matrix is

    res : Matrix(jm'range(1),v'range);
    x : constant Vector := Affine_Expand(c,b,v); --Combine(c,b,v);
    x1,x2 : Vector(1..n1);
   -- x1 : Vector(1..n1) := x(1..n1);
   -- x2 : Vector(1..n1) := x(n1+1..n2);
    eva : Matrix(jm'range(1),jm'range(2));

  begin
    for i in 1..n1 loop
      x1(i) := x(i);
      x2(i) := x(i+n1);
    end loop;
    for j in jm'range(2) loop
      eva(1,j) := Eval(jm(1,j),x1);
      eva(2,j) := Eval(jm(2,j),x2);
    end loop;
    for j in v'range loop
      res(1,j) := Create(0.0);
      res(2,j) := Create(0.0);
      for k in eva'range(2) loop
        res(1,j) := res(1,j) + eva(1,k)*v(j)(k);
        res(2,j) := res(2,j) + eva(2,k)*v(j)(k+n1);
      end loop;
    end loop;
    return res;
  end Diff;

-- TARGET ROUTINES :

  procedure Silent_Diagonal_Continuation
              ( n : in natural; p,q : in Poly;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    n2 : constant natural := 2*n;
    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    jm : Eval_Jaco_Mat(1..2,1..n) := Create(n,p,q);
    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(n,n2,ep,eq,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(n,n2,jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Sil_Cont(sols,false,Create(1.0));
    Clear(ep); Clear(eq); Clear(jm);
  end Silent_Diagonal_Continuation;

  procedure Reporting_Diagonal_Continuation
              ( file : in file_type; n : in natural; p,q : in Poly;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    n2 : constant natural := 2*n;
    ep : Eval_Poly := Create(p);
    eq : Eval_Poly := Create(q);
    jm : Eval_Jaco_Mat(1..2,1..n) := Create(n,p,q);
    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(n,n2,ep,eq,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(n,n2,jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Rep_Cont(file,sols,false,Create(1.0));
    Clear(ep); Clear(eq); Clear(jm);
  end Reporting_Diagonal_Continuation;

-- MORE EFFICIENT VERSIONS :

  procedure Silent_Diagonal_Continuation
              ( n : in natural; p,q : in Eval_Poly; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    n2 : constant natural := 2*n;
    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(n,n2,p,q,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(n,n2,jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Sil_Cont(sols,false,Create(1.0));
  end Silent_Diagonal_Continuation;

  procedure Reporting_Diagonal_Continuation
              ( file : in file_type;
                n : in natural; p,q : in Eval_Poly; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    n2 : constant natural := 2*n;
    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(n,n2,p,q,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(n,n2,jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Rep_Cont(file,sols,false,Create(1.0));
  end Reporting_Diagonal_Continuation;

end P_Intrinsic_Diagonal_Continuation;
