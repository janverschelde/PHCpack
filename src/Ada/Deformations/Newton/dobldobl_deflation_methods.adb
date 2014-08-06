with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_QR_Least_Squares; use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;  use DoblDobl_Complex_Singular_Values;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_Multiple_Solutions;       use DoblDobl_Multiple_Solutions;
with DoblDobl_Complex_Newton_Steps;
with DoblDobl_Deflate_Singularities;
with DoblDobl_Deflation_Trees;
with DoblDobl_Deflation_Trees_io;       use DoblDobl_Deflation_Trees_io;
with Monomial_Hashing;                  use Monomial_Hashing;

package body DoblDobl_Deflation_Methods is

-- INTERACTIVE NEWTON with DEFLATION :

  procedure Display_Results_of_One_Newton_Step
              ( file : in file_type; z : in DoblDobl_Complex_Vectors.Vector;
                tol : double_float; err,rco,res : in double_double;
                rank,cnt : in natural32 ) is

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := z'length;
    s : Solution(n);

  begin
    put(file,"The new approximation at step ");
    put(file,cnt,1); put_line(file," : "); 
    s.v := z;
    put_vector(file,s);
    put(file,"error :"); put(file,err,3);
    put(file,"  rcond :"); put(file,rco,3);
    put(file,"  resid :"); put(file,res,3); new_line(file);
    put(file,"The numerical rank : "); put(file,rank,1);
    if rank > natural32(n) then
      put(file," > "); put(file,n,1); 
      put_line(file," number of variables, BUG !??");
    else
      if rank = natural32(n)
       then put(file," = ");
       else put(file," < ");
      end if;
      put(file,n,1);
      put(file,"  Corank is "); put(file,n-integer32(rank),1);
      put(file,"  w.r.t. tol :"); put(file,tol,3);  new_line(file);
    end if;
  end Display_Results_of_One_Newton_Step;

  procedure One_Symbolic_Newton_Step
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Newton_Steps;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is
    begin
      return Eval(ejm,x);
    end jmf;

    procedure Newton is new Reporting_Newton_Step_with_Singular_Values(f,jmf);

  begin
    Newton(file,natural32(ep'last),z,tol,err,rco,res,s,rank);
  end One_Symbolic_Newton_Step;

  procedure One_Symbolic_Newton_Step
              ( ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Newton_Steps;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is
    begin
      return Eval(ejm,x);
    end jmf;

    procedure Newton is new Silent_Newton_Step_with_Singular_Values(f,jmf);

  begin
    Newton(natural32(ep'last),z,tol,err,rco,res,s,rank);
  end One_Symbolic_Newton_Step;

  function Stack ( n : integer32; x : DoblDobl_Complex_VecVecs.VecVec )
                 return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the values of the vector x in one n-vector.

  -- REQUIRED : number of elements in x is n.

    res : DoblDobl_Complex_Vectors.Vector(1..n);
    cnt : integer32 := 0;

    use DoblDobl_Complex_Vectors;

  begin
    for i in x'range loop
      exit when (x(i) = null);
      for j in x(i)'range loop
        cnt := cnt + 1;
        res(cnt) := x(i)(j);
      end loop;
      exit when (cnt >= n);
    end loop;
    return res;
  end Stack;

  procedure Split ( z : in DoblDobl_Complex_Vectors.Vector;
                    x : in out DoblDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Splits the values of the vector z over the vectors in x.

  -- REQUIRED :
  --   The number of elements in x equals the number of elements in z.

    cnt : integer32 := 0;

    use DoblDobl_Complex_Vectors;

  begin
    for i in x'range loop
      exit when (x(i) = null);
      for j in x(i)'range loop
        cnt := cnt + 1;
        x(i)(j) := z(cnt);
      end loop;
      exit when (cnt >= z'last);
    end loop;
  end Split;

  procedure One_Algorithmic_Newton_Step
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                z : out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Evaluate_Deflation;
    use DoblDobl_Complex_Newton_Steps;

    function f ( xx : Vector ) return Vector is
    begin
      Split(xx,x);
      return Eval(ep,ejm,evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end f;

    function jmf ( xx : Vector ) return Matrix is
    begin
      Split(xx,x);
      return Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end jmf;

    procedure Newton is new Reporting_Newton_Step_with_Singular_Values(f,jmf);

  begin
    z := Stack(integer32(nv(k)),x);
    Newton(file,nq(k),z,tol,err,rco,res,s,rank);
    Split(z,x);
  end One_Algorithmic_Newton_Step;

  procedure One_Algorithmic_Newton_Step
              ( ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                z : out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Evaluate_Deflation;
    use DoblDobl_Complex_Newton_Steps;

    function f ( xx : Vector ) return Vector is
    begin
      Split(xx,x);
      return Eval(ep,ejm,evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end f;

    function jmf ( xx : Vector ) return Matrix is
    begin
      Split(xx,x);
      return Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end jmf;

    procedure Newton is new Silent_Newton_Step_with_Singular_Values(f,jmf);

  begin
    z := Stack(integer32(nv(k)),x);
    Newton(nq(k),z,tol,err,rco,res,s,rank);
    Split(z,x);
  end One_Algorithmic_Newton_Step;

  procedure Add_Multipliers
               ( file : in file_type; output : in boolean;
                 z : in out DoblDobl_Complex_Vectors.Link_to_Vector;
                 f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 m : in natural32; res : out double_double ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Deflate_Singularities;

    la : Vector(1..integer32(m));
    nz : Vector(1..z'last+integer32(m));
    zero : constant double_double := create(0.0);

  begin
    Multipliers(f,z.all,m,la,res);
    if output then
      put_line(file,"Values computed for the multipliers :");
      put_line(file,la);
      put(file,"Residual of multiplier computation : ");
      put(file,res,3); new_line(file);
    end if;
    nz(z'range) := z.all;
    if res < 0.1 then
      for i in la'range loop
        nz(z'last+i) := la(i);
      end loop;
    else
      for i in la'range loop
        nz(z'last+i) := DoblDobl_Complex_Numbers.Create(zero);
      end loop;
    end if;
    Clear(z);
    z := new Vector'(nz);
  end Add_Multipliers;

  procedure Add_Multipliers_for_Corank_One
               ( file : in file_type; output : in boolean;
                 z : in out DoblDobl_Complex_Vectors.Link_to_Vector;
                 f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 n,m : in natural32; res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Deflate_Singularities;

    la : Vector(1..integer32(m));
    nz : Vector(1..z'last+integer32(m));
    zero : constant double_double := create(0.0);

  begin
    Multipliers(f,z.all,n,m,la,res);
    if output then
      put_line(file,"Values computed for the multipliers :");
      put_line(file,la);
      put(file,"Residual of multiplier computation : ");
      put(file,res,3); new_line(file);
    end if;
    nz(z'range) := z.all;
    if res < 0.1 then
      for i in la'range loop
        nz(z'last+i) := la(i);
      end loop;
    else
      for i in la'range loop
        nz(z'last+i) := Create(zero);
      end loop;
    end if;
    Clear(z);
    z := new Vector'(nz);
  end Add_Multipliers_for_Corank_One;

  procedure Add_Multipliers
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Evaluate_Deflation;

    k1 : constant integer32 := k-1;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    t1 : Link_to_Eval_Tree := Create(natural32(k1));
    n : constant natural32 := nq(k1)+1;
    m : constant natural32 := R1(k);
    A1 : constant Matrix(1..integer32(nq(k1)),1..integer32(nv(k1)))
       := Eval(t1,nd,monkeys,k1,nv(0..k1),nq(0..k1),R1(1..k1),
               B(1..k1),h(1..k1),x(0..k1));
    AB1 : constant Matrix(1..integer32(nq(k1)),1..integer32(R1(k)))
        := A1*B(k).all;
    A,wrk : Matrix(1..integer32(nq(k1))+1,1..integer32(R1(k)));
    rhs : Vector(1..integer32(nq(k1))+1)
        := (1..integer32(nq(k1))+1 => Create(zero));
    lambda : Vector(1..integer32(R1(k)));
    qraux : Vector(A'range(2)) := (A'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(A'range(2))
         := (A'range(2) => 0);
    rsd,dum : Vector(A'range(1));
    info : integer32;

  begin
    for i in AB1'range(1) loop
      for j in AB1'range(2) loop
        A(i,j) := AB1(i,j);
      end loop;
    end loop;
    for j in 1..integer32(R1(k)) loop
      A(A'last(1),j) := h(k)(j);
    end loop;
    rhs(rhs'last) := Create(one);
    wrk := A;
    QRD(wrk,qraux,jpvt,false);
    QRLS(wrk,integer32(n),integer32(m),qraux,rhs,
         dum,dum,lambda,rsd,dum,110,info);
    dum := rhs - A*lambda;
    res := Max_Norm(dum);
    if res < 0.0 then
      x(k) := new DoblDobl_Complex_Vectors.Vector'(lambda);
    else
      lambda := (1..integer32(R1(k)) => Create(zero));
      x(k) := new DoblDobl_Complex_Vectors.Vector'(lambda);
    end if;
    DoblDobl_Evaluate_Deflation.Clear(t1);
  end Add_Multipliers;

  procedure Add_Multipliers
              ( file : in file_type; output : in boolean;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Evaluate_Deflation;

    k1 : constant integer32 := k-1;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    t1 : Link_to_Eval_Tree := Create(natural32(k1));
    n : constant integer32 := integer32(nq(k1)+1);
    m : constant natural32 := R1(k);
    A1 : constant Matrix(1..integer32(nq(k1)),1..integer32(nv(k1)))
       := Eval(t1,nd,monkeys,k1,nv(0..k1),nq(0..k1),R1(1..k1),
               B(1..k1),h(1..k1),x(0..k1));
    AB1 : constant Matrix(1..integer32(nq(k1)),1..integer32(R1(k)))
        := A1*B(k).all;
    A,wrk : Matrix(1..integer32(nq(k1))+1,1..integer32(R1(k)));
    rhs : Vector(1..integer32(nq(k1))+1)
        := (1..integer32(nq(k1))+1 => Create(zero));
    lambda : Vector(1..integer32(R1(k)));
    qraux : Vector(A'range(2)) := (A'range(2) => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(A'range(2))
         := (A'range(2) => 0);
    rsd,dum : Vector(A'range(1));
    info : integer32;

  begin
    for i in AB1'range(1) loop
      for j in AB1'range(2) loop
        A(i,j) := AB1(i,j);
      end loop;
    end loop;
    for j in 1..integer32(R1(k)) loop
      A(A'last(1),j) := h(k)(j);
    end loop;
    rhs(rhs'last) := Create(one);
    wrk := A;
    QRD(wrk,qraux,jpvt,false);
    QRLS(wrk,n,integer32(m),qraux,rhs,dum,dum,lambda,rsd,dum,110,info);
    dum := rhs - A*lambda;
    res := Max_Norm(dum);
    if output then
      put_line(file,"Values computed for the multipliers :");
      put_line(file,lambda);
      put(file,"Residual of multiplier computation : ");
      put(file,res,3); new_line(file);
    end if;
    if res < 0.1 then
      x(k) := new DoblDobl_Complex_Vectors.Vector'(lambda);
    else
      lambda := (1..integer32(R1(k)) => Create(zero));
      x(k) := new DoblDobl_Complex_Vectors.Vector'(lambda);
    end if;
    DoblDobl_Evaluate_Deflation.Clear(t1);
  end Add_Multipliers;

  procedure Add_Deflation_Data
              ( k : in integer32; m : in natural32;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                B : in out DoblDobl_Complex_VecMats.VecMat;
                h : in out DoblDobl_Complex_VecVecs.VecVec ) is

    hk : DoblDobl_Complex_Vectors.Vector(1..integer32(m));
    Bk : DoblDobl_Complex_Matrices.Matrix
           (1..integer32(nv(k-1)),1..integer32(m));

  begin
    R1(k) := m;
    nv(k) := nv(k-1) + m;
    nq(k) := 2*nq(k-1) + 1;
    for j in 1..integer32(m) loop
      hk(j) := DoblDobl_Random_Numbers.Random1;
    end loop;
    h(k) := new DoblDobl_Complex_Vectors.Vector'(hk);
    for j1 in Bk'range(1) loop
       for j2 in 1..integer32(m) loop
         Bk(j1,j2) := DoblDobl_Random_Numbers.Random1;
       end loop;
    end loop;
    B(k) := new DoblDobl_Complex_Matrices.Matrix'(Bk);
  end Add_Deflation_Data;

  procedure Interactive_Symbolic_Newton
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                tol : in double_float; rank : out natural32 ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,z'range) := Create(p);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    cnt : natural32 := 0;
   -- k : natural := 0;
    n : constant integer32 := p'length;
    q : constant integer32 := z'length;
    m : constant integer32 := Min0(n+1,q);
    s : Vector(1..m);

  begin
    loop
      cnt := cnt + 1;
      One_Symbolic_Newton_Step(file,ep,ejm,z,tol,err,rco,res,s,rank);
      put_line("The singular values : "); put_line(s);
      Display_Results_of_One_Newton_Step
        (standard_output,z,tol,err,rco,res,rank,cnt);
      Display_Results_of_One_Newton_Step(file,z,tol,err,rco,res,rank,cnt);
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Clear(ep); Clear(jm); Clear(ejm);
  end Interactive_Symbolic_Newton;

  procedure Interactive_Algorithmic_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tol : in double_float; rank : out natural32 ) is

    ans : character;
    cnt : natural32 := 0;
    m : constant integer32 := Min0(integer32(nq(k))+1,integer32(nv(k)));
    s : DoblDobl_Complex_Vectors.Vector(1..m);
    z : DoblDobl_Complex_Vectors.Vector(1..integer32(nv(k)));

  begin
    loop
      cnt := cnt + 1;
      One_Algorithmic_Newton_Step
        (file,f,jm,evt,nd,monkeys,k,nv,nq,R1,B,h,x,z,tol,err,rco,res,s,rank);
      put_line("The singular values : "); put_line(s);
      Display_Results_of_One_Newton_Step
        (standard_output,z,tol,err,rco,res,rank,cnt);
      Display_Results_of_One_Newton_Step(file,z,tol,err,rco,res,rank,cnt);
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Algorithmic_Newton;

  procedure Deflate
              ( p : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32 ) is

    use DoblDobl_Complex_Poly_Systems;
    dp : constant Poly_Sys := DoblDobl_Deflate_Singularities.Deflate(p.all,m);

  begin
    Clear(p);
    p := new Poly_Sys'(dp);
  end Deflate;

  procedure Deflate_Corank_One
              ( p : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;
    dp : constant Poly_Sys
       := DoblDobl_Deflate_Singularities.Deflate_Corank_One(p.all);

  begin
    Clear(p);
    p := new Poly_Sys'(dp);
  end Deflate_Corank_One;

  procedure Deflate_Corank_One
              ( p : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                k,nq,nv : natural32 ) is

    use DoblDobl_Complex_Poly_Systems;
    dp : constant Poly_Sys
       := DoblDobl_Deflate_Singularities.Deflate_Corank_One(p.all,k,nq,nv);

  begin
    put("k = "); put(k,1); put("  nq = "); put(nq,1);
    put("  nv = "); put(nv,1); new_line;
    put("dimensions before deflation are "); put(p'last,1); put(" by "); 
    put(Number_of_Unknowns(p(p'first)),1); new_line;
    Clear(p);
    p := new Poly_Sys'(dp);
    put("dimensions after deflation are "); put(p'last,1); put(" by "); 
    put(Number_of_Unknowns(p(p'first)),1); new_line;
  end Deflate_Corank_One;

  procedure Write_Results
               ( file : in file_type; i : natural32;
                 p,dp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 z : in DoblDobl_Complex_Vectors.Vector;
                 err,rco,res : in double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    n : constant integer32 := z'length;
    one : constant double_double := create(1.0);
    sols,sols_last : Solution_List;
    s : Solution(n);

  begin
    new_line(file);
    put(file,"DEFLATED SYSTEM #");
    put(file,i,1); put_line(file," :");
    Write_Deflated_System(file,p,dp);
    s.t := Create(one);
    s.m := 1;
    s.v := z;
    s.err := err; s.rco := rco; s.res := res;
    Append(sols,sols_last,s);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,"1  "); put(file,n,1); new_line(file);
    put(file,sols);
  end Write_Results;

  procedure Interactive_Symbolic_Deflation
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sol : in DoblDobl_Complex_Vectors.Vector;
                 tol : in double_float ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    k1,nq1,nv1 : natural32 := 0;  -- for corank 1 deflation
    wz : Link_to_Vector := new Vector'(sol);
    err,rco,res,rsd : double_double;
    wp : Link_to_Poly_Sys;
    m,k,rank : natural32;
    ans : character;

  begin
    wp := new Poly_Sys(p'range);
    Copy(p,wp.all);
    k := 0;
    loop
      Interactive_Symbolic_Newton(file,wp.all,wz.all,err,rco,res,tol,rank);
      put("Do you want to deflate ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        m := rank+1;
        put("The number of multipliers : "); put(m,1);
        if m = natural32(wz'last) then
          put_line(", corank 1.");
          if k1 = 0 then
            k1 := 1;
            nq1 := natural32(wp'last);
            nv1 := Number_of_Unknowns(wp(wp'first));
          else
            k1 := k1+1;
          end if;
         -- Deflate_Corank_One(wp);
          Deflate_Corank_One(wp,k1,nq1,nv1);
        else
          put(", corank "); put(wz'last-integer32(m)+1,1); put_line(".");
          Deflate(wp,m);
        end if;
        k := k+1;
        if k1 > 1
         then Add_Multiplier_Symbols(k,nv1);     
              Add_Multipliers_for_Corank_One
                (standard_output,true,wz,wp.all,nq1+1,nv1,rsd);
         else Add_Multiplier_Symbols(k,m);     
              Add_Multipliers(standard_output,true,wz,wp.all,m,rsd);
        end if;
      else
        put("Do you want to continue ? (y/n) ");
        Ask_Yes_or_No(ans);
      end if;
      exit when (ans /= 'y');
    end loop;
    Write_Results(file,1,p,wp.all,wz.all,err,rco,res);
  end Interactive_Symbolic_Deflation;

  procedure Interactive_Algorithmic_Deflation
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sol : in DoblDobl_Complex_Vectors.Vector;
                 tol : in double_float ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Jacobian_Trees;
    use DoblDobl_Evaluate_Deflation;

    max_order : constant integer32 := 3;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nv : Standard_Natural_Vectors.Vector(0..max_order);
    nq : Standard_Natural_Vectors.Vector(0..max_order);
    R1 : Standard_Natural_Vectors.Vector(1..max_order);
    B : DoblDobl_Complex_VecMats.VecMat(1..max_order);
    h : DoblDobl_Complex_VecVecs.VecVec(1..max_order);
    x : DoblDobl_Complex_VecVecs.VecVec(0..max_order);
    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..max_order);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
    err,rco,res,rsd : double_double;
    m,rank : natural32;
    k : integer32;
    ans : character;

  begin
    nv(0) := natural32(n);
    nq(0) := natural32(p'last);
    R1(1) := 0;
    put("creating a remember table for Jacobian matrices up to order ");
    put(max_order,1); put_line("...");
    Create_Remember_Derivatives(jm,max_order,nd);
    monkeys := Monomial_Keys(natural32(max_order),nv(0));
    k := 0;
    x(0) := new DoblDobl_Complex_Vectors.Vector'(sol);
    loop
      if k = 0 then
        Interactive_Symbolic_Newton(file,p,x(0).all,err,rco,res,tol,rank);
      else
        Interactive_Algorithmic_Newton
          (file,ep,ejm,evt,nd,monkeys,k,nv(0..k),nq(0..k),R1(1..k),
           B(1..k),h(1..k),x(0..k),err,rco,res,tol,rank);
      end if;
      put("Do you want to deflate ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        m := rank+1;
        put("The number of multipliers : "); put(m,1); new_line; 
        k := k+1;
        evt := Create(natural32(k));
        Add_Multiplier_Symbols(natural32(k),m);     
        Add_Deflation_Data(k,m,nv,nq,R1,B,h);
        Add_Multipliers(standard_output,true,ep,ejm,nd,monkeys,k,nv(0..k),
                        nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
      else
        put("Do you want to continue ? (y/n) ");
        Ask_Yes_or_No(ans);
      end if;
      exit when (ans /= 'y');
    end loop;
    Clear(ep); Clear(jm); Clear(ejm);
  end Interactive_Algorithmic_Deflation;

  procedure Interactive_Symbolic_Deflation
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 tol : in double_float ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    k1,nq1,nv1 : natural32 := 0; -- for corank 1 deflation
    n : constant natural32 := Number_of_Unknowns(p(p'first));
    z : Vector(1..integer32(n));
    wz : Link_to_Vector;
    err,rco,res,rsd : double_double;
    wp : Link_to_Poly_Sys;
    m,k,rank,len : natural32;
    tmp : Solution_List := sols;
    ans : character;

  begin
    if Is_Null(sols) then
      put("Give "); put(n,1);
      put_line(" complex numbers to start at :");
      get(z);
      len := 1;
    else
      z := Head_Of(sols).v;
      len := Length_Of(sols);
    end if;
    for i in 1..len loop
      wp := new Poly_Sys(p'range);
      Copy(p,wp.all);
      k := 0;
      wz := new Vector'(z);
      loop
        Interactive_Symbolic_Newton(file,wp.all,wz.all,err,rco,res,tol,rank);
        put("Do you want to deflate ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          m := rank+1;
          put("The number of multipliers : "); put(m,1);
          if m = natural32(wz'last) then
            put_line(", corank 1.");
            if k1 = 0 then       -- first time corank 1 deflation
              k1 := 1;
              nq1 := natural32(wp'last);
              nv1 := Number_of_Unknowns(wp(wp'first));
            else 
              k1 := k1+1;
            end if;
           -- Deflate_Corank_One(wp);
            Deflate_Corank_One(wp,k1,nq1,nv1);
          else
            put(", corank "); put(wz'last-integer32(m)+1,1); put_line(".");
            Deflate(wp,m);
          end if;
          k := k+1;
          if k1 > 1 then
            Add_Multiplier_Symbols(k,nv1);
            Add_Multipliers_for_Corank_One
              (standard_output,true,wz,wp.all,nq1+1,nv1,rsd);
          else
            Add_Multiplier_Symbols(k,m);     
            Add_Multipliers(standard_output,true,wz,wp.all,m,rsd);
          end if;
        else
          put("Do you want to continue on this root ? (y/n) ");
          Ask_Yes_or_No(ans);
        end if;
        exit when (ans /= 'y');
      end loop;
      Write_Results(file,i,p,wp.all,wz.all,err,rco,res);
      if i < len then
        put("Move on to root #"); put(i+1,1); put(" ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Tail_Of(tmp);
        z := Head_Of(tmp).v;
      end if;
    end loop;
  end Interactive_Symbolic_Deflation;

  procedure Interactive_Algorithmic_Deflation
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 tol : in double_float ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Jacobian_Trees;
    use DoblDobl_Evaluate_Deflation;

    max_order : constant integer32 := 3;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nv : Standard_Natural_Vectors.Vector(0..max_order);
    nq : Standard_Natural_Vectors.Vector(0..max_order);
    R1 : Standard_Natural_Vectors.Vector(1..max_order);
    B : DoblDobl_Complex_VecMats.VecMat(1..max_order);
    h : DoblDobl_Complex_VecVecs.VecVec(1..max_order);
    x : DoblDobl_Complex_VecVecs.VecVec(0..max_order);
    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..max_order);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
    z : DoblDobl_Complex_Vectors.Vector(1..n);
    err,rco,res,rsd : double_double;
    m,rank,len : natural32;
    k : integer32;
    tmp : Solution_List := sols;
    ans : character;

  begin
    nv(0) := natural32(n);
    nq(0) := natural32(p'last);
    R1(1) := 0;
    put("creating a remember table for Jacobian matrices up to order ");
    put(max_order,1); put_line("...");
    Create_Remember_Derivatives(jm,max_order,nd);
    monkeys := Monomial_Keys(natural32(max_order),nv(0));
    if Is_Null(sols) then
      put("Give "); put(n,1);
      put_line(" complex numbers to start at :");
      get(z);
      len := 1;
    else
      z := Head_Of(sols).v;
      len := Length_Of(sols);
    end if;
    for i in 1..len loop
      k := 0;
      x(0) := new DoblDobl_Complex_Vectors.Vector'(z);
      loop
        if k = 0 then
          Interactive_Symbolic_Newton(file,p,x(0).all,err,rco,res,tol,rank);
        else
          Interactive_Algorithmic_Newton
            (file,ep,ejm,evt,nd,monkeys,k,nv(0..k),nq(0..k),R1(1..k),
             B(1..k),h(1..k),x(0..k),err,rco,res,tol,rank);
        end if;
        put("Do you want to deflate ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          m := rank+1;
          put("The number of multipliers : "); put(m,1); new_line; 
          k := k+1;
         -- DoblDobl_Evaluate_Deflation.Clear(evt); OOPS!!
          evt := Create(natural32(k));
          Add_Multiplier_Symbols(natural32(k),m);     
          Add_Deflation_Data(k,m,nv,nq,R1,B,h);
          Add_Multipliers(standard_output,true,ep,ejm,nd,monkeys,k,nv(0..k),
                          nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
        else
          put("Do you want to continue on this root ? (y/n) ");
          Ask_Yes_or_No(ans);
        end if;
        exit when (ans /= 'y');
      end loop;
     -- Write_Results(file,i,p,wp.all,wz.all,err,rco,res);
      if i < len then
        put("Move on to root #"); put(i+1,1); put(" ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Tail_Of(tmp);
        z := Head_Of(tmp).v;
      end if;
    end loop;
    Clear(ep); Clear(jm); Clear(ejm);
  end Interactive_Algorithmic_Deflation;

-- NEWTON with DEFLATION and CLUSTERING : 

  procedure Deflate_Solution
              ( file : in file_type; m : in integer32; output : in boolean;
                df : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Deflate_Singularities;

    ns : Solution(ls.n+m);
    la : Vector(1..m);
    res : double_double;

  begin
    Multipliers(df,ls.v,natural32(m),la,res);
    if output then
      put_line(file,"Values computed for the multipliers :");
      put_line(file,la);
      put(file,"Residual of multiplier computation :");
      put(file,res,3); new_line(file);
    end if;
    ns.v(ls.v'range) := ls.v;
    for i in la'range loop
      ns.v(ls.v'last+i) := la(i);
    end loop;
    ns.t := ls.t; ns.m := ls.m;
    ns.err := ls.err; ns.rco := ls.rco; ns.res := ls.res;
    Clear(ls);
    ls := new Solution'(ns);
  end Deflate_Solution;

  procedure Apply_Newton_Step
              ( file : in file_type; output : in boolean; step : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 ) is

    sv : DoblDobl_Complex_Vectors.Vector(ls.v'range);

  begin
    if output then
      One_Symbolic_Newton_Step
        (file,f,jf,ls.v,tolrnk,ls.err,ls.rco,ls.res,sv,rank);
      Display_Results_of_One_Newton_Step
        (file,ls.v,tolrnk,ls.err,ls.rco,ls.res,rank,step);
    else
      One_Symbolic_Newton_Step(f,jf,ls.v,tolrnk,ls.err,ls.rco,ls.res,sv,rank);
    end if;
  end Apply_Newton_Step;

  procedure Apply_Newton_Step
              ( file : in file_type; output : in boolean; step : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 ) is

    m : constant integer32 := Min0(integer32(nq(k))+1,integer32(nv(k)));
    s : DoblDobl_Complex_Vectors.Vector(1..m);
    z : DoblDobl_Complex_Vectors.Vector(1..integer32(nv(k)));
    evt : DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree
        := DoblDobl_Evaluate_Deflation.Create(natural32(k));

  begin
    if output then
      One_Algorithmic_Newton_Step
        (file,f,jf,evt,nd,monkeys,k,nv(0..k),nq(0..k),R1(1..k),
         B(1..k),h(1..k),x(0..k),z,tolrnk,err,rco,res,s,rank);
      Display_Results_of_One_Newton_Step
        (file,z,tolrnk,err,rco,res,rank,step);
    else
      One_Algorithmic_Newton_Step
        (f,jf,evt,nd,monkeys,k,nv(0..k),nq(0..k),R1(1..k),
         B(1..k),h(1..k),x(0..k),z,tolrnk,err,rco,res,s,rank);
    end if;
    DoblDobl_Evaluate_Deflation.Clear(evt);
  end Apply_Newton_Step;

  procedure Apply_Newton
              ( nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 ) is

    sv : DoblDobl_Complex_Vectors.Vector(ls.v'range);

  begin
    for i in 1..nit loop
      One_Symbolic_Newton_Step(f,jf,ls.v,tolrnk,ls.err,ls.rco,ls.res,sv,rank);
    end loop;
  end Apply_Newton;

  procedure Apply_Newton
              ( file : in file_type; output : in boolean; nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 ) is

  begin
    for i in 1..nit loop
      Apply_Newton_Step(file,output,i,f,jf,ls,tolrnk,rank);
    end loop;
  end Apply_Newton;

  procedure Apply_Newton
              ( nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 ) is

    m : constant integer32 := Min0(integer32(nq(k))+1,integer32(nv(k)));
    s : DoblDobl_Complex_Vectors.Vector(1..m);
    z : DoblDobl_Complex_Vectors.Vector(1..integer32(nv(k)));
    evt : DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;

  begin
    for i in 1..nit loop
      evt := DoblDobl_Evaluate_Deflation.Create(natural32(k));
      One_Algorithmic_Newton_Step
        (f,jf,evt,nd,monkeys,k,nv(0..k),nq(0..k),R1(1..k),
         B(1..k),h(1..k),x(0..k),z,tolrnk,err,rco,res,s,rank);
      DoblDobl_Evaluate_Deflation.Clear(evt);
    end loop;
  end Apply_Newton;

  procedure Apply_Newton
              ( file : in file_type; output : in boolean; nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 ) is

  begin
    for i in 1..nit loop
      Apply_Newton_Step
        (file,output,i,f,jf,nd,monkeys,k,nv,nq,R1,B,h,x,err,rco,res,
         tolrnk,rank);
    end loop;
  end Apply_Newton;

  procedure Write_Conclusion
              ( file : in file_type; fail : in boolean; i : in natural32 ) is
  begin
    put(file,"Convergence criteria");
    if not fail
     then put(file," are reached for solution ");
     else put(file," NOT reached for solution ");
    end if;
    put(file,i,1); put_line(file,".");
  end Write_Conclusion;

  procedure Symbolic_Deflation_and_Clustering
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Deflation_Trees;

    k,rank,m,nit : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    fail : boolean;
    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Node(ne,nv) := Create_Root(p);
    lt : Link_to_Node := new Node'(t);
    nd : Link_to_Node;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      if output
       then put(file,"SOLUTION #"); put(file,i,1); put_line(file," :");
      end if;
      nd := lt; k := 0; nit := 0; fail := true;
      loop
        nit := nit + 1;
        exit when (nit > maxitr);
        Apply_Newton_Step(file,output,nit,nd.f,nd.jf,ls,tolrnk,rank);
        if rank < natural32(ls.n) and k < maxdef then
          k := k + 1;
          if output then
            put(file,"Deflation #"); put(file,k,1);
            put(file," for solution "); put(file,i,1); new_line(file);
          end if;
          m := rank + 1;
          if nd.children(integer32(m)) = null then
            if output then Add_Multiplier_Symbols(k,m); end if;
            Create_Child(nd.all,integer32(m));
          end if;
          nd := nd.children(integer32(m));
          Deflate_Solution(file,integer32(m),output,nd.s,ls);
        else
          if (ls.err <= tolerr) or (ls.res <= tolres) then
            fail := false; exit;
          end if;
        end if;
      end loop;
      if output then Write_Conclusion(file,fail,i); end if;
      Append(nd.sols,nd.last,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    Compute_Multiplicities(lt.all,1.0E-8,natural32(nv));
    Write(file,name,lt.all);
    Clear(lt);
  end Symbolic_Deflation_and_Clustering;

  function Minimum ( a,b : natural32 ) return natural32 is
  begin
    if a < b 
     then return a;
     else return b;
    end if;
  end Minimum;

  procedure Algorithmic_Deflation_and_Clustering
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Jacobian_Trees;
    use DoblDobl_Evaluate_Deflation;

    max_order : constant natural32 := 3;
    order : constant integer32 := integer32(Minimum(max_order,maxdef));
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    B : DoblDobl_Complex_VecMats.VecMat(1..order);
    h : DoblDobl_Complex_VecVecs.VecVec(1..order);
    x : DoblDobl_Complex_VecVecs.VecVec(0..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    ep : Eval_Poly_Sys(p'range) := Create(p);
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
    tmp : Solution_list := sols;
    ls : Link_to_Solution;
    rank,m : natural32;
    k : integer32;
    done : boolean := true;
    rsd : double_double;

  begin
    if not Is_Null(sols) then
      nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
      Create_Remember_Derivatives(jm,order,nd);
      monkeys := Monomial_Keys(natural32(order),nv(0));
      for i in 1..Length_Of(sols) loop
        ls := Head_Of(tmp); k := 0;
        loop
          if k = 0 then
            Apply_Newton(maxitr,ep,ejm,ls,tolrnk,rank);
            done := (rank = natural32(ls.n));
            if not done then
              x(0) := new DoblDobl_Complex_Vectors.Vector'(ls.v);
            end if;
          else
            Apply_Newton(maxitr,ep,ejm,nd,monkeys,k,
                         nv,nq,R1,B,h,x,ls.err,ls.rco,ls.res,tolrnk,rank);
            ls.v := x(0).all;
            Set_Head(tmp,ls);
            done := (rank = nv(k));
          end if;
          exit when done or (k = order);
          k := k + 1; m := rank + 1;
          Add_Multiplier_Symbols(natural32(k),m);
          Add_Deflation_Data(k,m,nv,nq,R1,B,h);
          Add_Multipliers(ep,ejm,nd,monkeys,k,nv(0..k),
                          nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
          exit when (rsd > 0.1);  -- deflation fails
        end loop;
        DoblDobl_Complex_Vectors.Clear(x(0));
        for i in 1..k loop
          DoblDobl_Complex_Matrices.Clear(B(i));
          DoblDobl_Complex_Vectors.Clear(h(i));
          DoblDobl_Complex_Vectors.Clear(x(i));
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      DoblDobl_Complex_Jaco_Matrices.Clear(ejm);
      DoblDobl_Complex_Jaco_Matrices.Clear(jm); 
      DoblDobl_Complex_Poly_SysFun.Clear(ep);
      Standard_Natural64_VecVecs.Clear(monkeys);
     -- DoblDobl_Jacobian_Trees.Clear(nd);
      Compute_Multiplicities(sols,1.0E-8,natural32(n));
      Remove_Duplicates(sols,1.0E-8,natural32(n));
    end if;
  end Algorithmic_Deflation_and_Clustering;

  procedure Algorithmic_Deflation_and_Clustering
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Jacobian_Trees;
    use DoblDobl_Evaluate_Deflation;

    max_order : constant natural32 := 3;
    order : constant integer32 := integer32(Minimum(max_order,maxdef));
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    B : DoblDobl_Complex_VecMats.VecMat(1..order);
    h : DoblDobl_Complex_VecVecs.VecVec(1..order);
    x : DoblDobl_Complex_VecVecs.VecVec(0..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    ep : Eval_Poly_Sys(p'range) := Create(p);
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
    tmp : Solution_list := sols;
    ls : Link_to_Solution;
    rank,m : natural32;
    k : integer32;
    done : boolean := true;
    rsd : double_double;

  begin
    if not Is_Null(sols) then
      nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
      Create_Remember_Derivatives(jm,order,nd);
      monkeys := Monomial_Keys(natural32(order),nv(0));
      for i in 1..Length_Of(sols) loop
        ls := Head_Of(tmp); k := 0;
        if output
         then put(file,"SOLUTION #"); put(file,i,1); put_line(file," :");
        end if;
        loop
          if k = 0 then
            Apply_Newton(file,output,maxitr,ep,ejm,ls,tolrnk,rank);
            done := (rank = natural32(ls.n));
            if not done then
              x(0) := new DoblDobl_Complex_Vectors.Vector'(ls.v);
            end if;
          else
            Apply_Newton(file,output,maxitr,ep,ejm,nd,monkeys,k,
                         nv,nq,R1,B,h,x,ls.err,ls.rco,ls.res,tolrnk,rank);
            ls.v := x(0).all;
            Set_Head(tmp,ls);
            done := (rank = nv(k));
          end if;
          exit when done or (k = order);
          k := k + 1; m := rank + 1;
          Add_Multiplier_Symbols(natural32(k),m);
          Add_Deflation_Data(k,m,nv,nq,R1,B,h);
          Add_Multipliers(file,output,ep,ejm,nd,monkeys,k,nv(0..k),
                          nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
          exit when (rsd > 0.1);  -- deflation fails
        end loop;
        DoblDobl_Complex_Vectors.Clear(x(0));
        for i in 1..k loop
          DoblDobl_Complex_Matrices.Clear(B(i));
          DoblDobl_Complex_Vectors.Clear(h(i));
          DoblDobl_Complex_Vectors.Clear(x(i));
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      DoblDobl_Complex_Jaco_Matrices.Clear(ejm);
      DoblDobl_Complex_Jaco_Matrices.Clear(jm); 
      DoblDobl_Complex_Poly_SysFun.Clear(ep);
      Standard_Natural64_VecVecs.Clear(monkeys);
     -- DoblDobl_Jacobian_Trees.Clear(nd);
      Compute_Multiplicities(sols,1.0E-8,natural32(n));
      Remove_Duplicates(sols,1.0E-8,natural32(n));
    end if;
  end Algorithmic_Deflation_and_Clustering;

end DoblDobl_Deflation_Methods;
