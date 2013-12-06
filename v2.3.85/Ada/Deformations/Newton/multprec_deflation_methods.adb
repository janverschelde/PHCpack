with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Matrices;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Complex_QR_Least_Squares; use Multprec_Complex_QR_Least_Squares;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;
with Multprec_Complex_Newton_Steps;
with Multprec_Deflate_Singularities;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;
with Multprec_Deflation_Trees_io;       use Multprec_Deflation_Trees_io;
with Monomial_Hashing;                  use Monomial_Hashing;

package body Multprec_Deflation_Methods is

-- INTERACTIVE NEWTON with DEFLATION :

  procedure Display_Results_of_One_Newton_Step
              ( file : in file_type;
                z : in Multprec_Complex_Vectors.Vector; tol : in double_float;
                err,rco,res : in Floating_Number; rank,cnt : in natural32 ) is

    use Multprec_Complex_Solutions;

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
    if integer32(rank) > n then
      put(file," > "); put(file,n,1); 
      put_line(file," number of variables, BUG !??");
    else
      if integer32(rank) = n
       then put(file," = ");
       else put(file," < ");
      end if;
      put(file,n,1);
      put(file,"  Corank is "); put(file,natural32(n)-rank,1); 
      put(file,"  w.r.t. tol :"); put(file,tol,3); new_line(file);
    end if;
  end Display_Results_of_One_Newton_Step;

  procedure One_Symbolic_Newton_Step
              ( file : in file_type;
                ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Newton_Steps;

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
  exception
    when others => put_line("exception raised in One_Symbolic_Newton_Step");
                   raise;
  end One_Symbolic_Newton_Step;

  procedure One_Symbolic_Newton_Step
              ( ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Newton_Steps;

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
  exception
    when others => put_line("exception raised in One_Symbolic_Newton_Step");
                   raise;
  end One_Symbolic_Newton_Step;

  function Stack ( n : natural32; x : Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the values of the vector x in one n-vector.

  -- REQUIRED : number of elements in x is n.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;

    res : Multprec_Complex_Vectors.Vector(1..integer32(n));
    cnt : integer32 := 0;

  begin
    for i in x'range loop
      exit when (x(i) = null);
      for j in x(i)'range loop
        cnt := cnt + 1;
        Copy(x(i)(j),res(cnt));
      end loop;
      exit when (cnt >= integer32(n));
    end loop;
    return res;
  end Stack;

  procedure Split ( z : in Multprec_Complex_Vectors.Vector;
                    x : in out Multprec_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Splits the values of the vector z over the vectors in x.

  -- REQUIRED :
  --   The number of elements in x equals the number of elements in z.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;

    cnt : integer32 := 0;

  begin
    for i in x'range loop
      exit when (x(i) = null);
      for j in x(i)'range loop
        cnt := cnt + 1;
        Copy(z(cnt),x(i)(j));
      end loop;
      exit when (cnt >= z'last);
    end loop;
  end Split;

  procedure One_Algorithmic_Newton_Step
              ( file : in file_type;
                ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Evaluate_Deflation;
    use Multprec_Complex_Newton_Steps;

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
    Newton(file,nq(k),z,tol,err,rco,res,s,rank);
    Split(z,x);
  end One_Algorithmic_Newton_Step;

  procedure One_Algorithmic_Newton_Step
              ( ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 ) is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Evaluate_Deflation;
    use Multprec_Complex_Newton_Steps;

    function f ( xx : Vector ) return Vector is
    begin
      return Eval(ep,xx);
    end f;

    function jmf ( xx : Vector ) return Matrix is
    begin
      Split(xx,x);
      return Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end jmf;

    procedure Newton is new Silent_Newton_Step_with_Singular_Values(f,jmf);

  begin
    Newton(natural32(ep'last),z,tol,err,rco,res,s,rank);
    Split(z,x);
  end One_Algorithmic_Newton_Step;

  procedure Interactive_Symbolic_Newton
              ( file : in file_type;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number;
                tol : in double_float; rank : out natural32 ) is

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;

    ep : constant Eval_Poly_Sys(p'range) := Create(p);
   -- wp : constant Link_to_Poly_Sys := new Poly_Sys'(p);
    jm : constant Jaco_Mat(p'range,z'range) := Create(p);
    ejm : constant Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    cnt : natural32 := 0;
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
  end Interactive_Symbolic_Newton;

  procedure Interactive_Algorithmic_Newton
              ( file : in file_type;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                err,rco,res : out Floating_Number;
                tol : in double_float; rank : out natural32 ) is

    ans : character;
    cnt : natural32 := 0;
    m : constant integer32 := Min0(integer32(nq(k))+1,integer32(nv(k)));
    s : Multprec_Complex_Vectors.Vector(1..m);
    z : Multprec_Complex_Vectors.Vector(1..integer32(nv(k)));

  begin
    loop
      cnt := cnt + 1;
      z := Stack(nv(k),x);
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
               ( p : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                 m,size : in natural32 ) is

  -- DESCRIPTION :
  --   Replaces the system in p with the deflated system,
  --   added with m multipliers.

    use Multprec_Complex_Poly_Systems;
    dp : constant Poly_Sys
       := Multprec_Deflate_Singularities.Deflate(p.all,m,size);

  begin
    Clear(p);
    p := new Poly_Sys'(dp);
  end Deflate;

  procedure Deflate_Corank_One
               ( p : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                 size : in natural32 ) is

  -- DESCRIPTION :
  --   Replaces the system in p with the deflated system,
  --   for the corank one case.

    use Multprec_Complex_Poly_Systems;
    dp : constant Poly_Sys
       := Multprec_Deflate_Singularities.Deflate_Corank_One(p.all,size);

  begin
    Clear(p);
    p := new Poly_Sys'(dp);
  end Deflate_Corank_One;

  procedure Add_Multipliers
               ( z : in out Multprec_Complex_Vectors.Link_to_Vector;
                 f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 m : in natural32 ) is

  -- DESCRIPTION :
  --   Adds values for the multipliers to the current root z.

  -- ON ENTRY :
  --   z         current approximation for the root;
  --   f         polynomial system after deflation with m multipliers;
  --   m         number of multipliers to be added to z.

  -- ON RETURN :
  --   z         vector extended with m values for the multipliers.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Deflate_Singularities;

    la : Vector(1..integer32(m));
    res : Floating_Number;
    nz : Vector(1..z'last+integer32(m));

  begin
    Multipliers(f,z.all,m,la,res);
    put_line("Values computed for the multipliers :");
    put_line(la);
    put("Residual of multiplier computation : ");
    put(res,3); new_line;
    Clear(res);
    for i in z'range loop
      Copy(z(i),nz(i));
    end loop;
    for i in 1..integer32(m) loop
      Copy(la(i),nz(i+z'last));
      Clear(la(i));
    end loop;
    Clear(z);
    z := new Vector'(nz);
  end Add_Multipliers;

  procedure Add_Multipliers
              ( file : in file_type;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Evaluate_Deflation;

    k1 : constant integer32 := k-1;
    t1 : constant Link_to_Eval_Tree := Create(natural32(k1));
    n : constant natural32 := nq(k1)+1;
    m : constant natural32 := R1(k);
    A1 : Matrix(1..integer32(nq(k1)),1..integer32(nv(k1)))
       := Eval(t1,nd,monkeys,k1,nv(0..k1),nq(0..k1),R1(1..k1),
               B(1..k1),h(1..k1),x(0..k1));
    AB1 : Matrix(1..integer32(nq(k1)),1..integer32(R1(k))) := A1*B(k).all;
    A,wrk : Matrix(1..integer32(nq(k1))+1,1..integer32(R1(k)));
    rhs : Vector(1..integer32(nq(k1))+1);
    lambda : Vector(1..integer32(R1(k)));
    qraux : Vector(A'range(2)) := (A'range(2) => Create(integer(0)));
    jpvt : Standard_Integer_Vectors.Vector(A'range(2))
         := (A'range(2) => 0);
    rsd,dum : Vector(A'range(1));
    info : integer32;
    res : Floating_Number;

  begin
    for i in AB1'range(1) loop
      for j in AB1'range(2) loop
        Copy(AB1(i,j),A(i,j));
      end loop;
    end loop;
    for j in 1..integer32(R1(k)) loop
      Copy(h(k)(j),A(A'last(1),j));
    end loop;
    for i in 1..integer32(nq(k1)) loop
      rhs(i) := Create(integer(0));
    end loop;
    rhs(rhs'last) := Create(integer(1));
    Copy(A,wrk);
    QRD(wrk,qraux,jpvt,false);
    QRLS(wrk,integer32(n),integer32(n),integer32(m),qraux,rhs,dum,dum,
         lambda,rsd,dum,110,info);
    dum := A*lambda;
    Min(dum);
    Add(dum,rhs);
    res := Max_Norm(dum);
    Clear(A1); Clear(AB1); Clear(A); Clear(wrk);
    Multprec_Complex_Vectors.Clear(qraux);
    Multprec_Complex_Vectors.Clear(rsd);
    Multprec_Complex_Vectors.Clear(dum);
    put_line(file,"Values computed for the multipliers :");
    put_line(file,lambda);
    put(file,"Residual of multiplier computation : ");
    put(file,res,3); new_line(file);
    Clear(res);
    x(k) := new Multprec_Complex_Vectors.Vector'(lambda);
  end Add_Multipliers;

  procedure Write_Results
               ( file : in file_type; i : natural32;
                 p,dp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 z : in Multprec_Complex_Vectors.Vector;
                 err,rco,res : in Floating_Number ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    n : constant integer32 := z'length;
    sols,sols_last : Solution_List;
    s : Solution(n);

  begin
    new_line(file);
    put(file,"DEFLATED SYSTEM #");
    put(file,i,1); put_line(file," :");
    Write_Deflated_System(file,p,dp);
    s.t := Create(integer(1));
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
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 size : in natural32; tol : in double_float ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    z : Vector(1..n);
    wz : Link_to_Vector;
    err,rco,res : Floating_Number;
    wp : Link_to_Poly_Sys;
    m,k,rank,len : natural32;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ans : character;

  begin
    if Is_Null(sols) then
      put("Give "); put(n,1);
      put_line(" complex numbers to start at :");
      get(z);
      len := 1;
    else
      declare
        stz : constant Standard_Complex_Vectors.Vector(1..n)
            := Head_Of(sols).v;
      begin
        z := Create(stz);
        Set_Size(z,size);
      end;
      len := Standard_Complex_Solutions.Length_Of(sols);
    end if;
    for i in 1..len loop
      wp := new Poly_Sys(p'range);
      for i in p'range loop
        Copy(p(i),wp(i));
      end loop;
      k := 0;
      wz := new Vector'(z);
      loop
        Interactive_Symbolic_Newton(file,wp.all,wz.all,err,rco,res,tol,rank);
        put("Do you want to deflate ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          m := rank + 1;
          put("The number of multipliers : "); put(m,1);
          if integer32(m) = wz'last then
            put_line(", corank 1.");
            Deflate_Corank_One(wp,size);
          else
            put_line(", corank "); put(wz'last-integer32(m)+1,1);
            put_line(".");
            Deflate(wp,m,size);
          end if;
          k := k+1;
          Add_Multiplier_Symbols(k,m);     
          Add_Multipliers(wz,wp.all,m);
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
        tmp := Standard_Complex_Solutions.Tail_Of(tmp);
        declare
          stz : constant Standard_Complex_Vectors.Vector(1..n)
              := Head_Of(tmp).v;
        begin
          z := Create(stz);
          Set_Size(z,size);
        end;
      end if;
    end loop;
  end Interactive_Symbolic_Deflation;

  procedure Interactive_Algorithmic_Deflation
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 size : in natural32;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 tol : in double_float ) is

    use Multprec_Complex_Numbers;
    use Multprec_Random_Numbers;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Multprec_Jacobian_Trees;
    use Multprec_Evaluate_Deflation;

    max_order : constant integer32 := 3;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nv : Standard_Natural_Vectors.Vector(0..max_order);
    nq : Standard_Natural_Vectors.Vector(0..max_order);
    R1 : Standard_Natural_Vectors.Vector(1..max_order);
    B : Multprec_Complex_VecMats.VecMat(1..max_order);
    h : Multprec_Complex_VecVecs.VecVec(1..max_order);
    x : Multprec_Complex_VecVecs.VecVec(0..max_order);
    ep : constant Eval_Poly_Sys(p'range) := Create(p);
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    ejm : constant Eval_Jaco_Mat(p'range,1..n) := Create(jm);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..max_order);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
    z : Multprec_Complex_Vectors.Vector(1..n);
    err,rco,res : Floating_Number;
    m,rank,len : natural32;
    k : integer32;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ans : character;

  begin
    nv(0) := natural32(n);
    nq(0) := natural32(p'last);
    R1(1) := 0;
    put("creating a remember table for Jacobian matrices up to order ");
    put(max_order,1); put_line("...");
    Create_Remember_Derivatives(jm,max_order,nd);
    monkeys := Monomial_Keys(natural32(max_order),nv(0));
   -- put("creating a remember table for deflation matrices up to order ");
   -- put(max_order,1); put_line("...");
   -- evt := Create(max_order);
    if Is_Null(sols) then
      put("Give "); put(n,1);
      put_line(" complex numbers to start at :");
      get(z);
      len := 1;
    else
      declare
        stz : constant Standard_Complex_Vectors.Vector(1..n)
            := Head_Of(sols).v;
      begin
        z := Create(stz);
        Set_Size(z,size);
      end;
      len := Standard_Complex_Solutions.Length_Of(sols);
    end if;
    for i in 1..integer32(len) loop
      k := 0;
      x(0) := new Multprec_Complex_Vectors.Vector'(z);
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
         -- Standard_Evaluate_Deflation.Clear(evt); OOPS!!
          evt := Create(natural32(k));
          Add_Multiplier_Symbols(natural32(k),m);     
          R1(k) := m;
          nv(k) := nv(k-1) + m;
          nq(k) := 2*nq(k-1) + 1;
          declare
            hk : Multprec_Complex_Vectors.Vector(1..integer32(m));
            Bk : Multprec_Complex_Matrices.Matrix
                   (1..integer32(nv(k-1)),1..integer32(m));
          begin
            for j in 1..integer32(m) loop
              hk(j) := Random(size);
            end loop;
            h(k) := new Multprec_Complex_Vectors.Vector'(hk);
            for j1 in Bk'range(1) loop
              for j2 in 1..integer32(m) loop
                Bk(j1,j2) := Random(size);
              end loop;
            end loop;
            B(k) := new Multprec_Complex_Matrices.Matrix'(Bk);
          end;
          Add_Multipliers(standard_output,ep,ejm,nd,monkeys,k,nv(0..k),
                          nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k));
        else
          put("Do you want to continue on this root ? (y/n) ");
          Ask_Yes_or_No(ans);
        end if;
        exit when (ans /= 'y');
      end loop;
     -- Write_Results(file,i,p,wp.all,wz.all,err,rco,res);
      if i < integer32(len) then
        put("Move on to root #"); put(i+1,1); put(" ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Standard_Complex_Solutions.Tail_Of(tmp);
        declare
          stz : constant Standard_Complex_Vectors.Vector(1..n)
              := Head_Of(tmp).v;
        begin
          z := Create(stz);
          Set_Size(z,size);
        end;
      end if;
    end loop;
  end Interactive_Algorithmic_Deflation;

-- NEWTON with DEFLATION and CLUSTERING :

  function Equal ( n : natural32; tol : double_float;
                   s1,s2 : Multprec_Complex_Vectors.Vector )
                 return boolean is

    use Multprec_Complex_Numbers;

  begin
    for i in 1..integer32(n) loop
      declare
        tmp : Floating_Number := AbsVal(s1(i)-s2(i));
      begin
        if tmp > tol
         then Clear(tmp); return false;
        end if;
        Clear(tmp);
      end;
    end loop;
    return true;
  end Equal;

  procedure Set_Multiplicity
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                s : in Multprec_Complex_Solutions.Solution;
                tol : in double_float; n,m : in natural32 ) is

    use Multprec_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,s.v) then
        ls.m := integer32(m);
        Set_Head(tmp,ls);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Set_Multiplicity;

  function Number_of_Occurrences 
              ( sols : Multprec_Complex_Solutions.Solution_List;
                s : Multprec_Complex_Solutions.Solution;
                tol : in double_float; n : in natural32 ) return natural32 is

    use Multprec_Complex_Solutions;

    res : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,s.v)
       then res := res + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Number_of_Occurrences;

  function Is_In ( sols : Multprec_Complex_Solutions.Solution_List;
                   v : Multprec_Complex_Vectors.Vector;
                   tol : double_float; n : natural32 ) return boolean is

    use Multprec_Complex_Solutions;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,v)
       then return true;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Is_In;

  procedure Remove_Duplicates
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                tol : in double_float; n : in natural32 ) is

    use Multprec_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if not Is_In(res,ls.v,tol,n)
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(sols);
    sols := res;
  end Remove_Duplicates;

  procedure Compute_Multiplicities
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                tol : in double_float; n : in natural32 ) is

    use Multprec_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    m : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      m := Number_of_Occurrences(sols,ls.all,tol,n);
      Set_Multiplicity(sols,ls.all,tol,n,m);
      tmp := Tail_Of(tmp);
    end loop;
  end Compute_Multiplicities;

  procedure Compute_Multiplicities
              ( nd : in out Multprec_Deflation_Trees.Node;
                tol : in double_float; n : in natural32 ) is

    use Multprec_Complex_Solutions;
    use Multprec_Deflation_Trees;

  begin
    if not Is_Null(nd.sols)
     then Compute_Multiplicities(nd.sols,tol,n);
          Remove_Duplicates(nd.sols,tol,n);
    end if;
    for i in nd.children'range loop
      if nd.children(i) /= null
       then Compute_Multiplicities(nd.children(i).all,tol,n);
      end if;
    end loop;
  end Compute_Multiplicities;

  procedure Deflate_Solution
              ( file : in file_type; m : in natural32; output : in boolean;
                df : in Multprec_Complex_Poly_Systems.Poly_Sys;
                ls : in out Multprec_Complex_Solutions.Link_to_Solution ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_Solutions;
    use Multprec_Deflate_Singularities;

    ns : Solution(ls.n+integer32(m));
    la : Vector(1..integer32(m));
    res : Floating_Number;

  begin
    Multipliers(df,ls.v,m,la,res);
    if output then
      put_line(file,"Values computed for the multipliers :");
      put_line(file,la);
      put(file,"Residual of multiplier computation :");
      put(file,res,3); new_line(file);
    end if;
    Clear(res);
    for i in ls.v'range loop
       Copy(ls.v(i),ns.v(i));
    end loop;
    for i in la'range loop
      Copy(la(i),ns.v(ls.v'last+i));
      Clear(la(i));
    end loop;
    Copy(ls.t,ns.t);
    ns.m := ls.m;
    Copy(ls.err,ns.err);
    Copy(ls.rco,ns.rco);
    Copy(ls.res,ns.res);
    Clear(ls);
    ls := new Solution'(ns);
  end Deflate_Solution;

  procedure Apply_Newton
              ( file : in file_type; output : in boolean; cnt : in natural32;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Multprec_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 ) is

    sv : Multprec_Complex_Vectors.Vector(ls.v'range);

  begin
    if output then
      One_Symbolic_Newton_Step
        (file,f,jf,ls.v,tolrnk,ls.err,ls.rco,ls.res,sv,rank);
      Display_Results_of_One_Newton_Step
        (file,ls.v,tolrnk,ls.err,ls.rco,ls.res,rank,cnt);
    else
      One_Symbolic_Newton_Step(f,jf,ls.v,tolrnk,ls.err,ls.rco,ls.res,sv,rank);
    end if;
    Multprec_Complex_Vectors.Clear(sv);
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
              ( file : in file_type; outfilename : in string;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Multprec_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef,size : in natural32;
	        tolerr,tolres,tolrnk : in double_float ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Solutions;
    use Multprec_Deflation_Trees;
    use Multprec_Deflation_Methods;

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
        Apply_Newton(file,output,nit,nd.f,nd.jf,ls,tolrnk,rank);
        if integer32(rank) < ls.n and k < maxdef then
          k := k + 1;
          if output
           then put(file,"Deflation #"); put(file,k,1);
                put(file," for solution "); put(file,i,1); new_line(file);
          end if;
          m := rank + 1;
          if nd.children(integer32(m)) = null then
            if output then Add_Multiplier_Symbols(k,m); end if;
            Create_Child(nd.all,integer32(m),size);
          end if;
          nd := nd.children(integer32(m));
          Deflate_Solution(file,m,output,nd.s,ls);
        else
          if (ls.err < tolerr) or (ls.res < tolres) then
            fail := false; exit;
          end if;
        end if;
      end loop;
      if output then Write_Conclusion(file,fail,i); end if;
      Append(nd.sols,nd.last,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    Compute_Multiplicities(lt.all,1.0E-8,natural32(nv));
    Write(file,outfilename,lt.all);
    Clear(lt);
  end Symbolic_Deflation_and_Clustering;

end Multprec_Deflation_Methods;
