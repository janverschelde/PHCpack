with text_io;                           use text_io;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;          use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;          use DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;         use DoblDobl_Complex_Matrices;
with DoblDobl_Random_Numbers;           use DoblDobl_Random_Numbers;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with DoblDobl_Complex_Singular_Values;  use DoblDobl_Complex_Singular_Values;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;      use DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;    use DoblDobl_Complex_Jaco_Matrices;
with Projective_Transformations;        use Projective_Transformations;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;
with DoblDobl_Plane_Representations;    use DoblDobl_Plane_Representations;
with DoblDobl_Point_Coordinates;        use DoblDobl_Point_Coordinates;
with DoblDobl_Moving_Planes;            use DoblDobl_Moving_Planes;
with DoblDobl_Intrinsic_Solutions;      use DoblDobl_Intrinsic_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with DoblDobl_Intrinsic_Newton;         use DoblDobl_Intrinsic_Newton;

procedure ts_iddnewt is

-- DESCRIPTION :
--   Test on Newton's method with intrinsic coordinates,
--   in double double precision.  

  procedure Perturb ( x : in out Vector; tol : in double_double ) is

  -- DESCRIPTION :
  --   Perturbs every entry in x with magnitude equal to tol.

    dx : constant Vector(x'first..x'last) := Random_Vector(x'first,x'last);

  begin
    for i in x'range loop
      x(i) := x(i) + tol*dx(i);
    end loop;
  end Perturb;

  procedure Call_Affine_Newton
              ( f : in Poly_Sys; p : in Matrix; x : in out Vector;
                method : in character;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco : double_double;
    nit : natural32;

  begin
    case method is
      when '1' =>
        Affine_LU_Newton(Standard_Output,f,p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '2' =>
        Affine_LU_Newton(Standard_Output,f,p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      when '3' =>
        Affine_QR_Newton(Standard_Output,f,p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '4' =>
        Affine_SV_Newton(Standard_Output,f,p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,sv,fail);
      when others => null;
    end case;
  end Call_Affine_Newton;

  procedure Call_Generic_Affine_Newton
              ( f : in Poly_Sys; p : in Matrix; x : in out Vector;
                method : in character;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco : double_double;
    nit : natural32;
    ef : Eval_Poly_Sys(f'range) := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);

    function Eval ( x : Vector ) return Vector is
    begin
      return Eval(ef,x);
    end Eval;
    function Diff ( x : Vector ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff;
    procedure R_Affine_LU_Newton is
      new Reporting_Affine_LU_Newton(Eval,Diff);
    procedure R_Affine_LU_RCO_Newton is
      new Reporting_Affine_LU_RCO_Newton(Eval,Diff);
    procedure R_Affine_QR_Newton is
      new Reporting_Affine_QR_Newton(Eval,Diff);
    procedure R_Affine_SV_Newton is
      new Reporting_Affine_SV_Newton(Eval,Diff);

  begin
    case method is
      when '1' =>
        R_Affine_LU_Newton(Standard_Output,natural32(f'last),p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '2' =>
        R_Affine_LU_RCO_Newton(Standard_Output,natural32(f'last),p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      when '3' =>
        R_Affine_QR_Newton(Standard_Output,natural32(f'last),p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '4' =>
        R_Affine_SV_Newton(Standard_Output,natural32(f'last),p,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,sv,fail);
      when others => null;
    end case;
    Clear(ef); Clear(jm); Clear(jf);
  end Call_Generic_Affine_Newton;

 -- procedure Test_Projective_Eval
 --             ( f : in Poly_Sys; p : in Matrix; x : in Vector ) is
 --
 --   z : constant Vector := Projective_Expand(x,p);
 --   y : constant Vector := Eval(f,z);
 --
 -- begin
 --   put_line("The vector in extrinsic coordinates : "); put_line(z);
 --   put_line("Residual vector after evaluation : "); put_line(y);
 -- end Test_Projective_Eval;

  procedure Call_Projective_Newton
              ( f : in Poly_Sys; p : in Matrix; x : in out Vector;
                method : in character;
                epsax,epsrx,epsaf,epsrf : in double_double;
                incax,incrx,resaf,resrf : out double_double;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco,nrm : double_double;
    k : integer32;
    nit : natural32;
   -- z : constant Vector := Projective_Expand(x,p);
   -- y : constant Vector := Eval(f,z);

  begin
   -- put_line("The solution before scaling : "); put_line(x);
   -- Test_Projective_Eval(f,p,x);
    Max_Norm(x,k,nrm);
    put("k = "); put(k,1); put(" with norm = "); put(nrm,3); new_line;
    Scale(x,k);
   -- put_line("The solution after scaling : "); put_line(x);
   -- Test_Projective_Eval(f,p,x);
    case method is
      when '1' =>
        Projective_LU_Newton(Standard_Output,f,p,x,natural32(k),
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '2' =>
        Projective_LU_Newton(Standard_Output,f,p,x,natural32(k),
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      when '3' =>
        Projective_QR_Newton(Standard_Output,f,p,x,natural32(k),
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,fail);
      when '4' =>
        Projective_SV_Newton(Standard_Output,f,p,x,natural32(k),
          epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,sv,fail);
      when others => null;
    end case;
  end Call_Projective_Newton;

  function Make_Projective ( p : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Adds a row of ones to the matrix p.

    res : Matrix(p'first(1)..p'last(1)+1,p'range(2));
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    for j in p'range(2) loop
      for i in p'range(1) loop
        res(i,j) := p(i,j);
      end loop;
      res(res'last(1),j) := Create(zero);
    end loop;
    res(res'last(1),0) := Create(one);
    return res;
  end Make_Projective;

  procedure Newton ( f : in Poly_Sys; sols : in Solution_List;
                     p : in Matrix; method : in character ) is

    tmp : Solution_List := sols;
    epsax : constant double_double := create(1.0E-30);
    epsrx : constant double_double := create(1.0E-30);
    epsaf : constant double_double := create(1.0E-30);
    epsrf : constant double_double := create(1.0E-30);
    incax,incrx,resaf,resrf : double_double;
    cnt : natural32 := 0;
    max : constant natural32 := 5;
    fail : boolean;
    nbfail : natural32 := 0;
    ans : character;
    x : Vector(1..p'last(2));
    z : Vector(0..p'last(2));
    pf : Poly_Sys(f'range);
    pp : Matrix(p'first(1)..p'last(1)+1,p'range(2));
    projective,generic_newton : boolean;
    prbtol : constant double_double := create(0.0001);

  begin
    new_line;
    put("Use projective coordinates ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      pf := Projective_Transformation(f);
      pp := Make_Projective(p);
      new_line;
      put_line("The homogenized system : "); put(pf);
      projective := true; generic_newton := false;
    else
      put("Use the generic version of affine Newton's method ? (y/n) ");
      Ask_Yes_or_No(ans);
      projective := false; generic_newton := (ans = 'y');
    end if;
    new_line;
    while not Is_Null(tmp) loop
      x := Head_Of(tmp).v;
      cnt := cnt + 1;
      put("refining solution "); put(cnt,1); put_line(" ...");
      Perturb(x,prbtol);
      if not projective then
        if generic_newton then
          Call_Generic_Affine_Newton
            (f,p,x,method,epsax,epsrx,epsaf,epsrf,
             incax,incrx,resaf,resrf,max,fail);
        else
          Call_Affine_Newton(f,p,x,method,epsax,epsrx,epsaf,epsrf,
                             incax,incrx,resaf,resrf,max,fail);
        end if;
      else
        z := Projective_Coordinates(x);
        Call_Projective_Newton(pf,pp,z,method,epsax,epsrx,epsaf,epsrf,
                               incax,incrx,resaf,resrf,max,fail);
      end if;
      if fail then
        put_line("FAILED TO CONVERGE!!!");
        nbfail := nbfail + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put("#failures : "); put(nbfail,1); new_line;
  end Newton;

  procedure Call_Newton
              ( f : in Poly_Sys; sols : in Solution_List;
                p : in Matrix; method : in character ) is

  -- DESCRIPTION :
  --   If LU is used, then the system f is first made square
  --   before calling Newton's method.

    use_lu : constant boolean := (method = '1' or method = '2');
 
  begin
    if use_lu then
      declare
        q : constant integer32 := p'last(2);
        sf : constant Poly_Sys := Make_Square(f,natural32(q));
      begin
        Newton(sf,sols,p,method);
      end;
    else
      Newton(f,sols,p,method);
    end if;
  end Call_Newton;

  procedure Test_Newton ( n,d,k : in natural32;
                          ep : in Poly_Sys; esols : in Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    p : constant Poly_Sys := Remove_Embedding1(ep,d);
    s : VecVec(1..integer32(d)) := Slices(ep,d);
    eqs : constant Matrix(1..integer32(d),0..integer32(n))
        := Equations_to_Matrix(s,integer32(n));
    gen : constant Matrix(1..integer32(n),0..integer32(k))
        := Generators(eqs);
    pla : constant Matrix(1..integer32(n),0..integer32(k))
        := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    method : character;

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
   -- put_line("The coefficients of the slices :"); put(eqs,3);
   -- put_line("The parametric representation of the plane : ");
   -- put(pla,3);
    new_line;
    put_line("MENU to test intrinsic Newton's method :");
    put_line("  1. plain LU factorization;");
    put_line("  2. LU factorization with estimate of condition;");
    put_line("  3. least squares using QR factorization;");
    put_line("  4. singular value decomposition;");
    put("Type 1, 2, 3, or 4 to select : ");
    Ask_Alternative(method,"1234");
    Call_Newton(p,isols,pla,method);
    DoblDobl_Complex_VecVecs.Clear(s);
  end Test_Newton;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : natural32 := 0;

  begin
    DoblDobl_Read_Embedding(ep,sols,d);
    n := natural32(ep'last);
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component       : "); put(k,1); new_line;
    Test_Newton(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Newton's method in intrinsic coordinates...");
  Main;
end ts_iddnewt;
