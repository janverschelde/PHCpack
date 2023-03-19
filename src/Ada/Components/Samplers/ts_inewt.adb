with text_io;                           use text_io;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Projective_Transformations;        use Projective_Transformations;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;           use Continuation_Parameters;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;

procedure ts_inewt is

-- DESCRIPTION :
--   Interactive development of intrinsic sampling.

 -- procedure Set_Continuation_Parameter
 --             ( sa : in out Solu_Info_Array; vt : in Complex_Number ) is
 -- begin
 --   for i in sa'range loop
 --     sa(i).sol.t := vt;
 --   end loop;
 -- end Set_Continuation_Parameter;

 -- function Equal ( s1,s2 : Solu_Info; tol : double_float ) return boolean is
 -- begin
 --   for i in s1.sol.v'range loop
 --     if AbsVal(s1.sol.v(i) - s2.sol.v(i)) > tol
 --      then return false;
 --     end if;
 --   end loop;
 --   return true;
 -- end Equal;

  procedure Perturb ( x : in out Vector; tol : in double_float ) is

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
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco : double_float;
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
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco : double_float;
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
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                max : in natural32; fail : out boolean ) is

    n : constant integer32 := f'last;
    q : constant integer32 := p'last(2);
    mm : constant integer32 := Min0(n+1,q);
    sv : Vector(1..mm);
    rco,nrm : double_float;
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

  begin
    for j in p'range(2) loop
      for i in p'range(1) loop
        res(i,j) := p(i,j);
      end loop;
      res(res'last(1),j) := Create(0.0);
    end loop;
    res(res'last(1),0) := Create(1.0);
    return res;
  end Make_Projective;

  procedure Newton ( f : in Poly_Sys; sols : in Solution_List;
                     p : in Matrix; method : in character ) is

    tmp : Solution_List := sols;
    epsax : constant double_float := 1.0E-13;
    epsrx : constant double_float := 1.0E-13;
    epsaf : constant double_float := 1.0E-13;
    epsrf : constant double_float := 1.0E-13;
    incax,incrx,resaf,resrf : double_float;
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
      Perturb(x,0.001);
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

 -- function Special_Plane ( n,k : natural ) return Matrix is
 --
  -- DESCRIPTION :
  --   Special plane for cyclic 9-roots.
 --
 --   res : Matrix(1..n,0..k);
 --
 -- begin
 --   for i in 1..n loop
 --     res(i,0) := Random1;
 --   end loop;
 --   for i in 1..k loop
 --     for j in 1..k loop
 --       res(i,j) := Create(0.0);
 --     end loop;
 --     res(i,i) := Create(1.0);
 --   end loop;
 --   for i in k+1..n loop
 --     for j in 1..k loop
 --       res(i,j) := Random1;
 --     end loop;
 --   end loop;
 --   return res;
 -- end Special_Plane;

 -- function Random_Plane ( p : Matrix; n,k : natural ) return Matrix is
 --
  -- DESCRIPTION :
  --   Returns a matrix representing a k-plane in n-space,
  --   changing only last column of p.
 --
 --  res : Matrix(1..n,0..k);
 --
 -- begin
 --   for i in 1..n loop
 --     for j in 0..k-1 loop
 --       res(i,j) := p(i,j);
 --     end loop;
 --     res(i,k) := Random1;
 --   end loop;
 --   return Standard_Plane_Representations.Orthogonalize(res);
 -- end Random_Plane;

 -- function Random_Offset ( p : Matrix; n,k : natural ) return Matrix is
 --
  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.
 --
 --   res : Matrix(1..n,0..k) := p;
 --
 -- begin
 --   for i in 1..n loop
 --     res(i,0) := Random1;
 --   end loop;
 --   return res;
 -- end Random_Offset;

 -- function Random_Offset1 ( p : Matrix; n,k,i : natural ) return Matrix is
 --
  -- DESCRIPTION :
  --   Returns the k-plane with same directions, but the i-th
  --   component of the new offset is a random number.
 --
 --   res : Matrix(1..n,0..k) := p;
 --
 -- begin
 --   res(i,0) := Random1;
 --   return res;
 -- end Random_Offset1;

 -- function Random_Change ( p : Matrix; n,k,i,j : natural ) return Matrix is
 --
  -- DESCRIPTION :
  --   Returns the k-plane with same directions, but the (i,j)-th
  --   component of the new plane is a random number.
 --
 --   res : Matrix(1..n,0..k) := p;
 --
 -- begin
 --   res(i,j) := Random1;
 --   return res;
 -- end Random_Change;

  function Average_Iterations ( s : Solu_Info_Array ) return natural32 is

  -- DESCRIPTION :
  --   Returns the average number of iterations.

    sum : natural32 := 0;

  begin
    for i in s'range loop
      sum := sum + s(i).niter;
    end loop;
    return sum/natural32(s'last);
  end Average_Iterations;

  procedure Call_Tracker ( f : in Poly_Sys; sols : in out Solution_List;
                           p : in Matrix; method : in character ) is

    sf : constant Poly_Sys := Make_Square(f,natural32(p'last(2))); -- for LU
    sef : Eval_Poly_Sys := Create(sf);
    sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    ef : Eval_Poly_Sys := Create(f);                    -- for QR
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    sa : Solu_Info_Array(1..integer32(Length_Of(sols)));
    rg,sn,rl,cm,cl,fl : natural32;
    ans : character;
    oc : natural32;
    output : boolean;
    gamma : constant Complex_Number := Random1;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,gamma,t);
    end Path;
    procedure R_LU_Cont is new Reporting_Affine_LU_Continue(Path);
    procedure S_LU_Cont is new Silent_Affine_LU_Continue(Path);
    procedure R_QR_Cont is new Reporting_QR_Continue(Path);
    procedure S_QR_Cont is new Silent_QR_Continue(Path);

  begin
    new_line;
    Driver_for_Process_io(Standard_Output,oc);
    output := (oc > 0);
    Continuation_Parameters.Tune(2);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    loop
      tp := Random_Plane(p'last(1),p'last(2));
     -- tp := Random_Offset(p,p'last(1),p'last(2));
     -- tp := Special_Plane(p'last(1),p'last(2));
      sa := Deep_Create(sols);
      if method = '5' then
        if output
         then R_LU_Cont(Standard_Output,sef,sjf,sa,pp,cp);
         else S_LU_Cont(sef,sjf,sa,pp,cp);
        end if;
        LU_Validate(Standard_Output,sef,sjf,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      else
        if output
         then R_QR_Cont(Standard_Output,ef,jf,sa,pp,cp);
         else S_QR_Cont(ef,jf,sa,pp,cp);
        end if;
        SV_Validate(Standard_Output,ef,jf,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      end if;
      put("Average #iterations along a path : ");
      put(Average_Iterations(sa),1); put_line(".");
      Clear(sa);
      put("Do you want another run ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
    Standard_Complex_Poly_SysFun.Clear(sef);
    Standard_Complex_Jaco_Matrices.Clear(sjm);
    Standard_Complex_Jaco_Matrices.Clear(sjf);
  end Call_Tracker;

  procedure Call_Generic_Tracker
              ( f : in Poly_Sys; sols : in out Solution_List;
                p : in Matrix; method : in character ) is

    sf : constant Poly_Sys := Make_Square(f,natural32(p'last(2))); -- for LU
    sef : Eval_Poly_Sys := Create(sf);
    sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    ef : Eval_Poly_Sys := Create(f);                    -- for QR
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    sa : Solu_Info_Array(1..integer32(Length_Of(sols)));
    rg,sn,rl,cm,cl,fl : natural32;
    ne_lu : constant natural32 := natural32(sf'last);  -- #equations for LU
    ne_qr : constant natural32 := natural32(f'last);   -- #equations for QR
    nv : constant natural32 := natural32(p'last(1));   -- #variables in system
    ans : character;
    oc : natural32;
    output : boolean;
    gamma : constant Complex_Number := Random1;

    function Eval_for_LU ( x : Vector ) return Vector is
    begin
      return Eval(sef,x);
    end Eval_for_LU;
    function Eval_for_QR ( x : Vector ) return Vector is
    begin
      return Eval(ef,x);
    end Eval_for_QR;

    function Diff_for_LU ( x : Vector ) return Matrix is
    begin
      return Eval(sjf,x);
    end Diff_for_LU;
    function Diff_for_QR ( x : Vector ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff_for_QR;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,gamma,t);
    end Path;

    procedure R_LU_Cont is
      new G_Reporting_LU_Continue(Eval_for_LU,Diff_for_LU,Path);
    procedure S_LU_Cont is
      new G_Silent_LU_Continue(Eval_for_LU,Diff_for_LU,Path);
    procedure R_QR_Cont is
      new G_Reporting_QR_Continue(Eval_for_QR,Diff_for_QR,Path);
    procedure S_QR_Cont is
      new G_Silent_QR_Continue(Eval_for_QR,Diff_for_QR,Path);

    procedure R_LU_Vali is new Reporting_LU_Validate(Eval_for_LU,Diff_for_LU);
    procedure R_SV_Vali is new Reporting_SV_Validate(Eval_for_QR,Diff_for_QR);

  begin
    new_line;
    Driver_for_Process_io(Standard_Output,oc);
    output := (oc > 0);
    Continuation_Parameters.Tune(2);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    loop
      tp := Random_Plane(p'last(1),p'last(2));
     -- tp := Random_Offset(p,p'last(1),p'last(2));
     -- tp := Special_Plane(p'last(1),p'last(2));
      sa := Deep_Create(sols);
      if method = '5' then
        if output
         then R_LU_Cont(Standard_Output,ne_lu,nv,sa,pp,cp);
         else S_LU_Cont(ne_lu,nv,sa,pp,cp);
        end if;
        R_LU_Vali(Standard_Output,ne_lu,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      else
        if output
         then R_QR_Cont(Standard_Output,ne_qr,nv,sa,pp,cp);
         else S_QR_Cont(ne_qr,nv,sa,pp,cp);
        end if;
        R_SV_Vali(Standard_Output,ne_qr,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      end if;
      put("Average #iterations along a path : ");
      put(Average_Iterations(sa),1); put_line(".");
      Clear(sa);
      put("Do you want another run ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
    Standard_Complex_Poly_SysFun.Clear(sef);
    Standard_Complex_Jaco_Matrices.Clear(sjm);
    Standard_Complex_Jaco_Matrices.Clear(sjf);
  end Call_Generic_Tracker;

  procedure Write_Diagnostics ( rg,sn,cl,fl : in natural32 ) is
  begin
    put("  #regu : "); put(rg,1);
    put("  #sing : "); put(sn,1);
    put("  #clus : "); put(cl,1);
    put("  #fail : "); put(fl,1); new_line;
  end Write_Diagnostics;

  procedure Massive_Test ( f : in Poly_Sys; sols : in out Solution_List;
                           p : in Matrix; method : in character ) is

    sf : constant Poly_Sys := Make_Square(f,natural32(p'last(2))); -- for LU
    sef : Eval_Poly_Sys := Create(sf);
    sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    ef : Eval_Poly_Sys := Create(f);                    -- for QR
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    sa : Solu_Info_Array(1..integer32(Length_Of(sols)));
    rg,sn,rl,cm,cl,fl,N,bad : natural32 := 0;
    gamma : constant Complex_Number := Random1;
    timer : Timing_Widget;

    function Path ( t : Complex_Number ) return Matrix is
    begin
     -- return Rotating_Plane(p,1,2,t); --Moving_Plane(p,tp,t);
      return Moving_Plane(p,tp,gamma,t);
    end Path;
    procedure S_LU_Cont is new Silent_Affine_LU_Continue(Path);
    procedure S_QR_Cont is new Silent_QR_Continue(Path);

  begin
    new_line;
    put("Give the number of tests : "); get(N);
    new_line;
    put("Running "); put(N,1); put_line(" tests ...");
    new_line;
    Continuation_Parameters.Tune(0);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    for i in 1..N loop
      sa := Deep_Create(sols);
      tp := Random_Plane(p'last(1),p'last(2));
     -- tp := Random_Offset(p,p'last(1),p'last(2));
      tstart(timer);
      if method = '7' then
        S_LU_Cont(sef,sjf,sa,pp,cp);
        LU_Validate(sef,sjf,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      else
        S_QR_Cont(ef,jf,sa,pp,cp);
        SV_Validate(ef,jf,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      end if;
      tstop(timer);
      put(i,4); put(" : ");
      Write_Diagnostics(rg,sn,cl,fl);
      if cl > 0 or fl > 0
       then bad := bad + 1;
      end if;
      put("Average #iterations along a path : ");
      put(Average_Iterations(sa),1); put_line(".");
      put("  Elapsed user time : ");
      print_time(standard_output,Elapsed_User_Time(timer)); put_line(".");
      Clear(sa); 
    end loop;
    put("Ran "); put(N,1); put(" random tests and found ");
    put(bad,1); put_line(" bad cases.");
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Poly_SysFun.Clear(sef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
    Standard_Complex_Jaco_Matrices.Clear(sjm);
    Standard_Complex_Jaco_Matrices.Clear(sjf);
  end Massive_Test;

  procedure Generic_Massive_Test
               ( f : in Poly_Sys; sols : in out Solution_List;
                 p : in Matrix; method : in character ) is

    sf : constant Poly_Sys := Make_Square(f,natural32(p'last(2))); -- for LU
    sef : Eval_Poly_Sys := Create(sf);
    sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    ef : Eval_Poly_Sys := Create(f);                    -- for QR
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    ne_lu : constant natural32 := natural32(sf'last);  -- #equations for LU
    ne_qr : constant natural32 := natural32(f'last);   -- #equations for QR
    nv : constant natural32 := natural32(p'last(1));   -- #variables in system
    sa : Solu_Info_Array(1..integer32(Length_Of(sols)));
    rg,sn,rl,cm,cl,fl,N,bad : natural32 := 0;
    gamma : constant Complex_Number := Random1;
    timer : Timing_Widget;

    function Eval_for_LU ( x : Vector ) return Vector is
    begin
      return Eval(sef,x);
    end Eval_for_LU;
    function Eval_for_QR ( x : Vector ) return Vector is
    begin
      return Eval(ef,x);
    end Eval_for_QR;

    function Diff_for_LU ( x : Vector ) return Matrix is
    begin
      return Eval(sjf,x);
    end Diff_for_LU;
    function Diff_for_QR ( x : Vector ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff_for_QR;

    function Path ( t : Complex_Number ) return Matrix is
    begin
     -- return Rotating_Plane(p,1,2,t); --Moving_Plane(p,tp,t);
      return Moving_Plane(p,tp,gamma,t);
    end Path;

    procedure S_LU_Cont is
      new G_Silent_LU_Continue(Eval_for_LU,Diff_for_LU,Path);
    procedure S_QR_Cont is
      new G_Silent_QR_Continue(Eval_for_QR,Diff_for_QR,Path);

    procedure S_LU_Vali is new Silent_LU_Validate(Eval_for_LU,Diff_for_LU);
    procedure S_SV_Vali is new Silent_SV_Validate(Eval_for_QR,Diff_for_QR);

  begin
    new_line;
    put("Give the number of tests : "); get(N);
    new_line;
    put("Running "); put(N,1); put_line(" tests ...");
    new_line;
    Continuation_Parameters.Tune(0);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    for i in 1..N loop
      sa := Deep_Create(sols);
      tp := Random_Plane(p'last(1),p'last(2));
     -- tp := Random_Offset(p,p'last(1),p'last(2));
      tstart(timer);
      if method = '7' then
        S_LU_Cont(ne_lu,nv,sa,pp,cp);
        S_LU_Vali(ne_lu,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      else
        S_QR_Cont(ne_qr,nv,sa,pp,cp);
        S_SV_Vali(ne_qr,tp,sa,ecp,rg,sn,rl,cm,cl,fl);
      end if;
      tstop(timer);
      put(i,4); put(" : ");
      Write_Diagnostics(rg,sn,cl,fl);
      if cl > 0 or fl > 0
       then bad := bad + 1;
      end if;
      put("Average #iterations along a path : ");
      put(Average_Iterations(sa),1); put(".");
      put("  Elapsed user time : ");
      print_time(standard_output,Elapsed_User_Time(timer)); put_line(".");
      Clear(sa); 
    end loop;
    put("Ran "); put(N,1); put(" random tests and found ");
    put(bad,1); put_line(" bad cases.");
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
    Standard_Complex_Poly_SysFun.Clear(sef);
    Standard_Complex_Jaco_Matrices.Clear(sjm);
    Standard_Complex_Jaco_Matrices.Clear(sjf);
  end Generic_Massive_Test;

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
    ans,method : character;

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
   -- put_line("The coefficients of the slices :"); put(eqs,3);
   -- put_line("The parametric representation of the plane : ");
   -- put(pla,3);
    new_line;
    put_line("MENU to test intrinsic operations :");
    put_line("  1. Newton with plain LU factorization;");
    put_line("  2.             LU factorization with estimate of condition;");
    put_line("  3.             least squares using QR factorization;");
    put_line("  4.             singular value decomposition;");
    put_line("  5. path tracking with intrinsic Newton using LU;");
    put_line("  6.                                     using QR + SVD.");
    put_line("  7. massive testing on path trackers using LU;");
    put_line("  8.                                  using QR + SVD.");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to select : ");
    Ask_Alternative(method,"12345678");
    case method is
      when '1' | '2' | '3' | '4' => Call_Newton(p,isols,pla,method);
      when '5' | '6' | '7' | '8' =>
        new_line;
        put("Run the generic version of the trackers ? (y/n) ");
        Ask_Yes_or_No(ans);
        if method = '5' or method = '6' then
          if ans = 'y'
           then Call_Generic_Tracker(p,isols,pla,method);
           else Call_Tracker(p,isols,pla,method);
          end if;
        else
          if ans = 'y'
           then Generic_Massive_Test(p,isols,pla,method);
           else Massive_Test(p,isols,pla,method);
          end if;
        end if;
      when others => null;
    end case;
    Standard_Complex_VecVecs.Clear(s);
  end Test_Newton;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : natural32 := 0;

  begin
    Standard_Read_Embedding(ep,sols,d);
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
end ts_inewt;
