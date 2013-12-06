with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
--with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
--with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
--with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
--with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
--with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
--with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
--with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
--with Projective_Transformations;        use Projective_Transformations;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
--with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;
--with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
--with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
--with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
--with Standard_Intrinsic_Trackers;       use Standard_Intrinsic_Trackers;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;

procedure ts_ifactor is

-- DESCRIPTION :
--   Interactive development of factorization in intrinsic coordinates.

  procedure Set_Continuation_Parameter
              ( sa : in out Solu_Info_Array; vt : in Complex_Number ) is
  begin
    for i in sa'range loop
      sa(i).sol.t := vt;
    end loop;
  end Set_Continuation_Parameter;

  function Equal ( s1,s2 : Solu_Info; tol : double_float ) return boolean is
  begin
    for i in s1.sol.v'range loop
      if AbsVal(s1.sol.v(i) - s2.sol.v(i)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Find_Index ( sa : Solu_Info_Array;
                        s : Solu_Info ) return natural is

    tol : constant double_float := 1.0E-8;

  begin
    for i in sa'range loop
      if Equal(sa(i),s,tol)
       then return i;
      end if;
    end loop;
    return 0;
  end Find_Index;

  procedure Write_Diagnostics ( rg,sn,cl,fl : in natural ) is
  begin
    put("  #regu : "); put(rg,1);
    put("  #sing : "); put(sn,1);
    put("  #clus : "); put(cl,1);
    put("  #fail : "); put(fl,1); new_line;
  end Write_Diagnostics;

  procedure Monodromy_Permutation ( s1,s2 : in Solu_Info_Array )  is

    ind : natural;
    nothing_found : boolean := true;

  begin
    for i in s1'range loop
      ind := Find_Index(s1,s2(i));
      if ind /= i
       then put("  "); put(i,1); put(" --> ");
            put(ind,1); new_line;
            nothing_found := false;
      end if;
    end loop;
    if nothing_found
     then put_line("  no permutations found");
    end if;
  end Monodromy_Permutation;

  procedure Perturb ( x : in out Vector; tol : in double_float ) is

  -- DESCRIPTION :
  --   Perturbs every entry in x with magnitude equal to tol.

    dx : constant Vector(x'first..x'last) := Random_Vector(x'first,x'last);

  begin
    for i in x'range loop
      x(i) := x(i) + tol*dx(i);
    end loop;
  end Perturb;

  function Random_Plane ( n,k : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix representing a k-plane in n-space.

    res : Matrix(1..n,0..k);

  begin
    for i in 1..n loop
      for j in 0..k loop
        res(i,j) := Random1;
      end loop;
    end loop;
    return Standard_Plane_Representations.Orthogonalize(res);
  end Random_Plane;

  function Special_Plane ( n,k : natural ) return Matrix is

  -- DESCRIPTION :
  --   Special plane for cyclic 9-roots.

    res : Matrix(1..n,0..k);

  begin
    for i in 1..n loop
      res(i,0) := Random1;
    end loop;
    for i in 1..k loop
      for j in 1..k loop
        res(i,j) := Create(0.0);
      end loop;
      res(i,i) := Create(1.0);
    end loop;
    for i in k+1..n loop
      for j in 1..k loop
        res(i,j) := Random1;
      end loop;
    end loop;
    return res;
  end Special_Plane;

  function Random_Plane ( p : Matrix; n,k : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix representing a k-plane in n-space,
  --   changing only last column of p.

    res : Matrix(1..n,0..k);

  begin
    for i in 1..n loop
      for j in 0..k-1 loop
        res(i,j) := p(i,j);
      end loop;
      res(i,k) := Random1;
    end loop;
    return Standard_Plane_Representations.Orthogonalize(res);
  end Random_Plane;

  function Random_Offset ( p : Matrix; n,k : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.

    res : Matrix(1..n,0..k) := p;

  begin
    for i in 1..n loop
      res(i,0) := Random1;
    end loop;
    return res;
  end Random_Offset;

  function Random_Offset1 ( p : Matrix; n,k,i : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, but the i-th
  --   component of the new offset is a random number.

    res : Matrix(1..n,0..k) := p;

  begin
    res(i,0) := Random1;
    return res;
  end Random_Offset1;

  function Random_Change ( p : Matrix; n,k,i,j : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, but the (i,j)-th
  --   component of the new plane is a random number.

    res : Matrix(1..n,0..k) := p;

  begin
    res(i,j) := Random1;
    return res;
  end Random_Change;

  procedure Monodromy_Test1 ( f : in Poly_Sys; sols : in out Solution_List;
                              p : in Matrix; method : in character ) is

    ef : constant Eval_Poly_Sys := Create(f);
    jm : constant Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : constant Eval_Jaco_Mat := Create(jm);
    start,target,tp1,tp2 : Matrix(p'range(1),p'range(2));
    pp : Continuation_Parameters.Pred_Pars;
    cp,ecp : Continuation_Parameters.Corr_Pars;
    sa0,sa1 : Solu_Info_Array(1..Length_Of(sols));
    rg,sn,rl,cm,cl,fl,N,j1,j2 : natural := 0;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start,target,t);
    end Path;
    procedure S_LU_Cont is new Silent_Affine_LU_Continue(Path);
    procedure S_QR_Cont is new Silent_QR_Continue(Path);

  begin
    new_line;
    put("Give number of loops : "); get(N);
    Continuation_Parameters.Tune(0);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    sa0 := Deep_Create(sols);
    put_line("running the monodromy...");
    j1 := 1; j2 := 2;
    for i in 1..N loop
      sa1 := Deep_Create(sols);
      tp1 := Random_Offset1(p,p'last(1),p'last(2),j1);
      tp2 := Random_Offset1(p,p'last(1),p'last(2),j2);
      if method = '1'
       then start := p; target := tp1;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp1; target := tp2;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp2; target := p;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
       else start := p; target := tp1;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp1; target := tp2;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp2; target := p;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
      end if;
      Monodromy_Permutation(sa0,sa1);
      Clear(sa1);
      if j2 < p'last(1)
       then j2 := j2+1;
       elsif j1 < j2-1
           then j1 := j1+1;
           else j1 := 1;
                j2 := 2;
      end if;
    end loop;
  end Monodromy_Test1;

  procedure Monodromy_Test2 ( f : in Poly_Sys; sols : in out Solution_List;
                              p : in Matrix; method : in character ) is

    ef : constant Eval_Poly_Sys := Create(f);
    jm : constant Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : constant Eval_Jaco_Mat := Create(jm);
    start,target,tp1,tp2 : Matrix(p'range(1),p'range(2));
    pp : Continuation_Parameters.Pred_Pars;
    cp,ecp : Continuation_Parameters.Corr_Pars;
    sa0,sa1 : Solu_Info_Array(1..Length_Of(sols));
    rg,sn,rl,cm,cl,fl,N : natural := 0;
   -- i1,j1,i2,j2 : natural := 0;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start,target,t);
    end Path;
    procedure S_LU_Cont is new Silent_Affine_LU_Continue(Path);
    procedure S_QR_Cont is new Silent_QR_Continue(Path);

  begin
    new_line;
    put("Give number of loops : "); get(N);
    Continuation_Parameters.Tune(0);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    sa0 := Deep_Create(sols);
    put_line("running the monodromy...");
   -- j1 := 1; j2 := 2;
    for i in 1..N loop
     -- i1 := Random(1,p'last(1)); j1 := Random(0,p'last(2));
     -- i2 := Random(1,p'last(1)); j2 := Random(0,p'last(2));
     -- put("monodromy loop changing ");
     -- put("("); put(i1,1); put(","); put(j1,1); put(") and ");
     -- put("("); put(i2,1); put(","); put(j2,1); put_line(").");
      put("Executing loop "); put(i,1); put_line(" ...");
      sa1 := Deep_Create(sols);
     -- tp1 := Random_Change(p,p'last(1),p'last(2),i1,j1);
     -- tp2 := Random_Change(p,p'last(1),p'last(2),i2,j2);
      tp1 := Random_Plane(p'last(1),p'last(2));
      tp2 := Random_Plane(p'last(1),p'last(2));
      if method = '1'
       then start := p; target := tp1;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp1; target := tp2;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp2; target := p;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
       else start := p; target := tp1;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp1; target := tp2;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp2; target := p;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
      end if;
      Monodromy_Permutation(sa0,sa1);
      Clear(sa1);
    end loop;
  end Monodromy_Test2;

  procedure Monodromy_Test ( f : in Poly_Sys; sols : in out Solution_List;
                             p : in Matrix; method : in character ) is

    ef : constant Eval_Poly_Sys := Create(f);
    jm : constant Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : constant Eval_Jaco_Mat := Create(jm);
    start,target,tp : Matrix(p'range(1),p'range(2));
    pp : Continuation_Parameters.Pred_Pars;
    cp,ecp : Continuation_Parameters.Corr_Pars;
    sa0,sa1 : Solu_Info_Array(1..Length_Of(sols));
    rg,sn,rl,cm,cl,fl,N : natural := 0;
    gamma : Complex_Number;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start,target,gamma,t);
    end Path;
    procedure S_LU_Cont is new Silent_Affine_LU_Continue(Path);
    procedure S_QR_Cont is new Silent_QR_Continue(Path);

  begin
    new_line;
    put("Give number of loops : "); get(N);
    Continuation_Parameters.Tune(0);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Set_Continuation_Parameter(sols,Create(0.0));
    sa0 := Deep_Create(sols);
    put_line("running the monodromy...");
    for i in 1..N loop
      put("Executing loop "); put(i,1); put_line(" ...");
      sa1 := Deep_Create(sols);
      tp := Random_Plane(p'last(1),p'last(2));
      if method = '1'
       then start := p; target := tp; gamma := Random1;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp; target := p; gamma := Random1;
            S_LU_Cont(ef,jf,sa1,pp,cp);
            LU_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
       else start := p; target := tp; gamma := Random1;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
            Set_Continuation_Parameter(sa1,Create(0.0));
            start := tp; target := p; gamma := Random1;
            S_QR_Cont(ef,jf,sa1,pp,cp);
            SV_Validate(ef,jf,target,sa1,ecp,rg,sn,rl,cm,cl,fl);
            Write_Diagnostics(rg,sn,cl,fl);
      end if;
      Monodromy_Permutation(sa0,sa1);
      Clear(sa1);
    end loop;
  end Monodromy_Test;

  procedure Test_Newton ( n,d,k : in natural;
                          ep : in Poly_Sys; esols : in Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    p : constant Poly_Sys := Witness_Sets.Remove_Embedding1(ep,d);
    s : constant VecVec(1..d) := Slices(ep,d);
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
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
    put_line("MENU to test intrinsic operations :");
    put_line("  1. monodromy test with continuation using LU;");
    put_line("  2.                                  using QR + SVD.");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(method,"12");
    Monodromy_Test(p,isols,pla,method);
  end Test_Newton;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : natural;

  begin
    Standard_Read_Embedding(ep,sols,d);
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Test_Newton(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Newton's method in intrinsic coordinates...");
  Main;
end ts_ifactor;
