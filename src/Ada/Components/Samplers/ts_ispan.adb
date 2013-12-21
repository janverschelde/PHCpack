with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;    use Standard_Complex_Numbers_Polar;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;           use Continuation_Parameters;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Standard_Point_Coordinates;
with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
with Standard_Intrinsic_Trackers;       use Standard_Intrinsic_Trackers;
with Standard_Linear_Span;              use Standard_Linear_Span;

procedure ts_ispan is

-- DESCRIPTION :
--   Interactive development of computing the span of a component.

  function Random_Offset ( p : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.

    res : Matrix(p'range(1),p'range(2)) := p;

  begin
    for i in p'range(1) loop
      res(i,0) := Random1;
    end loop;
    return res;
  end Random_Offset;

  procedure Expand_and_Append ( first,last : in out Solution_List;
                                s : in Solution; p : in Matrix ) is

  -- DESCRIPTION :
  --   Expands the given solution s with respect to the plane p
  --   and appends the expanded solution to the list in first.

    esol : Solution(p'last(1));

  begin
    esol.t := s.t;
    esol.m := s.m;
    esol.v := Standard_Point_Coordinates.Affine_Expand(s.v,p);
    esol.err := s.err;
    esol.rco := s.rco;
    esol.res := s.res;
    Append(first,last,esol);
  end Expand_and_Append;

  procedure LU_Refine ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                        s : in out Solu_Info;
                        p : in Matrix; c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    nbit : natural32;
    fail : boolean;

  begin
    Affine_LU_Newton
      (f,jf,p,s.sol.v,c.epsax,c.epsrx,c.epsaf,c.epsrf,
       s.cora,s.corr,s.resa,s.resr,nbit,c.maxit,s.rcond,fail);
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
  end LU_Refine;

  procedure SV_Refine ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                        s : in out Solu_Info;
                        p : in Matrix; c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    mm : constant integer32 := Min0(f'last+1,p'last);
    sv : Vector(1..mm);
    nbit : natural32;
    fail : boolean;

  begin
    Affine_SV_Newton
      (f,jf,p,s.sol.v,c.epsax,c.epsrx,c.epsaf,c.epsrf,
       s.cora,s.corr,s.resa,s.resr,nbit,c.maxit,sv,fail);
    s.rcond := Radius(sv(s.sol.v'last)/sv(sv'first));
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
  end SV_Refine;

  function LU_Sample ( f : Eval_Poly_Sys; jf : Eval_Jaco_Mat;
                       sol : Solution; p : Matrix ) return Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

  -- REQUIRED :
  --   The system f defines a complete intersection.

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_Affine_LU_Track(Path);

  begin
    Continuation_Parameters.Tune(2);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(0.0);
      Cont(f,jf,s,pp,cp);
      LU_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end LU_Sample;

  function QR_Sample ( f : Eval_Poly_Sys; jf : Eval_Jaco_Mat;
                       sol : Solution; p : Matrix ) return Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_QR_Track(Path);

  begin
    Continuation_Parameters.Tune(2);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(0.0);
      Cont(f,jf,s,pp,cp);
      SV_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end QR_Sample;

  procedure Sample ( f : in Poly_Sys; s : in Solution;
                     p : in Matrix; method : in character;
		     samples : out Solution_List ) is
  begin
    if method = '1' then
      declare
        sf : Poly_Sys(1..p'last(2)) := Make_Square(f,natural32(p'last(2)));
        sef : Eval_Poly_Sys := Create(sf);
        sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
        sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
      begin
        samples := LU_Sample(sef,sjf,s,p);
        Clear(sf); Clear(sef); Clear(sjm); Clear(sjf); 
      end;
    else
      declare
        ef : Eval_Poly_Sys := Create(f);
        jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
        jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
      begin
        samples := QR_Sample(ef,jf,s,p);
        Clear(ef); Clear(jm); Clear(jf);
      end;
    end if;
  end Sample;

  procedure Evaluate_Samples_in_Span
              ( eqs : in Poly_Sys; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   To verify whether all sampled solutions satisfy the equations
  --   of the linear span, we evaluate the solution vectors.

    y : Vector(eqs'range);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    nrm : double_float;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      y := Eval(eqs,ls.v);
      nrm := Max_Norm(y);
      put("  value of span at sample ");
      put(cnt,1); put(" : "); put(nrm); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Evaluate_Samples_in_Span;

  procedure Determine_Linear_Span
              ( n : in integer32; p : in Poly_Sys; sols : in Solution_List ) is

    v : Matrix(1..n,1..n) := Create(sols);
    tol : constant double_float := 1.0E-8;
    r : natural32;

  begin
    Rank(v,tol,r);
    put("  The linear span has rank "); put(r,1); put_line(".");
    declare
      pv : constant Standard_Integer_Vectors.Vector(1..integer32(r))
         := Pivots(v,tol,r);
      point : constant Vector := Head_Of(sols).v;
      kr : constant Matrix(1..n-integer32(r),0..n)
         := Kernel(v,tol,r,pv,point);
      eqs : Poly_Sys(1..n-integer32(r)) := Equations(kr);
      eli : Poly_Sys(1..n-integer32(r)) := Eliminators(kr,pv);
      elp : Poly_Sys(p'range) := Eliminate_non_Pivots(p,pv,eli);
      flp : Poly_sys := Filter(elp,tol);
      esols : Solution_List := Eliminate_non_Pivots(sols,pv);
    begin
      put_line("The equations of the linear span : "); put(eqs);
      Evaluate_Samples_in_Span(eqs,sols);
      put_line("The equations to eliminate the nonpivot variables : ");
      put(eli);
      put_line("The system with nonpivots eliminated : "); put(flp);
      Evaluate_Samples_in_Span(flp,esols);
      Standard_Complex_Poly_Systems.Clear(eqs);
      Standard_Complex_Poly_Systems.Clear(eli);
      Standard_Complex_Poly_Systems.Clear(elp);
      Standard_Complex_Poly_Systems.Clear(flp);
      Standard_Complex_Solutions.Clear(esols);
    end;
  end Determine_Linear_Span;

  procedure Set_up_Sampling
              ( n,d,k : in integer32;
                ep : in Poly_Sys; esols : in Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the homotopy to call the intrinsic path trackers.

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    method : character;
    tmp,samples : Solution_List;
    ls : Link_to_Solution;

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
    new_line;
    put_line("MENU to test intrinsic operations :");
    put_line("  1. use LU for path tracking with intrinsic Newton;");
    put_line("  2. use QR for path tracking with intrinsic Newton;");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(method,"12");
    tmp := isols;
    for i in 1..Length_Of(isols) loop
      ls := Head_Of(tmp);
      put("Sampling from point "); put(i,1); put_line(" ...");
      Sample(p,ls.all,pla,method,samples);
     -- put_line("The list of samples : ");
     -- put(Standard_Output,Length_Of(samples),Head_Of(samples).n,samples);
      Determine_Linear_Span(n,p,samples);
      tmp := Tail_Of(tmp);
      Deep_Clear(samples);
    end loop;
    Standard_Complex_VecVecs.Clear(s);
    Standard_Complex_Solutions.Clear(isols);
  end Set_up_Sampling;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : integer32 := 0;

  begin
    Standard_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Set_up_Sampling(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Computing span of component in intrinsic coordinates...");
  Main;
end ts_ispan;
