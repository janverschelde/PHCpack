with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;    use Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;    use DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;    use QuadDobl_Complex_Numbers_Polar;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with DoblDobl_Random_Numbers;           use DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Continuation_Data;
with DoblDobl_Continuation_Data;
with QuadDobl_Continuation_Data;
with Continuation_Parameters;           use Continuation_Parameters;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Standard_Point_Coordinates;
with DoblDobl_Plane_Representations;    use DoblDobl_Plane_Representations;
with DoblDobl_Moving_Planes;            use DoblDobl_Moving_Planes;
with DoblDobl_Intrinsic_Solutions;      use DoblDobl_Intrinsic_Solutions;
with DoblDobl_Point_Coordinates;
with QuadDobl_Plane_Representations;    use QuadDobl_Plane_Representations;
with QuadDobl_Moving_Planes;            use QuadDobl_Moving_Planes;
with QuadDobl_Intrinsic_Solutions;      use QuadDobl_Intrinsic_Solutions;
with QuadDobl_Point_Coordinates;
with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
with Standard_Intrinsic_Trackers;
with DoblDobl_Intrinsic_Newton;         use DoblDobl_Intrinsic_Newton;
with DoblDobl_Intrinsic_Trackers;
with QuadDobl_Intrinsic_Newton;         use QuadDobl_Intrinsic_Newton;
with QuadDobl_Intrinsic_Trackers;
with Standard_Linear_Span;              use Standard_Linear_Span;
with DoblDobl_Linear_Span;              use DoblDobl_Linear_Span;
with QuadDobl_Linear_Span;              use QuadDobl_Linear_Span;

procedure ts_iwset is

-- DESCRIPTION :
--   Interactive development of computing witness set using its linear span.

  function Coefficients
             ( h : Standard_Complex_VecVecs.VecVec; n : integer32 )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    use Standard_Complex_VecVecs;
    use Standard_Complex_Matrices;

    res : Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( h : DoblDobl_Complex_VecVecs.VecVec; n : integer32 )
             return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Matrices;

    res : Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( h : QuadDobl_Complex_VecVecs.VecVec; n : integer32 )
             return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Matrices;

    res : Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Random_Offset
             ( p : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.

    use Standard_Complex_Matrices;

    res : Matrix(p'range(1),p'range(2)) := p;

  begin
    for i in p'range(1) loop
      res(i,0) := Random1;
    end loop;
    return res;
  end Random_Offset;

  function Random_Offset
             ( p : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.

    use DoblDobl_Complex_Matrices;

    res : Matrix(p'range(1),p'range(2)) := p;

  begin
    for i in p'range(1) loop
      res(i,0) := Random1;
    end loop;
    return res;
  end Random_Offset;

  function Random_Offset
             ( p : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the k-plane with same directions, with random offset.

    use QuadDobl_Complex_Matrices;

    res : Matrix(p'range(1),p'range(2)) := p;

  begin
    for i in p'range(1) loop
      res(i,0) := Random1;
    end loop;
    return res;
  end Random_Offset;

  procedure Expand_and_Append
              ( first,last : in out Standard_Complex_Solutions.Solution_List;
                s : in Standard_Complex_Solutions.Solution;
                p : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Expands the given solution s with respect to the plane p
  --   and appends the expanded solution to the list in first.

    use Standard_Complex_Solutions;

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

  procedure Expand_and_Append
              ( first,last : in out DoblDobl_Complex_Solutions.Solution_List;
                s : in DoblDobl_Complex_Solutions.Solution;
                p : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Expands the given solution s with respect to the plane p
  --   and appends the expanded solution to the list in first.

    use DoblDobl_Complex_Solutions;

    esol : Solution(p'last(1));

  begin
    esol.t := s.t;
    esol.m := s.m;
    esol.v := DoblDobl_Point_Coordinates.Affine_Expand(s.v,p);
    esol.err := s.err;
    esol.rco := s.rco;
    esol.res := s.res;
    Append(first,last,esol);
  end Expand_and_Append;

  procedure Expand_and_Append
              ( first,last : in out QuadDobl_Complex_Solutions.Solution_List;
                s : in QuadDobl_Complex_Solutions.Solution;
                p : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Expands the given solution s with respect to the plane p
  --   and appends the expanded solution to the list in first.

    use QuadDobl_Complex_Solutions;

    esol : Solution(p'last(1));

  begin
    esol.t := s.t;
    esol.m := s.m;
    esol.v := QuadDobl_Point_Coordinates.Affine_Expand(s.v,p);
    esol.err := s.err;
    esol.rco := s.rco;
    esol.res := s.res;
    Append(first,last,esol);
  end Expand_and_Append;

  procedure LU_Refine
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in out Standard_Continuation_Data.Solu_Info;
                p : in Standard_Complex_Matrices.Matrix;
                c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    nbit : natural32;
    fail : boolean;

  begin
    Affine_LU_Newton(f,jf,p,s.sol.v,c.epsax,c.epsrx,c.epsaf,c.epsrf,
                     s.cora,s.corr,s.resa,s.resr,nbit,c.maxit,s.rcond,fail);
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
  end LU_Refine;

  procedure LU_Refine
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in out DoblDobl_Continuation_Data.Solu_Info;
                p : in DoblDobl_Complex_Matrices.Matrix;
                c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    nbit : natural32;
    fail : boolean;
    scora,scorr,sresa,sresr,srcond : double_double;

  begin
    Affine_LU_Newton
      (f,jf,p,s.sol.v,
       create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
       scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
    s.cora := to_double(scora);
    s.corr := to_double(scorr);
    s.resa := to_double(sresa);
    s.resr := to_double(sresr);
    s.rcond := to_double(srcond);
    s.sol.err := scora;
    s.sol.rco := srcond;
    s.sol.res := sresa;
  end LU_Refine;

  procedure LU_Refine
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in out QuadDobl_Continuation_Data.Solu_Info;
                p : in QuadDobl_Complex_Matrices.Matrix;
                c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    nbit : natural32;
    fail : boolean;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    Affine_LU_Newton
      (f,jf,p,s.sol.v,
       create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
       scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
    s.cora := to_double(scora);
    s.corr := to_double(scorr);
    s.resa := to_double(sresa);
    s.resr := to_double(sresr);
    s.rcond := to_double(srcond);
    s.sol.err := scora;
    s.sol.rco := srcond;
    s.sol.res := sresa;
  end LU_Refine;

  procedure SV_Refine
             ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               s : in out Standard_Continuation_Data.Solu_Info;
               p : in Standard_Complex_Matrices.Matrix; c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    use Standard_Complex_Numbers;

    mm : constant integer32 := Min0(f'last+1,p'last);
    sv : Standard_Complex_Vectors.Vector(1..mm);
    nbit : natural32;
    fail : boolean;

  begin
    Affine_SV_Newton(f,jf,p,s.sol.v,c.epsax,c.epsrx,c.epsaf,c.epsrf,
                     s.cora,s.corr,s.resa,s.resr,nbit,c.maxit,sv,fail);
    s.rcond := Radius(sv(s.sol.v'last)/sv(sv'first));
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
  end SV_Refine;

  procedure SV_Refine
             ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               s : in out DoblDobl_Continuation_Data.Solu_Info;
               p : in DoblDobl_Complex_Matrices.Matrix; c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    use DoblDobl_Complex_Numbers;

    mm : constant integer32 := Min0(f'last+1,p'last);
    sv : DoblDobl_Complex_Vectors.Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    scora,scorr,sresa,sresr,srcond : double_double;

  begin
    Affine_SV_Newton
      (f,jf,p,s.sol.v,
       create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
       scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
    srcond := Radius(sv(s.sol.v'last)/sv(sv'first));
    s.cora := to_double(scora);
    s.corr := to_double(scorr);
    s.resa := to_double(sresa);
    s.resr := to_double(sresr);
    s.sol.err := scora;
    s.sol.rco := srcond;
    s.sol.res := sresa;
  end SV_Refine;

  procedure SV_Refine
             ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               s : in out QuadDobl_Continuation_Data.Solu_Info;
               p : in QuadDobl_Complex_Matrices.Matrix; c : in Corr_Pars ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton steps to refine the root.

    use QuadDobl_Complex_Numbers;

    mm : constant integer32 := Min0(f'last+1,p'last);
    sv : QuadDobl_Complex_Vectors.Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    Affine_SV_Newton
      (f,jf,p,s.sol.v,
       create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
       scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
    srcond := Radius(sv(s.sol.v'last)/sv(sv'first));
    s.cora := to_double(scora);
    s.corr := to_double(scorr);
    s.resa := to_double(sresa);
    s.resr := to_double(sresr);
    s.sol.err := scora;
    s.sol.rco := srcond;
    s.sol.res := sresa;
  end SV_Refine;

  function LU_Sample
             ( f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : Standard_Complex_Solutions.Solution;
               p : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

  -- REQUIRED :
  --   The system f defines a complete intersection.

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Intrinsic_Trackers;

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

  function LU_Sample
             ( f : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : DoblDobl_Complex_Solutions.Solution;
               p : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

  -- REQUIRED :
  --   The system f defines a complete intersection.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Intrinsic_Trackers;

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;
    zero : constant double_double := create(0.0);

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_Affine_LU_Track(Path);

  begin
    Continuation_Parameters.Tune(2,32);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(zero);
      Cont(f,jf,s,pp,cp);
      LU_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end LU_Sample;

  function LU_Sample
             ( f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : QuadDobl_Complex_Solutions.Solution;
               p : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

  -- REQUIRED :
  --   The system f defines a complete intersection.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Intrinsic_Trackers;

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;
    zero : constant quad_double := create(0.0);

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_Affine_LU_Track(Path);

  begin
    Continuation_Parameters.Tune(2,64);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(zero);
      Cont(f,jf,s,pp,cp);
      LU_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end LU_Sample;

  function QR_Sample
             ( f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : Standard_Complex_Solutions.Solution;
               p : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Intrinsic_Trackers;

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

  function QR_Sample
             ( f : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : DoblDobl_Complex_Solutions.Solution;
               p : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Intrinsic_Trackers;

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;
    zero : constant double_double := create(0.0);

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_QR_Track(Path);

  begin
    Continuation_Parameters.Tune(2,32);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(zero);
      Cont(f,jf,s,pp,cp);
      SV_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end QR_Sample;

  function QR_Sample
             ( f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
               jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
               sol : QuadDobl_Complex_Solutions.Solution;
               p : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns n+1 samples starting at the given solution.  The first
  --   sample in the list on return in the expanded given solution sol.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Intrinsic_Trackers;

    res,res_last : Solution_List;
    tp : Matrix(p'range(1),p'range(2));
    pp : Pred_Pars;
    cp,ecp : Corr_Pars;
    s : Solu_Info;
    zero : constant quad_double := create(0.0);

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(p,tp,t);
    end Path;
    procedure Cont is new Silent_QR_Track(Path);

  begin
    Continuation_Parameters.Tune(2,64);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    ecp := Continuation_Parameters.Create_End_Game;
    Expand_and_Append(res,res_last,sol,p);
    for i in p'range(1) loop
      tp := Random_Offset(p);
      s := Deep_Create(sol);
      s.sol.t := Create(zero);
      Cont(f,jf,s,pp,cp);
      SV_Refine(f,jf,s,tp,ecp);
      Expand_and_Append(res,res_last,s.sol.all,tp);
      Clear(s);
    end loop;
    return res;
  end QR_Sample;

  function Sample
             ( f : Standard_Complex_Poly_Systems.Poly_Sys;
               s : Standard_Complex_Solutions.Solution;
               p : Standard_Complex_Matrices.Matrix; method : character )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

    res : Solution_List;

  begin
    if method = '1' then
      declare
        sf : Poly_Sys(1..p'last(2)) := Make_Square(f,natural32(p'last(2)));
        sef : Eval_Poly_Sys := Create(sf);
        sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
        sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
      begin
        res := LU_Sample(sef,sjf,s,p);
        Clear(sf); Clear(sef); Clear(sjm); Clear(sjf); 
      end;
    else
      declare
        ef : Eval_Poly_Sys := Create(f);
        jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
        jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
      begin
        res := QR_Sample(ef,jf,s,p);
        Clear(ef); Clear(jm); Clear(jf);
      end;
    end if;
    return res;
  end Sample;

  function Sample
             ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : DoblDobl_Complex_Solutions.Solution;
               p : DoblDobl_Complex_Matrices.Matrix; method : character )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    res : Solution_List;

  begin
    if method = '1' then
      declare
        sf : Poly_Sys(1..p'last(2)) := Make_Square(f,natural32(p'last(2)));
        sef : Eval_Poly_Sys := Create(sf);
        sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
        sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
      begin
        res := LU_Sample(sef,sjf,s,p);
        Clear(sf); Clear(sef); Clear(sjm); Clear(sjf); 
      end;
    else
      declare
        ef : Eval_Poly_Sys := Create(f);
        jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
        jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
      begin
        res := QR_Sample(ef,jf,s,p);
        Clear(ef); Clear(jm); Clear(jf);
      end;
    end if;
    return res;
  end Sample;

  function Sample
             ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : QuadDobl_Complex_Solutions.Solution;
               p : QuadDobl_Complex_Matrices.Matrix; method : character )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    res : Solution_List;

  begin
    if method = '1' then
      declare
        sf : Poly_Sys(1..p'last(2)) := Make_Square(f,natural32(p'last(2)));
        sef : Eval_Poly_Sys := Create(sf);
        sjm : Jaco_Mat(sf'range,p'range(1)) := Create(sf);
        sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
      begin
        res := LU_Sample(sef,sjf,s,p);
        Clear(sf); Clear(sef); Clear(sjm); Clear(sjf); 
      end;
    else
      declare
        ef : Eval_Poly_Sys := Create(f);
        jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
        jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
      begin
        res := QR_Sample(ef,jf,s,p);
        Clear(ef); Clear(jm); Clear(jf);
      end;
    end if;
    return res;
  end Sample;

  procedure Evaluate_in_Span
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Vectors.Vector; ind : in natural32;
                nrm : out double_float ) is

  -- DESCRIPTION :
  --   To verify whether a solution satisfies the equations
  --   of the linear span, we evaluate its vector.

    use Standard_Complex_Poly_SysFun;

    y : constant Standard_Complex_Vectors.Vector(eqs'range) := Eval(eqs,v);

  begin
    nrm := Max_Norm(y);
    put(file,"  value of span at solution ");
    put(file,ind,1); put(file," : "); put(file,nrm); new_line(file);
  end Evaluate_in_Span;

  procedure Evaluate_in_Span
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Vectors.Vector; ind : in natural32;
                nrm : out double_double ) is

  -- DESCRIPTION :
  --   To verify whether a solution satisfies the equations
  --   of the linear span, we evaluate its vector.

    use DoblDobl_Complex_Poly_SysFun;

    y : constant DoblDobl_Complex_Vectors.Vector(eqs'range) := Eval(eqs,v);

  begin
    nrm := Max_Norm(y);
    put(file,"  value of span at solution ");
    put(file,ind,1); put(file," : "); put(file,nrm); new_line(file);
  end Evaluate_in_Span;

  procedure Evaluate_in_Span
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Vectors.Vector; ind : in natural32;
                nrm : out quad_double ) is

  -- DESCRIPTION :
  --   To verify whether a solution satisfies the equations
  --   of the linear span, we evaluate its vector.

    use QuadDobl_Complex_Poly_SysFun;

    y : constant QuadDobl_Complex_Vectors.Vector(eqs'range) := Eval(eqs,v);

  begin
    nrm := Max_Norm(y);
    put(file,"  value of span at solution ");
    put(file,ind,1); put(file," : "); put(file,nrm); new_line(file);
  end Evaluate_in_Span;

  procedure Evaluate_Samples_in_Span
              ( file : in file_type; n : in integer32;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   To verify whether all sampled solutions satisfy the equations
  --   of the linear span, we evaluate the solution vectors.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    nrm : double_float;
    cnt : natural32 := 0;
    sum : double_float := 0.0;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      Evaluate_in_Span(file,eqs,ls.v(1..n),cnt,nrm);
      sum := sum + nrm;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"sum of norms of values : "); put(file,sum); new_line(file);
  end Evaluate_Samples_in_Span;

  procedure Evaluate_Samples_in_Span
              ( file : in file_type; n : in integer32;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   To verify whether all sampled solutions satisfy the equations
  --   of the linear span, we evaluate the solution vectors.

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    nrm : double_double;
    cnt : natural32 := 0;
    sum : double_double := create(0.0);

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      Evaluate_in_Span(file,eqs,ls.v(1..n),cnt,nrm);
      sum := sum + nrm;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"sum of norms of values : "); put(file,sum); new_line(file);
  end Evaluate_Samples_in_Span;

  procedure Evaluate_Samples_in_Span
              ( file : in file_type; n : in integer32;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   To verify whether all sampled solutions satisfy the equations
  --   of the linear span, we evaluate the solution vectors.

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    nrm : quad_double;
    cnt : natural32 := 0;
    sum : quad_double := create(0.0);

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      Evaluate_in_Span(file,eqs,ls.v(1..n),cnt,nrm);
      sum := sum + nrm;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"sum of norms of values : "); put(file,sum); new_line(file);
  end Evaluate_Samples_in_Span;

  procedure Determine_Linear_Span
              ( file : in file_type; n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                samples : in Standard_Complex_Solutions.Solution_List;
                tol : in double_float;
                v : in Standard_Complex_Matrices.Matrix; r : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                i : in natural32 ) is

    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Complex_Poly_Systems;

    pv : constant Standard_Integer_Vectors.Vector(1..integer32(r))
       := Pivots(v,tol,r);
    point : constant Standard_Complex_Vectors.Vector := Head_Of(samples).v;
    kr : constant Matrix(1..n-integer32(r),0..n) := Kernel(v,tol,r,pv,point);
    eqs : Poly_Sys(1..n-integer32(r)) := Equations(kr);
    eli : Poly_Sys(1..n-integer32(r)) := Eliminators(kr,pv);
    elp : Poly_Sys(p'range) := Eliminate_non_Pivots(p,pv,eli);
    flp : constant Poly_sys := Filter(elp,tol);
    esols : Solution_List := Eliminate_non_Pivots(samples,pv);
    tmp : Solution_List := Tail_Of(sols);
    ls : Link_to_Solution;
    ind : natural32 := i;
    nrm : double_float;
    cnt : natural32 := 1;

  begin
    put_line(file,"The equations of the linear span : "); put(file,eqs);
    Evaluate_Samples_in_Span(file,n,eqs,samples);
    put_line(file,"The equations to eliminate the nonpivot variables : ");
    put(file,eli);
    put_line(file,"The system with nonpivots eliminated : "); put(file,flp);
    Evaluate_Samples_in_Span(file,n,flp,esols);
    put_line(file,"Verifying the other witness points ...");
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ind := ind + 1;
      if ls.m = 0 then
        Evaluate_in_Span(file,eqs,ls.v(1..n),ind,nrm);
        if nrm < tol then
          ls.m := integer32(i);
          Set_Head(tmp,ls);
          cnt := cnt + 1;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"Span of point "); put(file,i,1); 
    put(file," contains "); put(file,cnt,1); put_line(file," points.");
    Standard_Complex_Solutions.Clear(esols);
    Standard_Complex_Poly_Systems.Clear(eqs);
    Standard_Complex_Poly_Systems.Clear(eli);
    Standard_Complex_Poly_Systems.Clear(elp);
  end Determine_Linear_Span;

  procedure Determine_Linear_Span
              ( file : in file_type; n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                samples : in DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float;
                v : in DoblDobl_Complex_Matrices.Matrix; r : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                i : in natural32 ) is

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Poly_Systems;

    pv : constant Standard_Integer_Vectors.Vector(1..integer32(r))
       := Pivots(v,tol,r);
    point : constant DoblDobl_Complex_Vectors.Vector := Head_Of(samples).v;
    kr : constant Matrix(1..n-integer32(r),0..n) := Kernel(v,tol,r,pv,point);
    eqs : Poly_Sys(1..n-integer32(r)) := Equations(kr);
    eli : Poly_Sys(1..n-integer32(r)) := Eliminators(kr,pv);
    elp : Poly_Sys(p'range) := Eliminate_non_Pivots(p,pv,eli);
    flp : constant Poly_sys := Filter(elp,tol);
    esols : Solution_List := Eliminate_non_Pivots(samples,pv);
    tmp : Solution_List := Tail_Of(sols);
    ls : Link_to_Solution;
    ind : natural32 := i;
    nrm : double_double;
    cnt : natural32 := 1;

  begin
    put_line(file,"The equations of the linear span : "); put(file,eqs);
    Evaluate_Samples_in_Span(file,n,eqs,samples);
    put_line(file,"The equations to eliminate the nonpivot variables : ");
    put(file,eli);
    put_line(file,"The system with nonpivots eliminated : "); put(file,flp);
    Evaluate_Samples_in_Span(file,n,flp,esols);
    put_line(file,"Verifying the other witness points ...");
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ind := ind + 1;
      if ls.m = 0 then
        Evaluate_in_Span(file,eqs,ls.v(1..n),ind,nrm);
        if nrm < tol then
          ls.m := integer32(i);
          Set_Head(tmp,ls);
          cnt := cnt + 1;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"Span of point "); put(file,i,1); 
    put(file," contains "); put(file,cnt,1); put_line(file," points.");
    DoblDobl_Complex_Solutions.Clear(esols);
    DoblDobl_Complex_Poly_Systems.Clear(eqs);
    DoblDobl_Complex_Poly_Systems.Clear(eli);
    DoblDobl_Complex_Poly_Systems.Clear(elp);
  end Determine_Linear_Span;

  procedure Determine_Linear_Span
              ( file : in file_type; n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                samples : in QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float;
                v : in QuadDobl_Complex_Matrices.Matrix; r : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                i : in natural32 ) is

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Poly_Systems;

    pv : constant Standard_Integer_Vectors.Vector(1..integer32(r))
       := Pivots(v,tol,r);
    point : constant QuadDobl_Complex_Vectors.Vector := Head_Of(samples).v;
    kr : constant Matrix(1..n-integer32(r),0..n) := Kernel(v,tol,r,pv,point);
    eqs : Poly_Sys(1..n-integer32(r)) := Equations(kr);
    eli : Poly_Sys(1..n-integer32(r)) := Eliminators(kr,pv);
    elp : Poly_Sys(p'range) := Eliminate_non_Pivots(p,pv,eli);
    flp : constant Poly_sys := Filter(elp,tol);
    esols : Solution_List := Eliminate_non_Pivots(samples,pv);
    tmp : Solution_List := Tail_Of(sols);
    ls : Link_to_Solution;
    ind : natural32 := i;
    nrm : quad_double;
    cnt : natural32 := 1;

  begin
    put_line(file,"The equations of the linear span : "); put(file,eqs);
    Evaluate_Samples_in_Span(file,n,eqs,samples);
    put_line(file,"The equations to eliminate the nonpivot variables : ");
    put(file,eli);
    put_line(file,"The system with nonpivots eliminated : "); put(file,flp);
    Evaluate_Samples_in_Span(file,n,flp,esols);
    put_line(file,"Verifying the other witness points ...");
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ind := ind + 1;
      if ls.m = 0 then
        Evaluate_in_Span(file,eqs,ls.v(1..n),ind,nrm);
        if nrm < tol then
          ls.m := integer32(i);
          Set_Head(tmp,ls);
          cnt := cnt + 1;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"Span of point "); put(file,i,1); 
    put(file," contains "); put(file,cnt,1); put_line(file," points.");
    QuadDobl_Complex_Solutions.Clear(esols);
    QuadDobl_Complex_Poly_Systems.Clear(eqs);
    QuadDobl_Complex_Poly_Systems.Clear(eli);
    QuadDobl_Complex_Poly_Systems.Clear(elp);
  end Determine_Linear_Span;

  procedure Span_of_Component
              ( file : in file_type; n : in integer32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                p : in Standard_Complex_Matrices.Matrix;
                s : in Standard_Complex_Solutions.Solution;
                method : in character; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                i : in natural32 ) is

  -- DESCRIPTION :
  --   Determines the linear span of the first element in sols
  --   (which happens to be the i-th witness point) and scans
  --   the list sols to see which points belong to the linear span.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space;
  --   f        given polynomial system;
  --   p        plane intersection solution component;
  --   s        current solution to determine whether zero;
  --   method   1 for LU, 2 for QU;
  --   tol      tolerance to decide whether two floats are same;
  --   sols     solution list in extrinsic format;
  --   i        index of the current solution;

  -- ON RETURN :
  --   sols     updated m field connects points in same linear span.

    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;

    samples : Solution_List := Sample(f,s,p,method);
    v : Matrix(1..n,1..n) := Create(samples);
    r : natural32;

  begin
    Rank(v,tol,r);
    put(file,"The linear span has rank ");
    put(file,r,1); put_line(file,".");
    Determine_Linear_Span(file,n,f,samples,tol,v,r,sols,i);
    Deep_Clear(samples);
  end Span_of_Component;

  procedure Span_of_Component
              ( file : in file_type; n : in integer32;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                p : in DoblDobl_Complex_Matrices.Matrix;
                s : in DoblDobl_Complex_Solutions.Solution;
                method : in character; tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                i : in natural32 ) is

  -- DESCRIPTION :
  --   Determines the linear span of the first element in sols
  --   (which happens to be the i-th witness point) and scans
  --   the list sols to see which points belong to the linear span.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space;
  --   f        given polynomial system;
  --   p        plane intersection solution component;
  --   s        current solution to determine whether zero;
  --   method   1 for LU, 2 for QU;
  --   tol      tolerance to decide whether two floats are same;
  --   sols     solution list in extrinsic format;
  --   i        index of the current solution;

  -- ON RETURN :
  --   sols     updated m field connects points in same linear span.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;

    samples : Solution_List := Sample(f,s,p,method);
    v : Matrix(1..n,1..n) := Create(samples);
    r : natural32;

  begin
    Rank(v,tol,r);
    put(file,"The linear span has rank ");
    put(file,r,1); put_line(file,".");
    Determine_Linear_Span(file,n,f,samples,tol,v,r,sols,i);
    Deep_Clear(samples);
  end Span_of_Component;

  procedure Span_of_Component
              ( file : in file_type; n : in integer32;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                p : in QuadDobl_Complex_Matrices.Matrix;
                s : in QuadDobl_Complex_Solutions.Solution;
                method : in character; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                i : in natural32 ) is

  -- DESCRIPTION :
  --   Determines the linear span of the first element in sols
  --   (which happens to be the i-th witness point) and scans
  --   the list sols to see which points belong to the linear span.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space;
  --   f        given polynomial system;
  --   p        plane intersection solution component;
  --   s        current solution to determine whether zero;
  --   method   1 for LU, 2 for QU;
  --   tol      tolerance to decide whether two floats are same;
  --   sols     solution list in extrinsic format;
  --   i        index of the current solution;

  -- ON RETURN :
  --   sols     updated m field connects points in same linear span.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;

    samples : Solution_List := Sample(f,s,p,method);
    v : Matrix(1..n,1..n) := Create(samples);
    r : natural32;

  begin
    Rank(v,tol,r);
    put(file,"The linear span has rank ");
    put(file,r,1); put_line(file,".");
    Determine_Linear_Span(file,n,f,samples,tol,v,r,sols,i);
    Deep_Clear(samples);
  end Span_of_Component;

  procedure Initialize_Flag
              ( sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   All solutions in the list will have the m flag set to zero.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ptr : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      ptr.m := 0;
      Set_Head(tmp,ptr);
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Flag;

  procedure Initialize_Flag
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   All solutions in the list will have the m flag set to zero.

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ptr : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      ptr.m := 0;
      Set_Head(tmp,ptr);
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Flag;

  procedure Initialize_Flag
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   All solutions in the list will have the m flag set to zero.

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ptr : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ptr := Head_Of(tmp);
      ptr.m := 0;
      Set_Head(tmp,ptr);
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Flag;

  procedure Classify_to_Span
              ( file : in file_type; n,d,k : in integer32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys; 
                esols : in out Standard_Complex_Solutions.Solution_List;
                method : in character ) is

  -- DESCRIPTION :
  --   Classifies the witness points according to their linear span.

  -- ON ENTRY :
  --   file     for the results;
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system;
  --   esols    solutions to ep;
  --   method   1 for LU; 2 for QR.

    use Standard_Complex_VecVecs;
    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Complex_Poly_Systems;

    p : constant Poly_Sys := Witness_Sets.Remove_Embedding1(ep,natural32(d));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Coefficients(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    tol : constant double_float := 1.0E-8;
    isols_ptr : Solution_List := isols;
    esols_ptr : Solution_List := esols;
    i_ls,e_ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    put_line(file,"The original polynomial system :");
    put(file,p);
    Initialize_Flag(esols);
    for i in 1..Length_Of(isols) loop
      i_ls := Head_Of(isols_ptr);
      e_ls := Head_Of(esols_ptr);
      if e_ls.m = 0 then
        new_line(file);
        put(file,"Sampling from point ");
        put(file,i,1); put_line(file," ...");
        Span_of_Component(file,n,p,pla,i_ls.all,method,tol,esols_ptr,i);
        e_ls.m := integer32(i);
        Set_Head(esols_ptr,e_ls);
        cnt := cnt + 1;
      else
        put(file,"Point "); put(file,i,1);
        put(file," belongs to linear span of point ");
        put(file,e_ls.m,1); put_line(file,".");
      end if;
      isols_ptr := Tail_Of(isols_ptr);
      esols_ptr := Tail_Of(esols_ptr);
    end loop;
    put(file,"Number of different linear spans : ");
    put(file,cnt,1); put_line(file,".");
    Standard_Complex_VecVecs.Clear(s);
    Standard_Complex_Solutions.Clear(isols);
  end Classify_to_Span;

  procedure Classify_to_Span
              ( file : in file_type; n,d,k : in integer32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
                esols : in out DoblDobl_Complex_Solutions.Solution_List;
                method : in character ) is

  -- DESCRIPTION :
  --   Classifies the witness points according to their linear span.

  -- ON ENTRY :
  --   file     for the results;
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system;
  --   esols    solutions to ep;
  --   method   1 for LU; 2 for QR.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Poly_Systems;

    p : constant Poly_Sys := Witness_Sets.Remove_Embedding1(ep,natural32(d));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Coefficients(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    tol : constant double_float := 1.0E-16;
    isols_ptr : Solution_List := isols;
    esols_ptr : Solution_List := esols;
    i_ls,e_ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    put_line(file,"The original polynomial system :");
    put(file,p);
    Initialize_Flag(esols);
    for i in 1..Length_Of(isols) loop
      i_ls := Head_Of(isols_ptr);
      e_ls := Head_Of(esols_ptr);
      if e_ls.m = 0 then
        new_line(file);
        put(file,"Sampling from point ");
        put(file,i,1); put_line(file," ...");
        Span_of_Component(file,n,p,pla,i_ls.all,method,tol,esols_ptr,i);
        e_ls.m := integer32(i);
        Set_Head(esols_ptr,e_ls);
        cnt := cnt + 1;
      else
        put(file,"Point "); put(file,i,1);
        put(file," belongs to linear span of point ");
        put(file,e_ls.m,1); put_line(file,".");
      end if;
      isols_ptr := Tail_Of(isols_ptr);
      esols_ptr := Tail_Of(esols_ptr);
    end loop;
    put(file,"Number of different linear spans : ");
    put(file,cnt,1); put_line(file,".");
    DoblDobl_Complex_VecVecs.Clear(s);
    DoblDobl_Complex_Solutions.Clear(isols);
  end Classify_to_Span;

  procedure Classify_to_Span
              ( file : in file_type; n,d,k : in integer32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
                esols : in out QuadDobl_Complex_Solutions.Solution_List;
                method : in character ) is

  -- DESCRIPTION :
  --   Classifies the witness points according to their linear span.

  -- ON ENTRY :
  --   file     for the results;
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system;
  --   esols    solutions to ep;
  --   method   1 for LU; 2 for QR.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Poly_Systems;

    p : constant Poly_Sys := Witness_Sets.Remove_Embedding1(ep,natural32(d));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Coefficients(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    tol : constant double_float := 1.0E-32;
    isols_ptr : Solution_List := isols;
    esols_ptr : Solution_List := esols;
    i_ls,e_ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    put_line(file,"The original polynomial system :");
    put(file,p);
    Initialize_Flag(esols);
    for i in 1..Length_Of(isols) loop
      i_ls := Head_Of(isols_ptr);
      e_ls := Head_Of(esols_ptr);
      if e_ls.m = 0 then
        new_line(file);
        put(file,"Sampling from point ");
        put(file,i,1); put_line(file," ...");
        Span_of_Component(file,n,p,pla,i_ls.all,method,tol,esols_ptr,i);
        e_ls.m := integer32(i);
        Set_Head(esols_ptr,e_ls);
        cnt := cnt + 1;
      else
        put(file,"Point "); put(file,i,1);
        put(file," belongs to linear span of point ");
        put(file,e_ls.m,1); put_line(file,".");
      end if;
      isols_ptr := Tail_Of(isols_ptr);
      esols_ptr := Tail_Of(esols_ptr);
    end loop;
    put(file,"Number of different linear spans : ");
    put(file,cnt,1); put_line(file,".");
    QuadDobl_Complex_VecVecs.Clear(s);
    QuadDobl_Complex_Solutions.Clear(isols);
  end Classify_to_Span;

  procedure Standard_Test is

    ep : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    file : file_type;
    n,d,k : integer32 := 0;
    method : character;

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
    new_line;
    put_line("MENU to test intrinsic operations :");
    put_line("  1. use LU for path tracking with intrinsic Newton;");
    put_line("  2. use QR for path tracking with intrinsic Newton;");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(method,"12");
    new_line;
    put_line("Reading name of file to write results.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See output file for results...");
    new_line;
    Classify_to_Span(file,n,d,k,ep.all,sols,method);
  end Standard_Test;

  procedure DoblDobl_Test is

    ep : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    file : file_type;
    n,d,k : integer32 := 0;
    method : character;

  begin
    DoblDobl_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    new_line;
    put_line("MENU to test intrinsic operations :");
    put_line("  1. use LU for path tracking with intrinsic Newton;");
    put_line("  2. use QR for path tracking with intrinsic Newton;");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(method,"12");
    new_line;
    put_line("Reading name of file to write results.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See output file for results...");
    new_line;
    Classify_to_Span(file,n,d,k,ep.all,sols,method);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

    ep : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    file : file_type;
    n,d,k : integer32 := 0;
    method : character;

  begin
    QuadDobl_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    new_line;
    put_line("MENU to test intrinsic operations :");
    put_line("  1. use LU for path tracking with intrinsic Newton;");
    put_line("  2. use QR for path tracking with intrinsic Newton;");
    put("Type 1 or 2 to select : ");
    Ask_Alternative(method,"12");
    new_line;
    put_line("Reading name of file to write results.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See output file for results...");
    new_line;
    Classify_to_Span(file,n,d,k,ep.all,sols,method);
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Classifying witness set according to its linear span ...");
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012"); 
    case ans is
      when '1' => Standard_Test;
      when '2' => DoblDobl_Test;
      when '3' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_iwset;
