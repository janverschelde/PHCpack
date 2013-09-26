with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Hypersurface_Points;               use Hypersurface_Points;

procedure ts_hyppts is

  procedure Random_Poly ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Interactive generation of a random polynomial.

    d,m : natural32 := 0;

  begin
    put_line("Generating a random multivariate complex polynomial...");
    n := 0;
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
    p := Random_Sparse_Poly(n,d,m,0);
    put_line("The random polynomial : ");
    put_line(p);
  end Random_Poly;

  function Lead_Degrees ( p : Poly ) return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the degrees of the leading monomial of p.

    res : Standard_Natural_Vectors.Link_to_Vector;

    procedure Lead_Term ( t : in Term; continue : out boolean ) is
    begin
      res := Standard_Natural_Vectors.Link_to_Vector(t.dg);
      continue := false;
    end Lead_Term;
    procedure Lead_Terms is new Visiting_Iterator(Lead_Term);

  begin
    Lead_Terms(p);
    return res.all;
  end Lead_Degrees;

  procedure Test_Degree_Start_Hypersurface ( p : in Poly ) is

  -- DESCRIPTION :
  --   Creates a degree start hypersurface from the lead degrees of p.

    deg : constant Standard_Natural_Vectors.Vector := Lead_Degrees(p);
    q : Poly;
    eq : Eval_Poly;
    d : constant integer32 := integer32(Standard_Natural_Vectors.Sum(deg));
    b : constant Vector(deg'range) := (deg'range => Create(0.0));
    v : constant Vector(deg'range) := Random_Vector(deg'first,deg'last);
    t : Vector(1..d);
    maxres : double_float;

  begin
    put("The degrees : "); put(deg); new_line;
    Degree_Start_Hypersurface(deg,natural32(d),v,q,t);
    put("The degree start hypersurface : ");
    put(q); new_line;
    eq := Create(q);
    maxres := Maximal_Affine_Residual(eq,b,v,t);
    put("The maximal residual is "); put(maxres,3); put_line(".");
  end Test_Degree_Start_Hypersurface;

  procedure Test_Scaling ( p : in Eval_Poly; b,v,roots : in Vector;
                           z : in out VecVec; maxres : out double_float ) is

    x : Vector(b'range);

   begin
    put_line("The roots before projective scaling :");
    maxres := Maximal_Projective_Residual(p,b,v,roots,z);
    put("The maximal residual is "); put(maxres,3); new_line;
    for i in roots'range loop
      Scale(b,v,roots(i),x,z(i).all);
    end loop;
    put_line("The roots after projective scaling :");
    maxres := Maximal_Projective_Residual(p,b,v,roots,z);
    put("The maximal residual is "); put(maxres,3); new_line;
  end Test_Scaling;

  procedure Test_Affine_Correct_for_Line
              ( n : in integer32; p : in Poly; ep : in Eval_Poly;
                b,v : in out Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests the corrector, perturbing b and v slightly and correcting the
  --   roots on the surface intersected with the perturbed line.

    dp : Poly_Sys(1..n);
    sp : Eval_Poly_Sys(0..n);
    nit : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      sp(i) := Create(dp(i));
    end loop;
    for i in b'range loop
      b(i) := b(i) + Random1*0.001;
      v(i) := v(i) + Random1*0.001;
    end loop;
    for i in roots'range loop
      put("Correcting root "); put(i,1); put_line(" :");
      Affine_Correct_for_Line
        (Standard_Output,sp,b,v,roots(i),10,nit,1.0E-14,fail);
    end loop;
  end Test_Affine_Correct_for_Line;

  procedure Test_Projective_Correct_for_Line
              ( n : in integer32; p : in Poly; ep : in Eval_Poly;
                b,v : in out Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests the corrector, perturbing b and v slightly and correcting the
  --   roots on the surface intersected with the perturbed line.

    dp : Poly_Sys(1..n);
    sp : Eval_Poly_Sys(0..n);
    nit : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      sp(i) := Create(dp(i));
    end loop;
    for i in b'range loop
      b(i) := b(i) + Random1*0.001;
      v(i) := v(i) + Random1*0.001;
    end loop;
    for i in roots'range loop
      put("Correcting root "); put(i,1); put_line(" :");
      declare
        z : constant Vector(b'range) := (b'range => Create(1.0));
      begin
        Projective_Correct_for_Line
          (Standard_Output,sp,b,v,roots(i),z,10,nit,1.0E-14,fail);
      end;
    end loop;
  end Test_Projective_Correct_for_Line;

  procedure Test_Affine_Track_Moving_Line
              ( n : in integer32; p : in Poly; ep : in Eval_Poly;
                b0,v0 : in out Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests tracking the path when the line through the surface is moved.

    dp : Poly_Sys(1..n);
    sp : Eval_Poly_Sys(0..n);
    b1 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    nbsteps : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      sp(i) := Create(dp(i));
    end loop;
    for i in roots'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Affine_Track_Moving_Line(sp,b0,v0,b1,v1,roots(i),nbsteps,fail);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
    b0 := b1;
    v0 := v1;
  end Test_Affine_Track_Moving_Line;

  procedure Test_Projective_Track_Moving_Line
              ( n : in integer32; p : in Poly; ep : in Eval_Poly;
                b0,v0 : in out Vector; roots : in out Vector;
                z : in out VecVec; maxres : out double_float ) is

  -- DESCRIPTION :
  --   Tests tracking the path when the line through the surface is moved.

    dp : Poly_Sys(1..n);
    sp : Eval_Poly_Sys(0..n);
    b1 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    nbsteps : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      sp(i) := Create(dp(i));
    end loop;
    for i in roots'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Projective_Track_Moving_Line
        (sp,b0,v0,b1,v1,roots(i),z(i).all,nbsteps,fail);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
    b0 := b1;
    v0 := v1;
    maxres := Maximal_Projective_Residual(ep,b0,v0,roots,z);
  end Test_Projective_Track_Moving_Line;

  procedure Test_Affine_Correct_for_Surface
              ( n : in integer32;
                p : in Poly; ep : in Eval_Poly;
                q : in Poly; eq : in Eval_Poly;
                b,v : in Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests the corrector, after moving the start hypersurface q slightly,
  --   the line b + t*v intersecting the surface remains fixed.
  --   Everything happens in affine coordinates.

    dp,dq : Poly_Sys(1..n);
    sp,sq : Eval_Poly_Sys(0..n);
    nit : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    sq(0) := eq;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      dq(i) := Diff(q,i);
      sp(i) := Create(dp(i));
      sq(i) := Create(dq(i));
    end loop;
    for i in roots'range loop
      put("Correcting root "); put(i,1); put_line(" :");
      Affine_Correct_for_Surface
        (Standard_Output,sp,sq,Create(0.9),b,v,roots(i),10,nit,1.0E-14,fail);
    end loop;
  end Test_Affine_Correct_for_Surface;

  procedure Test_Projective_Correct_for_Surface
              ( n : in integer32;
                p : in Poly; ep : in Eval_Poly;
                q : in Poly; eq : in Eval_Poly;
                b,v : in Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests the corrector, after moving the start hypersurface q slightly,
  --   the line b + t*v intersecting the surface remains fixed.
  --   Everything happens in projective coordinates.

    dp,dq : Poly_Sys(1..n);
    sp,sq : Eval_Poly_Sys(0..n);
    nit : natural32 := 0;
    fail : boolean;

  begin
    sp(0) := ep;
    sq(0) := eq;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      dq(i) := Diff(q,i);
      sp(i) := Create(dp(i));
      sq(i) := Create(dq(i));
    end loop;
    for i in roots'range loop
      put("Correcting root "); put(i,1); put_line(" :");
      declare
        z : constant Vector(b'range) := (b'range => Create(1.0));
      begin
        Projective_Correct_for_Surface
          (Standard_Output,sp,sq,Create(0.9),b,v,
           roots(i),z,10,nit,1.0E-14,fail);
      end;
    end loop;
  end Test_Projective_Correct_for_Surface;

  procedure Test_Affine_Track_Moving_Surface
              ( n : in integer32;
                p : in Poly; ep : in Eval_Poly;
                q : in Poly; eq : in Eval_Poly;
                b,v : in Vector; roots : in out Vector ) is

  -- DESCRIPTION :
  --   Tests the tracking of the paths when moving q to p,
  --   the line b + t*v intersecting the surface remains fixed.

    dp,dq : Poly_Sys(1..n);
    sp,sq : Eval_Poly_Sys(0..n);
    nbsteps : natural32;
    fail : boolean;

  begin
    sp(0) := ep;
    sq(0) := eq;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      dq(i) := Diff(q,i);
      sp(i) := Create(dp(i));
      sq(i) := Create(dq(i));
    end loop;
    for i in roots'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Affine_Track_Moving_Surface(sp,sq,b,v,roots(i),nbsteps,fail);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
  end Test_Affine_Track_Moving_Surface;

  procedure Test_Projective_Track_Moving_Surface
              ( n : in integer32;
                p : in Poly; ep : in Eval_Poly;
                q : in Poly; eq : in Eval_Poly; b,v : in Vector;
                roots : in out Vector; z : in out VecVec;
                maxres : out double_float ) is

  -- DESCRIPTION :
  --   Tests the tracking of the paths when moving q to p,
  --   the line b + t*v intersecting the surface remains fixed.

    dp,dq : Poly_Sys(1..n);
    sp,sq : Eval_Poly_Sys(0..n);
    nbsteps : natural32;
    fail : boolean;

  begin
    sp(0) := ep;
    sq(0) := eq;
    for i in 1..n loop
      dp(i) := Diff(p,i);
      dq(i) := Diff(q,i);
      sp(i) := Create(dp(i));
      sq(i) := Create(dq(i));
    end loop;
    for i in roots'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Projective_Track_Moving_Surface
        (sp,sq,b,v,roots(i),z(i).all,nbsteps,fail);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
    maxres := Maximal_Projective_Residual(ep,b,v,roots,z);
  end Test_Projective_Track_Moving_Surface;

  function Affine_Solve ( n : integer32; p : Poly ) return double_float is

  -- DESCRIPTION :
  --   Solves the polynomial p in n variables.  On return is the maximal
  --   residual of the polynomial at the computed roots.

    ep : constant Eval_Poly := Create(p);
   -- t : Complex_Number := Random1;
    b : Vector(1..n) := Random_Vector(1,n);
    v : Vector(1..n) := Random_Vector(1,n);
    q : Poly;
    eq : Eval_Poly;
    d : constant integer32 := Degree(p);
   -- s : Eval_Poly_Sys(0..n);
    roots : Vector(1..d);
    maxres : double_float;

  begin
    Affine_Start_Hypersurface(natural32(n),natural32(d),1,q,b,v,roots);
   -- put("The start hypersurface : "); put(q); new_line;
    eq := Create(q);
   -- put_line("The roots of the start hypersurface :");
   -- maxres := Maximal_Residual(eq,b,v,roots);
   -- put("The maximal residual is "); put(maxres,3); new_line;
    Test_Affine_Track_Moving_Line(n,q,eq,b,v,roots);
    put_line("The roots at the new line : ");
    maxres := Maximal_Affine_Residual(eq,b,v,roots);
    put("The maximal residual is "); put(maxres,3); new_line;
    Test_Affine_Track_Moving_Surface(n,p,ep,q,eq,b,v,roots);
    maxres := Maximal_Affine_Residual(ep,b,v,roots);
   -- put("The maximal residual is "); put(maxres,3); new_line;
    return maxres;
   -- Test_Affine_Correct_for_Line(n,q,eq,b,v,roots);
   -- Test_Affine_Correct_for_Surface(n,p,ep,q,eq,b,v,roots);
  end Affine_Solve;

  function Affine_Degree_Solve
              ( n : integer32; p : Poly ) return double_float is

  -- DESCRIPTION :
  --   Solves the polynomial p in n variables.  On return is the maximal
  --   residual of the polynomial at the computed roots.  The start
  --   hypersurface is based on the leading degrees and the line through
  --   the surface passes through the origin.

    ep : constant Eval_Poly := Create(p);
    eq : Eval_Poly;
    deg : constant Standard_Natural_Vectors.Vector := Lead_Degrees(p);
    q : Poly;
    d : constant integer32 := integer32(Standard_Natural_Vectors.Sum(deg));
    b : constant Vector(deg'range) := (deg'range => Create(0.0));
    v : constant Vector(deg'range) := Random_Vector(deg'first,deg'last);
    t : Vector(1..d);
    maxres : double_float;

  begin
    put("The degrees : "); put(deg); new_line;
    Degree_Start_Hypersurface(deg,natural32(d),v,q,t);
    put("The start hypersurface : "); put(q); new_line;
    eq := Create(q);
    put_line("The roots of the start hypersurface :");
    maxres := Maximal_Affine_Residual(eq,b,v,t);
    put("The maximal residual is "); put(maxres,3); new_line;
    Test_Affine_Track_Moving_Surface(n,p,ep,q,eq,b,v,t);
    maxres := Maximal_Affine_Residual(ep,b,v,t);
    put("The maximal residual is "); put(maxres,3); new_line;
    return maxres;
  end Affine_Degree_Solve;

  procedure Solve_Multihomogeneous_Start_Hypersurface1
                  ( n : in integer32; p : in Poly; q : out Poly;
                    b,v,roots : out Vector; z : out VecVec;
                    maxres : out double_float ) is

  -- DESCRIPTION :
  --   Constructs and solves a multi-homogeneous start hypersurface from
  --   the leading term of the polynomial p.

    deg : constant Standard_Natural_Vectors.Vector(1..n) := Lead_Degrees(p);
   -- t : Complex_Number := Random1;
    b0,v0 : Vector(1..n);
    bb,vv : VecVec(1..n);
    eq : Eval_Poly;
    ind : integer32 := 0;
    dq : Poly_Sys(1..n);
    sq : Eval_Poly_Sys(0..n);
    b1 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    nbsteps : natural32;
    fail : boolean;

  begin
    Multihomogenization_Symbols(natural32(n));
    put("Degrees of the lead term : "); put(deg); new_line;
    for i in deg'range loop  -- set up the product of start polynomials
      if deg(i) > 0 then
        declare
          qi,mhqi : Poly;
          rts : Vector(1..integer32(deg(i)));
        begin
          Affine_Start_Hypersurface
            (natural32(n),deg(i),natural32(i),qi,b0,v0,rts);
          put_line("b0 : "); put_line(b0);
          put_line("v0 : "); put_line(v0);
          bb(i) := new Vector'(b0);
          vv(i) := new vector'(v0);
          mhqi := Multihomogenize(natural32(n),qi);
          put("The "); put(n,1);
          put_line("-homogeneous start polynomial :");
          put_line(mhqi);
          for j in rts'range loop
            ind := ind + 1;
            roots(ind) := rts(j);
            z(ind) := new Vector'(1..n => Create(1.0));
          end loop;
          if i = deg'first
           then q := mhqi;
           else Mul(q,mhqi);
          end if;
        end;
      end if;
    end loop;
    ind := 0;
    eq := Create(q);
    sq(0) := eq;
    for i in 1..n loop
      dq(i) := Diff(q,i);
      sq(i) := Create(dq(i));
    end loop;
    for i in deg'range loop  -- perform continuation to general line
      if deg(i) > 0 then
        for j in 1..deg(i) loop
          ind := ind + 1;
          put("Tracking path "); put(ind,1); put(" ... ");
          Projective_Track_Moving_Line
            (sq,bb(i).all,vv(i).all,b1,v1,roots(ind),z(ind).all,nbsteps,fail);
          put("done in "); put(nbsteps,1); put_line(" steps.");
        end loop;
      end if;
    end loop;
    maxres := Maximal_Projective_Residual(eq,b1,v1,roots,z);
    b := b1; v := v1;
  end Solve_Multihomogeneous_Start_Hypersurface1; 

  procedure Solve_Multihomogeneous_Start_Hypersurface
                  ( n : in integer32; p : in Poly; q : out Poly;
                    b,v,roots : out Vector; z : out VecVec;
                    maxres : out double_float ) is

    deg : constant Standard_Natural_Vectors.Vector(1..n) := Lead_Degrees(p);
    q1 : Poly;
    eq1,eq : Eval_Poly;
    d : constant natural32 := Standard_Natural_Vectors.Sum(deg);
    b0 : constant Vector(1..n) := (1..n => Create(0.0));
    sq : Eval_Poly_Sys(0..n);
    dq : Poly_Sys(1..n);
    nbsteps : natural32;
    fail : boolean;

  begin
    Multihomogenization_Symbols(natural32(n));
    put("Degrees of the lead term : "); put(deg); new_line;
    v := Random_Vector(1,n);
    Degree_Start_Hypersurface(deg,d,v,q1,roots);
    eq1 := Create(q1);
    put_line("Residuals of affine roots :");
    maxres := Maximal_Affine_Residual(eq1,b0,v,roots);
    for i in z'range loop
      z(i) := new Vector'(1..n => Create(1.0));
    end loop;
    q := Multihomogenize(natural32(n),q1);
    eq := Create(q);
    put_line("Residuals of projective roots :");
    maxres := Maximal_Projective_Residual(eq,b0,v,roots,z);
    sq(0) := eq;
    for i in 1..n loop
      dq(i) := Diff(q,i);
      sq(i) := Create(dq(i));
    end loop;
    b := Random_Vector(1,n);
    for i in roots'range loop  -- perform continuation to general line
      put("Tracking path "); put(i,1); put(" ... ");
      Projective_Track_Moving_Line
        (sq,b0,v,b,v,roots(i),z(i).all,nbsteps,fail);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
    maxres := Maximal_Projective_Residual(eq,b,v,roots,z);
  end Solve_Multihomogeneous_Start_Hypersurface;

  function Projective_Solve ( n : integer32; p : Poly ) return double_float is

   -- deg : constant Standard_Natural_Vectors.Vector(1..n) := Lead_Degrees(p);
   -- t : Complex_Number := Random1;
    b,v : Vector(1..n);
    q : Poly;
   -- mhq,mhp : Poly;
   -- eq,ep : Eval_Poly;
    d : constant integer32 := Degree(p);
   -- s : Eval_Poly_Sys(0..n);
    roots : Vector(1..d);
    z : VecVec(1..d);
    maxres : double_float := 1.0;
   -- ind : natural := 0;

  begin
    Solve_Multihomogeneous_Start_Hypersurface(n,p,q,b,v,roots,z,maxres);
    put("The start polynomial : "); put(q); new_line;
   -- Test_Scaling(eq,b,v,roots,z,maxres);
   -- Test_Projective_Track_Moving_Line(n,mhq,eq,b,v,roots,z,maxres);
   -- Test_Projective_Correct_for_Line(n,q,eq,b,v,roots);
   -- mhp := Multihomogenize(n,p);
   -- ep := Create(mhp);
   -- Test_Projective_Track_Moving_Surface(n,mhp,ep,mhq,eq,b,v,roots,z,maxres);
   -- Test_Projective_Correct_for_Surface(n,p,ep,q,eq,b,v,roots);
    return maxres;
  end Projective_Solve;

  procedure Random_Test ( nb : in natural32 ) is

    n,d,m : integer32 := 0;
    p : Poly;
    tol : constant double_float := 1.0E-8;
    res : double_float;

  begin
    put_line("Generation of random polynomials...");
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
    for i in 1..nb loop
      p := Random_Sparse_Poly(natural32(n),natural32(d),natural32(m),0);
      put_line(p);
      res := Affine_Degree_Solve(n,p);
     -- res := Affine_Solve(n,p);
     -- res := Projective_Solve(n,p);
      put("Maximal residual "); put(res,3);
      if res > tol then
        put(" > "); put(tol,3); put("  Failure at case ");
        put(i,1); put_line(".");
        return;
      else
        put(" <= "); put(tol,3); put_line("  Success.");
      end if;
    end loop;
    put("Tested "); put(nb,1); put_line(" cases successfully.");
  end Random_Test;

  procedure Main is

    n : natural32 := 0;
   -- p : Poly;

  begin
    new_line;
    put_line("Tracing points on a hypersurface.");
    new_line;
   -- Random_Poly(n,p);
   -- Test_Degree_Start_Hypersurface(p);
    put("Give the number of random tests : "); get(n);
    Random_Test(n);
  end Main;

begin
  Main;
end ts_hyppts;
