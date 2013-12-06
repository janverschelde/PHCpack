with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Continuation_Parameters;
with Standard_IncFix_Continuation;      use Standard_IncFix_Continuation;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets;                      use Witness_Sets;
with Witness_Sets_io;                   use Witness_Sets_io;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Standard_Intrinsic_Newton;

procedure ts_endgm is

-- DESCRIPTION :
--   Tracking paths with intrinsic coordinates, using the standard
--   path tracking routines in PHCpack.

--  procedure Set_Continuation_Parameter
--              ( sa : in out Solu_Info_Array; vt : in Complex_Number ) is
--  begin
--    for i in sa'range loop
--      sa(i).sol.t := vt;
--    end loop;
--  end Set_Continuation_Parameter;

  function Random_Plane ( n,k : integer32 ) return Matrix is

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

  function Special_Plane
             ( n,k : integer32; ind : Standard_Integer_Vectors.Vector )
             return Matrix is

  -- DESCRIPTION :
  --   Returns a k-plane in n-space in special position,
  --   according to the indices in the vector ind.
  --   If ind(i) = j, for some i, then x(j) will get a fixed value,
  --   provided j > 0.  If j = 0, then the offset vector is 0.

    res : Matrix(1..n,0..k) := Random_Plane(n,k);

  begin
    for j in ind'range loop
      if ind(j) = 0 then
        for i in 1..n loop
          res(i,0) := Create(0.0);
        end loop;
      else
        for i in 1..k loop
          res(ind(j),i) := Create(0.0);
        end loop;
      end if;
    end loop;
    return res;
  end Special_Plane;

  function Interactive_Special_Plane ( n,k : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   Allows the user to define a special k-plane in n-space,
  --   fixing several coordinates.

    res : Matrix(1..n,0..k);
    fix_nb : integer32 := 0;

  begin
    put("Give the number of fixed coordinates : ");
    get(fix_nb);
    if fix_nb = 0 then
      res := Random_Plane(n,k);
    else
      declare
        fix_ind : Standard_Integer_Vectors.Vector(1..fix_nb);
      begin
        put("Give "); put(fix_nb,1); put(" indices : ");
        get(fix_ind);
        put("Generating "); put(k,1);
        put("-plane with fixed coordinates");
        put(fix_ind); new_line;
        res := Special_Plane(n,k,fix_ind);
      end;
    end if;
    return res;
  end Interactive_Special_Plane;

  function Middle_Plane ( start_plane,target_plane : Matrix;
                          t : Complex_Number ) return Matrix is

  -- DESCRIPTION :
  --   Returns (1-t)*start_plane + t*target_plane.

    res : Matrix(start_plane'range(1),start_plane'range(2));
    one_min_t : constant Complex_Number := Create(1.0) - t;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := one_min_t*start_plane(i,j) + t*target_plane(i,j);
      end loop;
    end loop;
    return res;
  end Middle_Plane;

  procedure Call_Tracker ( f : in Poly_Sys; sols : in out Solution_List;
                           p : in Matrix ) is

    ef : constant Eval_Poly_Sys := Create(f);
    jm : constant Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : constant Eval_Jaco_Mat := Create(jm);
    sp,target,tp : Matrix(p'range(1),p'range(2));
    oc : natural32;
    output : boolean;
    wrk : Solution_List;
    ans : character;

    function H ( x : Vector; t : Complex_Number ) return Vector is

      mp : constant Matrix(p'range(1),p'range(2)) := Moving_Plane(sp,tp,t);
      z : constant Vector(p'range(1)) := Affine_Expand(x,mp);

    begin
      return Eval(ef,z);
    end H;

    function dH1 ( x : Vector; t : Complex_Number ) return Matrix is

      mp : constant Matrix(p'range(1),p'range(2)) := Moving_Plane(sp,tp,t);
      z : constant Vector(p'range(1)) := Affine_Expand(x,mp);
 
    begin
      return Standard_Intrinsic_Newton.Affine_Eval(jf,mp,z);
    end dH1;

    function dH2 ( x : Vector; t : Complex_Number ) return Vector is

      z1 : constant Vector(p'range(1)) := Affine_Expand(x,sp);
      z2 : constant Vector(p'range(1)) := Affine_Expand(x,tp);

    begin
      return (Eval(ef,z2) - Eval(ef,z1));
    end dH2;

    procedure R_Cont is new Reporting_Continue(Max_Norm,H,dH2,dH1);
    procedure S_Cont is new Silent_Continue(Max_Norm,H,dH2,dH1);

  begin
    new_line;
    Driver_for_Process_io(Standard_Output,oc);
    output := (oc > 0);
    Continuation_Parameters.Tune(2);
    Set_Continuation_Parameter(sols,Create(0.0));
    loop
      new_line;
      target := Interactive_Special_Plane(p'last(1),p'last(2));
      sp := p;
      tp := Middle_Plane(p,target,Create(0.9));
      Copy(sols,wrk);
      if output
       then R_Cont(Standard_Output,wrk,false);
       else S_Cont(wrk,false);
      end if;
      Set_Continuation_Parameter(wrk,Create(0.0));
      sp := tp;
      tp := target;
      if output
       then R_Cont(Standard_Output,wrk,false);
       else S_Cont(wrk,false);
      end if;
      put("Do you want another run ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Call_Tracker;

  procedure Test_Tracker ( n,d,k : in integer32;
                           ep : in Poly_Sys; esols : in Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    s : constant VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
    Call_Tracker(p,isols,pla);
  end Test_Tracker;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : integer32;

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
    Test_Tracker(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Path tracking with end game in intrinsic coordinates...");
  Main;
end ts_endgm;
