with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;           use Continuation_Parameters;
with Standard_IncFix_Continuation;      use Standard_IncFix_Continuation;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Standard_Intrinsic_Newton;
--with Standard_Intrinsic_Trackers;       use Standard_Intrinsic_Trackers;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;

procedure ts_itrack is

-- DESCRIPTION :
--   Tracking paths with intrinsic coordinates, using the standard
--   path tracking routines in PHCpack.

 -- procedure Set_Continuation_Parameter
 --             ( sa : in out Solu_Info_Array; vt : in Complex_Number ) is
 -- begin
 --   for i in sa'range loop
 --     sa(i).sol.t := vt;
 --   end loop;
 -- end Set_Continuation_Parameter;

  procedure Call_Extrinsic_Tracker
             ( f : in Poly_Sys; sols : in out Solution_List;
               p : in Matrix ) is

  -- DESCRIPTION :
  --   This extrinsic tracker expands the solutions each time,
  --   which voids the possible advantage from savings in the
  --   reduction of the dimension.

    ef : Eval_Poly_Sys := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    oc : natural32;
    output : boolean;
    wrk : Solution_List;
    ans : character;

    function H ( x : Vector; t : Complex_Number ) return Vector is

      mp : constant Matrix(p'range(1),p'range(2)) := Moving_Plane(p,tp,t);
      z : constant Vector(p'range(1)) := Affine_Expand(x,mp);

    begin
      return Eval(ef,z);
    end H;

    function dH1 ( x : Vector; t : Complex_Number ) return Matrix is

      mp : constant Matrix(p'range(1),p'range(2)) := Moving_Plane(p,tp,t);
      z : constant Vector(p'range(1)) := Affine_Expand(x,mp);
 
    begin
      return Standard_Intrinsic_Newton.Affine_Eval(jf,mp,z);
    end dH1;

    function dH2 ( x : Vector; t : Complex_Number ) return Vector is

      z1 : constant Vector(p'range(1)) := Affine_Expand(x,p);
      z2 : constant Vector(p'range(1)) := Affine_Expand(x,tp);

    begin
      return (Eval(ef,z2) - Eval(ef,z1));
    end dH2;

    procedure R_Cont is new Reporting_Continue(Max_Norm,H,dH2,dH1);
    procedure S_Cont is new Silent_Continue(Max_Norm,H,dH2,dH1);

  begin
    new_line;
    Driver_for_Continuation_Parameters;
    new_line;
    Driver_for_Process_io(Standard_Output,oc);
    output := (oc > 0);
   -- Continuation_Parameters.Tune(2);
    Set_Continuation_Parameter(sols,Create(0.0));
    loop
      tp := Random_Plane(p'last(1),p'last(2));
      Copy(sols,wrk);
      if output
       then R_Cont(Standard_Output,wrk,false);
       else S_Cont(wrk,false);
      end if;
      put("Do you want another run ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Call_Extrinsic_Tracker;

  procedure Call_Intrinsic_Tracker
             ( f : in Poly_Sys; sols : in out Solution_List;
               p : in Matrix ) is

  -- DESCRIPTION :
  --   This extrinsic tracker expands the solutions each time,
  --   which voids the possible advantage from savings in the
  --   reduction of the dimension.

    ef : Eval_Poly_Sys := Create(f);
    jm : Jaco_Mat(f'range,p'range(1)) := Create(f);
    jf : Eval_Jaco_Mat := Create(jm);
    tp : Matrix(p'range(1),p'range(2));
    oc : natural32;
    output : boolean;
    wrk : Solu_Info_Array(1..integer32(Length_Of(sols)));
    ans : character;
    pp : Pred_Pars;
    cp : Corr_Pars;

    function MP ( t : Complex_Number ) return Matrix is

      res : constant Matrix(p'range(1),p'range(2)) := Moving_Plane(p,tp,t);

    begin
      return res;
    end MP;

    procedure R_Cont is new Reporting_Affine_LU_Continue(MP);
    procedure S_Cont is new Silent_Affine_LU_Continue(MP);

  begin
    new_line;
    Driver_for_Continuation_Parameters;
    new_line;
    Driver_for_Process_io(Standard_Output,oc);
    output := (oc > 0);
   -- Continuation_Parameters.Tune(2);
    Set_Continuation_Parameter(sols,Create(0.0));
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    loop
      tp := Random_Plane(p'last(1),p'last(2));
      wrk := Deep_Create(sols);
      if output
       then R_Cont(Standard_Output,ef,jf,wrk,pp,cp);
       else S_Cont(ef,jf,wrk,pp,cp);
      end if;
      Clear(wrk);
      put("Do you want another run ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Poly_SysFun.Clear(ef);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Call_Intrinsic_Tracker;

  procedure Test_Affine_Tracker
              ( n,d,k : in integer32;
                ep : in Poly_Sys; esols : in Solution_List ) is

  -- DESCRIPTION :
  --   Test the path tracking in affine intrinsic coordinates.

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system;
  --   esols    solutions in extrinsic coordinates.

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    ans : character;

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
    new_line;
    put("Calling intrinsic or extrinsic tracker ? (i/e) ");
    Ask_Alternative(ans,"ie");
    if ans = 'i'
     then Call_Intrinsic_Tracker(p,isols,pla);
     else Call_Extrinsic_Tracker(p,isols,pla);
    end if;
    Standard_Complex_VecVecs.Clear(s);
  end Test_Affine_Tracker;

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
    Test_Affine_Tracker(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Path tracking in intrinsic coordinates...");
  Main;
end ts_itrack;
