with text_io;                           use text_io;
with timing_package;                    use timing_package;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;                  use String_Splitters;
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
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;      use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
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
with Standard_Continuation_Data;
with Standard_Continuation_Data_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Continuation_Data;
with DoblDobl_Continuation_Data_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Continuation_Data;
with QuadDobl_Continuation_Data_io;
with Continuation_Parameters;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;
with Standard_Plane_Representations;
with Standard_Moving_Planes;
with Standard_Rescaling_Coordinates;
with Standard_Intrinsic_Solutions;
with DoblDobl_Plane_Representations;
with DoblDobl_Moving_Planes;
with DoblDobl_Rescaling_Coordinates;
with QuadDobl_Plane_Representations;
with QuadDobl_Moving_Planes;
with QuadDobl_Rescaling_Coordinates;
with DoblDobl_Intrinsic_Solutions;
with QuadDobl_Intrinsic_Solutions;
with Standard_Intrinsic_Newton;
with Standard_Intrinsic_Trackers;
with Standard_Intrinsic_Continuation;
with DoblDobl_Intrinsic_Newton;
with DoblDobl_Intrinsic_Trackers;
with DoblDobl_Intrinsic_Continuation;
with QuadDobl_Intrinsic_Newton;
with QuadDobl_Intrinsic_Trackers;
with QuadDobl_Intrinsic_Continuation;

procedure ts_ilocal is

-- DESCRIPTION :
--   Interactive development of local intrinsic coordinates
--   for numerically stable witness set representations.

  procedure Check_Orthogonality
              ( p : in Standard_Complex_Matrices.Matrix;
                d : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows all inner products of the vector d with the vectors that
  --   span the plane p, in standard double precision.

    v : Standard_Complex_Vectors.Vector(p'range(1));
    c : Standard_Complex_Numbers.Complex_Number;

  begin
    put_line("Checking orthogonality of direction vector ...");
    for j in 1..p'last(2) loop
      for i in v'range loop
        v(i) := p(i,j);
      end loop;
      c := Conjugated_Inner_Product(d,v);
      put("  * with vector "); put(j,1); put(" : ");
      put(c); new_line;
    end loop;
  end Check_Orthogonality;

  procedure Check_Orthogonality
              ( p : in DoblDobl_Complex_Matrices.Matrix;
                d : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows all inner products of the vector d with the vectors that
  --   span the plane p, in double double precision.

    v : DoblDobl_Complex_Vectors.Vector(p'range(1));
    c : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    put_line("Checking orthogonality of direction vector ...");
    for j in 1..p'last(2) loop
      for i in v'range loop
        v(i) := p(i,j);
      end loop;
      c := Conjugated_Inner_Product(d,v);
      put("  * with vector "); put(j,1); put(" : ");
      put(c); new_line;
    end loop;
  end Check_Orthogonality;

  procedure Check_Orthogonality
              ( p : in QuadDobl_Complex_Matrices.Matrix;
                d : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows all inner products of the vector d with the vectors that
  --   span the plane p, in quad double precision.

    v : QuadDobl_Complex_Vectors.Vector(p'range(1));
    c : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    put_line("Checking orthogonality of direction vector ...");
    for j in 1..p'last(2) loop
      for i in v'range loop
        v(i) := p(i,j);
      end loop;
      c := Conjugated_Inner_Product(d,v);
      put("  * with vector "); put(j,1); put(" : ");
      put(c); new_line;
    end loop;
  end Check_Orthogonality;

  procedure Evaluate_Intrinsic_Solutions
              ( p : in Standard_Complex_Matrices.Matrix;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates f at the intrinsic coordinates for the solutions in sols,
  --   after expanding using the plane p.

  -- ON ENTRY :
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   s        intrinsic coordinates of the solutions.

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;
    use Standard_Intrinsic_Solutions;

    n : constant integer32 := p'last(1);
    ef : Eval_Poly_Sys(f'range) := Create(f);
    y : Standard_Complex_Vectors.Vector(f'range);
    tmp : Solution_List := s;
    ls : Link_to_Solution;
    res : double_float := 0.0;

  begin
    for i in 1..Length_Of(s) loop
      put("Value at solution "); put(i,1); put(" :");
      ls := Head_Of(tmp);
      declare
        es : constant Solution(n) := Expand(ls.all,p);
        rs : double_float;
      begin
        y := Eval(ef,es.v);
        rs := Max_Norm(y);
        put(rs); new_line;
        res := res + rs;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put("Sum of all residuals :"); put(res); new_line;
    Standard_Complex_Poly_SysFun.Clear(ef);
  end Evaluate_Intrinsic_Solutions;

  procedure Evaluate_Intrinsic_Solutions
              ( p : in DoblDobl_Complex_Matrices.Matrix;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates f at the intrinsic coordinates for the solutions in sols,
  --   after expanding using the plane p.

  -- ON ENTRY :
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   s        intrinsic coordinates of the solutions.

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Intrinsic_Solutions;

    n : constant integer32 := p'last(1);
    ef : Eval_Poly_Sys(f'range) := Create(f);
    y : DoblDobl_Complex_Vectors.Vector(f'range);
    tmp : Solution_List := s;
    ls : Link_to_Solution;
    res : double_double := create(0.0);

  begin
    for i in 1..Length_Of(s) loop
      put("Value at solution "); put(i,1); put(" :");
      ls := Head_Of(tmp);
      declare
        es : constant Solution(n) := Expand(ls.all,p);
        rs : double_double;
      begin
        y := Eval(ef,es.v);
        rs := Max_Norm(y);
        put(rs); new_line;
        res := res + rs;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put("Sum of all residuals :"); put(res); new_line;
    DoblDobl_Complex_Poly_SysFun.Clear(ef);
  end Evaluate_Intrinsic_Solutions;

  procedure Evaluate_Intrinsic_Solutions
              ( p : in QuadDobl_Complex_Matrices.Matrix;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates f at the intrinsic coordinates for the solutions in sols,
  --   after expanding using the plane p.

  -- ON ENTRY :
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   s        intrinsic coordinates of the solutions.

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Intrinsic_Solutions;

    n : constant integer32 := p'last(1);
    ef : Eval_Poly_Sys(f'range) := Create(f);
    y : QuadDobl_Complex_Vectors.Vector(f'range);
    tmp : Solution_List := s;
    ls : Link_to_Solution;
    res : quad_double := create(0.0);

  begin
    for i in 1..Length_Of(s) loop
      put("Value at solution "); put(i,1); put(" :");
      ls := Head_Of(tmp);
      declare
        es : constant Solution(n) := Expand(ls.all,p);
        rs : quad_double;
      begin
        y := Eval(ef,es.v);
        rs := Max_Norm(y);
        put(rs); new_line;
        res := res + rs;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put("Sum of all residuals :"); put(res); new_line;
    QuadDobl_Complex_Poly_SysFun.Clear(ef);
  end Evaluate_Intrinsic_Solutions;

  function Evaluate_Prediction
              ( f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                x : Standard_Complex_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the 2-norm of f evaluated at x.
  --   The prediction will turn right if the number on return is
  --   of the same magnitude of the step size.

    y : constant Standard_Complex_Vectors.Vector
      := Standard_Complex_Poly_SysFun.Eval(f,x);
    res : constant double_float := Norm2(y);

  begin
    return res;
  end Evaluate_Prediction;

  function Evaluate_Prediction
              ( f : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                x : DoblDobl_Complex_Vectors.Vector ) return double_double is

  -- DESCRIPTION :
  --   Returns the 2-norm of f evaluated at x.
  --   The prediction will turn right if the number on return is
  --   of the same magnitude of the step size.

    y : constant DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Complex_Poly_SysFun.Eval(f,x);
    res : constant double_double := Norm2(y);

  begin
    return res;
  end Evaluate_Prediction;

  function Evaluate_Prediction
              ( f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                x : QuadDobl_Complex_Vectors.Vector ) return quad_double is

  -- DESCRIPTION :
  --   Returns the 2-norm of f evaluated at x.
  --   The prediction will turn right if the number on return is
  --   of the same magnitude of the step size.

    y : constant QuadDobl_Complex_Vectors.Vector
      := QuadDobl_Complex_Poly_SysFun.Eval(f,x);
    res : constant quad_double := Norm2(y);

  begin
    return res;
  end Evaluate_Prediction;

  procedure Show_Value_of_Prediction
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                p : in Standard_Complex_Matrices.Matrix;
                h : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates f at the offset of the k-plane p and also writes the
  --   current step size h for comparison.

    x : Standard_Complex_Vectors.Vector(p'range(1));
    v : double_float;

  begin
    for i in p'range(1) loop
      x(i) := p(i,0);
    end loop;
    v := Evaluate_Prediction(f,x);
    put("Norm of prediction : "); put(v);
    put("  step : "); put(h,3); new_line;
  end Show_Value_of_Prediction;

  procedure Show_Value_of_Prediction
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                p : in DoblDobl_Complex_Matrices.Matrix;
                h : in double_double ) is

  -- DESCRIPTION :
  --   Evaluates f at the offset of the k-plane p and also writes the
  --   current step size h for comparison.

    x : DoblDobl_Complex_Vectors.Vector(p'range(1));
    v : double_double;

  begin
    for i in p'range(1) loop
      x(i) := p(i,0);
    end loop;
    v := Evaluate_Prediction(f,x);
    put("Norm of prediction : "); put(v);
    put("  step : "); put(h,3); new_line;
  end Show_Value_of_Prediction;

  procedure Show_Value_of_Prediction
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                p : in QuadDobl_Complex_Matrices.Matrix;
                h : in quad_double ) is

  -- DESCRIPTION :
  --   Evaluates f at the offset of the k-plane p and also writes the
  --   current step size h for comparison.

    x : QuadDobl_Complex_Vectors.Vector(p'range(1));
    v : quad_double;

  begin
    for i in p'range(1) loop
      x(i) := p(i,0);
    end loop;
    v := Evaluate_Prediction(f,x);
    put("Norm of prediction : "); put(v);
    put("  step : "); put(h,3); new_line;
  end Show_Value_of_Prediction;

  procedure Global_versus_Local_Coordinates
              ( n,k : in integer32; h : in double_float;
                p : in Standard_Complex_Matrices.Matrix;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                esols,isols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions;
  --   isols    intrinsic coordinates of the solutions.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Rescaling_Coordinates;
    use Standard_Intrinsic_Newton;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.Random_Plane(n,k);
    step : constant Complex_Number := Create(h);
    movpla : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.Moving_Plane(p,target,step);
    newpla : Matrix(1..n,0..k) := target;
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp,itmp : Solution_List;
    els,ils : Link_to_Solution;
    epsax : constant double_float := 1.0e-12;
    epsrx : constant double_float := 0.0;
    epsaf : constant double_float := 1.0e-12;
    epsrf : constant double_float := 0.0;
    incax,incrx,extres,intres,resrf,rco : double_float;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    b : Vector(1..n);
    better_cnt,worse_cnt : natural32 := 0;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,1.0E-12,fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    itmp := isols; etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      ils := Head_Of(itmp); els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        x : Vector(1..ils.n) := ils.v;
      begin
        Affine_LU_Newton(Standard_Output,eq,jf,movpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,extres,resrf,nit,max,rco,fail);
      end;
     -- put_line("els.v : "); put_line(els.v);
     -- put_line("ils.v : "); put_line(ils.v);
      put("perturbing solution "); put(i,1);
      put_line(" using local coordinates ...");
      declare
        x : Vector(1..k) := (1..k => Create(0.0));
       -- c : constant Vector(1..n) := Linear_Offset_Shift(els.v,b,step);
        to_b : constant Vector(1..n) := b - els.v(1..n);
        c : Vector(1..n) := Complement_of_Projection(target,to_b);
      begin
        Normalize(c);
        c := step*c;
       -- put_line("the new offset vector :"); put_line(c);
        for i in c'range loop
          newpla(i,0) := els.v(i) + c(i);
        end loop;
        Affine_LU_Newton(Standard_Output,eq,jf,newpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,intres,resrf,nit,max,rco,fail);
      end;
      if intres < extres then
        put_line("Improved residual by local coordinates, ok.");
        better_cnt := better_cnt + 1;
      else
        put_line("Local coordinates did NOT improve residual!");
        worse_cnt := worse_cnt + 1;
      end if;
      itmp := Tail_Of(itmp); etmp := Tail_Of(etmp);
    end loop;
    put("# improvements : "); put(better_cnt,1); new_line;
    put("# worse cases  : "); put(worse_cnt,1); new_line;
    Standard_Complex_Poly_Systems.Clear(q);
    Standard_Complex_Poly_SysFun.Clear(eq);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Global_versus_Local_Coordinates;

  procedure Global_versus_Local_Coordinates
              ( n,k : in integer32; h : in double_double;
                p : in DoblDobl_Complex_Matrices.Matrix;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols,isols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions;
  --   isols    intrinsic coordinates of the solutions.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Rescaling_Coordinates;
    use DoblDobl_Intrinsic_Newton;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.Random_Plane(n,k);
    step : constant Complex_Number := Create(h);
    movpla : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.Moving_Plane(p,target,step);
    newpla : Matrix(1..n,0..k) := target;
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp,itmp : Solution_List;
    els,ils : Link_to_Solution;
    epsax : constant double_double := create(1.0e-12);
    epsrx : constant double_double := create(0.0);
    epsaf : constant double_double := create(1.0e-12);
    epsrf : constant double_double := create(0.0);
    incax,incrx,extres,intres,resrf,rco : double_double;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    b : Vector(1..n);
    better_cnt,worse_cnt : natural32 := 0;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    itmp := isols; etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      ils := Head_Of(itmp); els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        x : Vector(1..ils.n) := ils.v;
      begin
        Affine_LU_Newton(Standard_Output,eq,jf,movpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,extres,resrf,nit,max,rco,fail);
      end;
     -- put_line("els.v : "); put_line(els.v);
     -- put_line("ils.v : "); put_line(ils.v);
      put("perturbing solution "); put(i,1);
      put_line(" using local coordinates ...");
      declare
        zero : constant double_double := create(0.0);
        x : Vector(1..k) := (1..k => Create(zero));
       -- c : constant Vector(1..n) := Linear_Offset_Shift(els.v,b,step);
        to_b : constant Vector(1..n) := b - els.v(1..n);
        c : Vector(1..n) := Complement_of_Projection(target,to_b);
      begin
        Normalize(c);
        c := step*c;
       -- put_line("the new offset vector :"); put_line(c);
        for i in c'range loop
          newpla(i,0) := els.v(i) + c(i);
        end loop;
        Affine_LU_Newton(Standard_Output,eq,jf,newpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,intres,resrf,nit,max,rco,fail);
      end;
      if intres < extres then
        put_line("Improved residual by local coordinates, ok.");
        better_cnt := better_cnt + 1;
      else
        put_line("Local coordinates did NOT improve residual!");
        worse_cnt := worse_cnt + 1;
      end if;
      itmp := Tail_Of(itmp); etmp := Tail_Of(etmp);
    end loop;
    put("# improvements : "); put(better_cnt,1); new_line;
    put("# worse cases  : "); put(worse_cnt,1); new_line;
    DoblDobl_Complex_Poly_Systems.Clear(q);
    DoblDobl_Complex_Poly_SysFun.Clear(eq);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Global_versus_Local_Coordinates;

  procedure Global_versus_Local_Coordinates
              ( n,k : in integer32; h : in quad_double;
                p : in QuadDobl_Complex_Matrices.Matrix;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols,isols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions;
  --   isols    intrinsic coordinates of the solutions.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Rescaling_Coordinates;
    use QuadDobl_Intrinsic_Newton;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.Random_Plane(n,k);
    step : constant Complex_Number := Create(h);
    movpla : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.Moving_Plane(p,target,step);
    newpla : Matrix(1..n,0..k) := target;
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp,itmp : Solution_List;
    els,ils : Link_to_Solution;
    epsax : constant quad_double := create(1.0e-12);
    epsrx : constant quad_double := create(0.0);
    epsaf : constant quad_double := create(1.0e-12);
    epsrf : constant quad_double := create(0.0);
    incax,incrx,extres,intres,resrf,rco : quad_double;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    b : Vector(1..n);
    better_cnt,worse_cnt : natural32 := 0;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    itmp := isols; etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      ils := Head_Of(itmp); els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        x : Vector(1..ils.n) := ils.v;
      begin
        Affine_LU_Newton(Standard_Output,eq,jf,movpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,extres,resrf,nit,max,rco,fail);
      end;
     -- put_line("els.v : "); put_line(els.v);
     -- put_line("ils.v : "); put_line(ils.v);
      put("perturbing solution "); put(i,1);
      put_line(" using local coordinates ...");
      declare
        zero : constant quad_double := create(0.0);
        x : Vector(1..k) := (1..k => Create(zero));
       -- c : constant Vector(1..n) := Linear_Offset_Shift(els.v,b,step);
        to_b : constant Vector(1..n) := b - els.v(1..n);
        c : Vector(1..n) := Complement_of_Projection(target,to_b);
      begin
        Normalize(c);
        c := step*c;
       -- put_line("the new offset vector :"); put_line(c);
        for i in c'range loop
          newpla(i,0) := els.v(i) + c(i);
        end loop;
        Affine_LU_Newton(Standard_Output,eq,jf,newpla,x,
          epsax,epsrx,epsaf,epsrf,incax,incrx,intres,resrf,nit,max,rco,fail);
      end;
      if intres < extres then
        put_line("Improved residual by local coordinates, ok.");
        better_cnt := better_cnt + 1;
      else
        put_line("Local coordinates did NOT improve residual!");
        worse_cnt := worse_cnt + 1;
      end if;
      itmp := Tail_Of(itmp); etmp := Tail_Of(etmp);
    end loop;
    put("# improvements : "); put(better_cnt,1); new_line;
    put("# worse cases  : "); put(worse_cnt,1); new_line;
    QuadDobl_Complex_Poly_Systems.Clear(q);
    QuadDobl_Complex_Poly_SysFun.Clear(eq);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Global_versus_Local_Coordinates;

  procedure Track_one_Path
              ( n,k : in integer32; h : in double_float;
                f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                b : in Standard_Complex_Vectors.Vector;
                p : in Standard_Complex_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tracks one path using local coordinates in standard double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   f        polynomial system made square, of k equations;
  --   jf       Jacobian matrix of the system f.
  --   b        offset point of the target plane p;
  --   p        parametric representation of the target k-plane;
  --   x        start solution in extrinsic coordinates.

  -- ON RETURN :
  --   x        solution at the end of the path.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Rescaling_Coordinates;
    use Standard_Intrinsic_Newton;

    step : constant Complex_Number := Create(h);
    epsax : constant double_float := 1.0e-06;
    epsrx : constant double_float := 1.0e-06;
    epsaf : constant double_float := 1.0e-10;
    epsrf : constant double_float := 1.0e-10;
    incax,incrx,resaf,resrf,rco : double_float;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    ans : character;
    cnt : natural32 := 0;
    q : Matrix(1..n,0..k) := p;
    to_b : constant Vector(1..n) := b - x(1..n);
    c : constant Vector(1..n) := Complement_of_Projection(p,to_b);
    d : Vector(1..n) := c;
    nrm_b : constant double_float := Norm2(b);
    nrm_c : constant double_float := Norm2(c);
    z : Vector(1..k);

  begin
    put("Norm of to_b : "); put(nrm_b,3);
    put("  Norm of projection : "); put(nrm_c,3); new_line;
    Check_Orthogonality(p,c);
    Normalize(d);
    Check_Orthogonality(p,d);
    loop
      cnt := cnt + 1;
      put("Step "); put(cnt,1);
      put_line(", the solution vector :"); put_line(x);
      put("Distance to target : "); put(Distance(p,x)); new_line;
      z := (1..k => Create(0.0));
      for i in c'range loop
        q(i,0) := x(i) + step*d(i);
      end loop;
      Show_Value_of_Prediction(f,q,h);
      Affine_LU_Newton(Standard_Output,f,jf,q,z,
        epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      if fail
       then put_line("Reported failure to meet accuracy requirements!");
       else put_line("Met the accuracy requirements, no failure reported.");
      end if;
      put("Next step ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      for i in q'range(1) loop
        x(i) := q(i,0);
        for j in z'range loop
          x(i) := x(i) + z(j)*q(i,j);
        end loop;
      end loop;
    end loop;
  end Track_one_Path;

  procedure Track_one_Path
              ( n,k : in integer32; h : in double_double;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                b : in DoblDobl_Complex_Vectors.Vector;
                p : in DoblDobl_Complex_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tracks one path using local coordinates in double double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   f        polynomial system made square, of k equations;
  --   jf       Jacobian matrix of the system f.
  --   b        offset point of the target plane p;
  --   p        parametric representation of the target k-plane;
  --   x        start solution in extrinsic coordinates.

  -- ON RETURN :
  --   x        solution at the end of the path.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Rescaling_Coordinates;
    use DoblDobl_Intrinsic_Newton;

    step : constant Complex_Number := Create(h);
    epsax : constant double_double := create(1.0e-06);
    epsrx : constant double_double := create(1.0e-06);
    epsaf : constant double_double := create(1.0e-10);
    epsrf : constant double_double := create(1.0e-10);
    incax,incrx,resaf,resrf,rco : double_double;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    ans : character;
    cnt : natural32 := 0;
    q : Matrix(1..n,0..k) := p;
    to_b : constant Vector(1..n) := b - x(1..n);
    c : constant Vector(1..n) := Complement_of_Projection(p,to_b);
    d : Vector(1..n) := c;
    nrm_b : constant double_double := Norm2(b);
    nrm_c : constant double_double := Norm2(c);
    z : Vector(1..k);
    zero : constant double_double := create(0.0);

  begin
    put("Norm of to_b : "); put(nrm_b,3);
    put("  Norm of projection : "); put(nrm_c,3); new_line;
    Check_Orthogonality(p,c);
    Normalize(d);
    Check_Orthogonality(p,d);
    loop
      cnt := cnt + 1;
      put("Step "); put(cnt,1);
      put_line(", the solution vector :"); put_line(x);
      put("Distance to target : "); put(Distance(p,x)); new_line;
      z := (1..k => Create(zero));
      for i in c'range loop
        q(i,0) := x(i) + step*d(i);
      end loop;
      Show_Value_of_Prediction(f,q,h);
      Affine_LU_Newton(Standard_Output,f,jf,q,z,
        epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      if fail
       then put_line("Reported failure to meet accuracy requirements!");
       else put_line("Met the accuracy requirements, no failure reported.");
      end if;
      put("Next step ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      for i in q'range(1) loop
        x(i) := q(i,0);
        for j in z'range loop
          x(i) := x(i) + z(j)*q(i,j);
        end loop;
      end loop;
    end loop;
  end Track_one_Path;

  procedure Track_one_Path
              ( n,k : in integer32; h : in quad_double;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                b : in QuadDobl_Complex_Vectors.Vector;
                p : in QuadDobl_Complex_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tracks one path using local coordinates in double double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   f        polynomial system made square, of k equations;
  --   jf       Jacobian matrix of the system f.
  --   b        offset point of the target plane p;
  --   p        parametric representation of the target k-plane;
  --   x        start solution in extrinsic coordinates.

  -- ON RETURN :
  --   x        solution at the end of the path.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Rescaling_Coordinates;
    use QuadDobl_Intrinsic_Newton;

    step : constant Complex_Number := Create(h);
    epsax : constant quad_double := create(1.0e-06);
    epsrx : constant quad_double := create(1.0e-06);
    epsaf : constant quad_double := create(1.0e-10);
    epsrf : constant quad_double := create(1.0e-10);
    incax,incrx,resaf,resrf,rco : quad_double;
    max : constant natural32 := 3;
    nit : natural32 := 0;
    fail : boolean;
    ans : character;
    cnt : natural32 := 0;
    q : Matrix(1..n,0..k) := p;
    to_b : constant Vector(1..n) := b - x(1..n);
    c : constant Vector(1..n) := Complement_of_Projection(p,to_b);
    d : Vector(1..n) := c;
    nrm_b : constant quad_double := Norm2(b);
    nrm_c : constant quad_double := Norm2(c);
    z : Vector(1..k);
    zero : constant quad_double := create(0.0);

  begin
    put("Norm of to_b : "); put(nrm_b,3);
    put("  Norm of projection : "); put(nrm_c,3); new_line;
    Check_Orthogonality(p,c);
    Normalize(d);
    Check_Orthogonality(p,d);
    loop
      cnt := cnt + 1;
      put("Step "); put(cnt,1);
      put_line(", the solution vector :"); put_line(x);
      put("Distance to target : "); put(Distance(p,x)); new_line;
      z := (1..k => Create(zero));
      for i in c'range loop
        q(i,0) := x(i) + step*d(i);
      end loop;
      Show_Value_of_Prediction(f,q,h);
      Affine_LU_Newton(Standard_Output,f,jf,q,z,
        epsax,epsrx,epsaf,epsrf,incax,incrx,resaf,resrf,nit,max,rco,fail);
      if fail
       then put_line("Reported failure to meet accuracy requirements!");
       else put_line("Met the accuracy requirements, no failure reported.");
      end if;
      put("Next step ? (y/n) "); 
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      for i in q'range(1) loop
        x(i) := q(i,0);
        for j in z'range loop
          x(i) := x(i) + z(j)*q(i,j);
        end loop;
      end loop;
    end loop;
  end Track_one_Path;

  procedure Deform_using_Local_Coordinates
              ( n,k : in integer32; h : in double_float;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp : Solution_List;
    els : Link_to_Solution;
    b : Vector(1..n);
    ans : character;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,1.0E-12,fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      Track_one_Path(n,k,h,eq,jf,b,target,els.v);
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        etmp := Tail_Of(etmp);
      end if;
    end loop;
    Standard_Complex_Poly_Systems.Clear(q);
    Standard_Complex_Poly_SysFun.Clear(eq);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Deform_using_Local_Coordinates;

  procedure Deform_using_Local_Coordinates
              ( n,k : in integer32; h : in double_double;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp : Solution_List;
    els : Link_to_Solution;
    b : Vector(1..n);
    ans : character;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      Track_one_Path(n,k,h,eq,jf,b,target,els.v);
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        etmp := Tail_Of(etmp);
      end if;
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(q);
    DoblDobl_Complex_Poly_SysFun.Clear(eq);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Deform_using_Local_Coordinates;

  procedure Deform_using_Local_Coordinates
              ( n,k : in integer32; h : in quad_double;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   h        step size to move in the direction of the target;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    etmp : Solution_List;
    els : Link_to_Solution;
    b : Vector(1..n);
    ans : character;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    etmp := esols;
    for i in 1..n loop
      b(i) := target(i,0);
    end loop;
    for i in 1..Length_Of(esols) loop
      els := Head_Of(etmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      Track_one_Path(n,k,h,eq,jf,b,target,els.v);
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        etmp := Tail_Of(etmp);
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(q);
    QuadDobl_Complex_Poly_SysFun.Clear(eq);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Deform_using_Local_Coordinates;

  procedure Call_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Rescaling_Coordinates;
    use Standard_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    ans : character;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,1.0E-12,fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    tmp := esols;
    for i in 1..Length_Of(esols) loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if output then
          Reporting_Local_LU_Track(standard_output,eq,jf,target,s,pp,cp,fail);
        else
          Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
        end if;
      end;
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    Standard_Complex_Poly_Systems.Clear(q);
    Standard_Complex_Poly_SysFun.Clear(eq);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Rescaling_Coordinates;
    use DoblDobl_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    ans : character;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    tmp := esols;
    for i in 1..Length_Of(esols) loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if output then
          Reporting_Local_LU_Track(standard_output,eq,jf,target,s,pp,cp,fail);
        else
          Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
        end if;
      end;
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(q);
    DoblDobl_Complex_Poly_SysFun.Clear(eq);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   p        parametric representation of a k-plane;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Rescaling_Coordinates;
    use QuadDobl_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    ans : character;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    tmp := esols;
    for i in 1..Length_Of(esols) loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if output then
          Reporting_Local_LU_Track(standard_output,eq,jf,target,s,pp,cp,fail);
        else
          Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
        end if;
      end;
      if i < Length_Of(esols) then
        put("Go to next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(q);
    QuadDobl_Complex_Poly_SysFun.Clear(eq);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Call_Path_Trackers;

  procedure Run_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                start : in Standard_Complex_Matrices.Matrix;
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   calling the trackers one after the other on the solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   start    start plane that contains the solutions;
  --   esols    extrinsic coordinates of the solutions.

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Rescaling_Coordinates;
    use Standard_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;
    fail_cnt : natural32 := 0;
    totstep,totiter : natural32 := 0;
    deg : constant natural32 := Length_Of(esols);
    timer : Timing_Widget;
    ans : character;
    recenter : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,1.0E-12,fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    new_line;
    put("Using start plane to recenter ? (y/n) ");
    Ask_Yes_or_No(ans);
    recenter := (ans = 'y');
    tmp := esols;
    tstart(timer);
    for i in 1..deg loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if recenter then
          if output then
            Reporting_Recentered_LU_Track
              (standard_output,eq,jf,start,target,false,s,pp,cp,fail);
          else
            Silent_Recentered_LU_Track(eq,jf,start,target,false,s,pp,cp,fail);
          end if;
        else
          if output then
            Reporting_Local_LU_Track
              (standard_output,eq,jf,target,s,pp,cp,fail);
          else
            Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
          end if;
        end if;
        put(i,3); put(" #steps : "); put(s.nstep,1);
        put("  err : "); put(s.corr,3);
        put("  rco : "); put(s.rcond,3);
        put("  res : "); put(s.resa,3);
        if fail
         then put_line(" FAIL"); fail_cnt := fail_cnt + 1;
         else put_line(" ok");
        end if;
        totstep := totstep + s.nstep;
        totiter := totiter + s.niter;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    put("#paths : "); put(deg,1);
    put("  #failures : "); put(fail_cnt,1);
    put("  average #steps : "); put(totstep/deg,1);
    put("  average #iterations : "); put(totiter/deg,1); new_line;
    new_line;
    print_times(standard_output,timer,"local intrinsic sampling");
    Standard_Complex_Poly_Systems.Clear(q);
    Standard_Complex_Poly_SysFun.Clear(eq);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Run_Path_Trackers;

  procedure Run_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                start : in DoblDobl_Complex_Matrices.Matrix;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   calling the trackers one after the other on the solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   start    start plane that contains the solutions;
  --   esols    extrinsic coordinates of the solutions.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Rescaling_Coordinates;
    use DoblDobl_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;
    fail_cnt : natural32 := 0;
    totstep,totiter : natural32 := 0;
    deg : constant natural32 := Length_Of(esols);
    timer : Timing_Widget;
    ans : character;
    recenter : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    new_line;
    put("Using start plane to recenter ? (y/n) ");
    Ask_Yes_or_No(ans);
    recenter := (ans = 'y');
    tmp := esols;
    tstart(timer);
    for i in 1..deg loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if recenter then
          if output then
            Reporting_Recentered_LU_Track
              (standard_output,eq,jf,start,target,false,s,pp,cp,fail);
          else
            Silent_Recentered_LU_Track(eq,jf,start,target,false,s,pp,cp,fail);
          end if;
        else
          if output then
            Reporting_Local_LU_Track
              (standard_output,eq,jf,target,s,pp,cp,fail);
          else
            Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
          end if;
        end if;
        put(i,3); put(" #steps : "); put(s.nstep,1);
        put("  err : "); put(s.corr,3);
        put("  rco : "); put(s.rcond,3);
        put("  res : "); put(s.resa,3);
        if fail
         then put_line(" FAIL"); fail_cnt := fail_cnt + 1;
         else put_line(" ok");
        end if;
        totstep := totstep + s.nstep;
        totiter := totiter + s.niter;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    put("#paths : "); put(deg,1);
    put("  #failures : "); put(fail_cnt,1);
    put("  average #steps : "); put(totstep/deg,1);
    put("  average #iterations : "); put(totiter/deg,1); new_line;
    new_line;
    print_times(standard_output,timer,"local intrinsic sampling");
    DoblDobl_Complex_Poly_Systems.Clear(q);
    DoblDobl_Complex_Poly_SysFun.Clear(eq);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Run_Path_Trackers;

  procedure Run_Path_Trackers
              ( n,k : in integer32; output : in boolean;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                start : in QuadDobl_Complex_Matrices.Matrix;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   calling the trackers one after the other on the solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   start    start plane that contains the solutions;
  --   esols    extrinsic coordinates of the solutions.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Rescaling_Coordinates;
    use QuadDobl_Intrinsic_Trackers;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.Random_Plane(n,k);
    q : Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : Eval_Poly_Sys(1..k) := Create(q);
    jm : Jaco_Mat(1..k,1..n) := Create(q);
    jf : Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    tmp : Solution_List;
    ls : Link_to_Solution;
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    fail : boolean;
    fail_cnt : natural32 := 0;
    totstep,totiter : natural32 := 0;
    deg : constant natural32 := Length_Of(esols);
    timer : Timing_Widget;
    ans : character;
    recenter : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    new_line;
    put("Using start plane to recenter ? (y/n) ");
    Ask_Yes_or_No(ans);
    recenter := (ans = 'y');
    tmp := esols;
    tstart(timer);
    for i in 1..deg loop
      ls := Head_Of(tmp);
      put("perturbing solution "); put(i,1); put_line(" ...");
      declare
        s : Solu_Info := Shallow_Create(ls);
      begin
        if recenter then
          if output then
            Reporting_Recentered_LU_Track
              (standard_output,eq,jf,start,target,false,s,pp,cp,fail);
          else
            Silent_Recentered_LU_Track(eq,jf,start,target,false,s,pp,cp,fail);
          end if;
        else
          if output then
            Reporting_Local_LU_Track
              (standard_output,eq,jf,target,s,pp,cp,fail);
          else
            Silent_Local_LU_Track(eq,jf,target,s,pp,cp,fail);
          end if;
        end if;
        put(i,3); put(" #steps : "); put(s.nstep,1);
        put("  err : "); put(s.corr,3);
        put("  rco : "); put(s.rcond,3);
        put("  res : "); put(s.resa,3);
        if fail
         then put_line(" FAIL"); fail_cnt := fail_cnt + 1;
         else put_line(" ok");
        end if;
        totstep := totstep + s.nstep;
        totiter := totiter + s.niter;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    put("#paths : "); put(deg,1);
    put("  #failures : "); put(fail_cnt,1);
    put("  average #steps : "); put(totstep/deg,1);
    put("  average #iterations : "); put(totiter/deg,1); new_line;
    new_line;
    print_times(standard_output,timer,"local intrinsic sampling");
    QuadDobl_Complex_Poly_Systems.Clear(q);
    QuadDobl_Complex_Poly_SysFun.Clear(eq);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Run_Path_Trackers;

  procedure Local_Intrinsic_Continuation
              ( n,k : in integer32; output : in boolean;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                start : in Standard_Complex_Matrices.Matrix;
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   running intrinsic continuation on all solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Intrinsic_Continuation;
    use Standard_Continuation_Data;
    use Standard_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.One_Random_Direction(start);
          -- := Standard_Moving_Planes.Random_Plane(n,k);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
   -- totstep,totiter : natural := 0;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,1.0E-12,fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    declare
      file : file_type;
      name : Link_to_String;
    begin
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file,name);
      new_line;
      if not output then
        tstart(timer);
        Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      else
        tstart(timer);
        Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      end if;
      Write_Witness_Set_to_File
        (name.all,natural32(n),natural32(k),f,target,sa);
      new_line(file);
      Standard_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
    end;
  end Local_Intrinsic_Continuation;

  procedure Local_Intrinsic_Continuation
              ( n,k : in integer32; output : in boolean;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                start : in DoblDobl_Complex_Matrices.Matrix;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   running intrinsic continuation on all solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Intrinsic_Continuation;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.One_Random_Direction(start);
          -- := DoblDobl_Moving_Planes.Random_Plane(n,k);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
   -- totstep,totiter : natural := 0;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    declare
      file : file_type;
      name : Link_to_String;
    begin
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file,name);
      new_line;
      if not output then
        tstart(timer);
        Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      else
        tstart(timer);
        Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      end if;
      Write_Witness_Set_to_File
        (name.all,natural32(n),natural32(k),f,target,sa);
      new_line(file);
      DoblDobl_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
    end;
  end Local_Intrinsic_Continuation;

  procedure Local_Intrinsic_Continuation
              ( n,k : in integer32; output : in boolean;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                start : in QuadDobl_Complex_Matrices.Matrix;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   running intrinsic continuation on all solutions.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Intrinsic_Continuation;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Rescaling_Coordinates;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.One_Random_Direction(start);
          -- := QuadDobl_Moving_Planes.Random_Plane(n,k);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
   -- totstep,totiter : natural := 0;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);
    fail : boolean;

  begin
    new_line;
    put_line("The square system : "); put(q);
    put_line("moving to a random target plane : "); put(target,3);
    Check_Orthonormality(target,create(1.0E-12),fail);
    if fail
     then put_line("Target plane did not pass the orthonormality test!");
     else put_line("Target plane passed the orthonormality test.");
    end if;
    declare
      file : file_type;
      name : Link_to_String;
    begin
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file,name);
      new_line;
      if not output then
        tstart(timer);
        Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      else
        tstart(timer);
        Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
        tstop(timer);
      end if;
      Write_Witness_Set_to_File
        (name.all,natural32(n),natural32(k),f,target,sa);
      new_line(file);
      QuadDobl_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
    end;
  end Local_Intrinsic_Continuation;

  procedure Setup_Local_Coordinates
              ( n,d,k : in integer32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys; 
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the setup of working with local intrinsic coordinates,
  --   in standard double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables,
  --   d        dimension of the solution set;
  --   k        codimension of the solution set;
  --   ep       embedded polynomial system;
  --   esols    embedded extrinsic coordinates of the solutions.

    use Standard_Complex_Matrices;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Intrinsic_Solutions;
    use Standard_Plane_Representations;

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant natural32 := Number_of_Zero_Equations(p);
    pp : constant Poly_Sys := p(p'first..p'last-integer32(nbzero));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    step : double_float := 0.0;
    ans : character;
    output : boolean := true;
    oc : natural32 := 0;

  begin
    new_line;
    put_line("The original polynomial system : "); put(pp);
    put_line("The coefficients of the slices : "); put(eqs,3);
    put_line("The parametric representation of the plane : "); put(pla,3);
    Evaluate_Intrinsic_Solutions(pla,p,isols);
    new_line;
    put_line("MENU to test local intrinsic coordinates : ");
    put_line("  1. compare global versus local intrinsic coordinates;");
    put_line("  2. explore path tracking in local intrinsic coordinates;");
    put_line("  3. call trackers to run path after path;");
    put_line("  4. run trackers for an entire new witness set;");
    put_line("  5. call intrinsic continuation in local coordinates;");
    put("Type 1, 2, 3, 4, or 5 to choose : ");
    Ask_Alternative(ans,"12345");
    if ans = '3' or ans = '4' or ans = '5' then
      new_line;
      Driver_for_Continuation_Parameters;
      new_line;
      Driver_for_Process_io(Standard_Output,oc);
      output := (oc > 0);
    else
      new_line;
      put("Give the step size : "); get(step);
    end if;
    case ans is
      when '1' => Global_versus_Local_Coordinates(n,k,step,pla,p,esols,isols);
      when '2' => Deform_using_Local_Coordinates(n,k,step,p,esols);
      when '3' => Call_Path_Trackers(n,k,output,p,esols);
      when '4' => Run_Path_Trackers(n,k,output,p,pla,esols);
      when '5' => Local_Intrinsic_Continuation(n,k,output,pp,pla,esols);
      when others => null;
    end case;
    Standard_Complex_VecVecs.Clear(s);
    Standard_Complex_Solutions.Clear(isols);
  end Setup_Local_Coordinates;

  procedure Setup_Local_Coordinates
              ( n,d,k : in integer32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the setup of working with local intrinsic coordinates,
  --   in double double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables,
  --   d        dimension of the solution set;
  --   k        codimension of the solution set;
  --   ep       embedded polynomial system;
  --   esols    embedded extrinsic coordinates of the solutions.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Intrinsic_Solutions;
    use DoblDobl_Plane_Representations;

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant natural32 := Number_of_Zero_Equations(p);
    pp : constant Poly_Sys := p(p'first..p'last-integer32(nbzero));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    step : double_double := create(0.0);
    ans : character;
    output : boolean := true;
    oc : natural32 := 0;

  begin
    new_line;
    put_line("The original polynomial system : "); put(pp);
    put_line("The coefficients of the slices : "); put(eqs,3);
    put_line("The parametric representation of the plane : "); put(pla,3);
    Evaluate_Intrinsic_Solutions(pla,p,isols);
    new_line;
    put_line("MENU to test local intrinsic coordinates : ");
    put_line("  1. compare global versus local intrinsic coordinates;");
    put_line("  2. explore path tracking in local intrinsic coordinates;");
    put_line("  3. call trackers to run path after path;");
    put_line("  4. run trackers for an entire new witness set;");
    put_line("  5. call intrinsic continuation in local coordinates;");
    put("Type 1, 2, 3, 4, or 5 to choose : ");
    Ask_Alternative(ans,"12345");
    if ans = '3' or ans = '4' or ans = '5' then
      new_line;
      Driver_for_Continuation_Parameters;
      new_line;
      Driver_for_Process_io(standard_output,oc);
      output := (oc > 0);
    else
      new_line;
      put("Give the step size : "); get(step);
    end if;
    case ans is
      when '1' => Global_versus_Local_Coordinates(n,k,step,pla,p,esols,isols);
      when '2' => Deform_using_Local_Coordinates(n,k,step,p,esols);
      when '3' => Call_Path_Trackers(n,k,output,p,esols);
      when '4' => Run_Path_Trackers(n,k,output,p,pla,esols);
      when '5' => Local_Intrinsic_Continuation(n,k,output,pp,pla,esols);
      when others => null;
    end case;
    DoblDobl_Complex_VecVecs.Clear(s);
    DoblDobl_Complex_Solutions.Clear(isols);
  end Setup_Local_Coordinates;

  procedure Setup_Local_Coordinates
              ( n,d,k : in integer32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the setup of working with local intrinsic coordinates,
  --   in quad double precision.

  -- ON ENTRY :
  --   n        ambient dimension, number of original variables,
  --   d        dimension of the solution set;
  --   k        codimension of the solution set;
  --   ep       embedded polynomial system;
  --   esols    embedded extrinsic coordinates of the solutions.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Intrinsic_Solutions;
    use QuadDobl_Plane_Representations;

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant natural32 := Number_of_Zero_Equations(p);
    pp : constant Poly_Sys := p(p'first..p'last-integer32(nbzero));
    s : VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    isols : Solution_List := Project(esols,pla);
    step : quad_double := create(0.0);
    ans : character;
    output : boolean := true;
    oc : natural32 := 0;

  begin
    new_line;
    put_line("The original polynomial system : "); put(pp);
    put_line("The coefficients of the slices : "); put(eqs,3);
    put_line("The parametric representation of the plane : "); put(pla,3);
    Evaluate_Intrinsic_Solutions(pla,p,isols);
    new_line;
    put_line("MENU to test local intrinsic coordinates : ");
    put_line("  1. compare global versus local intrinsic coordinates;");
    put_line("  2. explore path tracking in local intrinsic coordinates;");
    put_line("  3. call trackers to run path after path;");
    put_line("  4. run trackers for an entire new witness set;");
    put_line("  5. call intrinsic continuation in local coordinates;");
    put("Type 1, 2, 3, 4, or 5 to choose : ");
    Ask_Alternative(ans,"12345");
    if ans = '3' or ans = '4' or ans = '5' then
      new_line;
      Driver_for_Continuation_Parameters;
      new_line;
      Driver_for_Process_io(standard_output,oc);
      output := (oc > 0);
    else
      new_line;
      put("Give the step size : "); get(step);
    end if;
    case ans is
      when '1' => Global_versus_Local_Coordinates(n,k,step,pla,p,esols,isols);
      when '2' => Deform_using_Local_Coordinates(n,k,step,p,esols);
      when '3' => Call_Path_Trackers(n,k,output,p,esols);
      when '4' => Run_Path_Trackers(n,k,output,p,pla,esols);
      when '5' => Local_Intrinsic_Continuation(n,k,output,pp,pla,esols);
      when others => null;
    end case;
    QuadDobl_Complex_VecVecs.Clear(s);
    QuadDobl_Complex_Solutions.Clear(isols);
  end Setup_Local_Coordinates;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a witness set
  --   and then calls the test on local intrinsic coordinates,
  --   in standard double precision.

    use Standard_Complex_Solutions;

    ep : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
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
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    Setup_Local_Coordinates(n,d,k,ep.all,sols);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a witness set
  --   and then calls the test on local intrinsic coordinates,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    ep : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    n,d,k : integer32 := 0;

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
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    Setup_Local_Coordinates(n,d,k,ep.all,sols);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a witness set
  --   and then calls the test on local intrinsic coordinates,
  --   in double double precision.

    use QuadDobl_Complex_Solutions;

    ep : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    n,d,k : integer32 := 0;

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
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    Setup_Local_Coordinates(n,d,k,ep.all,sols);
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Path tracking in local intrinsic coordinates ...");
    new_line;
    put_line("MENU to set the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_ilocal;
