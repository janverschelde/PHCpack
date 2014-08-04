with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Complex_VecLists;
with DoblDobl_Complex_VecLists;
with QuadDobl_Complex_VecLists;
with Multprec_Complex_VecLists;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_to_Real_Poly;     use Standard_Complex_to_Real_Poly;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Floating_Poly_Systems;
with Standard_Complex_Poly_Systems; 
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_to_Multprec_Convertors;   use Standard_to_Multprec_Convertors;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;
with Standard_Complex_Newton_Steps;
with DoblDobl_Complex_Newton_Steps;
with QuadDobl_Complex_Newton_Steps;
with Multprec_Complex_Newton_Steps;
with Standard_Aitken_Extrapolation;
with Multprec_Aitken_Extrapolation;

procedure ts_newton is

-- DESCRIPTION :
--   Interactive testing of Newton with SVD and Aitken extrapolation.

  procedure Write_Residuals
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in Standard_Complex_VecLists.List ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_VecLists;

    tmp : List := z;
    cnt : natural32 := 0;
    lv : Standard_Complex_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      declare
        y : Standard_Complex_Vectors.Vector(lv'range) := Eval(p,lv.all);
        r : double_float := Standard_Complex_Norms_equals.Norm2(y);
      begin
        put(cnt,2); put(" : "); put(r); new_line;
      end;
      tmp := Tail_Of(tmp);
      cnt := cnt + 1;
    end loop;
  end Write_Residuals;

  procedure Write_Residuals
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in DoblDobl_Complex_VecLists.List ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_VecLists;

    tmp : List := z;
    cnt : natural32 := 0;
    lv : DoblDobl_Complex_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      declare
        y : DoblDobl_Complex_Vectors.Vector(lv'range) := Eval(p,lv.all);
        r : double_double := DoblDobl_Complex_Vector_Norms.Norm2(y);
      begin
        put(cnt,2); put(" : "); put(r); new_line;
      end;
      tmp := Tail_Of(tmp);
      cnt := cnt + 1;
    end loop;
  end Write_Residuals;

  procedure Write_Residuals
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in Multprec_Complex_VecLists.List ) is

    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_VecLists;

    tmp : List := z;
    cnt : natural32 := 0;
    lv : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      declare
        y : Multprec_Complex_Vectors.Vector(lv'range) := Eval(p,lv.all);
        r : Floating_Number := Multprec_Complex_Norms_equals.Norm2(y);
      begin
        put(cnt,2); put(" : "); put(r); new_line;
        Multprec_Complex_Vectors.Clear(y); Clear(r);
      end;
      tmp := Tail_Of(tmp);
      cnt := cnt + 1;
    end loop;
  end Write_Residuals;

  procedure Standard_Apply_Aitken
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in Standard_Complex_VecLists.List;
                e : out Standard_Complex_Vectors.Vector ) is

    ez : Standard_Complex_VecLists.List
       := Standard_Aitken_Extrapolation.Extrapolate(standard_output,z);
    tmp : Standard_Complex_VecLists.List;

    use Standard_Complex_VecLists;

  begin
    put_line("Residuals at the original sequence : ");   
    Write_Residuals(p,z);
    put_line("Residuals at the new sequence : ");
    Write_Residuals(p,ez);
    tmp := ez;
    while not Is_Null(Tail_Of(tmp)) loop
      tmp := Tail_Of(tmp);
    end loop;
    e := Head_Of(tmp).all;
  end Standard_Apply_Aitken;

  procedure Multprec_Apply_Aitken
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in Multprec_Complex_VecLists.List;
                e : out Multprec_Complex_Vectors.Vector ) is

    ez : Multprec_Complex_VecLists.List
       := Multprec_Aitken_Extrapolation.Extrapolate(standard_output,z);
    tmp : Multprec_Complex_VecLists.List;

    use Multprec_Complex_VecLists;

  begin
    put_line("Residuals at the original sequence : ");   
    Write_Residuals(p,z);
    put_line("Residuals at the new sequence : ");
    Write_Residuals(p,ez);
    tmp := ez;
    while not Is_Null(Tail_Of(tmp)) loop
      tmp := Tail_Of(tmp);
    end loop;
    e := Head_Of(tmp).all;
  end Multprec_Apply_Aitken;

  procedure Write_Standard_Solution
              ( z : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector.

    use Standard_Complex_Solutions;
    s : Solution(z'length);

  begin
    s.v := z;
    put_vector(s);
  end Write_Standard_Solution;

  procedure Write_DoblDobl_Solution
              ( z : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector.

    use DoblDobl_Complex_Solutions;
    s : Solution(z'length);

  begin
    s.v := z;
    put_vector(s);
  end Write_DoblDobl_Solution;

  procedure Write_QuadDobl_Solution
              ( z : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector.

    use QuadDobl_Complex_Solutions;
    s : Solution(z'length);

  begin
    s.v := z;
    put_vector(s);
  end Write_QuadDobl_Solution;

  procedure Write_Multprec_Solution
              ( z : in Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector.

    use Multprec_Complex_Solutions;
    s : Solution(z'length);

  begin
    s.v := z;
    put_vector(s);
  end Write_Multprec_Solution;

  procedure Write_Standard_Solution
              ( z : in Standard_Complex_Vectors.Vector;
                err,rco,res : in double_float; cnt,rank : in natural32 ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector z,
  --   but then uses the values of err, rco, res for the diagnostics.

  begin
    put("The new approximation at step "); put(cnt,1); put_line(" : "); 
    Write_Standard_Solution(z);
    put("error :"); put(err,3);
    put("  rcond :"); put(rco,3);
    put("  resid :"); put(res,3); new_line;
    put("The numerical rank : "); put(rank,1); new_line;
  end Write_Standard_Solution;

  procedure Write_DoblDobl_Solution
              ( z : in DoblDobl_Complex_Vectors.Vector;
                err,rco,res : in double_double; cnt,rank : in natural32 ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector z,
  --   but then uses the values of err, rco, res for the diagnostics.

  begin
    put("The new approximation at step "); put(cnt,1); put_line(" : "); 
    Write_DoblDobl_Solution(z);
    put("error :"); put(err,3);
    put("  rcond :"); put(rco,3);
    put("  resid :"); put(res,3); new_line;
    put("The numerical rank : "); put(rank,1); new_line;
  end Write_DoblDobl_Solution;

  procedure Write_QuadDobl_Solution
              ( z : in QuadDobl_Complex_Vectors.Vector;
                err,rco,res : in quad_double; cnt,rank : in natural32 ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector z,
  --   but then uses the values of err, rco, res for the diagnostics.

  begin
    put("The new approximation at step "); put(cnt,1); put_line(" : "); 
    Write_QuadDobl_Solution(z);
    put("error :"); put(err,3);
    put("  rcond :"); put(rco,3);
    put("  resid :"); put(res,3); new_line;
    put("The numerical rank : "); put(rank,1); new_line;
  end Write_QuadDobl_Solution;

  procedure Write_Multprec_Solution
              ( z : in Multprec_Complex_Vectors.Vector;
                err,rco,res : in Floating_Number; cnt,rank : in natural32 ) is

  -- DESCRIPTION :
  --   Uses the symbol table to write the solution vector z,
  --   but then uses the values of err, rco, res for the diagnostics.

  begin
    put("The new approximation at step "); put(cnt,1); put_line(" : "); 
    Write_Multprec_Solution(z);
    put("error :"); put(err,3);
    put("  rcond :"); put(rco,3);
    put("  resid :"); put(res,3); new_line;
    put("The numerical rank : "); put(rank,1); new_line;
  end Write_Multprec_Solution;

  procedure Call_Standard_Newton
              ( file : in file_type; tol : in double_float;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Complex_VecLists;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,z'range) := Create(p);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    err,rco,res : double_float;
    rank : natural32;
    cnt : natural32 := 0;
    ez : Standard_Complex_Vectors.Vector(z'range);
    sz,sz_last : List;
    order : natural32 := 1;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is

      res : Matrix(jm'range(1),jm'range(2)) := Eval(ejm,x);

    begin
      return res;
    end jmf;

    procedure Newton is
      new Standard_Complex_Newton_Steps.Reporting_Newton_Step(f,jmf);

  begin
    Append(sz,sz_last,z);
    loop
      cnt := cnt + 1;
      Newton(file,natural32(p'last),z,tol,err,rco,res,rank);
      Write_Standard_Solution(z,err,rco,res,cnt,rank);
      Append(sz,sz_last,z);
      if cnt > 1 then
        put("Apply Aitken extrapolation ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          Standard_Apply_Aitken(p,sz,ez);
          put_line("The extrapolated root :");
          Write_Standard_Solution(ez);
        end if;
       -- order := order + 1;
      end if;
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Call_Standard_Newton;

  procedure Call_DoblDobl_Newton
              ( file : in file_type; tol : in double_float;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z, in double double precision.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_VecLists;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,z'range) := Create(p);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    err,rco,res : double_double;
    rank : natural32;
    cnt : natural32 := 0;
    sz,sz_last : List;
    order : natural32 := 1;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is

      res : Matrix(jm'range(1),jm'range(2)) := Eval(ejm,x);

    begin
      return res;
    end jmf;

    procedure Newton is
      new DoblDobl_Complex_Newton_Steps.Reporting_Newton_Step(f,jmf);

  begin
    Append(sz,sz_last,z);
    loop
      cnt := cnt + 1;
      Newton(file,natural32(p'last),z,tol,err,rco,res,rank);
      Write_DoblDobl_Solution(z,err,rco,res,cnt,rank);
      Append(sz,sz_last,z);
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Call_DoblDobl_Newton;

  procedure Call_QuadDobl_Newton
              ( file : in file_type; tol : in double_float;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z, in quad double precision.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_VecLists;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,z'range) := Create(p);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    err,rco,res : quad_double;
    rank : natural32;
    cnt : natural32 := 0;
    sz,sz_last : List;
    order : natural32 := 1;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is

      res : Matrix(jm'range(1),jm'range(2)) := Eval(ejm,x);

    begin
      return res;
    end jmf;

    procedure Newton is
      new QuadDobl_Complex_Newton_Steps.Reporting_Newton_Step(f,jmf);

  begin
    Append(sz,sz_last,z);
    loop
      cnt := cnt + 1;
      Newton(file,natural32(p'last),z,tol,err,rco,res,rank);
      Write_QuadDobl_Solution(z,err,rco,res,cnt,rank);
      Append(sz,sz_last,z);
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Call_QuadDobl_Newton;

  procedure Append_Copy
              ( first,last : in out Multprec_Complex_VecLists.List;
                z : in Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Appends a copy of z to the list pointed to by first.

    az : Multprec_Complex_Vectors.Vector(z'range);

  begin
    Multprec_Complex_Vectors.Copy(z,az);
    Multprec_Complex_VecLists.Append(first,last,az);
  end Append_Copy;

  procedure Call_Multprec_Newton
              ( file : in file_type; tol : in double_float;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z.

  -- NOTE :
  --   One could think that doubling the working precision at each
  --   Newton step is a good idea, but it makes the execution time
  --   to grow exponentially as well, so it is commented out.
  --   The iterations happen in fixed precision.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Solutions;
    use Multprec_Complex_VecLists;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,z'range) := Create(p);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    err,rco,res : Floating_Number;
    rank : natural32;
    cnt : natural32 := 0;
    ez : Multprec_Complex_Vectors.Vector(z'range);
    sz,sz_last : List;
    order : natural32 := 1;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(ep,x);
    end f;

    function jmf ( x : Vector ) return Matrix is

      res : Matrix(jm'range(1),jm'range(2)) := Eval(ejm,x);

    begin
      return res;
    end jmf;

    procedure Newton is
      new Multprec_Complex_Newton_Steps.Reporting_Newton_Step(f,jmf);

  begin
    Append_Copy(sz,sz_last,z);
    loop
      cnt := cnt + 1;
      Newton(file,natural32(p'last),z,tol,err,rco,res,rank);
      Write_Multprec_Solution(z,err,rco,res,cnt,rank);
      Append_Copy(sz,sz_last,z);
      if cnt > 1 then
        put("Apply Aitken extrapolation ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          Multprec_Apply_Aitken(p,sz,ez);
          put_line("The extrapolated root :");
          Write_Multprec_Solution(ez);
        end if;
      end if;
      put("Do you want another iteration ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Call_Multprec_Newton;

  procedure Standard_Newton
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                n : in natural32; tol : in double_float;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method in standard double precision.

    z : Standard_Complex_Vectors.Vector(1..integer32(n));
    tmp : Standard_Complex_Solutions.Solution_List;
    ans : character;

  begin
    if not Standard_Complex_Solutions.Is_Null(sols) then
      z := Standard_Complex_Solutions.Head_Of(sols).v;
      tmp := Standard_Complex_Solutions.Tail_Of(sols);
    else
      put("Give "); put(n,1);
      put_line(" complex numbers to start at:"); get(z);
    end if;
    loop
      Call_Standard_Newton(Standard_Output,tol,p,z);
      exit when Standard_Complex_Solutions.Is_Null(tmp);
      put("Move on to the next root ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      z := Standard_Complex_Solutions.Head_Of(tmp).v;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Standard_Newton;

  procedure DoblDobl_Newton
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                n : in natural32; tol : in double_float;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision.

    z : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    tmp : DoblDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      z := DoblDobl_Complex_Solutions.Head_Of(sols).v;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(sols);
    else
      put("Give "); put(n,1);
      put_line(" complex numbers to start at : "); 
      DoblDobl_Complex_Vectors_io.get(z);
    end if;
    loop
      Call_DoblDobl_Newton(Standard_Output,tol,p,z);
      exit when DoblDobl_Complex_Solutions.Is_Null(tmp);
      put("Move on to the next root ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      z := DoblDobl_Complex_Solutions.Head_Of(tmp).v;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end DoblDobl_Newton;

  procedure QuadDobl_Newton
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                n : in natural32; tol : in double_float;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision.

    z : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    tmp : QuadDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      z := QuadDobl_Complex_Solutions.Head_Of(sols).v;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(sols);
    else
      put("Give "); put(n,1);
      put_line(" complex numbers to start at : "); 
      QuadDobl_Complex_Vectors_io.get(z);
    end if;
    loop
      Call_QuadDobl_Newton(Standard_Output,tol,p,z);
      exit when QuadDobl_Complex_Solutions.Is_Null(tmp);
      put("Move on to the next root ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      z := QuadDobl_Complex_Solutions.Head_Of(tmp).v;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Newton;

  procedure Multprec_Newton
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                n : in natural32; tol : in double_float;
                sols : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method in multiprecision.

    z : Multprec_Complex_Vectors.Vector(1..integer32(n));
    tmp : Multprec_Complex_Solutions.Solution_List;
    ans : character;

  begin
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      z := Multprec_Complex_Solutions.Head_Of(sols).v;
      tmp := Multprec_Complex_Solutions.Tail_Of(sols);
    else
      put("Give "); put(n,1);
      put_line(" complex numbers to start at : ");
      Multprec_Complex_Vectors_io.get(z);
    end if;
    loop
      Call_Multprec_Newton(Standard_Output,tol,p,z);
      exit when Multprec_Complex_Solutions.Is_Null(tmp);
      put("Move on to the next root ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      z := Multprec_Complex_Solutions.Head_Of(tmp).v;
      tmp := Multprec_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Multprec_Newton;

  procedure Newton_with_Standard_Complex_Arithmetic is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    infile : file_type;
    ans : character;
    found : boolean;
    sols : Solution_List;
    tol : double_float := 1.0E-8;
    n : natural32 := 0;

  begin
    new_line;
    put("Is the polynomial system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      n := Number_of_Unknowns(lp(lp'first));
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if found then
        get(infile,sols);
        new_line;
        put("Read "); put(Length_Of(sols),1);
        put_line(" solutions from file.");
      end if;
    else
      new_line;
      put("Give the dimension of the system : "); get(n);
      Symbol_Table.Init(n);
      lp := new Poly_Sys(1..integer32(n));
      put("Give "); put(n,1); put_line(" polynomials :");
      get(lp.all); skip_line;
      found := false;
    end if;
    new_line;
    put("Tolerance to decide the numerical rank is ");
    put(tol,3); new_line;
    Standard_Newton(lp.all,n,tol,sols);
  end Newton_with_Standard_Complex_Arithmetic;

  procedure Newton_with_DoblDobl_Complex_Arithmetic is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    infile : file_type;
    ans : character;
    found : boolean;
    sols : Solution_List;
    tol : double_float := 1.0E-8;
    n : natural32 := 0;

  begin
    new_line;
    put("Is the polynomial system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      n := Number_of_Unknowns(lp(lp'first));
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if found then
        get(infile,sols);
        new_line;
        put("Read "); put(Length_Of(sols),1);
        put_line(" solutions from file.");
      end if;
    else
      new_line;
      put("Give the dimension of the system : "); get(n);
      Symbol_Table.Init(n);
      lp := new Poly_Sys(1..integer32(n));
      put("Give "); put(n,1); put_line(" polynomials :");
      get(standard_input,lp.all); skip_line;
      found := false;
    end if;
    new_line;
    put("Tolerance to decide the numerical rank is ");
    put(tol,3); new_line;
    DoblDobl_Newton(lp.all,n,tol,sols);
  end Newton_with_DoblDobl_Complex_Arithmetic;

  procedure Newton_with_QuadDobl_Complex_Arithmetic is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    infile : file_type;
    ans : character;
    found : boolean;
    sols : Solution_List;
    tol : double_float := 1.0E-8;
    n : natural32 := 0;

  begin
    new_line;
    put("Is the polynomial system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      n := Number_of_Unknowns(lp(lp'first));
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if found then
        get(infile,sols);
        new_line;
        put("Read "); put(Length_Of(sols),1);
        put_line(" solutions from file.");
      end if;
    else
      new_line;
      put("Give the dimension of the system : "); get(n);
      Symbol_Table.Init(n);
      lp := new Poly_Sys(1..integer32(n));
      put("Give "); put(n,1); put_line(" polynomials :");
      get(standard_input,lp.all); skip_line;
      found := false;
    end if;
    new_line;
    put("Tolerance to decide the numerical rank is ");
    put(tol,3); new_line;
    QuadDobl_Newton(lp.all,n,tol,sols);
  end Newton_with_QuadDobl_Complex_Arithmetic;

  procedure Newton_with_Multprec_Complex_Arithmetic is

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    infile : file_type;
    ans : character;
    found : boolean;
    sols : Solution_List;
    tol : double_float := 1.0E-8;
    n,dgts,size : natural32 := 0;

  begin
    new_line;
    put("Is the polynomial system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      n := Number_of_Unknowns(lp(lp'first));
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if found then
        get(infile,sols);
        new_line;
        put("Read "); put(Length_Of(sols),1);
        put_line(" solutions from file.");
      end if;
    else
      new_line;
      put("Give the dimension of the system : "); get(n);
      Symbol_Table.Init(n);
      lp := new Poly_Sys(1..integer32(n));
      put("Give "); put(n,1); put_line(" polynomials :");
      get(lp.all); skip_line;
      found := false;
    end if;
    new_line;
    put("Give the number of digits : ");
    get(dgts);
    size := Decimal_to_Size(dgts);
    Set_Size(lp.all,size);
    Set_Size(sols,size);
    new_line;
    put("Tolerance to decide the numerical rank is ");
    put(tol,3); new_line;
    Multprec_Newton(lp.all,n,tol,sols);
  end Newton_with_Multprec_Complex_Arithmetic;

  procedure Call_Newton
              ( p : in Standard_Floating_Poly_Systems.Poly_Sys;
                tol : in double_float ) is
  begin
    null;
  end Call_Newton;

  procedure Newton_in_Real_Arithmetic is

    clp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlp : Standard_Floating_Poly_Systems.Link_to_Poly_Sys;
    tol : double_float := 1.0E-8;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(clp);
    rlp := new Standard_Floating_Poly_Systems.Poly_Sys'
                 (Convert_Complex_to_Real(clp.all));
    put("Tolerance to decide the numerical rank is ");
    put(tol,3); new_line;
    Call_Newton(rlp.all,tol);
  end Newton_in_Real_Arithmetic;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing Newton's method ...");
    new_line;
    put_line("MENU of testing operations : ");
    put_line("  1. Newton with standard double complex arithmetic; or");
    put_line("  2. Newton with double double complex arithmetic; or");
    put_line("  3. Newton with quad double complex arithmetic; or");
    put_line("  4. Newton with multiprecision complex arithmetic.");
    put("Type 1, 2, 3, or 4 to make a choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Newton_with_Standard_Complex_Arithmetic;
      when '2' => Newton_with_DoblDobl_Complex_Arithmetic;
      when '3' => Newton_with_QuadDobl_Complex_Arithmetic;
      when '4' => Newton_with_Multprec_Complex_Arithmetic;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_newton;
