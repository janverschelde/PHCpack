with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Vector_Splitters;
with Standard_Complex_Linear_Solvers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;            use Standard_Circuit_Makers;

procedure ts_newcirc is

-- DESCRIPTION :
--   Test the development of Newton's method on coefficient circuits.

  procedure LU_Newton_Step
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err : out double_float;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one step with Newton's method on s, starting at v,
  --   using LU factorization to solve the linear system.

  -- REQUIRED : s.neq = s.dim, i.e.: the system is square,
  --   and xr'range = xi'range = v'range = ipvt'range.

  -- ON ENTRY :
  --   s        coefficient system of circuits;
  --   v        vector with approximate values for a solution to s;
  --   xr       work space allocated for the real parts of v;
  --   xi       work space allocated for the imaginary parts of v;
  --   ipvt     pivoting vector for the LU factorization;
  --   verbose  flag for extra output.

  -- ON RETURN :
  --   s        s.fx contains the update value to v, if info = 0,
  --            otherwise s.fx contains the function value of s at v,
  --            s.jm contains the LU factorization of the Jacobian
  --            at the given v;
  --   xr       real parts of the complex numbers in v;
  --   xi       imaginary parts of the complex numbers in v;
  --   ipvt     pivoting information of the LU factorization;
  --   info     if nonzero, then the Jacobian may be singular;
  --   res      residual, max norm of the function value at v;
  --   err      max norm of the update vector to v.

  begin
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.EvalDiff(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    if verbose then
      put_line("The function value : "); put_line(s.fx);
      put("The residual :"); put(res,3); new_line;
    end if;
    Standard_Complex_Linear_Solvers.lufac(s.jm,s.dim,ipvt,info);
    if info /= 0 then
      if verbose
       then put("info : "); put(info,1); put_line(" singular Jacobian?");
      end if;
    else
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      if verbose then
        put_line("The update : "); put_line(s.fx);
        put("Forward error :"); put(err,3); new_line;
      end if;
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Steps
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton steps on the system s starting at the vector v.
 
  -- REQUIRED : The system is square.

    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,err : double_float;
    ans : character;

  begin
    put_line("The vector v :"); put_line(v);
    loop
      LU_Newton_Step(s,v,xr,xi,ipvt,info,res,err,true);
      put_line("The vector v :"); put_line(v);
      put("Another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end LU_Newton_Steps;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions,
  --   makes a coefficient circuit and then does a Newton step.

    p : Link_to_Poly_Sys;
    sols,ptr : Solution_List;
    nbq,len,dim : integer32 := 0;
    s : Link_to_System;
    ls : Link_to_Solution;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols);
    nbq := p'last;
    len := integer32(Length_Of(sols));
    if len = 0 then
      put_line("No solutions found on file.");
    else
      dim := Head_Of(sols).n;
      new_line;
      put("Read system of "); put(nbq,1); put(" polynomials and ");
      put(len,1); put(" solutions in dimension "); put(dim,1); put_line(".");
      s := Make_Coefficient_System(p);
      ptr := sols;
      for k in 1..len loop
        put("Running Newton's method on solution ");
        put(k,1); put_line(" ...");
        ls := Head_Of(ptr);
        LU_Newton_Steps(s,ls.v);
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        ptr := Tail_Of(ptr);
      end loop;
    end if;
  end Main;

begin
  Main;
end ts_newcirc;
