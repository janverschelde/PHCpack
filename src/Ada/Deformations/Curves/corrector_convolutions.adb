with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;

package body Corrector_Convolutions is

  function Step_Coefficient
              ( c : Standard_Complex_Vectors.Vector; t : double_float )
              return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  function Step_Coefficient
              ( c : DoblDobl_Complex_Vectors.Vector; t : double_double )
              return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  function Step_Coefficient
              ( c : QuadDobl_Complex_Vectors.Vector; t : quad_double )
              return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_float ) is

    use Standard_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_double ) is

    use DoblDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in quad_double ) is

    use QuadDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                t : in double_float ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                t : in double_double ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                t : in quad_double ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    Standard_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    Standard_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    DoblDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    DOblDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true ) is
  begin
    fail := true;
    for k in 1..maxit loop
      LU_Newton_Step(file,hom,sol,dx,ipvt,info,verbose);
      nrmdx := Standard_Complex_Vector_Norms.Max_Norm(dx);
      if verbose then
        put(file,"after step "); put(k,1);
        put(file,", |dx| :"); put(file,nrmdx,3); new_line(file);
      end if;
      if nrmdx <= tol
       then nbrit := k; fail := false; exit;
      end if;
    end loop;
    nbrit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true ) is
  begin
    fail := true;
    for k in 1..maxit loop
      LU_Newton_Step(file,hom,sol,dx,ipvt,info,verbose);
      nrmdx := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      if verbose then
        put(file,"after step "); put(k,1);
        put(file,", |dx| : "); put(file,nrmdx,3); new_line(file);
      end if;
      if nrmdx <= tol
       then nbrit := k; fail := false; exit;
      end if;
    end loop;
    nbrit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true ) is
  begin
    fail := true;
    for k in 1..maxit loop
      LU_Newton_Step(file,hom,sol,dx,ipvt,info,verbose);
      nrmdx := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      if verbose then
        put(file,"after step "); put(k,1);
        put(file,", |dx| : "); put(file,nrmdx,3); new_line(file);
      end if;
      if nrmdx <= tol
       then nbrit := k; fail := false; exit;
      end if;
    end loop;
    nbrit := maxit;
  end LU_Newton_Steps;

end Corrector_Convolutions;
