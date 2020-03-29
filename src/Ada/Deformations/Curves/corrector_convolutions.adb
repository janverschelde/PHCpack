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

-- MANAGEMENT OF LEADING COEFFICIENTS :

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := Standard_Complex_Numbers.Create(0.0);
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := DoblDobl_Complex_Numbers.Create(integer(0));
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := QuadDobl_Complex_Numbers.Create(integer(0));
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_Vectors.Link_to_Vector;
                c : in Standard_Speelpenning_Convolutions.Link_to_Circuit ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_Vectors.Link_to_Vector;
                c : in  DoblDobl_Speelpenning_Convolutions.LInk_to_Circuit ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_Vectors.Link_to_Vector;
                c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Allocate_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    cff : Standard_Complex_VecVecs.VecVec(c'range);

    use Standard_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : Standard_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => Standard_Complex_Numbers.Create(0.0));
          cff(k) := new Standard_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new Standard_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Allocate_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : out DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    cff : DoblDobl_Complex_VecVecs.VecVec(c'range);

    use DoblDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : DoblDobl_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => DoblDobl_Complex_Numbers.Create(integer(0)));
          cff(k) := new DoblDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new DoblDobl_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Allocate_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : out QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    cff : QuadDobl_Complex_VecVecs.VecVec(c'range);

    use QuadDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : QuadDobl_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => QuadDobl_Complex_Numbers.Create(integer(0)));
          cff(k) := new QuadDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new QuadDobl_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : in Standard_Complex_VecVecs.Link_to_VecVec ) is

    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant Standard_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    use DoblDobl_Complex_Vectors,DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant DoblDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant QuadDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_VecVecs.Link_to_VecVec;
                c : in Standard_Speelpenning_Convolutions.Circuits ) is

    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant Standard_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits ) is

    use DoblDobl_Complex_Vectors,DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant DoblDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits ) is

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant QuadDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

-- EVALUATION OF STEP SIZE :

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

-- RUNNING NEWTON'S METHOD :

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
