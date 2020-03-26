with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy_Convolutions_io;
with Test_Predictor_Convolutions;

procedure ts_corcnv is

-- DESCRIPTION :
--   Development of the corrector convolutions.

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   in double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose;
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular.

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

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   in double double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose;
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular.

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

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   in quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose;
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular.

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

  procedure Standard_Run_Newton
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in double precision.

    ls : constant Standard_Complex_Solutions.Link_to_Solution
       := Standard_Complex_Solutions.Head_Of(sols);
    dx : Standard_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info : integer32;
    ans : character;

  begin
    loop
      LU_Newton_Step(standard_output,chom,ls.v,dx,ipvt,info);
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Run_Newton;

  procedure DoblDobl_Run_Newton
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in double double precision.

    ls : constant DoblDobl_Complex_Solutions.Link_to_Solution
       := DoblDobl_Complex_Solutions.Head_Of(sols);
    dx : DoblDobl_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info : integer32;
    ans : character;

  begin
    loop
      LU_Newton_Step(standard_output,chom,ls.v,dx,ipvt,info);
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Run_Newton;

  procedure QuadDobl_Run_Newton
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in quad double precision.

    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := QuadDobl_Complex_Solutions.Head_Of(sols);
    dx : QuadDobl_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info : integer32;
    ans : character;

  begin
    loop
      LU_Newton_Step(standard_output,chom,ls.v,dx,ipvt,info);
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Run_Newton;

  procedure Standard_Test_Correction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    Standard_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.Standard_Check_Solutions(cnvhom,sols);
    end if;
    Standard_Run_Newton(cnvhom,sols);
  end Standard_Test_Correction;

  procedure DoblDobl_Test_Correction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    DoblDobl_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.DoblDobl_Check_Solutions(cnvhom,sols);
    end if;
    DoblDobl_Run_Newton(cnvhom,sols);
  end DoblDobl_Test_Correction;

  procedure QuadDobl_Test_Correction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    QuadDobl_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.QuadDobl_Check_Solutions(cnvhom,sols);
    end if;
    QuadDobl_Run_Newton(cnvhom,sols);
  end QuadDobl_Test_Correction;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Test_Correction;
      when '1' => DoblDobl_Test_Correction;
      when '2' => QuadDobl_Test_Correction;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_corcnv;
