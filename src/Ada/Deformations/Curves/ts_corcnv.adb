with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
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
with Corrector_Convolutions;             use Corrector_Convolutions;

procedure ts_corcnv is

-- DESCRIPTION :
--   Development of the corrector convolutions.

  function Step_Coefficient
              ( c : Standard_Complex_Vectors.Vector; t : double_float )
              return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Given in c is the coefficient vector of a series
  --   and in t is some step size, t > 0.
  --   Returns the value of the series c evaluated at t.

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

  -- DESCRIPTION :
  --   Given in c is the coefficient vector of a series
  --   and in t is some step size, t > 0.
  --   Returns the value of the series c evaluated at t.

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

  -- DESCRIPTION :
  --   Given in c is the coefficient vector of a series
  --   and in t is some step size, t > 0.
  --   Returns the value of the series c evaluated at t.

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

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.

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

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.

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

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.

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

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.
 
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                t : in double_double ) is

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.
 
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                t : in quad_double ) is

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t.
 
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Standard_Run_Newton
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in double precision.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution := Head_Of(sols);
    dx : Standard_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info,nbrit : integer32;
    ans : character;
    maxit : constant integer32 := 4;
    tol : constant double_float := 1.0E-12;
    nrm : double_float;
    fail : boolean;

  begin
    while not Is_Null(tmp) loop
      ls := Head_of(tmp);
      LU_Newton_Steps
        (standard_output,chom,ls.v,maxit,nbrit,tol,nrm,dx,ipvt,info,fail);
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Run_Newton;

  procedure DoblDobl_Run_Newton
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in double double precision.

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution := Head_Of(sols);
    dx : DoblDobl_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info,nbrit : integer32;
    ans : character;
    maxit : constant integer32 := 4;
    tol : constant double_float := 1.0E-24;
    nrm : double_double;
    fail : boolean;

  begin
    while not Is_Null(tmp) loop
      ls := Head_of(tmp);
      LU_Newton_Steps
        (standard_output,chom,ls.v,maxit,nbrit,tol,nrm,dx,ipvt,info,fail);
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Run_Newton;

  procedure QuadDobl_Run_Newton
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy in chom,
  --   starting at the solutions in sols, in quad double precision.

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution := Head_Of(sols);
    dx : QuadDobl_Complex_Vectors.Vector(ls.v'range);
    ipvt : Standard_Integer_Vectors.Vector(dx'range);
    info,nbrit : integer32;
    ans : character;
    maxit : constant integer32 := 4;
    tol : constant double_float := 1.0E-48;
    nrm : quad_double;
    fail : boolean;

  begin
    while not Is_Null(tmp) loop
      ls := Head_of(tmp);
      LU_Newton_Steps
        (standard_output,chom,ls.v,maxit,nbrit,tol,nrm,dx,ipvt,info,fail);
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Run_Newton;

  procedure Standard_Test_Correction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;
    step : double_float := 0.0;

  begin
    Standard_Homotopy_Convolutions_io.get(2,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.Standard_Check_Solutions(cnvhom,sols);
    end if;
    put("Evaluate for a positive step ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the step : "); get(step);
      Step_Coefficient(cnvhom,step);
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
    step : double_float := 0.0;
    ddstep : double_double;

  begin
    DoblDobl_Homotopy_Convolutions_io.get(2,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.DoblDobl_Check_Solutions(cnvhom,sols);
    end if;
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Evaluate for a positive step ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the step : "); get(step);
      ddstep := Create(step);
      Step_Coefficient(cnvhom,ddstep);
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
    step : double_float := 0.0;
    qdstep : quad_double;

  begin
    QuadDobl_Homotopy_Convolutions_io.get(2,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.QuadDobl_Check_Solutions(cnvhom,sols);
    end if;
    put("Evaluate for a positive step ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the step : "); get(step);
      qdstep := Create(step);
      Step_Coefficient(cnvhom,qdstep);
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
