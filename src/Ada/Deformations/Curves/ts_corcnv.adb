with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
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

  procedure Newton_Step
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Does one Newton step in standard double precision.

    A : Standard_Complex_Matrices.Matrix(hom.crc'range,sol'range);
    y : Standard_Complex_Vectors.Vector(hom.crc'range);

  begin
    null;
  end Newton_Step;

  procedure Newton_Step
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Does one Newton step in double double precision.

    A : DoblDobl_Complex_Matrices.Matrix(hom.crc'range,sol'range);
    y : DoblDobl_Complex_Vectors.Vector(hom.crc'range);

  begin
    null;
  end Newton_Step;

  procedure Newton_Step
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Does one Newton step in quad double precision.

    A : QuadDobl_Complex_Matrices.Matrix(hom.crc'range,sol'range);
    y : QuadDobl_Complex_Vectors.Vector(hom.crc'range);

  begin
    null;
  end Newton_Step;

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
    Newton_Step(cnvhom,Standard_Complex_Solutions.Head_Of(sols).v);
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
    Newton_Step(cnvhom,DoblDobl_Complex_Solutions.Head_Of(sols).v);
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
    Newton_Step(cnvhom,QuadDobl_Complex_Solutions.Head_Of(sols).v);
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
