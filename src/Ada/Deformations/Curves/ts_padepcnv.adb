with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Homotopy;
with Solution_Drops;
with Standard_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Series_and_Homotopies;
with Test_Series_Predictors;
with Standard_Speelpenning_Convolutions;
with System_Convolution_Circuits;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

  procedure Standard_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    c : Standard_Speelpenning_Convolutions.Circuits(1..nq);
    deg : integer32 := 0;

  begin
    put("Give the degree of the series : "); get(deg);
    Complex_Series_and_Polynomials.Set_Degree(s,deg);
    c := System_Convolution_Circuits.Make_Convolution_Circuits(s);
    for k in c'range loop
      Standard_Integer_VecVecs_io.put_line(c(k).xps);
    end loop;
  end Standard_Test_Prediction;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq,idxpar : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      Standard_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        Standard_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end Standard_Main;

begin
  Standard_Main;
end ts_padepcnv;
