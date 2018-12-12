with text_io;                            use text_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Standard_SeriesPade_Tracker;

procedure ts_nxtpadsol is

-- DESCRIPTION :
--   Interactive test on the get_next() method to track paths
--   with a series Pade predictor.

  procedure Main is

  -- DESCRIPTION :
  --   Allows the user to tune the homotopy continuation parameters,
  --   prompts for a target, start system, and start solutions.

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin
    new_line;
    put_line("Tuning the homotopy continuation parameters ...");
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    Standard_SeriesPade_Tracker.Init(pars);
    new_line;
    put_line("The stored values of the parameters : ");
    new_line;
    Homotopy_Continuation_Parameters_io.put
      (Standard_SeriesPade_Tracker.Get_Parameters.all);
  end Main;

begin
  Main;
end ts_nxtpadsol;
