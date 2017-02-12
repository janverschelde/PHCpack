with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Path_Parameters;                   use Path_Parameters;

procedure ts_pathpars is

-- DESCRIPTION :
--   Test on the parameters for the path trackers in the Path library.

  function Prompt_for_Precision return integer32 is

  -- DESCRIPTION :
  --   Asks the user for the working precision.

    res : integer32 := 0;

  begin
    put_line("Reading the working precision ...");
    put("Enter 16, 32, or 64 : ");
    get(res);
    return res;
  end Prompt_for_Precision;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision, shows the default
  --   values for the parameters.

    prec : constant integer32 := Prompt_for_Precision;
    pars : Parameters := Default_Parameters(prec);

  begin
    Tune(pars);
  end Main;

begin
  Main;
end ts_pathpars;
