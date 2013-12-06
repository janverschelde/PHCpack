with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Homotopy_Evaluator_Packages;

procedure ts_homline is

-- DESCRIPTION :
--   Test the creation of an inline evaluator.

  lp,lq : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Interactive testing of the creation of an inline evaluator.");
  new_line;
  put_line("Reading the target system");
  get(lp);
  new_line;
  put_line("Reading the start system");
  get(lq);
  new_line;
  Homotopy_Evaluator_Packages.Create(lp.all,lq.all);
end ts_homline;
