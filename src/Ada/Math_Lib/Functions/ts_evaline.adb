with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Evaluator_Packages;

procedure ts_evaline is

-- DESCRIPTION :
--   Test the creation of an inline evaluator.

  lp : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Interactive testing of the creation of an inline evaluator.");
  new_line;
  get(lp);
  new_line;
  Standard_Evaluator_Packages.Create(lp.all);
end ts_evaline;
