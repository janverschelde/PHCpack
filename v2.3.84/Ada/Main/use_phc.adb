with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with PHCpack;

procedure use_phc is

  infile,outfile : file_type;        -- input and output file
  p,q : Link_to_Poly_Sys;            -- target and start system
  mixed_volume : natural32;          -- root count is mixed volume
  sols : Solution_List;              -- list of solutions

begin
  Open(infile,in_file,"test.in");
  get(infile,p);
  Create(outfile,out_file,"test.out");
  put(outfile,p.all);
  q := new Poly_Sys(p'range);
  PHCpack.Static_Lifting(outfile,p.all,mixed_volume,q.all,sols);
  PHCpack.Artificial_Parameter_Continuation(outfile,p.all,q.all,sols);
  PHCpack.Refine_Roots(outfile,p.all,sols);
end use_phc;
