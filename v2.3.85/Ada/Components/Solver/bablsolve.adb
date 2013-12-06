with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Drivers_to_Eqn_by_Eqn_Solvers;      use Drivers_to_Eqn_by_Eqn_Solvers;

procedure bablsolve ( p : in Poly_Sys ) is

  file : file_type;
  name : Link_to_String;

begin
  new_line;
  put_line("Calling the equation-by-equation solver ...");
  new_line;
  put_line("Reading the name of the output file ...");
  Read_Name_and_Create_File(file,name);
  Shuffle_Polynomials_and_Solve(file,name.all,p);
end bablsolve;
