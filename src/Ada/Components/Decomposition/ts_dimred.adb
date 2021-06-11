with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Substitutors;     use Standard_Complex_Substitutors;
with Witness_Sets;                      use Witness_Sets;

procedure ts_dimred is

-- DESCRIPTION :
--   This is a test on "dimension reduction", which stands for
--   the "reducing the cost of determining the dimension of a
--   solution set of polynomial systems".  Usuall we cut with
--   a generic hyperplane, but we may as well substitute one or
--   more variables with random constants, which corresponds to
--   taking very special slices.  Note that after that substitution
--   we need to take random combinations to form a square system.
--   This is a very primitive facility: the routine substitutes
--   the last variable and collapses thereafter.
--   Apply repeatedly to substitute more than one variable.

  function Collapse ( p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Adds a random times the last polynomial to the other
  --   polynomials in the system.

    res : Poly_Sys(p'first..p'last-1);

  begin
    for i in res'range loop
      res(i) := Random1*p(p'last);
      Add(res(i),p(i));
    end loop;
    return res;
  end Collapse;

  procedure Main is

    lp : Link_to_Poly_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Reading some polynomial system.");
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put(file,integer32(lp'length),1);
    new_line(file);
    put(file,lp.all);
    declare
      sp : Poly_Sys(lp'range) := Substitute(lp'last,Random1,lp.all);
    begin
      new_line(file);
      put(file,"The system after substituting variable ");
      put(file,lp'last,1); put_line(file," :");
      put_line(file,sp);
      new_line(file);
      put_line(file,"The collapsed system : ");
      put_line(file,Collapse(sp));
      put_line(file,"The embedded system : ");
      put_line(file,Embed(sp));
    end;
    close(file);
  end Main;

begin
  Main;
end ts_dimred;
