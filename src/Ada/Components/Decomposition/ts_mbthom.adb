with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Witness_Sets_io;                    use Witness_Sets_io;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;

procedure ts_mbthom is

-- DESCRIPTION :
--   Test on membership with homotopy.

  procedure Main is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural;
    restol : constant double_float := 1.0E-10;
    homtol : constant double_float := 1.0E-6;

  begin
    new_line;
    put_line("Membership test with homotopy :");
    put_line("  Input : embedded polynomial system with generic points, and");
    put_line("          list of test points.");
    put_line("  Output : decision whether test point lies on component.");
    Standard_Read_Embedding(lp,genpts,dim);
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
  end Main;

begin
  Main;
end ts_mbthom;
