with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;
with Witness_Sets_io;                    use Witness_Sets_io;
with Drivers_to_Factor_Components;       use Drivers_to_Factor_Components;

procedure ts_filfac is

-- DESCRIPTION :
--   Interactive development of the main filtering and factoring module.

  procedure Homotopy_Membership_Test is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
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
  end Homotopy_Membership_Test;

  procedure Driver_to_Monodromy_Breakup is

  -- DESCRIPTION :
  --   Reads the embedded polynomial system and calls the monodromy
  --   breakup algorithm.

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
  end Driver_to_Monodromy_Breakup;

  procedure Trace_Form_Interpolation is

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to certify monodromy breakup with by interpolation :");
    put_line
      ("  1. on given decomposition : use bootstrapping Newton to certify;");
    put_line
      ("  2.                        : use full trace form to certify;");
    put_line
      ("  3.                        : use Newton identities on trace form;");
    put_line
      ("  4.                        : use linear trace only to certify.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Call_Newton_Interpolate(file,lp.all,sols,dim);
      when '2' => Call_Trace_Interpolate(file,lp.all,sols,dim);
      when '3' => Call_Power_Trace_Interpolate(file,lp.all,sols,dim);
      when '4' => Call_Linear_Trace_Interpolate(file,lp.all,sols,dim);
      when others => null;
    end case;
  end Trace_Form_Interpolation;

  procedure Incremental_Interpolation is

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    ans : character;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to decompose with incremental use of interpolation :");
    put_line("  1. massive interpolate with standard arithmetic;");
    put_line("  2. incremental interpolate with standard arithmetic;");
    put_line("  3.  + determine span with standard arithmetic;");
    put_line("  4.  + use central projections;");
    put_line("  5. massive interpolate with multi-precision arithmetic;");
    put_line("  6. incremental interpolate with multi-precision arithmetic;");
    put_line("  7.  + determine span with multi-precision arithmetic;");
    put_line("  8.  + use central projections;");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Call_Standard_Interpolate(file,lp.all,sols,dim,0);
      when '2' => Call_Standard_Interpolate(file,lp.all,sols,dim,1);
      when '3' => Call_Standard_Interpolate(file,lp.all,sols,dim,2);
      when '4' => Call_Standard_Interpolate(file,lp.all,sols,dim,3);
      when '5' => Call_Multprec_Interpolate(file,lp.all,sols,dim,0);
      when '6' => Call_Multprec_Interpolate(file,lp.all,sols,dim,1);
      when '7' => Call_Multprec_Interpolate(file,lp.all,sols,dim,2);
      when '8' => Call_Multprec_Interpolate(file,lp.all,sols,dim,3);
      when others => null;
    end case;
  end Incremental_Interpolation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to filter junk and factor components :");
    put_line("  1. filter junk with homotopy membership test;");
    put_line("  2. apply monodromy breakup on filtered witness point set;");
    put_line("  3. given breakup partition, compute trace form of filter;");
    put_line("  4. perform tasks 1, 2, and 3 by incremental interpolation.");
    put("Type 1, 2, 3, or 4 to select a task : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Homotopy_Membership_Test;
      when '2' => Driver_to_Monodromy_Breakup;
      when '3' => Trace_Form_Interpolation;
      when '4' => Incremental_Interpolation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_filfac;
