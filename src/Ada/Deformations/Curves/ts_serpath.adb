with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials_io;  use Complex_Series_and_Polynomials_io;
with Series_and_Homotopies;
with Series_and_Trackers;
with Homotopy_Series_Readers;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;

procedure ts_serpath is

-- DESCRIPTION :
--   Developing path tracers with power series.

  procedure Set_Output
              ( file : in out file_type; verbose,tofile : out boolean ) is

  -- DESCRIPTION :
  --   Prompts the user if verbose or not, and if so, whether the output
  --   should be written to a file or not, which sets tofile to true.
  --   If tofile, then file is created, ready for output.

    ans : character;

  begin
    put("Verbose?  Want to see extra output ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
    tofile := (ans = 'y');
    if tofile then
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
    end if;
  end Set_Output;

  procedure Prompt_for_Degrees ( numdeg,dendeg : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the degrees of the numerator and denominator,
  --   returned respectively in numdeg and dendeg.

  begin
    new_line;
    put("Give the degree of numerator of the Pade approximants : ");
    numdeg := 0; get(numdeg);
    put("Give the degree of denominator of the Pade approximants : ");
    dendeg := 0; get(dendeg);
  end Prompt_for_Degrees;

  procedure Write_Timer
              ( file : in file_type;
                numdeg,dendeg,precision : in natural32;
                timer : in Timing_Widget ) is

  -- DESCRIPTION :
  --   Writes the times with as well the degrees of numerator
  --   and denominator of the Pade approximants.
  --   The precision is 0, 1, or 2, respectively
  --   for double, double double, or quad double precision.

    s : constant string
      := "[" & Characters_and_Numbers.Convert(integer32(numdeg))
       & "," & Characters_and_Numbers.Convert(integer32(dendeg))
       & "]-Tracking";
 
  begin
    new_line(file);
    case precision is
      when 0 => print_times(file,timer,s & " in double precision.");
      when 1 => print_times(file,timer,s & " in double double precision.");
      when 2 => print_times(file,timer,s & " in quad double precision.");
      when others => null;
    end case;
  end Write_Timer;

  procedure Standard_Test
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in standard double precision.

    use Standard_Complex_Solutions;

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
   -- numdeg,dendeg : integer32 := 0;
    timer : Timing_Widget;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
   -- Prompt_for_Degrees(numdeg,dendeg);
    Homotopy_Continuation_Parameters_io.Tune(p);
    Set_Output(file,verbose,tofile);
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,0,timer);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,0,timer);
    end if;
    Standard_Complex_Poly_Systems.Clear(h);
    Standard_CSeries_Poly_Systems.Clear(s);
  end Standard_Test;

  procedure DoblDobl_Test
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
   -- numdeg,dendeg : integer32 := 0;
    timer : Timing_Widget;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
   -- Prompt_for_Degrees(numdeg,dendeg);
    Homotopy_Continuation_Parameters_io.Tune(p);
    Set_Output(file,verbose,tofile);
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Continue to the next path ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,1,timer);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,1,timer);
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(h);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
   -- numdeg,dendeg : integer32 := 0;
    timer : Timing_Widget;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
   -- Prompt_for_Degrees(numdeg,dendeg);
    Homotopy_Continuation_Parameters_io.Tune(p);
    Set_Output(file,verbose,tofile);
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,2,timer);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,2,timer);
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(h);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
    new_line;
    Standard_Test(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
    new_line;
    DoblDobl_Test(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
    new_line;
    QuadDobl_Test(nbeq,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpath;
