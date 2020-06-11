with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Point_Lists;
with DoblDobl_Point_Lists;
with QuadDobl_Point_Lists;
with Standard_Select_Solutions;
with Standard_Condition_Report;
with DoblDobl_Condition_Report;
with QuadDobl_Condition_Report;
with Drivers_for_Failed_Paths;

package body Drivers_for_Condition_Tables is

  procedure Standard_Read_and_Compute_Condition_Tables is

    use Standard_Complex_Solutions;
    use Standard_Condition_Report;

    sols : Solution_List;
    n,m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Condition Tables for Solution Lists");
    new_line;
    Read(sols);
    m := Length_Of(sols);
    if m > 0
     then n := natural32(Head_Of(sols).n);
     else n := 0;
    end if;
    new_line;
    put("Read list of "); put(m,1);
    put(" solution vectors of length "); put(n,1); put_line(".");
    if m > 0 then
      put("Do you wish to see all diagnostics ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Diagnostics(sols);
      end if;
      Compute_Condition_Tables(sols);
    end if;
  end Standard_Read_and_Compute_Condition_Tables;

  procedure DoblDobl_Read_and_Compute_Condition_Tables is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Condition_Report;

    sols : Solution_List;
    n,m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Condition Tables for Solution Lists");
    new_line;
    Read(sols);
    m := Length_Of(sols);
    if m > 0
     then n := natural32(Head_Of(sols).n);
     else n := 0;
    end if;
    new_line;
    put("Read list of "); put(m,1);
    put(" solution vectors of length "); put(n,1); put_line(".");
    if m > 0 then
      put("Do you wish to see all diagnostics ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Diagnostics(sols);
      end if;
      Compute_Condition_Tables(sols);
    end if;
  end DoblDobl_Read_and_Compute_Condition_Tables;

  procedure QuadDobl_Read_and_Compute_Condition_Tables is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Condition_Report;

    sols : Solution_List;
    n,m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Condition Tables for Solution Lists");
    new_line;
    Read(sols);
    m := Length_Of(sols);
    if m > 0
     then n := natural32(Head_Of(sols).n);
     else n := 0;
    end if;
    new_line;
    put("Read list of "); put(m,1);
    put(" solution vectors of length "); put(n,1); put_line(".");
    if m > 0 then
      put("Do you wish to see all diagnostics ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Diagnostics(sols);
      end if;
      Compute_Condition_Tables(sols);
    end if;
  end QuadDobl_Read_and_Compute_Condition_Tables;

  procedure Interactive_Driver_to_Scan_Solution_Lists is

    use Standard_Point_Lists;
    use Standard_Select_Solutions;
    use Standard_Condition_Report;

    file : file_type;
    m,n : natural32;
    bannered : boolean;
    fail : boolean := true;
    tol_real : constant double_float := 1.0E-10;
    tol_clus : constant double_float := 1.0E-12;
    e,c,r : Standard_Natural_Vectors.Vector(0..15);
    nbreal,i_end,f_val : natural32;
    pl : Point_List;

  begin
    new_line;
    put_line("Scanning solution lists and computing condition tables.");
    new_line;
    loop
      put_line("Reading the name of the input file...");
      Read_Name_and_Open_File(file);
      Scan_Banner_Dimensions(file,m,n,bannered,fail);
      exit when not fail;
      close(file);
      put_line("Incorrect format or wrong file.  Please try again...");
    end loop;
    Scan_for_Condition_Tables
      (file,standard_output,bannered,false,m,n,tol_real,tol_clus,
       i_end,f_val,e,c,r,nbreal,pl);
    Write_Condition_Results
     (standard_output,i_end,f_val,e,c,r,nbreal,tol_real);
    Write_Cluster_Report
     (file,standard_output,bannered,false,pl,i_end-1,tol_clus);
  end Interactive_Driver_to_Scan_Solution_Lists;

  procedure Standard_Scan_Solution_Lists
              ( infilename,outfilename : in string ) is

    use Standard_Point_Lists;
    use Standard_Select_Solutions;
    use Standard_Condition_Report;

    infile,outfile : file_type;
    m,n : natural32 := 0;
    bannered : boolean;
    fail : boolean := true;
    retry : boolean := false;
    tol_real : constant double_float := 1.0E-10;
    tol_clus : constant double_float := 1.0E-12;
    e,c,r : Standard_Natural_Vectors.Vector(0..15);
    nbreal,i_end,f_val : natural32;
    pl : Point_List;

  begin
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      Scan_Banner_Dimensions(infile,m,n,bannered,fail);
      exit when not fail;
      close(infile);
      put_line("Incorrect format or wrong file.  Please try again...");
      retry := true;
    end loop;
    if outfilename = "" then
      Scan_for_Condition_Tables
        (infile,standard_output,bannered,false,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results
        (standard_output,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report
        (infile,standard_output,bannered,false,pl,i_end-1,tol_clus);
    else
      new_line;
      put_line("Creating file " & outfilename & "...");
      Create_Output_File(outfile,outfilename);
      Scan_for_Condition_Tables
        (infile,outfile,bannered,true,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results(outfile,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report(infile,outfile,bannered,true,pl,i_end-1,tol_clus);
      close(outfile);
    end if;
    close(infile);
  end Standard_Scan_Solution_Lists;

  procedure DoblDobl_Scan_Solution_Lists
              ( infilename,outfilename : in string ) is

    use DoblDobl_Point_Lists;
    use Standard_Select_Solutions;
    use DoblDobl_Condition_Report;

    infile,outfile : file_type;
    m,n : natural32 := 0;
    bannered : boolean;
    fail : boolean := true;
    retry : boolean := false;
    tol_real : constant double_float := 1.0E-10;
    tol_clus : constant double_float := 1.0E-12;
    e,c,r : Standard_Natural_Vectors.Vector(0..30);
    nbreal,i_end,f_val : natural32;
    pl : Point_List;

  begin
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      Scan_Banner_Dimensions(infile,m,n,bannered,fail);
      exit when not fail;
      close(infile);
      put_line("Incorrect format or wrong file.  Please try again...");
      retry := true;
    end loop;
    if outfilename = "" then
      Scan_for_Condition_Tables
        (infile,standard_output,bannered,false,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results
        (standard_output,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report
        (infile,standard_output,bannered,false,pl,i_end-1,tol_clus);
    else
      new_line;
      put_line("Creating file " & outfilename & "...");
      Create_Output_File(outfile,outfilename);
      Scan_for_Condition_Tables
        (infile,outfile,bannered,true,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results(outfile,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report(infile,outfile,bannered,true,pl,i_end-1,tol_clus);
      close(outfile);
    end if;
    close(infile);
  end DoblDobl_Scan_Solution_Lists;

  procedure QuadDobl_Scan_Solution_Lists
              ( infilename,outfilename : in string ) is

    use QuadDobl_Point_Lists;
    use Standard_Select_Solutions;
    use QuadDobl_Condition_Report;

    infile,outfile : file_type;
    m,n : natural32 := 0;
    bannered : boolean;
    fail : boolean := true;
    retry : boolean := false;
    tol_real : constant double_float := 1.0E-10;
    tol_clus : constant double_float := 1.0E-12;
    e,c,r : Standard_Natural_Vectors.Vector(0..60);
    nbreal,i_end,f_val : natural32;
    pl : Point_List;

  begin
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      Scan_Banner_Dimensions(infile,m,n,bannered,fail);
      exit when not fail;
      close(infile);
      put_line("Incorrect format or wrong file.  Please try again...");
      retry := true;
    end loop;
    if outfilename = "" then
      Scan_for_Condition_Tables
        (infile,standard_output,bannered,false,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results
        (standard_output,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report
        (infile,standard_output,bannered,false,pl,i_end-1,tol_clus);
    else
      new_line;
      put_line("Creating file " & outfilename & "...");
      Create_Output_File(outfile,outfilename);
      Scan_for_Condition_Tables
        (infile,outfile,bannered,true,m,n,tol_real,tol_clus,
         i_end,f_val,e,c,r,nbreal,pl);
      Write_Condition_Results(outfile,i_end,f_val,e,c,r,nbreal,tol_real);
      Write_Cluster_Report(infile,outfile,bannered,true,pl,i_end-1,tol_clus);
      close(outfile);
    end if;
    close(infile);
  end QuadDobl_Scan_Solution_Lists;

  procedure Main_Driver_to_Scan_Solution_Lists
              ( infilename,outfilename : in string;
                verbose : in integer32 := 0 ) is

    prc,ans : character;

    use Drivers_for_Failed_Paths;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1); put_line(",");
      put_line("in drivers_for_condition_tables."
             & "Main_Driver_to_Scan_Solution_Lists ...");
    end if;
    new_line;
    put_line("Scanning solution lists and computing condition tables.");
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    put("Output file of a path tracker ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Scan_Failed_Paths(infilename,outfilename);
        when '1' => DoblDobl_Scan_Failed_Paths(infilename,outfilename);
        when '2' => QuadDobl_Scan_Failed_Paths(infilename,outfilename);
        when others => null;
      end case;
    else
      new_line;
      case prc is
        when '0' => Standard_Scan_Solution_Lists(infilename,outfilename);
        when '1' => DoblDobl_Scan_Solution_Lists(infilename,outfilename);
        when '2' => QuadDobl_Scan_Solution_Lists(infilename,outfilename);
        when others => null;
      end case;
    end if;
  end Main_Driver_to_Scan_Solution_Lists;

end Drivers_for_Condition_Tables;
