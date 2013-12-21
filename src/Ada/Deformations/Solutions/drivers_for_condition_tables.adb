with Communications_with_User;           use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Point_Lists;               use Standard_Point_Lists;
with Standard_Condition_Report;          use Standard_Condition_Report;

package body Drivers_for_Condition_Tables is

  procedure Read_and_Compute_Condition_Tables is

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
  end Read_and_Compute_Condition_Tables;

  procedure Read_Dimensions
              ( file : in file_type; len,dim : out natural32;
                fail : out boolean ) is

  begin
    fail := true;
    len := 0; get(file,len);
    dim := 0; get(file,dim);
    fail := false;
  exception
    when others
      => put_line("Something bad happened while reading dimensions.");
         fail := true;
  end Read_Dimensions;

  procedure Prompt_to_Scan_Banner
              ( infile : in file_type; bannered,fail : out boolean ) is

    ans : character;
    found : boolean;

  begin
    new_line;
    put("Are the solutions preceeded by a system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      bannered := true;
      put_line("Scanning for THE SOLUTIONS banner...");
      File_Scanning.Scan_and_Skip(infile,"THE SOLUTIONS",found);
      if found then
        put_line("  found banner, ready to continue reading dimensions...");
        fail := false;
      else
        put_line("  did not find banner, format of file maybe wrong...");
        fail := true;
      end if;
    else
      bannered := false;
      fail := false;
    end if;
  end Prompt_to_Scan_Banner;

  procedure Scan_Banner_Dimensions
              ( infile : in file_type; len,dim : out natural32;
                bannered,fail : out boolean ) is
  begin
    Prompt_to_Scan_Banner(infile,bannered,fail);
    if not fail
     then Read_Dimensions(infile,len,dim,fail);
    end if;
  end Scan_Banner_Dimensions;

  procedure Interactive_Driver_to_Scan_Solution_Lists is

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

  procedure Main_Driver_to_Scan_Solution_Lists
              ( infilename,outfilename : in string ) is

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
    new_line;
    put_line("Scanning solution lists and computing condition tables.");
    new_line;
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
  end Main_Driver_to_Scan_Solution_Lists;

end Drivers_for_Condition_Tables;
