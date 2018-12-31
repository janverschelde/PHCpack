with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Tracked_Solutions_io;
with DoblDobl_Tracked_Solutions_io;
with QuadDobl_Tracked_Solutions_io;
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with DoblDobl_Solution_Filters;          use DoblDobl_Solution_Filters;
with QuadDobl_Solution_Filters;          use QuadDobl_Solution_Filters;

package body Drivers_for_Failed_Paths is

  procedure Prompt_for_File
              ( file : in out file_type; len : in natural32;
                tofile : out boolean ) is

    ans : character := 'n';
 
  begin
    put("Selected "); put(len,1); put_line(" failed solutions.");
    if len > 0 then
      new_line;
      put("Write start solutions corresponding to failed paths ");
      put("to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of an output file ...");
        Read_Name_and_Create_File(file);
      end if;
    end if;
    tofile := (ans = 'y');
  end Prompt_for_File;

  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put(file,len,natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put(file,len,natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( file : in file_type; len : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put(file,len,natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( file : in file_type; len : natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put(file,len,natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( file : in file_type; len : natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    put(file,natural32(p'last),p);
    new_line(file);
    put_line(file,
      "TITLE : start system with solutions corresponding to failed paths");
    new_line(file);
    put(file,len,natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end Write_to_File;

  procedure Write_to_File
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    tofile : boolean;
 
  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    tofile : boolean;
 
  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    tofile : boolean;

  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    tofile : boolean;

  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    tofile : boolean;
 
  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    file : file_type;
    len : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    tofile : boolean;
 
  begin
    Prompt_for_File(file,len,tofile);
    if tofile then
      Write_to_File(file,len,p,sols);
      close(file);
    end if;
  end Write_to_File;

  procedure Standard_Scan_Failed_Paths
              ( infilename,outfilename : in string ) is

    ans : character;
    laurent : boolean;
    infile,outfile : file_type;
    lps,lqs : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    llp,llq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    psols,qsols,failed : Standard_Complex_Solutions.Solution_List;
    fail : boolean := true;
    retry : boolean := false;
    tol : constant double_float := 1.0E-8;
    len : natural32;

  begin
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    laurent := (ans = 'y');
    new_line;
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      declare
      begin
        fail := false;
        if laurent then
          Standard_Tracked_Solutions_io.get(infile,llp,llq,psols,qsols,true);
        else
          Standard_Tracked_Solutions_io.get(infile,lps,lqs,psols,qsols,true);
        end if;
      exception
        when others => fail := true;
      end;
      if not fail then
        failed := Select_Failed_Solutions(psols,qsols,tol,true);
      else
        put_line("Incorrect format or wrong file.  Please try again...");
        retry := true;
      end if;
      close(infile);
      exit when not fail;
    end loop;
    len := Standard_Complex_Solutions.Length_Of(failed);
    if len > 0 then
      if outfilename = "" then
        if laurent
         then Write_to_File(llq.all,failed); -- prompt for file
         else Write_to_File(lqs.all,failed); -- prompt for file
        end if;
      else
        new_line;
        put_line("Creating file " & outfilename & "...");
        Create_Output_File(outfile,outfilename);
        new_line;
        put_line("See the file " & outfilename & " for results.");
        new_line;
        if laurent
         then Write_to_File(outfile,len,llq.all,failed);
         else Write_to_File(outfile,len,lqs.all,failed);
        end if;
      end if;
    end if;
  end Standard_Scan_Failed_Paths;

  procedure DoblDobl_Scan_Failed_Paths
              ( infilename,outfilename : in string ) is

    ans : character;
    laurent : boolean;
    infile,outfile : file_type;
    lps,lqs : Dobldobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    llp,llq : Dobldobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    psols,qsols,failed : Dobldobl_Complex_Solutions.Solution_List;
    fail : boolean := true;
    retry : boolean := false;
    tol : constant double_float := 1.0E-8;
    len : natural32;

  begin
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    laurent := (ans = 'y');
    new_line;
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      declare
      begin
        fail := false;
        if laurent then
          Dobldobl_Tracked_Solutions_io.get(infile,llp,llq,psols,qsols,true);
        else
          Dobldobl_Tracked_Solutions_io.get(infile,lps,lqs,psols,qsols,true);
        end if;
      exception
        when others => fail := true;
      end;
      if not fail then
        failed := Select_Failed_Solutions(psols,qsols,tol,true);
      else
        put_line("Incorrect format or wrong file.  Please try again...");
        retry := true;
      end if;
      close(infile);
      exit when not fail;
    end loop;
    len := Dobldobl_Complex_Solutions.Length_Of(failed);
    if len > 0 then
      if outfilename = "" then
        if laurent
         then Write_to_File(llq.all,failed); -- prompt for file
         else Write_to_File(lqs.all,failed); -- prompt for file
        end if;
      else
        new_line;
        put_line("Creating file " & outfilename & "...");
        Create_Output_File(outfile,outfilename);
        new_line;
        put_line("See the file " & outfilename & " for results.");
        new_line;
        if laurent
         then Write_to_File(outfile,len,llq.all,failed);
         else Write_to_File(outfile,len,lqs.all,failed);
        end if;
      end if;
    end if;
  end DoblDobl_Scan_Failed_Paths;

  procedure QuadDobl_Scan_Failed_Paths
              ( infilename,outfilename : in string ) is

    ans : character;
    laurent : boolean;
    infile,outfile : file_type;
    lps,lqs : Dobldobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    llp,llq : Dobldobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    psols,qsols,failed : Dobldobl_Complex_Solutions.Solution_List;
    fail : boolean := true;
    retry : boolean := false;
    tol : constant double_float := 1.0E-8;
    len : natural32;

  begin
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    laurent := (ans = 'y');
    new_line;
    loop
      if infilename = "" or retry then
        put_line("Reading the name of the input file...");
        Read_Name_and_Open_File(infile);
      else
        put_line("Opening the file " & infilename & "...");
        Open_Input_File(infile,infilename);
      end if;
      declare
      begin
        fail := false;
        if laurent then
          Dobldobl_Tracked_Solutions_io.get(infile,llp,llq,psols,qsols,true);
        else
          Dobldobl_Tracked_Solutions_io.get(infile,lps,lqs,psols,qsols,true);
        end if;
      exception
        when others => fail := true;
      end;
      if not fail then
        failed := Select_Failed_Solutions(psols,qsols,tol,true);
      else
        put_line("Incorrect format or wrong file.  Please try again...");
        retry := true;
      end if;
      close(infile);
      exit when not fail;
    end loop;
    len := Dobldobl_Complex_Solutions.Length_Of(failed);
    if len > 0 then
      if outfilename = "" then
        if laurent
         then Write_to_File(llq.all,failed); -- prompt for file
         else Write_to_File(lqs.all,failed); -- prompt for file
        end if;
      else
        new_line;
        put_line("Creating file " & outfilename & "...");
        Create_Output_File(outfile,outfilename);
        new_line;
        put_line("See the file " & outfilename & " for results.");
        new_line;
        if laurent
         then Write_to_File(outfile,len,llq.all,failed);
         else Write_to_File(outfile,len,lqs.all,failed);
        end if;
      end if;
    end if;
  end QuadDobl_Scan_Failed_Paths;

end Drivers_for_Failed_Paths;
