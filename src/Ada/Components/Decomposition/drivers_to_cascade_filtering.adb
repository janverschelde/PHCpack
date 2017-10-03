with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Prompt_for_Systems;
with Black_Box_Solvers;
with Witness_Sets_io;                    use Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Greeting_Banners;
with Write_Seed_Number;
with Cascade_Homotopies;                 use Cascade_Homotopies;

package body Drivers_to_Cascade_Filtering is

  procedure Standard_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Laur_Poly_Convertors;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      Standard_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         Standard_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end Standard_Witness_Generate;

  procedure DoblDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Laur_Poly_Convertors;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      DoblDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         DoblDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end DoblDobl_Witness_Generate;

  procedure QuadDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Laur_Poly_Convertors;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      QuadDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         QuadDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end QuadDobl_Witness_Generate;

  procedure Driver_to_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Witness_Generate(nt,inpname,outname);
      when '1' => DoblDobl_Witness_Generate(nt,inpname,outname);
      when '2' => QuadDobl_Witness_Generate(nt,inpname,outname);
      when others => null;
    end case;
  end Driver_to_Witness_Generate;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    null;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    null;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use Standard_Laur_Poly_Convertors;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use DoblDobl_Laur_Poly_Convertors;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use QuadDobl_Laur_Poly_Convertors;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end QuadDobl_Embed_and_Cascade;

  procedure Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    prc : constant character := Prompt_for_Precision;

  begin
    case prc is
      when '0' => Standard_Embed_and_Cascade(nt,inpname,outname);
      when '1' => DoblDobl_Embed_and_Cascade(nt,inpname,outname);
      when '2' => QuadDobl_Embed_and_Cascade(nt,inpname,outname);
      when others => null;
    end case;
  end Embed_and_Cascade;

end Drivers_to_Cascade_Filtering;
