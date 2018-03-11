with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Solution_Manipulators;
with DoblDobl_Solution_Manipulators;
with QuadDobl_Solution_Manipulators;
with Prompt_for_Systems;
with Black_Box_Solvers;
with Witness_Sets_io;                    use Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Greeting_Banners;
with Write_Seed_Number;
with Path_Counts_Table;
with Cascade_Homotopies;                 use Cascade_Homotopies;
with Cascade_Homotopy_Filters;           use Cascade_Homotopy_Filters;
with Running_Cascades;                   use Running_Cascades;

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
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,0,1.0E-8);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,0,1.0E-8);
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
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,0,1.0E-8);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,0,1.0E-8);
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
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,0,1.0E-8);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,0,1.0E-8);
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

  procedure Prompt_for_Top_Dimension
              ( nq,nv : in natural32; topdim,lowdim : out natural32 ) is
  begin
    if nq >= nv
     then lowdim := 0;  -- allow no negative values for lower bound
     else lowdim := nv-nq;
    end if;
    loop
      put("The number of equations : "); put(nq,1); new_line;
      put("The number of variables : "); put(nv,1); new_line;
      put("-> the default, largest top dimension is "); put(nv-1,1);
      put_line(" ...");
      put("Give the expected top dimension : ");
      Numbers_io.Read_Natural(topdim);
      exit when (topdim < nv) and (topdim >= lowdim); -- check bounds
      if topdim >= nv then
        put("Error: The top dimension cannot be larger than ");
        put(nv-1,1); put_line(".");
      end if;
      if topdim < lowdim then
        put("Error: The top dimension should be at least "); 
        put(lowdim,1); put_line(".");
      end if;
      put("Please enter a number between "); put(lowdim,1);
      put(" and "); put(nv-1,1); put_line(".");
      put("The default, largest top dimension is ");
      put(nv-1,1); put_line(".");
    end loop;
  end Prompt_for_Top_Dimension;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        Standard_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        Standard_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        Standard_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        Standard_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        DoblDobl_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        DoblDobl_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        QuadDobl_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        QuadDobl_Run_Cascade(nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Solutions;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

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
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all,true,true);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all,true,true);
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
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all,true,true);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all,true,true);
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
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all,true,true);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all,true,true);
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
