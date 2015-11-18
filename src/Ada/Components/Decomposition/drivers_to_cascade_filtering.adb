with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Splitters;        use DoblDobl_Solution_Splitters;
with QuadDobl_Solution_Splitters;        use QuadDobl_Solution_Splitters;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Scaling;                   use Standard_Scaling;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_BlackBox_Continuations;    use DoblDobl_BlackBox_Continuations;
with QuadDobl_BlackBox_Continuations;    use QuadDobl_BlackBox_Continuations;
with Black_Box_Solvers;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Standard_Irreducible_Decomp;
with Standard_Irreducible_Decomp_io;     use Standard_Irreducible_Decomp_io;
with Multprec_Irreducible_Decomp;
with Multprec_Irreducible_Decomp_io;     use Multprec_Irreducible_Decomp_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Homotopy_Cascade_Filter;            use Homotopy_Cascade_Filter;
with Drivers_to_Breakup_Components;      use Drivers_to_Breakup_Components;

package body Drivers_to_Cascade_Filtering is

  procedure Standard_Square_and_Embed is

    use Standard_Complex_Poly_Systems;

    lp,ep : Link_to_Poly_Sys;
    file : file_type;
    k : natural32;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Interactive_Square_and_Embed(file,lp.all,ep,k);
    new_line;
    put_line("See the output file for results...");
    new_line;
  end Standard_Square_and_Embed;

  procedure DoblDobl_Square_and_Embed is

    use DoblDobl_Complex_Poly_Systems;

    lp,ep : Link_to_Poly_Sys;
    file : file_type;
    k : natural32;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Interactive_Square_and_Embed(file,lp.all,ep,k);
    new_line;
    put_line("See the output file for results...");
    new_line;
  end DoblDobl_Square_and_Embed;

  procedure QuadDobl_Square_and_Embed is

    use QuadDobl_Complex_Poly_Systems;

    lp,ep : Link_to_Poly_Sys;
    file : file_type;
    k : natural32;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Interactive_Square_and_Embed(file,lp.all,ep,k);
    new_line;
    put_line("See the output file for results...");
    new_line;
  end QuadDobl_Square_and_Embed;

  procedure Driver_to_Square_and_Embed is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Square_and_Embed;
      when '1' => DoblDobl_Square_and_Embed;
      when '2' => QuadDobl_Square_and_Embed;
      when others => null;
    end case;
  end Driver_to_Square_and_Embed;

  procedure Standard_Remove_Embedding is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    file : file_type;
    k,ns : natural32 := 0;

  begin
    new_line;
    put_line("Removing an Embedding of a Polynomial System.");
    Standard_Read_Embedding(lp,sols,k,ns);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system...");
    new_line;
    declare
      rp : constant Poly_Sys := Remove_Embedding(lp.all,k,ns);
      nq : constant natural32 := natural32(rp'last);
      nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
      rsols : Solution_List;
    begin
      if k + ns > 0
       then rsols := Remove_Embedding(sols,k+ns);
       else rsols := sols;
      end if;
      put(file,nq,1);
      put(file," ");
      put(file,nv,1);
      new_line(file);
      put(file,rp);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
    end;
    Close(file);
  end Standard_Remove_Embedding;

  procedure DoblDobl_Remove_Embedding is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    file : file_type;
    k,ns : natural32 := 0;

  begin
    new_line;
    put_line("Removing an Embedding of a Polynomial System.");
    DoblDobl_Read_Embedding(lp,sols,k,ns);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system...");
    new_line;
    declare
      rp : constant Poly_Sys := Remove_Embedding(lp.all,k,ns);
      nq : constant natural32 := natural32(rp'last);
      nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
      rsols : Solution_List;
    begin
      if k + ns > 0
       then rsols := Remove_Embedding(sols,k+ns);
       else rsols := sols;
      end if;
      put(file,nq,1);
      put(file," ");
      put(file,nv,1);
      new_line(file);
      put(file,rp);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
    end;
    Close(file);
  end DoblDobl_Remove_Embedding;

  procedure QuadDobl_Remove_Embedding is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    file : file_type;
    k,ns : natural32 := 0;

  begin
    new_line;
    put_line("Removing an Embedding of a Polynomial System.");
    QuadDobl_Read_Embedding(lp,sols,k,ns);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system...");
    new_line;
    declare
      rp : constant Poly_Sys := Remove_Embedding(lp.all,k,ns);
      nq : constant natural32 := natural32(rp'last);
      nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
      rsols : Solution_List;
    begin
      if k + ns > 0
       then rsols := Remove_Embedding(sols,k+ns);
       else rsols := sols;
      end if;
      put(file,nq,1);
      put(file," ");
      put(file,nv,1);
      new_line(file);
      put(file,rp);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
    end;
    Close(file);
  end QuadDobl_Remove_Embedding;

  procedure Driver_to_Remove_Embedding is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Remove_Embedding;
      when others => null;
    end case;
  end Driver_to_Remove_Embedding;

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Down_Continuation
              ( file : in file_type;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    Black_Box_Polynomial_Continuation(file,target,embsys,sols,pocotime);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    Black_Box_Polynomial_Continuation(file,target,embsys,sols,pocotime);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    Black_Box_Polynomial_Continuation(file,target,embsys,sols,pocotime);
  end Down_Continuation;

  procedure Witness_Generate
              ( outfile,resfile : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32;
                zerotol : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32;
                zerotol : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32;
                zerotol : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  function Append_ck ( name : string; k : natural32 ) return string is

  -- DESCRIPTION :
  --   Appends "_swk" to the name.

    nbk : constant string := Convert(integer32(k));

  begin
     return name & "_sw" & nbk;
  end Append_ck;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32; zerotol : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    Write_Witness_Superset(name,ep,sols0,k);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32; zerotol : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    Write_Witness_Superset(name,ep,sols0,k);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32; zerotol : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-k;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(k));

  begin
    tstart(timer);
    embsys(integer32(k)) := new Poly_Sys'(ep);
    for i in 0..k-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,k-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(k),zerotol,sols0,sols1);
    Write_Witness_Superset(name,ep,sols0,k);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(k) loop
        Down_Continuation(outfile,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
  end Witness_Generate;

  procedure Standard_Witness_Generate is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;

  begin
    new_line;
    put_line("Reading the name of the file with the embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      Standard_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      Continuation_Parameters.Tune(0);
      Driver_for_Continuation_Parameters(outfile);
      new_line;
      put_line("See the input and output file for results...");
      new_line;
      close(infile);
     -- open(infile,out_file,name);
     -- Witness_Generate(outfile,infile,lp.all,sols,dim,1.0E-8);
      Witness_Generate(name,outfile,lp.all,sols,dim,1.0E-8);
    end;
  end Standard_Witness_Generate;

  procedure DoblDobl_Witness_Generate is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;

  begin
    new_line;
    put_line("Reading the name of the file with the embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      DoblDobl_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      Continuation_Parameters.Tune(0);
      Driver_for_Continuation_Parameters(outfile);
      new_line;
      put_line("See the input and output file for results...");
      new_line;
      close(infile);
     -- open(infile,out_file,name);
     -- Witness_Generate(outfile,infile,lp.all,sols,dim,1.0E-8);
      Witness_Generate(name,outfile,lp.all,sols,dim,1.0E-8);
    end;
  end DoblDobl_Witness_Generate;

  procedure QuadDobl_Witness_Generate is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;

  begin
    new_line;
    put_line("Reading the name of the file with the embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      QuadDobl_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      Continuation_Parameters.Tune(0);
      Driver_for_Continuation_Parameters(outfile);
      new_line;
      put_line("See the input and output file for results...");
      new_line;
      close(infile);
     -- open(infile,out_file,name);
     -- Witness_Generate(outfile,infile,lp.all,sols,dim,1.0E-8);
      Witness_Generate(name,outfile,lp.all,sols,dim,1.0E-8);
    end;
  end QuadDobl_Witness_Generate;

  procedure Driver_to_Witness_Generate is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Witness_Generate;
      when '1' => DoblDobl_Witness_Generate;
      when '2' => QuadDobl_Witness_Generate;
      when others => null;
    end case;
  end Driver_to_Witness_Generate;

-- OLD CODE FOR WITNESS GENERATE + CLASSIFY:

  procedure Timing_Summary ( file : in file_type;
                             roco,hoco,poco,total : in duration ) is

  -- DESCRIPTION :
  --   Displays the timing summary for the black-box solver.

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |         TIMING INFORMATION SUMMARY for Black-Box Solver           |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Black_Box_Solver
              ( file : in file_type;
                sys : in Standard_Complex_Poly_Systems.Poly_Sys;
                deg : in boolean;
                sols : out Standard_Complex_Solutions.Solution_List;
                rc : out natural32; totaltime : out duration ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    q : Poly_Sys(sys'range);
    roco,hoco,poco,total : duration;
    sols0 : Solution_List;

  begin
    tstart(timer);
    declare
      pp : Poly_Sys(sys'range);
    begin
      Copy(sys,pp);
      Black_Box_Root_Counting(file,0,pp,deg,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        Scale(pp);
        Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
      Clear(pp);
    end;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    new_line(file);
    print_times(file,timer,"Solving the polynomial system");
    new_line(file);
    Timing_Summary(file,roco,hoco,poco,total);
    totaltime := total;
  end Black_Box_Solver;

  procedure Write ( file : in file_type; arvv : in Array_of_VecVecs ) is

  -- DESCRIPTION :
  --   Writes the array of vectors of vectors in a formatted way on file.

  begin
    for i in arvv'range loop
      put(file,"Component "); put(file,i,1); put_line(file," of the array :");
      for j in arvv(i)'range loop
        put(file,"The "); put(file,j,1); put_line(file,"-th vector :");
        for k in arvv(i)(j)'range loop
          put(file,k,3); put(file," : "); put(file,arvv(i)(j)(k));
          new_line(file);
        end loop;
      end loop;
    end loop;
  end Write;

  procedure Write_Generate_Summary
              ( file : in file_type;
                npaths,nsols0,nsols1,ndiv : in Standard_Natural_Vectors.Vector;
                timings : in Array_of_Duration ) is

  -- DESCRIPTION :
  --   Writes the summary for the cascade of homotopies.

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |       TIMING INFORMATION SUMMARY for Cascade of Homotopies        |";
    b2 : constant string :=
     "  | level | #paths | #isols | #comps | #infty |     user cpu time     |";

    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    for i in reverse npaths'range loop
      put(file,"  | ");
      put(file,i,4); put(file,"  | ");
      put(file,npaths(i),5); put(file,"  | ");
      put(file,nsols1(i),5); put(file,"  | ");
      put(file,nsols0(i),5); put(file,"  | ");
      put(file,ndiv(i),5); put(file,"  |    ");
      print_hms(file,timings(i)); put_line(file,"     |");
    end loop;
    put_line(file,b0);
    put(file,"  | total | ");
    put(file,Standard_Natural_Vectors.Sum(npaths),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(nsols1),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(nsols0),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(ndiv),5); put(file,"  |    ");
    print_hms(file,total); put_line(file,"     |");
    put_line(file,b0);
  end Write_Generate_Summary;

  procedure Write ( file : in file_type; fp : in List; i : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the i-th component of every element in the list.

    tmp : List := fp;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      put(file,lpt(i),4); put(file," |");
      tmp := Tail_Of(tmp);
    end loop;
  end Write;

  procedure Write_Banner ( file : in file_type; 
                           m : in natural32; sep : character ) is
  begin
    for i in 1..m loop
      for j in 1..6 loop
        put(file,sep);
      end loop;
    end loop;
  end Write_Banner;

  procedure Write_Classify_Summary
               ( file : in file_type; fp : in List;
                 timings : in Array_of_Duration ) is

    n : integer32;
    m : constant natural32 := Length_Of(fp);
    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,"  -------------------------------------------------------");
    put_line(file,"  |  TIMING INFORMATION SUMMARY for Classifying Points  |");
    put_line(file,"  -------------------------------------------------------");
    if m /= 0 then
      n := Head_Of(fp)'last;
      put(file,"  | dimension "); put(file," |");
      Write(file,fp,n); 
      put_line(file,"  user cpu time  |");
      put(file,"  =============="); Write_Banner(file,m,'=');
      put_line(file,"==================");
      if n > 0 then
        put(file,"  | level ");
        put(file,n-1,1);  put(file,"    |");
        Write(file,fp,n-1);
        put(file," "); print_hms(file,timings(n-1));
        put_line(file,"  |");
        for i in reverse 0..n-2 loop
          put(file,"  |       ");
          put(file,i,1);  put(file,"    |");
          Write(file,fp,i);
          put(file," "); print_hms(file,timings(i));
          put_line(file,"  |");
        end loop;
      end if;
      put(file,"  --------------"); Write_Banner(file,m,'-');
      put_line(file,"------------------");
    end if;
    put(file,"  | total time : ");
    Write_Banner(file,m,' ');
    print_hms(file,total); put_line(file,"  |");
    put(file,"  --------------"); Write_Banner(file,m,'-');
    put_line(file,"------------------");
  end Write_Classify_Summary;

  procedure Driver_for_Cascade_Filter
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    nbequ : constant natural32 := natural32(p'length);
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));

  -- ROUTINES :

    function Max_Dim ( p : Poly_Sys ) return natural32 is

    -- DESCRIPTION :
    --   Returns either number of equations or the number of unknowns,
    --   depending on what is maximal.

    begin
      if nbequ >= nbunk
       then return nbequ;
       else return nbunk;
      end if;
    end Max_Dim;

    function Eliminate_Slices
                ( sqp : Standard_Complex_Poly_Systems.Poly_Sys;
                  m : natural32 )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    -- DESCRIPTION :
    --   Given is a polynomial system sqp that is the result of making
    --   p square with nbunk - nbequ additional slices.  The system on
    --   return will have min(nbunk-nbequ-2,m) slices and variables less.
    --   The last two slices remain.

      use Standard_Complex_Polynomials;

      res : Poly_Sys(sqp'range);
      cnt : natural32 := nbunk - nbequ;
      reslast : integer32 := sqp'last;

    begin
      Copy(sqp,res);
      if cnt > 2 then
        for i in 1..m loop
          declare
            elim : Poly_Sys(res'first..reslast-1)
                 := Eliminate_Slice
                      (res,natural32(reslast-2),natural32(reslast));
          begin
            for i in elim'range loop
              Copy(elim(i),res(i));
            end loop;
            Clear(res(reslast));
            Clear(elim);
          end;
          cnt := cnt - 1;
          reslast := reslast - 1;
          exit when cnt <= 2;
        end loop;
      end if;
      return res(res'first..reslast);
    end Eliminate_Slices;

    procedure Get_Multprec_System 
              ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                mpsys : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                size : in natural32 ) is

    -- DESCRIPTION :
    --   Depending on the answer of the user, the mpsys on return is read in
    --   from file or is the conversion of the given stsys.

      res : Multprec_Complex_Poly_Systems.Poly_Sys(stsys'range);
      ans : character;
      infile : file_type;
      m : natural32 := 0;
  
    begin
      put("Do you wish to read in the multi-precision coefficients? (y/n) "); 
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("Reading the polynomial system"
                  & " with multi-precision coefficients.");
        Read_Name_and_Open_File(infile);
        get(infile,m,res);
        Set_Size(res,size);
      else
        res := Convert(stsys);
        Set_Size(res,size);
      end if;
      put_line(file,"The system in multi-precision format : ");
      put_line(file,res);
      mpsys := new Multprec_Complex_Poly_Systems.Poly_Sys'(res);
    end Get_Multprec_System;

    procedure Cascade_Filter
                ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mpsys : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                  embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                  n,size,itp : natural32; skewproj,deg : in boolean ) is

    -- DESCRIPTION :
    --   Calls the routines to perform the filtering in the homotopy cascade.
  
    -- ON ENTRY :
    --   stsys    target system with standard complex coefficients;
    --   mpsys    system with multi-precision coefficients, null if size = 0;
    --   embp     sequence of embedded polynomial systems;
    --   n        original dimension;
    --   size     size of the numbers;
    --   itp      interpolator type;
    --   deg      if true, only degree based root counting, otherwise full;

      use Standard_Complex_Solutions;

      k : constant integer32 := embp'last;
      sols,sols0,sols1 : Solution_List;
      stansoco : Standard_Irreducible_Decomp.Solution_Components(0..k);
      multsoco : Multprec_Irreducible_Decomp.Solution_Components(0..k);
      tol : constant double_float := 10.0E-10;
      tol_sing : constant double_float := 10.0E-8;
      tol_eval : constant double_float := 10.0E-5;
      npa,ns0,ns1,div : Standard_Natural_Vectors.Vector(0..k) := (0..k => 0);
        -- npa = #paths,
        -- ns0 = #solutions with zz == 0
        -- ns0 = #solutions with zz /= 0
        -- dvi = #solutions diverged to infinity
      fp,fp_last : List;
      gentimes,clatimes : Array_of_Duration(0..k) := (0..k => 0.0);
      firstzero : boolean := false;  -- flag for first time zero added var

    begin
      if size = 0 
       then Standard_Initialize(stansoco,integer32(n),stsys);
       else Multprec_Initialize(multsoco,integer32(n),mpsys.all);
      end if;
      put(file,"EMBEDDED SYSTEM at the top dimension k = ");
      put(file,k,1); put_line(file," :");
      put_line(file,embp(k).all);
      Black_Box_Solver(file,embp(k).all,deg,sols,npa(k),gentimes(k));
      Filter_and_Split_Solutions(file,sols,integer32(n),k,tol,sols0,sols1);
      ns0(k) := Length_Of(sols0);
      ns1(k) := Length_Of(sols1);
      div(k) := npa(k) - ns0(k) - ns1(k);
      firstzero := (ns0(k) /= 0);
      if firstzero then
        if size = 0 then
          Standard_Update_Hypersurfaces
            (file,stansoco,n,natural32(k),natural32(k),itp,skewproj,
             embp(k).all,sols0,tol_sing,clatimes(k),fp,fp_last);
        else
          Multprec_Update_Hypersurfaces
            (file,multsoco,n,natural32(k),natural32(k),size,itp,skewproj,
             embp(k).all,mpsys.all,sols0,tol_sing,clatimes(k),fp,fp_last);
        end if;
      end if;
      new_line(file);
      if not Is_Null(sols1) then
        if k = 0 then
          if size = 0
           then Copy(sols1,stansoco(0).pts);
           else Copy(sols1,multsoco(0).pts);
          end if;
        else
          Copy(sols1,sols);
          if size = 0 then
            Standard_Cascade_Loop
              (file,n,natural32(k),itp,skewproj,embp,sols,stansoco,tol,tol_sing,
               tol_eval,npa,ns0,ns1,div,fp,fp_last,gentimes,clatimes);
          else
            Multprec_Cascade_Loop
              (file,n,natural32(k),size,itp,skewproj,embp,mpsys.all,sols,
               multsoco,tol,tol_sing,tol_eval,npa,ns0,ns1,div,
               fp,fp_last,gentimes,clatimes);
          end if;
        end if;
      end if;
      new_line(file);
      put_line(file,"THE SOLUTION COMPONENTS : ");
      if size = 0
       then put(file,stansoco);
       else put(file,multsoco);
      end if;
      Write_Generate_Summary(file,npa,ns0,ns1,div,gentimes);
      Write_Classify_Summary(file,fp,clatimes);
    end Cascade_Filter;

    procedure Main_Driver is

      n : natural32 := natural32(p'last);
      sqp : Poly_Sys(1..integer32(Max_Dim(p))) := Square(p);
      embp : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k);
      stsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
      mpsys : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
      ans : character;
      itp,size : natural32;
      deg,skewproj : boolean;
      timer : Timing_Widget;

    begin
      tstart(timer);
      if nbequ < nbunk then
        --embp := Slice_and_Embed(sqp,k,1);
        Copy(sqp(sqp'last-1),sqp(sqp'last));
        if nbequ < nbunk - 1 then
          declare
            esqp : constant Poly_Sys := Eliminate_Slices(sqp,natural32(k));
          begin
            stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(esqp);
          end;
        else
          stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(sqp);
        end if;
       else stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(sqp);
      end if;
      put_line(file,"The original system after squaring : ");
      put_line(file,stsys.all);
      embp := Slice_and_Embed(stsys.all,natural32(k));
      n := natural32(stsys'last);     -- WARNING : this could be bad....
      Breakup_Menu(itp,size,skewproj);
      if size > 0
       then Get_Multprec_System(stsys.all,mpsys,size);
      end if;
      put("Full root counting or only based on degrees? (f/d) ");
      Ask_Alternative(ans,"fd");
      deg := (ans = 'd');
      new_line;
      put_line("See the output file for results...");
      new_line;
      tstart(timer);
     -- put_line(file,"The original system (after squaring) : ");
     -- if nbequ > nbunk
     --  then Add_New_Symbols(k+nbequ-nbunk);
     --  else
      Add_Embed_Symbols(natural32(k));
     -- end if;
     -- put(file,stsys'last,stsys.all);
      Cascade_Filter(stsys.all,mpsys,embp,n,size,itp,skewproj,deg);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Numerical Irreducible Decomposition");
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Cascade_Filter;

  procedure Embed_and_Cascade is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp,ep : Link_to_Poly_Sys;
    file : file_type;
    timer : timing_widget;
    sols : Solution_List;
   -- sols0,sols1 : Solution_List;
    tol : constant double_float := 1.0E-8;
   -- n : natural;
    k,rc : natural32;
    outfilename : Link_to_String;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file,outfilename);
    new_line;
    Interactive_Square_and_Embed(file,lp.all,ep,k);
    new_line;
    put_line("Calling the blackbox solver on the top embedding.");
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    Black_Box_Solvers.Solve(file,ep.all,rc,sols);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
     -- n := Head_Of(sols).n - k;
     -- Filter_and_Split_Solutions(file,sols,n,k,tol,sols0,sols1);
      Witness_Generate(outfilename.all,file,ep.all,sols,k,tol);
    end if;
  end Embed_and_Cascade;

end Drivers_to_Cascade_Filtering;
