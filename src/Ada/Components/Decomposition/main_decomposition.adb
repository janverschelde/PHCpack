with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Prompt_for_Systems;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;         use Standard_Monomial_Maps_io;
with Standard_Monomial_Map_Filters;
with Black_Box_Binomial_Solvers;        use Black_Box_Binomial_Solvers;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;
with Drivers_to_Witness_Generate;       use Drivers_to_Witness_Generate;
with Add_and_Remove_Embedding;
with Extrinsic_Diagonal_Continuation;   use Extrinsic_Diagonal_Continuation;
with Drivers_to_Cascade_Filtering;      use Drivers_to_Cascade_Filtering;
with Extrinsic_Diagonal_Homotopies_io;
with Extrinsic_Diagonal_Solvers;        use Extrinsic_Diagonal_Solvers;
with Drivers_to_Intersect_Varieties;    use Drivers_to_Intersect_Varieties;
with Driver_to_Rank_Supports;

package body Main_Decomposition is

  procedure Read_Two_Witness_Sets
              ( lp1,lp2 : out Link_to_Poly_Sys;
                sols1,sols2 : out Solution_List;
                dim1,dim2 : out natural32; vrb : in integer32 := 0 ) is
  begin
    if vrb > 0
     then put_line("-> in main_decomposition.Read_Two_Witness_Sets 1 ...");
    end if;
    new_line;
    put_line("Reading the first embedded polynomial system.");
    Standard_Read_Embedding(lp1,sols1,dim1);
    new_line;
    put_line("Reading the second embedded polynomial system.");
    Standard_Read_Embedding(lp2,sols2,dim2);
  end Read_Two_Witness_Sets;

  procedure Read_Two_Witness_Sets
              ( lp1,lp2 : out Link_to_Poly_Sys;
                sols1,sols2 : out Solution_List;
                dim1,dim2 : out natural32;
                lsym1,lsym2 : out Symbol_Table.Link_to_Array_of_Symbols;
                vrb : in integer32 := 0 ) is
  begin
    if vrb > 0
     then put_line("-> in main_decomposition.Read_Two_Witness_Sets 2 ...");
    end if;
    new_line;
    put_line("Reading the first embedded polynomial system.");
    Standard_Read_Embedding(lp1,sols1,dim1);
    lsym1 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
    new_line;
    put_line("Reading the second embedded polynomial system.");
    Standard_Read_Embedding(lp2,sols2,dim2);
    lsym2 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
  end Read_Two_Witness_Sets;

  procedure Call_Extrinsic_Diagonal_Homotopies
              ( outfilename : in string; vrb : in integer32 := 0 ) is

    lp1,lp2 : Link_to_Poly_Sys;
    esols1,esols2,sols1,sols2 : Solution_List;
    dim1,dim2,oc : natural32;
    lsym1,lsym2 : Symbol_Table.Link_to_Array_of_Symbols;
    file : file_type;
    name : Link_to_String;
    report : boolean;

  begin
    if vrb > 0 then
      put("-> in main_decomposition.");
      put_line("Call_Extrinsic_Diagonal_Homotopies ...");
    end if;
    new_line;
    put_line("Executing diagonal homotopies extrinsically...");
    Read_Two_Witness_Sets(lp1,lp2,esols1,esols2,dim1,dim2,lsym1,lsym2,vrb-1);
    new_line;
    put("Symbols of first system : ");
    Extrinsic_Diagonal_Homotopies_io.Write(lsym1.all);
    put("Symbols of second system : ");
    Extrinsic_Diagonal_Homotopies_io.Write(lsym2.all);
    sols1 := Remove_Embedding(esols1,dim1);
    sols2 := Remove_Embedding(esols2,dim2);
    Create_Output_File(file,outfilename,name);
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("Running continuation, see the output file for results...");
    new_line;
    if dim1 >= dim2 then
      Extrinsic_Diagonal_Homotopy
        (file,name.all,report,lp1.all,lp2.all,dim1,dim2,sols1,sols2);
    else
      Extrinsic_Diagonal_Homotopy
        (file,name.all,report,lp2.all,lp1.all,dim2,dim1,sols2,sols1);
    end if;
  end Call_Extrinsic_Diagonal_Homotopies;

  procedure Call_Intrinsic_Diagonal_Homotopies
              ( outfilename : in string; vrb : in integer32 := 0 ) is

    lp1,lp2 : Link_to_Poly_Sys;
    sols1,sols2 : Solution_List;
    dim1,dim2,oc : natural32;
    file : file_type;
    name : Link_to_String;
    ans : character;
    generic_version : boolean;
    report : boolean;
    n : natural32;
    f : Link_to_Poly_Sys;
    p : Link_to_Matrix;
    s : Solution_List;

  begin
    if vrb > 0 then
      put("-> in main_decomposition.");
      put_line("-> Call_Intrinsic_Diagonal_Homotopies ...");
    end if;
    new_line;
    put_line("Executing diagonal homotopies intrinsically...");
    Read_Two_Witness_Sets(lp1,lp2,sols1,sols2,dim1,dim2,vrb-1);
    Create_Output_File(file,outfilename,name);
    new_line;
    put("Do you want the version using generics ? (y/n) ");
    Ask_Yes_or_No(ans);
    generic_version := (ans = 'y');
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("Running continuation, see the output file for results...");
    new_line;
    if generic_version then
      if dim1 >= dim2
       then A_Call_Generic_Diagonal_Homotopy
              (file,lp1.all,lp2.all,sols1,sols2,dim1,dim2);
       else A_Call_Generic_Diagonal_Homotopy
              (file,lp2.all,lp1.all,sols2,sols1,dim2,dim1);
      end if;
    else
      if dim1 >= dim2
       then Intrinsic_Diagonal_Homotopy
              (file,report,lp1.all,lp2.all,sols1,sols2,dim1,dim2,f,p,s);
       else Intrinsic_Diagonal_Homotopy
              (file,report,lp2.all,lp1.all,sols2,sols1,dim2,dim1,f,p,s);
      end if;
      n := natural32(lp1'last) - dim1;
      new_line(file);
      Write_Witness_Set(file,name.all,n,n-natural32(p'last(2)),f.all,s,p.all);
    end if;
  end Call_Intrinsic_Diagonal_Homotopies;

  procedure Call_Binomial_Solver
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    infile,file : file_type;
    lp : Link_to_Laur_Sys;
    sysonfile : boolean;
    ans : character;
    sols : Link_to_Array_of_Monomial_Map_Lists;
    fail : boolean;
    nv : natural32;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in maindeco.Call_Binomial_Solver ...");
    end if;
    Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(file,outfilename);
    nv := Number_of_Unknowns(lp(lp'first));
    put(file,lp'last,1); put(file," ");
    put(file,nv,1); new_line(file);
    put(file,lp.all);
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Black_Box_Binomial_Solver(standard_output,lp.all,sols,fail);
    else
      put("Is the solution set pure dimensional ? (y/n) ");
      Ask_Yes_or_No(ans);
      new_line;
      put_line("See the output file for results ...");
      new_line;
      tstart(timer);
      if ans = 'y'
       then Black_Box_Binomial_Solver(lp.all,true,sols,fail);
       else Black_Box_Binomial_Solver(lp.all,false,sols,fail);
      end if;
      tstop(timer);
    end if;
    if not fail and sols /= null then
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,sols.all);
    end if;
    new_line(file);
    print_times(file,timer,"the blackbox binomial system solver");
  end Call_Binomial_Solver;

  procedure Degrees_of_Monomial_Maps
              ( infilename : in string; vrb : in integer32 := 0 ) is

    file : file_type;
    lp : Link_to_Laur_Sys;
    maps : Monomial_Map_List;

  begin
    if vrb > 0 then
      put_line("-> in main_decomposition.");
      put_line("Degrees_of_Monomial_Maps ...");
    end if;
    if infilename = "" then
      new_line;
      Read_System_and_Maps(lp,maps);
    else
      Open_Input_File(file,infilename);
      Read_System_and_Maps(file,lp,maps);
    end if;
    new_line;
    declare
      td : constant natural32 := Top_Dimension(maps);
      c : constant Array_of_Monomial_Map_Lists(0..integer32(td))
        := Standard_Monomial_Map_Filters.Pure_Dimensional_Maps(maps);
    begin
      Show_Degrees(c);
    end;
  end Degrees_of_Monomial_Maps;

  procedure Transform_Positive_Corank 
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    infile,file : file_type;
    lp : Link_to_Laur_Sys;
    sysonfile : boolean;

  begin
    if vrb > 0
     then put_line("-> in main_decomposition.Transform_Positive_Corank ...");
    end if;
    Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(file,outfilename);
    Driver_to_Rank_Supports(file,lp.all);
  end Transform_Positive_Corank;

  procedure Run_Cascade_Filter
               ( infilename,outfilename : in string;
                 vrb : in integer32 := 0 ) is

    lp : Link_to_Poly_Sys;
    k : natural32 := 0;
    infile,file : file_type;
    sysonfile : boolean;

  begin
    if vrb > 0
     then put_line("-> in maindeco.Run_Cascade_Filter ...");
    end if;
    Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(file,outfilename);
    new_line;
    put("Give the top dimension : "); get(k); skip_line;
    Driver_for_Cascade_Filter(file,lp.all,integer32(k));
  end Run_Cascade_Filter;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    ans : character;

    use Add_and_Remove_Embedding; -- for 1 and 4

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in maindeco.Main ...");
    end if;
    new_line;
    put_line("MENU to generate and classify witness points :");
    put_line("  0. embed the system and run the cascade of homotopies;");
    put_line("  1. generate an embedding for the top dimensional component;");
    put_line("  2. given the top embedding, compute all witness points;");
    put_line("  3. perform steps 1, 2, and classify the witness points;");
    put_line("  4. remove the embedding from a polynomial system.");
    put_line("MENU to intersect positive dimensional solution components :");
    put_line("  5. build extrinsic diagonal cascade from 2 embedded systems;");
    put_line("  6. eliminate the diagonal, go to original at end of cascade;");
    put_line("  7. intersect varieties with extrinsic diagonal homotopies;");
    put_line("  8.                          intrinsic diagonal homotopies.");
    put_line("MENU for an irreducible decomposition of binomial systems :");
    put_line("  9. apply permanent factor & filter solver;");
    put_line("  A. compute degrees of a given list of monomial maps;");
    put_line("  B. transform a Laurent system with positive corank supports.");
    put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, or B to choose : ");
    Ask_Alternative(ans,"0123456789AB");
    case ans is
      when '0' => Embed_and_Cascade(nt,infilename,outfilename,verbose-1);
      when '1' => Driver_to_Square_and_Embed(infilename,outfilename);
      when '2' => Driver_to_Witness_Generate(nt,infilename,outfilename);
      when '3' => Run_Cascade_Filter(infilename,outfilename,verbose-1);
      when '4' => Driver_to_Remove_Embedding(infilename,outfilename);
      when '5' => Build_Diagonal_Cascade;
      when '6' => Collapse_Diagonal_System;
      when '7' => Call_Extrinsic_Diagonal_Homotopies(outfilename,verbose-1);
      when '8' => Call_Intrinsic_Diagonal_Homotopies(outfilename,verbose-1);
      when '9' => Call_Binomial_Solver(infilename,outfilename,verbose-1);
      when 'A' => Degrees_of_Monomial_Maps(infilename,verbose-1);
      when 'B' => Transform_Positive_Corank(infilename,outfilename,verbose-1);
      when others => null;
    end case;
  end Main;

end Main_Decomposition;
