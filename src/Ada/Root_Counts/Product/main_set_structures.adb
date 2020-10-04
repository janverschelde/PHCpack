with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Set_Structure,Set_Structure_io;
with Degree_Sets_Tables;                 use Degree_Sets_Tables;
with Standard_Linear_Product_System;
with Random_Product_Start_Systems;       use Random_Product_Start_Systems;
with Standard_Complex_Prod_Planes;

package body Main_Set_Structures is

  procedure Set_Structure_Info is

    i : array(1..18) of string(1..65);

  begin
    i( 1):="  A generalized Bezout  number  is  based  on  a  supporting  set";
    i( 2):="structure.   A  set  structure is a tuple of arrays of subsets of";
    i( 3):="unknowns.                                                        ";
    i( 4):="  The corresponding start system is a linear-product system:  the";
    i( 5):="i-th  equation  is  the  product  of linear equations with random";
    i( 6):="coefficient in the unknowns of the set of the  i-th  array.   The";
    i( 7):="number  of  factors  in  the product for the i-th equation of the";
    i( 8):="start system equals the number of subsets in the  i-th  array  of";
    i( 9):="the set structure.                                               ";
    i(10):="  A set structure is supporting for a polynomial system if  every";
    i(11):="monomial  in  the system also occurs in the corresponding linear-";
    i(12):="product start system.                                            ";
    i(13):="  Given a supporting set structure, the generalized Bezout number";
    i(14):="equals  the  number  of  solutions  of  the corresponding linear-";
    i(15):="product start system.   Before  the  construction  of  the  start";
    i(16):="system, a generalized Bezout number is first computed in a formal";
    i(17):="way as a generalized permanent of a degree matrix.   A  heuristic";
    i(18):="procedure is available for generating a supporting set structure.";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Set_Structure_Info;

  procedure Read_Set_Structure ( n : in natural32 ) is

    ns : Standard_Natural_Vectors.Vector(1..integer32(n));

  begin
    for i in ns'range loop
      put("  Give the number of sets for polynomial ");
      put(i,1); put(" : "); Read_Natural(ns(i));
    end loop;
    Set_Structure.Init(ns);
    put_line("Give the set structure :");
    Set_Structure_io.get;
  end Read_Set_Structure;

  procedure Write_Results ( file : in file_type; bb : in natural32 ) is
  begin
    new_line(file);
    put(file,"  generalized Bezout number is "); put(file,bb,1);
    new_line(file);
    put_line(file,"  based on the set structure :");
    Set_Structure_io.put(file);
  end Write_Results;

  procedure Save_Results ( q : in Prod_Sys; qsols : in Solution_List ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading file name to write start system.");
    Read_Name_and_Create_File(file);
    put_line(file,natural32(q'last),q);
    put_line(file,
             "TITLE : a random coefficient linear-product start system");
    if not Is_Null(qsols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS : ");
      new_line(file);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      Close(file);
    end if;
  end Save_Results;

  procedure Menu_Prompt ( choice : out character; bb : in natural32 ) is

    ans : character;

  begin
    new_line;
    put_line("MENU for generalized Bezout Numbers based on Set Structures :");
    put     ("  0. exit - current Bezout number is "); put(bb,1); new_line;
    put_line("  1. Apply heuristic constructor for set structure");
    put_line("  2. Evaluate your own set structure");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(ans,"012"); choice := ans;
  end Menu_Prompt;

  procedure Menu_Handler ( file : in file_type; choice : in character;
                           p : in Poly_Sys; bb : in out natural32 ) is
  begin
    case choice is
      when '1' =>  
        Random_Product_Start_Systems.Build_Set_Structure(p);
        bb := natural32(Matching_Permanent(Degree_Sets_Tables.Create));
      when '2' => 
        Read_Set_Structure(natural32(p'last));
        bb := natural32(Matching_Permanent(Degree_Sets_Tables.Create));
      when others => null;
    end case;
    Write_Results(Standard_Output,bb); Write_Results(file,bb);
  end Menu_Handler;

  procedure Compute_Bezout_Number
                ( file : in file_type;
                  p : in Poly_Sys; bb : in out natural32 ) is

    method : character;
    timer : timing_widget;

  begin
    new_line(file);
    put_line(file,"SET STRUCTURE ANALYSIS :");
    tstart(timer);
    loop
      Menu_Prompt(method,bb);
      exit when method = '0';
      Menu_Handler(file,method,p,bb);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"set structure analysis");
  end Compute_Bezout_Number;

  procedure Construct_Start_System
              ( file : in file_type; p : in Poly_Sys;
                bb : in natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

    n : constant natural32 := natural32(p'last);
    ans : character;
    timer : timing_widget;
    qq : Poly_Sys(p'range);
    rq : Prod_Sys(p'range);
    qqsols : Solution_List;
    nl : natural32;

  begin
    new_line;
    put("Do you want a start system based on the set structure ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("There are "); put(bb,1); put_line(" start solutions.");
      put_line("phc -q can compute start solutions later whenever needed.");
      put("Do you want to compute all solutions now ? (y/n) ");
      Ask_Yes_or_No(ans);
      tstart(timer);
      Standard_Linear_Product_System.Init(n);
      Build_Random_Product_System(n);
      rq := Standard_Complex_Prod_Planes.Create;
      qq := Standard_Linear_Product_System.Polynomial_System;
     -- Random_Product_System.Solve(qqsols,nl,lpos);
      if ans = 'y'
       then Standard_Linear_Product_System.Solve(qqsols,nl);
      end if;
      tstop(timer);
      Save_Results(rq,qqsols);
      q := qq; qsols := qqsols;
      new_line(file);
      put_line(file,"RANDOM LINEAR-PRODUCT START SYSTEM : ");
      put_line(file,natural32(rq'last),rq); -- put_line(file,qq);
      if not Is_Null(qqsols) then
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        new_line(file);
        put(file,Length_Of(qqsols),natural32(Head_Of(qqsols).n),qqsols);
        new_line(file);
        print_times(file,timer,"constructing and solving the start system");
      else
        new_line(file);
        print_times(file,timer,"constructing a linear-product start system");
      end if;
      Set_Structure.Clear;
      Standard_Linear_Product_System.Clear;
    else
      Set_Structure.Clear;
     -- Clear(lpos);
    end if;
  end Construct_Start_System;

  procedure Main ( file : in file_type; p : in Poly_Sys; 
                   b : in out natural32; -- lpos : in out List;
                   q : out Poly_Sys; qsols : out Solution_List ) is

    bb : natural32 := b;

  begin
    Compute_Bezout_Number(file,p,bb);
    if not Set_Structure.Empty then
      b := bb;
      Construct_Start_System(file,p,bb,q,qsols);
    end if;
  end Main;

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;
   -- lpos : List;
    b : natural32 := 0;

  begin
    get(lp);
    declare
      q : Poly_Sys(lp'range);
      qsols : Solution_List;
    begin
      put_line("Reading the output file.");
      Read_Name_and_Create_File(file);
     -- Main_Set_Structures.Main(file,lp.all,b,lpos,q,qsols);
      Main_Set_Structures.Main(file,lp.all,b,q,qsols);
    end;
  end Main;

end Main_Set_Structures;
