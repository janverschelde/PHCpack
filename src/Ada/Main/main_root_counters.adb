with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Driver_for_Own_Start_System;
with Main_m_Homogenization;
with Main_Multi_Homogenization;
with Main_Set_Structures;
with Driver_for_Symmetric_Set_Structure; use Driver_for_Symmetric_Set_Structure;
with Drivers_for_Implicit_Lifting;       use Drivers_for_Implicit_Lifting;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Drivers_for_Dynamic_Lifting;        use Drivers_for_Dynamic_Lifting;
with Drivers_for_Symmetric_Lifting;      use Drivers_for_Symmetric_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;
with Drivers_for_DEMiCs_Algorithm;       use Drivers_for_DEMiCs_Algorithm;
with Write_Seed_Number;
with Greeting_Banners;
with Bye_Bye_Message;

package body Main_Root_Counters is

  procedure High_Total_Degree ( file : in file_type; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Displays a warning about a total degree, higher than what can be
  --   represented exactly by integers supported by hardware.

    totdeg : Natural_Number;

    use Standard_Complex_Polynomials;

  begin
    put("The total degree : ");
    for i in p'range loop
      put(Degree(p(i)),1);    put(file,Degree(p(i)),1);
      exit when i = p'last;
      put("*");               put(file,"*");
    end loop;
    totdeg := Total_Degree(p);
    put(file," = "); put(file,totdeg); new_line(file);
    put(" = "); put(totdeg); new_line;
    new_line;
    put_line("  this is larger than the largest machine integer...");
    Clear(totdeg);
  end High_Total_Degree;

  procedure Display_General_Menu
               ( rc : in natural64; own : in boolean;
                 method : in character ) is

  -- DESCRIPTION :
  --   Shows the menu to count roots and construct start systems
  --   for general polynomial systems.

    m : array(0..11) of string(1..66);

  begin
    new_line;
    put_line("MENU with ROOT COUNTS and Methods to Construct START SYSTEMS :");
    put("  0. exit - current root count is ");
    case method is
      when '0' => put("based on total degree : ");
      when '1' => put("based on multi-homogenization : ");
      when '2' => put("based on partitioned linear-product : ");
      when '3' => put("based on set structure : ");
      when '4' => put("based on symmetric set structure : ");
      when '5' => put("based on Bezout and BKK Bound : ");
      when '6' => put("based on static mixed-volume computation : ");
      when '7' => put("based on dynamic mixed-volume computation : ");
      when '8' => put("based on symmetric mixed-volume computation : ");
      when '9' => put("based on the mixed volume (via MixedVol) : ");
      when 'A' => put("based on the mixed volume (via DEMiCs) : ");
      when others => put("based on your start system : ");
    end case;
    put(rc,1); new_line;
    m(0):="PRODUCT HOMOTOPIES based on DEGREES ------------------------------";
    m(1):="  1. multi-homogeneous Bezout number               (one partition)";
    m(2):="  2. partitioned linear-product Bezout number    (many partitions)";
    m(3):="  3. general linear-product Bezout number          (set structure)";
    m(4):="  4. symmetric general linear-product Bezout number (group action)";
    m(5):="POLYHEDRAL HOMOTOPIES based on NEWTON POLYTOPES ------------------";
    m(6):="  5. combination between Bezout and BKK Bound   (implicit lifting)";
    m(7):="  6. mixed-volume computation                     (static lifting)";
    m(8):="  7. incremental mixed-volume computation        (dynamic lifting)";
    m(9):="  8. symmetric mixed-volume computation        (symmetric lifting)";
   m(10):="  9. using MixedVol Algorithm to compute the mixed volume fast (!)";
   m(11):="  A. compute the mixed volume with DEMiCs to see the cells quickly";
    for i in m'range loop
      put_line(m(i));
    end loop;
    if own then
      put_line
        ("START SYSTEM DEFINED BY USER -------------------------------------");
      put_line("  B. you can give your own start system");
    else
      put_line
        ("------------------------------------------------------------------");
    end if;
  end Display_General_Menu;

  procedure Display_Laurent_Menu
              ( rc : in natural32; method : in character ) is

  -- DESCRIPTION :
  --   Shows the menu to count roots and construct start systems
  --   for Laurent polynomial systems.

    m : array(0..5) of string(1..66);

  begin
    new_line;
    put_line("MENU for MIXED VOLUMES and POLYHEDRAL HOMOTOPIES :");
    put("  0. exit - current root count is ");
    put_line
        ("------------------------------------------------------------------");
    case method is
      when '1' => put("based on the mixed volume (via DEMiCs) : ");
      when '2' => put("based on the mixed volume (via MixedVol) : ");
      when '3' => put("based on Bezout and BKK Bound : ");
      when '4' => put("based on static mixed-volume computation : ");
      when '5' => put("based on dynamic mixed-volume computation : ");
      when '6' => put("based on symmetric mixed-volume computation : ");
      when others =>  null;
    end case;
    put(rc,1); new_line;
    m(0):="  1. compute the mixed volume with DEMiCs to see the cells quickly";
    m(1):="  2. using MixedVol Algorithm to compute the mixed volume fast (!)";
    m(2):="  3. combination between Bezout and BKK Bound   (implicit lifting)";
    m(3):="  4. mixed-volume computation                     (static lifting)";
    m(4):="  5. incremental mixed-volume computation        (dynamic lifting)";
    m(5):="  6. symmetric mixed-volume computation        (symmetric lifting)";
    for i in m'range loop
      put_line(m(i));
    end loop;
    put_line
        ("------------------------------------------------------------------");
  end Display_Laurent_Menu;

  procedure Display_Info ( method : in character ) is

  -- DESCRIPTION :
  --   Displays the information that corresponds with the current method.

  begin
    new_line;
    case method is
      when '0' => Total_Degree_Info;
      when '1' => Main_m_Homogenization.m_Homogenization_Info;
      when '2' => Main_Multi_Homogenization.Multi_Homogenization_Info;
      when '3' => Main_Set_Structures.Set_Structure_Info;
      when '4' => Symmetric_Set_Structure_Info;
      when '5' => Implicit_Lifting_Info;
      when '6' => Static_Lifting_Info;
      when '7' => Dynamic_Lifting_Info;
      when '8' => Symmetric_Lifting_Info;
      when '9' => MixedVol_Algorithm_Info;
      when 'A' => DEMiCs_Algorithm_Info;
      when others => put_line("No information available.");
    end case;
    new_line;
  end Display_Info;

  procedure Write_Start_System
               ( q : in Poly_Sys; qsols : in Solution_List ) is

  -- DESCRIPTION :
  --   This procedures asks the user if the start system needs
  --   to be written on a separate file.  This routine is needed
  --   only in case of a total degree start system.

    qfile : file_type;
    ans : character;

  begin
    put("Do you wish start system and solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading file name to write start system.");
      Read_Name_and_Create_File(qfile);
      put_line(qfile,q);
      new_line(qfile);
      put_line(qfile,"TITLE : start system based on total degree");
      new_line(qfile);
      if not Is_Null(qsols) then
        put_line(qfile,"THE SOLUTIONS : ");
        new_line(qfile);
        put(qfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      end if;
      Close(qfile);
    end if;
  end Write_Start_System;

  procedure Driver_for_Total_Degree_Start_System 
              ( p : in Poly_Sys; rc : in natural32;
                q : out Poly_Sys; qsols : out Solution_List ) is

  -- DESCRIPTION :
  --   Offers the user the choice whether or not to compute all
  --   start solutions, which may be a high number.

    cq : Standard_Complex_Vectors.Vector(q'range);
    dp : Standard_Natural_Vectors.Vector(p'range);
    ans : character;

  begin
    q := Total_Degree_Start_Systems.Start_System(p);
    put("There are "); put(rc,1); put_line(" start solutions.");
    put_line("phc -q can compute start solutions later whenever needed.");
    put("Do you want to compute all solutions now ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      dp := Total_Degree_Start_Systems.Degrees(p);
      cq := Coefficients(q);
      qsols := Solve(dp,cq);
    end if;
    Write_Start_System(q,qsols);
  end Driver_for_Total_Degree_Start_System;

  function To32Bit ( n : natural64 ) return natural32 is

  -- DESCRIPTION :
  --   Graceful conversion to 32-bit natural number.
  --   If n is out of the 32-bit range, 0 is returned.

  begin
    if n > natural64(natural'last)
     then return 0;
     else return natural32(n);
    end if;
  end To32Bit;

  procedure Apply_General_Method
              ( file : in file_type; nt : in integer32;
                p,q : in out Poly_Sys;
                qsols : in out Solution_List; method : in character;
                irc : in out natural64; lpos : out List ) is

  -- DESCRIPTION :
  --   Applies the root count that corresponds with the current method,
  --   for general polynomial systems.

    rc : natural32 := To32Bit(irc);
    rq : Prod_Sys(p'range);
    qsols0 : Solution_List;
    smv,tmv : natural32;

  begin
    case method is
      when '0' => Driver_for_Total_Degree_Start_System(p,rc,q,qsols);
      when '1' => Main_m_Homogenization.Main(file,p,irc,q,qsols);
      when '2' => Main_Multi_Homogenization.Main(file,p,rc,q,rq,qsols);
      when '3' => Main_Set_Structures.Main(file,p,rc,q,qsols);
      when '4' => Driver_for_Symmetric_Random_Product_Systems
                    (file,p,q,qsols,rc,lpos);
      when '5' => Driver_for_Mixture_Bezout_BKK(file,p,false,q,qsols,rc);
      when '6' => Driver_for_Mixed_Volume_Computation
                    (file,nt,p,false,q,qsols,qsols0,rc,smv,tmv);
      when '7' => Driver_for_Dynamic_Mixed_Volume_Computation
                    (file,p,false,q,qsols,rc);
      when '8' => Driver_for_Symmetric_Mixed_Volume_Computation
                    (file,p,false,q,qsols,rc);
      when '9' => Driver_for_MixedVol_Algorithm
                    (file,nt,p,false,false,q,qsols,qsols0,rc,smv,tmv);
      when 'A' => Driver_for_DEMiCs_Algorithm
                    (file,nt,p,q,qsols,qsols0,rc,smv,tmv);
      when 'B' => Driver_for_Own_Start_System(file,p,q,qsols);
                  rc := Length_of(qsols);
      when others => null;
    end case;
    irc := natural64(rc);
  end Apply_General_Method;

  procedure Apply_Laurent_Method
              ( file : in file_type; nt : in integer32;
                p,q : in out Laur_Sys;
                qsols : in out Solution_List; method : in character;
                rc : in out natural32 ) is

  -- DESCRIPTION :
  --   Applies the root count that corresponds with the current method,
  --   for Laurent polynomial systems.

    qsols0 : Solution_List;
    smv,tmv : natural32;

  begin
    case method is
      when '1' => Driver_for_DEMiCs_Algorithm
                    (file,nt,p,q,qsols,qsols0,rc,smv,tmv);
      when '2' => Driver_for_MixedVol_Algorithm
                    (file,nt,p,false,false,q,qsols,qsols0,rc,smv,tmv);
      when '3' => Driver_for_Mixture_Bezout_BKK(file,p,false,q,qsols,rc);
      when '4' => Driver_for_Mixed_Volume_Computation
                    (file,nt,p,false,q,qsols,qsols0,rc,smv,tmv);
      when '5' => Driver_for_Dynamic_Mixed_Volume_Computation
                    (file,p,false,q,qsols,rc);
      when '6' => Driver_for_Symmetric_Mixed_Volume_Computation
                    (file,p,false,q,qsols,rc);
      when others => null;
    end case;
  end Apply_Laurent_Method;

  procedure Polynomial_Main
               ( file : in file_type; nt : in natural32;
                 p,q : in out Poly_Sys; own : in boolean;
                 qsols : in out Solution_List; roco : out natural32;
                 verbose : in integer32 := 0 ) is

    rc : natural64;
    lpos : List;
    choice : string(1..2) := "  ";
    ans : character := 'y';
    method : character := '0';
    noqsols : natural32 := 0;
    first : boolean := true; -- total degree only when first time

  begin
    if verbose > 0
     then put_line("-> in main_root_counters.Polynomial_Main ...");
    end if;
    declare
    begin
      rc := Total_Degree(p);
    exception
      when others => rc := 0;
    end;
    new_line(file); put_line(file,"ROOT COUNTS :"); new_line(file);
    put(file,"total degree : ");
    if rc > 0
     then put(file,rc,1); -- put(rc,1);
     else High_Total_Degree(file,p);
    end if;
    new_line(file);
    loop
      Display_General_Menu(rc,own,method);
      if own then
        put("Type a number between 0 and 9, or A, B" 
            & " preceded by i for info : ");
        Ask_Alternative(choice,"0123456789AB",'i');
      else
        put("Type a number between 0 and 9, or A" 
            & " preceded by i for info : ");
        Ask_Alternative(choice,"0123456789A",'i');
      end if;
      if choice(1) = 'i' then
        method := choice(2);
        Display_Info(method);
        put("Do you want to apply this root count ? (y/n) ");
        Ask_Yes_or_No(ans);
      else
        method := choice(1);
      end if;
      exit when not first and (method = '0');
      first := false;
      if ans = 'y' then
        Apply_General_Method(file,integer32(nt),p,q,qsols,method,rc,lpos);
        noqsols := Length_Of(qsols);
        if method /= '0' then
          new_line;
          put("The current root count equals "); put(rc,1); put_line(".");
          if noqsols /= 0 then
            put("The number of start solutions equals ");
            put(noqsols,1); put_line(".");
          end if;
          put("Do you want to perform more root counting ? (y/n) ");
          Ask_Yes_or_No(ans);
        else
          ans := 'n';
        end if;
      else
        ans := 'y';
      end if;
      exit when ans /= 'y';
    end loop;
    roco := natural32(rc);
    Clear(lpos);
  end Polynomial_Main;

  procedure Laurent_Main
               ( file : in file_type; nt : in natural32;
                 p,q : in out Laur_Sys; qsols : in out Solution_List;
                 roco : out natural32; verbose : in integer32 := 0 ) is

    rc : natural32 := 0;
    choice : character := '0';

  begin
    if verbose > 0
     then put_line("-> in main_root_counters.Laurent_Main ...");
    end if;
    loop
      Display_Laurent_Menu(rc,choice);
      put("Type a number in the range from 0 to 6 : ");
      Ask_Alternative(choice,"0123456");
      exit when (choice = '0');
      Apply_Laurent_Method(file,integer32(nt),p,q,qsols,choice,rc);
    end loop;
    roco := rc;
  end Laurent_Main;
 
  procedure Read_System
              ( filename : in string; lp : out Link_to_Laur_Sys ) is

    file : file_type;

  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      Close(file);
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

  function Is_Square ( p : Laur_Sys ) return boolean is

    use Standard_Complex_Laurentials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if nv = nq then
      return true;
    else
      new_line;
      put("The number of equations : "); put(nq,1); new_line;
      put("The number of variables : "); put(nv,1); new_line;
      put_line("The system is not square!");
      new_line;
      return false;
    end if;
  end Is_Square;

  procedure Count_Roots
              ( nt : in natural32; outfilename : in string;
                p : in out Laur_Sys; v : in integer32 := 0 ) is

    outft : file_type;
    q : Laur_Sys(p'range);
    qsols : Solution_List;
    rc : natural32;

  begin
    if v > 0
     then put_line("-> in main_root_counters.Count_Roots 1 ...");
    end if;
    Create_Output_File(outft,outfilename);
    put(outft,p'last,1);
    new_line(outft);
    put(outft,p);
    Laurent_Main(outft,nt,p,q,qsols,rc,v-1);
    new_line(outft);
    put_line(outft,Bye_Bye_Message);
    Write_Seed_Number(outft);
    put_line(outft,Greeting_Banners.Version);
    Close(outft);
  end Count_Roots;

  procedure Count_Roots
              ( nt : in natural32; outfilename : in string;
                p : in out Poly_Sys; v : in integer32 := 0 ) is

    outft : file_type;
    q : Poly_Sys(p'range);
    qsols : Solution_List;
    rc : natural32;

  begin
    if v > 0
     then put_line("-> in main_root_counters.Count_Roots 2 ...");
    end if;
    Create_Output_File(outft,outfilename);
    put(outft,p'last,1);
    new_line(outft);
    put(outft,p);
    Polynomial_Main(outft,nt,p,q,false,qsols,rc,v-1);
    new_line(outft);
    put_line(outft,Bye_Bye_Message);
    Write_Seed_Number(outft);
    put_line(outft,Greeting_Banners.Version);
    Close(outft);
  end Count_Roots;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    lp : Link_to_Laur_Sys := null;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_root_counters.Main ...");
    end if;
    Read_System(infilename,lp);
    if lp = null
     then new_line; get(lp);
    end if;
    if Is_Square(lp.all) then
      if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
        Count_Roots(nt,outfilename,lp.all,verbose-1);
      else
        declare
          use Standard_Laur_Poly_Convertors;
          q : Poly_Sys(lp'range)
            := Positive_Laurent_Polynomial_System(lp.all);
        begin
          Count_Roots(nt,outfilename,q,verbose-1);
        end;
      end if;
    end if;
  end Main;

end Main_Root_Counters;
