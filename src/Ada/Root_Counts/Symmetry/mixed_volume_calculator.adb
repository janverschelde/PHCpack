with Communications_with_User;           use Communications_with_User;
with Time_Stamps;                        use Time_Stamps;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Polynomial_Convertors;     use QuadDobl_Polynomial_Convertors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Drivers_for_Implicit_Lifting;       use Drivers_for_Implicit_Lifting;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Drivers_for_Dynamic_Lifting;        use Drivers_for_Dynamic_Lifting;
with Drivers_for_Symmetric_Lifting;      use Drivers_for_Symmetric_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;
with Drivers_for_DEMiCs_Algorithm;       use Drivers_for_DEMiCs_Algorithm;
with Write_Seed_Number;
with Greeting_Banners;

package body Mixed_Volume_Calculator is

  procedure Read_System
              ( filename : in string;
                lq : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    file : file_type;

  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lq);
      Close(file);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lq := null; return;
  end Read_System;

  function Prompt_for_Lifting return natural32 is

    choice : string(1..2) := "  ";
    ans : character;

  begin
    loop
      new_line;
      put_line("MENU with available Lifting Strategies (0 is default) :");
      put_line("  0. Static lifting     : lift points and prune lower hull.");
      put_line("  1. Implicit lifting   : based on recursive formula.");
      put_line("  2. Dynamic lifting    : incrementally add the points.");
      put_line
          ("  3. Symmetric lifting  : points in same orbit get same lifting.");
      put_line("  4. MixedVol Algorithm : a faster mixed volume computation.");
      put_line
          ("  5. DEMiCs Algorithm   : dynamic enumeration for mixed cells.");
      put("Type 0, 1, 2, 3, 4, or 5 to select,"
                & " eventually preceded by i for info : ");
      Ask_Alternative(choice,"012345",'i');
      exit when choice(1) /= 'i';
      new_line;
      case choice(2) is
        when '0'
          => Static_Lifting_Info; new_line;
             put("Do you want to apply static lifting ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '0';
             end if;
        when '1'
          => Implicit_Lifting_Info; new_line;
             put("Do you want to apply implicit lifting ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '1';
             end if;
        when '2'
          => Dynamic_Lifting_Info; new_line;
             put("Do you want to apply dynamic lifting ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '2';
             end if;
        when '3'
          => Symmetric_Lifting_Info; new_line;
             put("Do you want to apply implicit lifting ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '3';
             end if;
        when '4'
          => MixedVol_Algorithm_Info; new_line;
             put("Do you want to run the MixedVol Algorithm ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '4';
             end if;
        when '5'
          => DEMiCs_Algorithm_Info; new_line;
             put("Do you want to run the DEMiCs Algorithm ? (y/n) ");
             Ask_Yes_or_No(ans);
             if ans = 'y'
              then choice(1) := '5';
             end if;
        when others
          => put_line("No information available.");
      end case;
      exit when choice(1) /= 'i';
    end loop;
    case choice(1) is
      when '1'    => return 1;
      when '2'    => return 2;
      when '3'    => return 3;
      when '4'    => return 4;
      when '5'    => return 5;
      when others => return 0;
    end case;
  end Prompt_for_Lifting;

  procedure Call_MixedVol
               ( file : in file_type; nt : in natural32;
                 lq : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 v : in integer32 := 0 ) is

    d_qq : Standard_Complex_Laur_Systems.Laur_Sys(lq'range);
    d_qsols,d_qsols0 : Standard_Complex_Solutions.Solution_List;
    dd_p,dd_qq : DoblDobl_Complex_Laur_Systems.Laur_Sys(lq'range);
    dd_qsols,dd_qsols0 : DoblDobl_Complex_Solutions.Solution_List;
    qd_p,qd_qq : QuadDobl_Complex_Laur_Systems.Laur_Sys(lq'range);
    qd_qsols,qd_qsols0 : QuadDobl_Complex_Solutions.Solution_List;
    mv,smv,tmv : natural32;
    ans : character;
    nostart : boolean;
  
  begin
    if v > 0
     then put_line("-> in mixed_volume_calculator.Call_MixedVol ...");
    end if;
    new_line;
    put_line("MENU for precision to compute random coefficient system :");
    put_line("  0. do not solve a random coefficient system;");
    put_line("  1. run polyhedral homotopies in standard double precision;");
    put_line("  2. run polyhedral homotopies in double double precision;");
    put_line("  3. run polyhedral homotopies in quad double precision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    nostart := (ans = '0');
    case ans is
      when '0' | '1' =>
        Driver_for_MixedVol_Algorithm
          (file,integer32(nt),lq.all,true,nostart,
           d_qq,d_qsols,d_qsols0,mv,smv,tmv,verbose=>v-1);
      when '2' =>
        dd_p := Standard_Laur_Sys_to_DoblDobl_Complex(lq.all);
        Driver_for_MixedVol_Algorithm
          (file,integer32(nt),dd_p,true,nostart,
           dd_qq,dd_qsols,dd_qsols0,mv,smv,tmv,verbose=>v-1);
      when '3' =>
        qd_p := Standard_Laur_Sys_to_QuadDobl_Complex(lq.all);
        Driver_for_MixedVol_Algorithm
          (file,integer32(nt),qd_p,true,nostart,
           qd_qq,qd_qsols,qd_qsols0,mv,smv,tmv,verbose=>v-1);
      when others => null;
    end case;
  end Call_MixedVol;

  procedure Lift_Set_and_Run
              ( nt : in natural32; outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                lq : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                v : in integer32 := 0 ) is

    outft : file_type;
    pq : Poly_Sys(lq'range);
    qq : Standard_Complex_Laur_Systems.Laur_Sys(lq'range);
    qsols,qsols0 : Standard_Complex_Solutions.Solution_List;
    mv,smv,tmv : natural32;
    strategy : natural32;
    deflate : boolean := false;
    ended_moment : Ada.Calendar.Time;

  begin
    if v > 0
     then put_line("-> in mixed_volume_calculator.Lift_Set_and_Run ...");
    end if;
    Create_Output_File(outft,outfilename);
    put(outft,natural32(lq'last),lq.all);
    strategy := Prompt_for_Lifting;
    new_line(outft);
    case strategy is
      when 0 =>
        put_line(outft,"Static lifting applies lift-and-prune algorithms");
        Driver_for_Mixed_Volume_Computation
          (outft,integer32(nt),lq.all,true,qq,qsols,qsols0,mv,smv,tmv,
           vrblvl=>v-1);
      when 1 =>
        put_line(outft,
          "Implicit lifting applies a recursive formula for the mixed volume");
        Driver_for_Mixture_Bezout_BKK(outft,lq.all,true,qq,qsols,mv);
      when 2 =>
        put_line
          (outft,"Dynamic lifting to place points in a regular subdivsion");
        Driver_for_Dynamic_Mixed_Volume_Computation
          (outft,lq.all,true,qq,qsols,mv);
      when 3 =>
        put_line(outft,"Symmetric lifting to exploit permutation symmetry"); 
        Driver_for_Symmetric_Mixed_Volume_Computation
          (outft,lq.all,true,qq,qsols,mv);
      when 4 =>
        put_line(outft,"MixedVol Algorithm to compute the mixed volume");
        Call_MixedVol(outft,nt,lq,v=>v-1);
      when 5 =>
        put_line(outft,
          "DEMiCs Algorithm applies dynamic enumeration for all mixed cells");
        Driver_for_DEMiCs_Algorithm(outft,integer32(nt),lq.all,vrblvl=>v-1);
      when others => null;
    end case;
    if Standard_Complex_Solutions.Length_Of(qsols) > 0 then
      pq := Laurent_to_Polynomial_System(qq);
      declare
        epsxa,epsfa : constant double_float := 1.0E-8;
        tolsing : constant double_float := 1.0E-8;
        nb : natural32 := 0;
      begin
        new_line(outft);
        Reporting_Root_Refiner
          (outft,pq,qsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
      end;
    end if;
    new_line(outft);
    ended_moment := Ada.Calendar.Clock;
    put(outft,"phc -m ran from "); Write_Time_Stamp(outft,start_moment);
    put(outft," till "); Write_Time_Stamp(outft,ended_moment);
    put_line(outft,".");
    Write_Elapsed_Time(outft,start_moment,ended_moment);
    Write_Seed_Number(outft);
    put_line(outft,Greeting_Banners.Version);
    Close(outft);
  end Lift_Set_and_Run;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Standard_Complex_Laur_Systems;
    lq : Link_to_Laur_Sys := null;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mixed_volume_calculator.Main ...");
    end if;
    Read_System(infilename,lq);
    if lq = null
     then new_line; get(lq);
    end if;
    Lift_Set_and_Run(nt,outfilename,start_moment,lq,verbose-1);
  end Main;

end Mixed_Volume_Calculator;
