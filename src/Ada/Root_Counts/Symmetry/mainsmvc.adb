with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Drivers_for_Implicit_Lifting;       use Drivers_for_Implicit_Lifting;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Drivers_for_Dynamic_Lifting;        use Drivers_for_Dynamic_Lifting;
with Drivers_for_Symmetric_Lifting;      use Drivers_for_Symmetric_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;

procedure mainsmvc ( nt : in natural32; infilename,outfilename : in string ) is

  procedure Read_System
              ( filename : in string; lq : in out Link_to_Laur_Sys ) is

  -- DESCRIPTION :
  --   Attempts to open the file with the given name for reading of
  --   a polynomial system.
  
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

  function Lifting_Strategy return natural32 is

  -- DESCRIPTION :
  --   Prompts the user to select a lifting strategy and returns this
  --   choice as a natural number between 0 and 4.

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
      put("Type 0, 1, 2, 3, or 4 to select,"
                & " eventually preceded by i for info : ");
      Ask_Alternative(choice,"01234",'i');
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
        when others => put_line("No information available.");
      end case;
      exit when choice(1) /= 'i';
    end loop;
    case choice(1) is
      when '1'    => return 1;
      when '2'    => return 2;
      when '3'    => return 3;
      when '4'    => return 4;
      when others => return 0;
    end case;
  end Lifting_Strategy;

  procedure Ask_and_Dispatch_Lifting_Strategy ( lq : in Link_to_Laur_Sys ) is

    outft : file_type;
    pq : Poly_Sys(lq'range);
    qq : Laur_Sys(lq'range);
    qsols,qsols0 : Solution_List;
    mv,smv,tmv : natural32;
    strategy : natural32;
    deflate : boolean := false;

  begin
    Create_Output_File(outft,outfilename);
    put(outft,natural32(lq'last),lq.all);
    strategy := Lifting_Strategy;
    new_line(outft);
    case strategy is
      when 0 => put_line(outft,"STATIC LIFTING");
                Driver_for_Mixed_Volume_Computation
                  (outft,integer32(nt),lq.all,true,qq,qsols,qsols0,mv,smv,tmv);
      when 1 => put_line(outft,"IMPLICIT LIFTING");
                Driver_for_Mixture_Bezout_BKK(outft,lq.all,true,qq,qsols,mv);
      when 2 => put_line(outft,"DYNAMIC LIFTING");
                Driver_for_Dynamic_Mixed_Volume_Computation
                  (outft,lq.all,true,qq,qsols,mv);
      when 3 => put_line(outft,"SYMMETRIC LIFTING"); 
                Driver_for_Symmetric_Mixed_Volume_Computation
                  (outft,lq.all,true,qq,qsols,mv);
      when 4 => put_line(outft,"MixedVol Algorithm to compute mixed volume");
                Driver_for_MixedVol_Algorithm
                  (outft,integer32(nt),lq.all,true,qq,qsols,qsols0,mv,smv,tmv);
      when others => null;
    end case;
    if Length_Of(qsols) > 0 then
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
    Close(outft);
  end Ask_and_Dispatch_Lifting_Strategy;

  procedure Main is

    lq : Link_to_Laur_Sys := null;

  begin
    Read_System(infilename,lq);
    if lq = null
     then new_line; get(lq);
    end if;
    Ask_and_Dispatch_Lifting_Strategy(lq);
  end Main;

begin
  Main;
end mainsmvc;
