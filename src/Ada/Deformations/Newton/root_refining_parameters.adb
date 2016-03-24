with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Numbers_io;                         use Numbers_io;

package body Root_Refining_Parameters is

  procedure Standard_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean ) is
  begin
    epsxa := 1.0E-8;        -- precision for correction on x
    epsfa := 1.0E-8;        -- precision for residual 
    tolsing := 1.0E-8;      -- tolerance on inverse condition numbers
    maxit := 3;             -- maximal number of Newton iterations
    deflate := true;        -- if deflation is wanted
    wout := false;          -- if intermediate output is wanted
  end Standard_Default_Root_Refining_Parameters;

  procedure DoblDobl_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean ) is
  begin
    epsxa := 1.0E-12;       -- precision for correction on x
    epsfa := 1.0E-12;       -- precision for residual 
    tolsing := 1.0E-8;     -- tolerance on inverse condition numbers
    maxit := 3;             -- maximal number of Newton iterations
    deflate := true;        -- if deflation is wanted
    wout := false;          -- if intermediate output is wanted
  end DoblDobl_Default_Root_Refining_Parameters;

  procedure QuadDobl_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean ) is
  begin
    epsxa := 1.0E-16;       -- precision for correction on x
    epsfa := 1.0E-16;       -- precision for residual 
    tolsing := 1.0E-8;     -- tolerance on inverse condition numbers
    maxit := 3;             -- maximal number of Newton iterations
    deflate := true;        -- if deflation is wanted
    wout := false;          -- if intermediate output is wanted
  end QuadDobl_Default_Root_Refining_Parameters;

  procedure Multprec_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out Floating_Number;
                maxit,deci : out natural32; deflate,wout : out boolean ) is
  begin
    epsxa := Create(1.0E-8);       -- precision for correction on x
    epsfa := Create(1.0E-8);       -- precision for residual
    tolsing := Create(1.0E-8);     -- tolerance on inverse condition numbers
    maxit := 3;                    -- maximal number of Newton iterations
    deci := 16;                    -- number of decimal places
    deflate := false;              -- if deflation is wanted
    wout := false;                 -- if intermediate output is wanted
  end Multprec_Default_Root_Refining_Parameters;

  procedure Standard_Put_Root_Refining_Parameters
              ( file : in file_type; epsxa,epsfa,tolsing : in double_float;
                maxit : in natural32; deflate,wout : in boolean ) is
  begin
    put(file,"  1. output during the iterations    : ");
    if wout
     then put(file," yes"); new_line(file);
     else put(file," no"); new_line(file);
    end if;
    put(file,"  2. deflate isolated singularities  : ");
    if deflate
     then put(file," yes"); new_line(file);
     else put(file," no"); new_line(file);
    end if;
    put(file,"  3. tolerance for error on the root : ");
    put(file,epsxa,2,3,3); new_line(file);
    put(file,"  4. tolerance for the residual      : ");
    put(file,epsfa,2,3,3); new_line(file);
    put(file,"  5. tolerance for singular roots    : ");
    put(file,tolsing,2,3,3); new_line(file);
    put(file,"  6. maximum number of iterations    : ");
    put(file,maxit,2); new_line(file);
  end Standard_Put_Root_Refining_Parameters;

  procedure Multprec_Put_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in Floating_Number;
                maxit,deci : in natural32; deflate,wout : in boolean ) is
  begin
    put(file,"  1. output during the iterations    : ");
    if wout
     then put(file," yes"); new_line(file);
     else put(file," no"); new_line(file);
    end if;
    put(file,"  2. deflate isolated singularities  : ");
    if deflate
     then put(file," yes"); new_line(file);
     else put(file," no"); new_line(file);
    end if;
    put(file,"  3. tolerance for error on the root : ");
    put(file,epsxa,2,3,3); new_line(file);
    put(file,"  4. tolerance for the residual      : ");
    put(file,epsfa,2,3,3); new_line(file);
    put(file,"  5. tolerance for singular roots    : ");
    put(file,tolsing,2,3,3); new_line(file);
    put(file,"  6. maximum number of iterations    : ");
    put(file,maxit,2); new_line(file);
    put(file,"  7. number of decimal places        : ");
    put(file,deci,2); new_line(file);
  end Multprec_Put_Root_Refining_Parameters;

  procedure Standard_Menu_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in out double_float;
                maxit : in out natural32; deflate,wout : in out boolean ) is

    ans : character;

  begin
    new_line;
    loop
      put_line("MENU with current Settings for the Root Refiner :");
      Standard_Put_Root_Refining_Parameters
        (Standard_Output,epsxa,epsfa,tolsing,maxit,deflate,wout);
      put("Type 1, 2, 3, 4, 5, or 6 to change, type 0 to exit : ");
      Ask_Alternative(ans,"0123456");
      exit when ans = '0';
      case ans is
        when '1' => put("Do you want output during the iterations ? (y/n) ");
                    Ask_Yes_or_No(ans); wout := (ans = 'y');
        when '2' => put("Are isolated singularities to be deflated ? (y/n) ");
                    Ask_Yes_or_No(ans); deflate := (ans = 'y');
        when '3' => put("Give new tolerance for error on the root : ");
                    Read_Double_Float(epsxa);
        when '4' => put("Give new tolerance for residual : ");
                    Read_Double_Float(epsfa);
        when '5' => put("Give new tolerance for singular roots : ");
                    Read_Double_Float(tolsing);
        when '6' => put("Give new maximum number of iterations : ");
                    Read_Natural(maxit);
        when others => null;
      end case;
    end loop;
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS : ");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
  end Standard_Menu_Root_Refining_Parameters;

  procedure Multprec_Menu_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in out Floating_Number;
                maxit,deci : in out natural32;
                deflate,wout : in out boolean ) is

    ans : character;

  begin
    new_line;
    loop
      put_line("MENU with current Settings for the Root Refiner :");
      Multprec_Put_Root_Refining_Parameters
        (Standard_Output,epsxa,epsfa,tolsing,maxit,deci,deflate,wout);
      put("Type 1, 2, 3, 4, 5, 6, or 7 to change, type 0 to exit : ");
      Ask_Alternative(ans,"01234567");
      exit when ans = '0';
      case ans is
        when '1' => put("Do you want output during the iterations ? (y/n) ");
                    Ask_Yes_or_No(ans); wout := (ans = 'y');
        when '2' => put("Are isolated singularities to be deflated ? (y/n) ");
                    Ask_Yes_or_No(ans); deflate := (ans = 'y');
        when '3' => put("Give new tolerance for error on the root : ");
                    get(epsxa);
        when '4' => put("Give new tolerance for residual : ");
                    get(epsfa);
        when '5' => put("Give new tolerance for singular roots : ");
                    get(tolsing);
        when '6' => put("Give new maximum number of iterations : ");
                    Read_Natural(maxit);
        when '7' => put("Give new number of decimal places : ");
                    Read_Natural(deci);
        when others => null;
      end case;
    end loop;
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS : ");
    Multprec_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deci,deflate,wout);
  end Multprec_Menu_Root_Refining_Parameters;

end Root_Refining_Parameters;
