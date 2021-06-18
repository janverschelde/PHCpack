with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Run_Power_Series_Methods;           use Run_Power_Series_Methods;
with Main_Laurent_Series_Newton;         use Main_Laurent_Series_Newton;
with Regular_Newton_Puiseux;
with Series_Path_Trackers;
with Interactive_Pade_Trackers;
with Track_Path_Convolutions;
with Multitasked_Path_Convolutions;
with Newton_Fabry_on_Homotopy;

package body Main_Pade_Trackers is

  procedure Run_Regular_Newton_Puiseux ( valprc : in character ) is
  begin
    put_line("Using as lifting the powers of the first variable,");
    put_line("assuming coefficients are sufficiently generic ...");
    case valprc is
      when '1' =>
        put_line("The working precision is double precision");
        Regular_Newton_Puiseux.Standard_Main;
      when '2' =>
        put_line("The working precision is double double precision.");
        Regular_Newton_Puiseux.DoblDobl_Main;
      when '4' =>
        put_line("The working precision is quad double precision.");
        Regular_Newton_Puiseux.QuadDobl_Main;
      when others => null;
    end case;
  end Run_Regular_Newton_Puiseux;

  procedure Run_Power_Series_Newton
              ( infilename,outfilename : in string;
                valprc : in character; vrb : in integer32 := 0 ) is

    ans : character;

  begin
    new_line;
    put("Start Newton's method at a constant term ? (y/n) ");
    Ask_Yes_or_No(ans);
    case valprc is
      when '1' =>
        put_line("The working precision is double precision");
        if ans = 'y'
         then Standard_Main_at_Constant(infilename,outfilename,vrb-1);
         else Standard_Main_at_Series(infilename,outfilename,vrb-1);
        end if;
      when '2' =>
        put_line("The working precision is double double precision.");
        if ans = 'y'
         then DoblDobl_Main_at_Constant(infilename,outfilename,vrb-1);
         else DoblDobl_Main_at_Series(infilename,outfilename,vrb-1);
        end if;
      when '3' =>
        put_line("The working precision is triple double precision.");
        if ans = 'y'
         then TripDobl_Main_at_Constant(infilename,outfilename,vrb-1);
        end if;
      when '4' =>
        put_line("The working precision is quad double precision.");
        if ans = 'y'
         then QuadDobl_Main_at_Constant(infilename,outfilename,vrb-1);
         else QuadDobl_Main_at_Series(infilename,outfilename,vrb-1);
        end if;
      when '5' =>
        put_line("The working precision is penta double precision.");
        if ans = 'y'
         then PentDobl_Main_at_Constant(infilename,outfilename,vrb-1);
        end if;
      when '6' =>
        put_line("The working precision is octo double precision.");
        if ans = 'y'
         then OctoDobl_Main_at_Constant(infilename,outfilename,vrb-1);
        end if;
      when '7' =>
        put_line("The working precision is deca double precision.");
        if ans = 'y'
         then DecaDobl_Main_at_Constant(infilename,outfilename,vrb-1);
        end if;
      when others => null;
    end case;
  end Run_Power_Series_Newton;

  procedure Run_Path_Trackers
              ( valprc : in character; vrb : in integer32 := 0 ) is

    ans : character;

  begin
    new_line;
    put("Step-by-step interactive execution ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case valprc is
        when '1' => Interactive_Pade_Trackers.Standard_Main(vrb-1);
        when '2' => Interactive_Pade_Trackers.DoblDobl_Main(vrb-1);
        when '4' => Interactive_Pade_Trackers.QuadDobl_Main(vrb-1);
        when others => null;
      end case;
    else
      case valprc is
        when '1' => Series_Path_Trackers.Standard_Main(vrb-1);
        when '2' => Series_Path_Trackers.DoblDobl_Main(vrb-1);
        when '4' => Series_Path_Trackers.QuadDobl_Main(vrb-1);
        when others => null;
      end case;
    end if;
  end Run_Path_Trackers;

  procedure Run_Path_Convolution_Trackers
              ( nbtasks : in natural32; valprc : in character;
                vrb : in integer32 := 0 ) is

    ans : character;

  begin
    if nbtasks > 0 then
      case valprc is
        when '1' => Multitasked_Path_Convolutions.Standard_Main(nbtasks,vrb-1);
        when '2' => Multitasked_Path_Convolutions.DoblDobl_Main(nbtasks,vrb-1);
        when '4' => Multitasked_Path_Convolutions.QuadDobl_Main(nbtasks,vrb-1);
        when others => null;
      end case;
    else
      new_line;
      put("Apply multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'n' then
        case valprc is
          when '1' => Track_Path_Convolutions.Standard_Main(vrb-1);
          when '2' => Track_Path_Convolutions.DoblDobl_Main(vrb-1);
          when '4' => Track_Path_Convolutions.QuadDobl_Main(vrb-1);
          when others => null;
        end case;
      else
        case valprc is
          when '1' => Multitasked_Path_Convolutions.Standard_Main(0,vrb-1);
          when '2' => Multitasked_Path_Convolutions.DoblDobl_Main(0,vrb-1);
          when '4' => Multitasked_Path_Convolutions.QuadDobl_Main(0,vrb-1);
          when others => null;
        end case;
      end if;
    end if;
  end Run_Path_Convolution_Trackers;

  function Prompt_for_Precision_Level return character is

    res : character;

  begin
    new_line;
    put_line("MENU for the precision :");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(res,"1234567");
    return res;
  end Prompt_for_Precision_Level;

  procedure Run_Regular_Newton_Puiseux is

    prc : constant character := Prompt_for_Precision_Level;

  begin
    Run_Regular_Newton_Puiseux(prc);
  end Run_Regular_Newton_Puiseux;

  procedure Run_Power_Series_Newton
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    prc : constant character := Prompt_for_Precision_Level;

  begin
    Run_Power_Series_Newton(infilename,outfilename,prc,vrb-1);
  end Run_Power_Series_Newton;

  procedure Run_Path_Trackers ( vrb : in integer32 := 0 ) is

    prc : constant character := Prompt_for_Precision_Level;

  begin
    Run_Path_Trackers(prc,vrb-1);
  end Run_Path_Trackers;

  procedure Run_Path_Convolution_Trackers
              ( nbt : in natural32; vrb : in integer32 := 0 ) is

    prc : constant character := Prompt_for_Precision_Level;

  begin
    Run_Path_Convolution_Trackers(nbt,prc,vrb-1);
  end Run_Path_Convolution_Trackers;

  procedure Run_Newton_Fabry
              ( nbtasks : in natural32; prc : in character;
                vrblvl : in integer32 := 0 ) is

    valprc : character := prc;

  begin
    if vrblvl > 0 then
      put_line("-> in main_pade_trackers.Run_Newton_Fabry ...");
    end if;
    if valprc = '0'
     then valprc := Newton_Fabry_on_Homotopy.Prompt_for_Precision;
    end if;
    Newton_Fabry_on_Homotopy.Run_Newton_Fabry(nbtasks,valprc,vrblvl);
  end Run_Newton_Fabry;

  function Prompt_for_Method return character is

    res : character;

  begin
    new_line;
    put_line("MENU for power series methods :");
    put_line("  1. apply polyhedral methods for tropisms");
    put_line("  2. run Newton's method starting at a series or a point");
    put_line("  3. track paths with Pade approximants as predictor");
    put_line("  4. run a faster Newton-Fabry-Pade-Hesse path tracker");
    put_line("  5. compute the convergence radius of a solution series");
    put("Type 1, 2, 3, 4, or 5 to select the method : ");
    Ask_Alternative(res,"12345");
    return res;
  end Prompt_for_Method;

  procedure Nonzero_Precision_Main
              ( infilename,outfilename : in string;
                nbtasks : in natural32; valprc : in character;
                vrb : in integer32 := 0 ) is

    method,ans : character;

  begin 
    if vrb > 0 then
      put_line("-> in main_pade_trackers.Nonzero_Precision_Main ...");
    end if;
    method := Prompt_for_Method;
    case method is
      when '1' => Run_Regular_Newton_Puiseux(valprc);
      when '2' => 
        new_line;
        put("Run on Laurent series ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          Run_Laurent_Series_Newton(infilename,outfilename,vrb-1);
        else
          Run_Power_Series_Newton(infilename,outfilename,valprc,vrb-1);
        end if;
      when '3' => Run_Path_Trackers(valprc,vrb-1);
      when '4' => Run_Path_Convolution_Trackers(nbtasks,valprc,vrb-1);
      when '5' => Run_Newton_Fabry(nbtasks,valprc,vrb-1);
      when others => null;
    end case;
  end Nonzero_Precision_Main;

  procedure Zero_Precision_Main
              ( infilename,outfilename : in string;
                nbtasks : in natural32; vrb : in integer32 := 0 ) is

    method,ans : character;

  begin 
    if vrb > 0 then
      put_line("-> in main_pade_trackers.Nonzero_Precision_Main ...");
    end if;
    method := Prompt_for_Method;
    case method is
      when '1' => Run_Regular_Newton_Puiseux;
      when '2' =>
        new_line;
        put("Run on Laurent series ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          Run_Laurent_Series_Newton(infilename,outfilename,vrb-1);
        else
          Run_Power_Series_Newton(infilename,outfilename,vrb-1);
        end if;
      when '3' => Run_Path_Trackers(vrb-1);
      when '4' => Run_Path_Convolution_Trackers(nbtasks,vrb-1);
      when '5' => Run_Newton_Fabry(nbtasks,'0',vrb-1);
      when others => null;
    end case;
  end Zero_Precision_Main;

  procedure Main ( infilename,outfilename : in string;
                   nbtasks : in natural32; precision : in character;
                   verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put_line("-> in main_pade_trackers.Main ...");
    end if;
    if precision = '0' then
      Zero_Precision_Main(infilename,outfilename,nbtasks,verbose-1);
    else
      Nonzero_Precision_Main
        (infilename,outfilename,nbtasks,precision,verbose-1);
    end if;
  end Main;

end Main_Pade_Trackers;
