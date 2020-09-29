with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Run_Power_Series_Methods;          use Run_Power_Series_Methods;
with Regular_Newton_Puiseux;
with Series_Path_Trackers;
with Interactive_Pade_Trackers;
with Track_Path_Convolutions;
with Multitasked_Path_Convolutions;

procedure mainseries ( nt : in natural32; precision : in character;
                       infilename,outfilename : in string;
                       verbose : in integer32 := 0 ) is

  procedure Nonzero_Precision_Main
              ( valprc : in character; vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   This main procedure is called when the precision is determined,
  --   either interactively or via the value of precision.
  --   The value for the precision in valprc is not '0', it is
  --   '1' for double, '2' for double double, or '4' for quad double.
  --   The procedure can be called immediately if the precision is
  --   set at the command line.

    ans : character;

  begin 
    if vrb > 0
     then put_line("-> in mainseries.Nonzero_Precision_Main ...");
    end if;
    new_line;
    put_line("MENU for power series methods :");
    put_line("  1. apply polyhedral methods for tropisms;");
    put_line("  2. run Newton's method starting at a series or a point;");
    put_line("  3. track paths with Pade approximants as predictor;");
    put_line("  4. run a faster Newton-Fabry-Pade-Hesse path tracker.");
    put("Type 1, 2, 3, or 4 to select the method : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' =>
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
      when '2' =>
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
          when '4' =>
            put_line("The working precision is quad double precision.");
            if ans = 'y'
             then QuadDobl_Main_at_Constant(infilename,outfilename,vrb-1);
             else QuadDobl_Main_at_Series(infilename,outfilename,vrb-1);
            end if;
          when others => null;
        end case;
      when '3' =>
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
      when '4' =>
        if nt > 0 then
          case valprc is
            when '1' => Multitasked_Path_Convolutions.Standard_Main(nt,vrb-1);
            when '2' => Multitasked_Path_Convolutions.DoblDobl_Main(nt,vrb-1);
            when '4' => Multitasked_Path_Convolutions.QuadDobl_Main(nt,vrb-1);
            when others => null;
          end case;
        else
          new_line;
          put("Apply multitasking ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            case valprc is
              when '1' => Multitasked_Path_Convolutions.Standard_Main(0,vrb-1);
              when '2' => Multitasked_Path_Convolutions.DoblDobl_Main(0,vrb-1);
              when '4' => Multitasked_Path_Convolutions.QuadDobl_Main(0,vrb-1);
              when others => null;
            end case;
          else
            case valprc is
              when '1' => Track_Path_Convolutions.Standard_Main(vrb-1);
              when '2' => Track_Path_Convolutions.DoblDobl_Main(vrb-1);
              when '4' => Track_Path_Convolutions.QuadDobl_Main(vrb-1);
              when others => null;
            end case;
          end if;
        end if;
      when others => null;
    end case;
  end Nonzero_Precision_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the working precision
  --   if precision /= '0' and then launches Nonzero_Precision_Main.

    prc : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainseries.Main ...");
    end if;
    if precision /= '0' then
      Nonzero_Precision_Main(precision,verbose-1);
    else
      new_line;
      put_line("MENU to select the working precision :");
      put_line("  0. standard double precision;");
      put_line("  1. double double precision;");
      put_line("  2. quad double precision.");
      put("Type 0, 1, or 2 to select the working precision : ");
      Ask_Alternative(prc,"012");
      case prc is
        when '0' => Nonzero_Precision_Main('1',verbose-1);
        when '1' => Nonzero_Precision_Main('2',verbose-1);
        when '2' => Nonzero_Precision_Main('4',verbose-1);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end Mainseries;
