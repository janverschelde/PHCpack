with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Algorithmic_DiffEval_Trackers;      use Algorithmic_DiffEval_Trackers;

procedure ts_ademanypaths is

-- DESCRIPTION :
--   Procedure to test the development of Newton's method and path trackers
--   using algorithmic differentiation for the evaluation.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precsion and then calls the test.

    tst,prc : character;

  begin
    new_line;
    put_line("MENU to test the algorithmic differentiation :");
    put_line("  0. running Newton's method;");
    put_line("  1. tracking one path;");
    put_line("  2. tracking many paths;");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(tst,"012");
    new_line;
    put_line("MENU to test the ADE code to track many paths :");
    put_line("  0. in standard double precision;");
    put_line("  1. in double double precision;");
    put_line("  2. in quad double precision.");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(prc,"012");
    case tst is
      when '0' =>
        case prc is
          when '0' => Standard_Newton;
          when '1' => DoblDobl_Newton;
          when '2' => QuadDobl_Newton;
          when others => null;
        end case;
      when '1' =>
        case prc is
          when '0' => Standard_Track_one_Path;
          when '1' => DoblDobl_Track_one_Path;
          when '2' => QuadDobl_Track_one_Path;
          when others => null;
        end case;
      when '2' =>
        case prc is
          when '0' => Standard_Track_many_Paths;
          when '1' => DoblDobl_Track_many_Paths;
          when '2' => QuadDobl_Track_many_Paths;
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_ademanypaths;
