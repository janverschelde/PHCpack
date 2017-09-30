with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Add_and_Remove_Embedding;           use Add_and_Remove_Embedding;
with Drivers_to_Cascade_Filtering;       use Drivers_to_Cascade_Filtering;

procedure ts_cascade is

-- DESCRIPTION :
--   A cascade of homotopies is used to compute witness points on
--   all positive dimensional solution components.

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Witness Generate: a cascade of homotopies for witness points");
    new_line;
    put_line("Choose one of the following :");
    put_line("  1. generate an embedding for the top dimensional component;");
    put_line("  2. given the top embedding, compute all witness points.");
    put("Type 1 or 2 to choose : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Driver_to_Square_and_Embed("",""); -- empty file names
      when '2' => Driver_to_Witness_Generate(0,"","");
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cascade;
