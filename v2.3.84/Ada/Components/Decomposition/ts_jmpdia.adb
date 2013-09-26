with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Jumpstart_Diagonal_Homotopies;     use Jumpstart_Diagonal_Homotopies;

procedure ts_jmpdia is

-- DESCRIPTION :
--   Interactive development of jumpstarting a diagonal homotopy.

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to jumpstart homotopies to intersect algebraic sets :");
    put_line("  1. first extrinsic homotopy at the start of the cascade;");
    put_line("  2. given witness set at one level, go down the cascade.");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Jumpstart_Diagonal_Homotopy;
     else Jumpstart_Cascade_Homotopy;
    end if;
  end Main;

begin
  Main;
end ts_jmpdia;
