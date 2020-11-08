with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Poly_Randomizers;  use Standard_Complex_Poly_Randomizers;
with Standard_Complex_Laur_Randomizers;  use Standard_Complex_Laur_Randomizers;
with Main_Poly_Continuation;             use Main_Poly_Continuation;

package body Drivers_for_Coefficient_Systems is

  function Show_Menu_and_Prompt return character is

  -- DESCRIPTION :
  --   Shows the user the menu and returns the character
  --   corresponding to the user's choice.

    res : character;

  begin
    new_line;
    put_line("MENU for Polyhedral Continuation : ");
    put_line("  0. No polyhedral continuation, leave the menu.");
    put_line("  1. Create and solve random coefficient system.");
    put_line("  2. Solve given system by polyhedral continuation.");
    put("Type 0,1, or 2 to choose : ");
    Ask_Alternative(res,"012");
    return res;
  end Show_Menu_and_Prompt;

  procedure Driver_for_Coefficient_System
              ( file : in file_type; p : in Poly_Sys; k : in natural32;
                byebye : in boolean;
                q : out Poly_Sys; qfile,solsfile : in out file_type;
                tosolve,ranstart,contrep : out boolean ) is

    ans : constant character := Show_Menu_and_Prompt;
    oc : natural32;
    qq : Poly_Sys(p'range);

  begin
    tosolve := (ans /= '0');
    ranstart := (ans = '1');
    if ans /= '0' then
      if ans = '1' then
        put_line("Reading a file name to write random coefficient system.");
        Read_Name_and_Create_File(qfile);
        qq := Complex_Randomize1(p); q := qq;
        if k = 0 then
          new_line(file);
          put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
          new_line(file);
          put_line(file,qq);
          put_line(qfile,qq);
        end if;
      else
        q := p;
        put_line("Reading a name of a file to write start solutions on.");
        Read_Name_and_Create_File(solsfile);
      end if;
      new_line;
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      contrep := (oc /= 0);
    end if;
    if byebye then
      new_line;
      put_line("No more input expected.  See output file for results.");
      new_line;
    end if;
  end Driver_for_Coefficient_System;

  procedure Driver_for_Coefficient_System
              ( file : in file_type; p : in Laur_Sys; k : in natural32;
                byebye : in boolean;
                q : out Laur_Sys; qfile,solsfile : in out file_type;
                tosolve,ranstart,contrep : out boolean ) is

    ans : constant character := Show_Menu_and_Prompt;
    oc : natural32;
    qq : Laur_Sys(p'range);

  begin
    tosolve := (ans /= '0');
    ranstart := (ans = '1');
    if ans /= '0' then
      if ans = '1' then
        put_line("Reading a file name to write random coefficient system.");
        Read_Name_and_Create_File(qfile);
        qq := Complex_Randomize1(p); q := qq;
        if k = 0 then
          new_line(file);
          put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
          new_line(file);
          put_line(file,qq);
          put_line(qfile,qq);
        end if;
      else
        q := p;
        put_line("Reading a name of a file to write start solutions on.");
        Read_Name_and_Create_File(solsfile);
      end if;
      new_line;
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      contrep := (oc /= 0);
    end if;
    if byebye then
      new_line;
      put_line("No more input expected.  See output file for results.");
      new_line;
    end if;
  end Driver_for_Coefficient_System;

end Drivers_for_Coefficient_Systems;
