with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Brackets,Brackets_io;             use Brackets,Brackets_io;
with Bracket_Monomials;                use Bracket_Monomials;
with Bracket_Monomials_io;             use Bracket_Monomials_io;

procedure ts_brackmons is

  procedure Test_Input_Output is

    bm : Bracket_Monomial;

  begin
    new_line;
    put_line("Give a bracket monomial, terminate with semicolon ; ");
    get(bm);
    put("Your bracket monomial : "); put(bm);
  end Test_Input_Output;

  procedure Test_Multiplication is

    d,m : integer32 := 0;
    bm : Bracket_Monomial;

  begin
    put("Give the number of entries in the brackets : "); get(d);
    put("Give the number of brackets : "); get(m);
    for i in 1..m loop
      declare
        b : Bracket(1..d);
      begin
        put("Give "); put(d,1); put(" numbers for the ");
        put(i,1); put("th bracket : "); get(b);
        Multiply(bm,b);
      end;
    end loop;
    put_line("The bracket monomial : "); put(bm); new_line;
  end Test_Multiplication;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for testing the operations on bracket monomials :");
    put_line("  1. test input and output;");
    put_line("  2. test multiplication of monomials.");
    put("Type 1 or 2 to make a choice : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_Input_Output;
     else Test_Multiplication;
    end if;
  end Main;

begin
  Main;
end ts_brackmons;
