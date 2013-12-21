with text_io,integer_io;                  use text_io,integer_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;        use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;     use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;       use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Permutations;                        use Permutations;
with Permute_Operations;                  use Permute_Operations;

procedure ts_reorder is

-- DESCRIPTION :
--   Test on reordering the variables in a polynomial happens
--   in two stages :
--     1) permuting the variables in the polynomial;
--     2) permuting the symbols in the table.

  procedure Write_Symbol_Table is

  -- DESCRPITION :
  --   Writes the content of the symbol table in the current order.
  begin
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); put(" ");
    end loop;
    new_line;
  end Write_Symbol_Table;

  procedure Write_Symbol_Order ( perm : in Permutation ) is

  -- DESCRIPTION :
  --   Writes the symbols according to the given permutation.

  begin
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(perm(i))); put(" ");
    end loop;
  end Write_Symbol_Order;

  procedure Permute_Symbol_Table ( perm : in Permutation ) is

  -- DESCRIPTION :
  --   Permutes the symbols in the symbol table.

    nb : constant natural := Symbol_Table.Number;
    use Symbol_Table;
    syms : array(1..nb) of Symbol;

  begin
    for i in syms'range loop
      syms(i) := Symbol_Table.Get(i);
    end loop;
    for i in syms'range loop
      Symbol_Table.Replace(i,syms(perm(i)));
    end loop;
  end Permute_Symbol_Table;

  procedure Read_Permutation ( perm : out Permutation ) is

    nb : constant natural := Symbol_Table.Number;
    use Symbol_Table;
    sb : Symbol;
    error : boolean;

  begin
    for i in 1..nb loop
      loop
        put("Give symbol at place "); put(i,1); put(" : ");
        Symbol_Table_io.Get(sb);
        put("  Your symbol : "); Symbol_Table_io.put(sb);
        perm(i) := Symbol_Table.Get(sb);     
        error := false;
        if perm(i) = 0
         then put_line(" did not occur before!  Please try again...");
              error := true;
         else for j in 1..i-1 loop
                if perm(j) = perm(i)
                 then error := true;
                      put(" occurs already at new place ");
                      put(j,1); put_line("!  Please try again...");
                end if;
                exit when error;
              end loop;
              if not error
               then put(" was at position "); put(perm(i),1); new_line;
              end if;
        end if;
        exit when not error;
      end loop;
    end loop;
  end Read_Permutation;

  procedure Change_Order ( p : in Poly ) is

  -- DESCRIPTION :
  --   Reorders the symbols in the polynomial.

    nbunk : constant natural := Number_of_Unknowns(p);
    perm : Standard_Integer_Vectors.Vector(1..nbunk);
    pp : Poly;
    ans : character;

  begin
    put("Number of unknowns  : "); put(nbunk,1); new_line; 
    put("The order of the unknowns : "); Write_Symbol_Table;
    loop
      Read_Permutation(Permutation(perm));
     -- put("Give a permutation : "); get(perm);
      put("Your permutation : "); put(perm); new_line;
      put("The new order of symbols : ");
      Write_Symbol_Order(Permutation(perm));
      new_line;
      put("Happy with this new order ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'y');
    end loop;
    pp := Permutation(perm)*p;
    put("The permuted polynomial : "); put(pp); new_line;
    Permute_Symbol_Table(Permutation(perm));
    put("The permuted order : "); Write_Symbol_Table;
    put("The permuted polynomial : "); put(pp); new_line;
  end Change_Order;

  procedure Change_Order ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Reorders the symbols in the polynomial.

    nbunk : constant natural := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..nbunk);
    pp : Poly_Sys(p'range);
    ans : character;

  begin
    put("Number of unknowns  : "); put(nbunk,1); new_line; 
    put("The order of the unknowns : "); Write_Symbol_Table;
    loop
      Read_Permutation(Permutation(perm));
     -- put("Give a permutation : "); get(perm);
      put("The permutation : "); put(perm); new_line;
      put("The new order of symbols : ");
      Write_Symbol_Order(Permutation(perm));
      new_line;
      put("Happy with this new order ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'y');
    end loop;
    pp := p*Permutation(perm);
    put_line("The permuted polynomial system : "); 
    put(pp);
    Permute_Symbol_Table(Permutation(perm));
    put("The permuted order : "); Write_Symbol_Table;
    put_line("The permuted polynomial system : ");
    put(pp);
  end Change_Order;

  procedure Main is

    n : natural;
    p : Poly;
    lp : Link_to_Poly_Sys;
    ans : character;

  begin
    new_line;
    put_line("Reorder the unknowns in a polynomial.");
    new_line;
    put("Work on polynomial or on polynomial system ? (p/s) ");
    Ask_Alternative(ans,"ps");
    if ans = 's'
     then get(lp);
          Change_Order(lp.all);
     else put("Give the number of variables : "); get(n);
          Symbol_Table.Init(n);
          put("Give a polynomial : "); get(p);
          put("your polynomial : "); put(p); new_line;
          Change_Order(p);
    end if;
  end Main;

begin
  Main;
end ts_reorder;
