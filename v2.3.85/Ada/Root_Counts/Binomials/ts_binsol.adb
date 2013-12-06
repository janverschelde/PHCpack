with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Affine_Binomial_Iterator;
with Standard_Affine_Binomials;
with Standard_Permanent_Factors;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;         use Standard_Monomial_Maps_io;

procedure ts_binsol is

-- DESCRIPTION :
--   Interactive development of computing all affine solution sets
--   of a binomial system.

  procedure Subset_Iterator ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Calls an iterator to enumerate all subsets of cardinality
  --   at most m of a set of n variables.
  --   This code was developed as warm up for the actual iterator.

    cnt : integer32 := 0;
    choice : Standard_Integer_Vectors.Vector(1..n);
    count : integer32;

  begin
    put_line("enumeration of subsets with iterator : ");
    Affine_Binomial_Iterator.Initialize_Subsets(n,m);
    loop
      Affine_Binomial_Iterator.Next_Subset(choice,count);
      exit when (count < 0);
      cnt := cnt + 1;
      put(cnt,3);
      put(" : choice : "); put(choice);
      put(" : "); put(count,1); new_line;
    end loop;
    Affine_Binomial_Iterator.Clear_Subsets;
  end Subset_Iterator;

  procedure Call_Subset_Iterator is

  -- DESCRIPTION :
  --   Prompts the user for the number n of variables 
  --   and the cardinality m of the subsets.  Then all subsets 
  --   of cardinality at most m of a set of n variables are shown.

    n,m : integer32 := 0;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the cardinality of the subset : "); get(m);
    Subset_Iterator(n,m);
  end Call_Subset_Iterator;

  procedure Frequency_Table ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   For every variable in the system, we count the number of monomials
  --   the variable occurs with a positive exponent.  If the variable
  --   occurs anywhere with a negative exponent, then the count is -1.
  --   The current implementation uses "extract_two_terms"
  --   which is not the fastest way.

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    f : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    t1,t2 : Term;
    fail : boolean;

  begin
    for i in p'range loop
      Standard_Affine_Binomials.Extract_Two_Terms(p(i),t1,t2,fail);
      if fail then
        put("Polynomial "); put(i,1); put_line(" is not binomial!");
        exit;
      else
        for k in 1..n loop
          if f(k) >= 0 then  -- once negative remains negative
            if t1.dg(k) < 0 or t2.dg(k) < 0 then
              f(k) := -1;
            else
              if t1.dg(k) > 0 
               then f(k) := f(k) + 1;
              end if;
              if t2.dg(k) > 0 
               then f(k) := f(k) + 1;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    put("The frequency table : "); put(f); new_line;
    put("The order of the symbols :");
    for i in 1..n loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(i));
      begin
        put(" "); Symbol_Table_io.put(sb);
      end;
    end loop;
    new_line;
  end Frequency_Table;

  procedure Compute_Frequency_Table is

  -- DESCRIPTION :
  --   Prompts the user for a binomial system and then computes the
  --   frequency table of the occurring variables for use in a greedy
  --   algorithm to compute affine solution sets.

    p : Link_to_Laur_Sys;

  begin
    put_line("Reading a binomial system ...");
    get(p);
    Frequency_Table(p.all);
  end Compute_Frequency_Table;

  function Common_Factor
             ( t1,t2 : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the common factor of the terms t1 and t2.

    res : Standard_Integer_Vectors.Vector(t1'range);

  begin
    put("  t1 : "); put(t1); 
    put("  t2 : "); put(t2); new_line; 
    for i in res'range loop
      if t1(i) > 0 and t2(i) > 0 then
        if t1(i) <= t2(i)
         then res(i) := t1(i);
         else res(i) := t2(i);
        end if;
      else
        res(i) := 0;
      end if;
    end loop;
    return res;
  end Common_Factor;

  function Is_Zero ( v : Standard_Integer_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the vector v consists entirely of zeroes.

  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  procedure Monomial_Common_Factors ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Scans the polynomials of p and extracts the common monomial factors
  --   of the polynomials in p.

    t1,t2 : Term;
    fail : boolean;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    f : Standard_Integer_Vectors.Vector(1..n);
    cnt : integer32 := 0;

  begin
    new_line;
    for i in p'range loop
      Standard_Affine_Binomials.Extract_Two_Terms(p(i),t1,t2,fail);
      if not fail then 
        f := Common_Factor(t1.dg.all,t2.dg.all);
       -- put(i,1); put(" : "); put(f); new_line;
        if not Is_Zero(f) then
          put("polynomial "); put(i,1);
          put(" has common factor : "); put(f); new_line;
          cnt := cnt + 1;
        end if;
      end if;
    end loop;
    put("Found "); put(cnt,1); put_line(" common monomial factors.");
  end Monomial_Common_Factors;

  procedure Find_Monomial_Common_Factors is

  -- DESCRIPTION :
  --   Prompts the user for a binomial system 
  --   and then looks for monomial common factors.

    p : Link_to_Laur_Sys;

  begin
    put_line("Reading a binomial system ...");
    get(p);
    Monomial_Common_Factors(p.all);
  end Find_Monomial_Common_Factors;

  procedure Write_to_File ( maps : in Monomial_Map_List ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and writes the maps to file.

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of an output file to write the maps ...");
    Read_Name_and_Create_File(file);
    Insert_Parameter_Symbols(natural32(Head_Of(maps).d));
    put(file,maps);
    close(file);
  end Write_to_File;

  procedure Recursive_Affine_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a binomial system 
  --   and then enumerates candidate subsets for affine solution sets
  --   with a recursive procedure.

    p : Link_to_Laur_Sys;
    sols : Monomial_Map_List;
    fail : boolean;
    use Standard_Permanent_Factors;

  begin
    put_line("Reading a binomial system ...");
    get(p);
    Interactive_Affine_Solutions_with_Recursion(p.all,sols,fail);
    if Length_Of(sols) > 0 then
      put("Found "); put(Length_Of(sols),1); put_line(" solution maps.");
      Write_to_File(sols);
    end if;
  end Recursive_Affine_Solver;

  procedure Iterative_Affine_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a binomial system 
  --   and then enumerates candidate subsets for affine solution sets
  --   with an iterator.

    p : Link_to_Laur_Sys;
    sols : Monomial_Map_List;
    fail : boolean;
    use Standard_Permanent_Factors;

  begin
    put_line("Reading a binomial system ...");
    get(p);
    Interactive_Affine_Solutions_with_Iterator(p.all,sols,fail);
    if Length_Of(sols) > 0 then
      put("Found "); put(Length_Of(sols),1); put_line(" solution maps.");
      Write_to_File(sols);
    end if;
  end Iterative_Affine_Solver;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Computing affine solution sets of binomial systems.");
    new_line;
    put_line("MENU for testing the operations on affine solution sets :");
    put_line("  1. enumerate all subsets of a set of variables;");
    put_line("  2. compute frequency table of the occurring variables;");
    put_line("  3. find common factors in a binomial system.");
    put_line("  4. enumerate candidate subsets for affine sets recursively;");
    put_line("  5. compute affine solution sets with iterator.");
    put("Type 1, 2, 3, 4, or 5 to select : "); Ask_Alternative(ans,"12345");
    new_line;
    case ans is
      when '1' => Call_Subset_Iterator;
      when '2' => Compute_Frequency_Table;
      when '3' => Find_Monomial_Common_Factors;
      when '4' => Recursive_Affine_Solver;
      when '5' => Iterative_Affine_Solver;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_binsol;
