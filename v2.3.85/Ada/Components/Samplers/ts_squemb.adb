with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Communications_with_User;          use Communications_with_User;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;

procedure ts_squemb is

-- DESCRIPTION :
--   A polynomial system with any number of equations and unknowns
--   can be turned into a square system:
--     1) by adding random slices in case it is underdetermined;
--     2) by adding slack variables in case it is overdetermined. 
--   To find generic points on all k-dimensional solution components
--   of the system, we add k extra hyperplanes to the square system.
--   For underdetermined systems, the embedding clears a number of
--   hyperplanes that were added in the squaring process. 

  procedure Write_Symbol_Table is
  begin
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); put(" ");
    end loop;
    new_line;
  end Write_Symbol_Table;

  function Maximum ( a,b : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  procedure Ask_to_Save ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   The user is asked whether the system p should be save to a file.

    ans : character;
    outfile : file_type;

  begin
    put("Do you want to save the system to a file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(outfile);
      put_line(outfile,p);
      close(outfile);
    end if;
  end Ask_to_Save;

  procedure Square_and_Embed ( p : in Poly_Sys; k : in natural32 ) is

  -- DESCRIPTION :
  --   The system p is made square and k extra hyperplanes are added to
  --   compute generic points on all k-dimensional solution components.

    nbequ : constant natural32 := natural32(p'length);
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : Poly_Sys(1..integer32(max)) := Square(p);
    ep : Poly_Sys(1..integer32(max+k));

  begin
    if nbequ = nbunk then
      put_line("The system is square.");
    else
      put_line("The system is not square.");
      put("Number of unknowns  : "); put(nbunk,1); new_line; 
      put_line("The unknowns :"); Write_Symbol_Table;
      put("Number of equations : "); put(nbequ,1); new_line; 
      if nbequ > nbunk
       then Add_Slack_Symbols(nbequ-nbunk);
      end if;
      put_line("The system after squaring : "); put_line(sp);
      Ask_to_Save(sp);
    end if;
    if k > 0 then
      Add_Embed_Symbols(k);
      if nbunk > nbequ then
        for i in nbunk-k+1..nbunk loop
          exit when (i = nbequ);  -- do not wipe out original eqs
          Clear(sp(integer32(i)));
        end loop;
      end if;
      ep := Slice_and_Embed(sp,k);
      put_line("The system after embedding : "); put_line(ep);
      Ask_to_Save(ep);
    end if;
  end Square_and_Embed;

  procedure Main is

    lp : Link_to_Poly_Sys;
    k : natural32 := 0;

  begin
    new_line;
    put_line("Turning a polynomial system into a square system,");
    put_line("by adding random slices or slack variables.");
    new_line;
    get(lp);
    new_line;
    put("Give the suspected top dimension of solution components : ");
    get(k);
    Square_and_Embed(lp.all,k);
  end Main;

begin
  Main;
end ts_squemb;
