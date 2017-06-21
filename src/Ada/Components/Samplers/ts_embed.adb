with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;       use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;       use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;    use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with Witness_Sets,Witness_Sets_io;        use Witness_Sets,Witness_Sets_io;

procedure ts_embed is

-- DESCRIPTION :
--   Test on slicing and embedding polynomial systems.

  procedure Write_Symbol_Table is

  begin
    put("The symbol table:");
    for i in 1..Symbol_Table.Number loop
      put(" ");
      Symbol_Table_io.put(Symbol_Table.Get(i));
    end loop;
    new_line;
  end Write_Symbol_Table;

  procedure Square_Slice_and_Embed
              ( file : in file_type; p : in Poly_Sys; k : in natural32 ) is

  -- NOTE :
  --   This procedure does not add dummy equations in case p
  --   is underdetermined.

    use Standard_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'last));

  begin
   -- put("Number of equations : "); put(nq,1); new_line;
   -- put("Number of variables : "); put(nv,1); new_line;
    if nv < nq then
      declare
        sp : constant Poly_Sys := Square(p);
        ep : Poly_Sys(sp'first..sp'last+integer32(k));
      begin
        Add_Slack_Symbols(nq-nv);
       -- put_line(file,sp);
        Add_Embed_Symbols(k);
        ep := Slice_and_Embed(sp,k);
        put_line(file,ep);
      end;
    elsif nv = nq then
      Add_Embed_Symbols(k);
      declare
        ep : constant Poly_Sys(p'first..p'last+integer32(k))
           := Slice_and_Embed(p,k);
      begin
        put_line(file,ep);
      end;
    else
      if nv - nq > k then
        Add_Embed_Symbols(nv-nq-k);
        declare
          ep : constant Poly_Sys(p'first..p'last+integer32(k))
             := Slice_and_embed(p,k);
        begin
          put_line(file,ep);
        end;
      end if;
    end if;
  end Square_Slice_and_Embed;

  procedure Square_Slice_and_Embed
              ( file : in file_type; p : in Laur_Sys; k : in natural32 ) is

  -- NOTE :
  --   This procedure does not add dummy equations in case p
  --   is underdetermined.

    use Standard_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'last));

  begin
   -- put("Number of equations : "); put(nq,1); new_line;
   -- put("Number of variables : "); put(nv,1); new_line;
    if nv < nq then
      declare
        sp : constant Laur_Sys := Square(p);
        ep : Laur_Sys(sp'first..sp'last+integer32(k));
      begin
        Add_Slack_Symbols(nq-nv);
       -- put_line(file,sp);
        Add_Embed_Symbols(k);
        ep := Slice_and_Embed(sp,k);
        put_line(file,ep);
      end;
    elsif nv = nq then
      Add_Embed_Symbols(k);
      declare
        ep : constant Laur_Sys(p'first..p'last+integer32(k))
           := Slice_and_Embed(p,k);
      begin
        put_line(file,ep);
      end;
    else
      if nv - nq > k then
        Add_Embed_Symbols(nv-nq-k);
        declare
          ep : constant Laur_Sys(p'first..p'last+integer32(k))
             := Slice_and_embed(p,k);
        begin
          put_line(file,ep);
        end;
      end if;
    end if;
  end Square_Slice_and_Embed;

  function Remove_Last_Variables
             ( p : Poly_Sys; n : natural32 ) return Poly_Sys is

  -- DESCRIPTION :
  --   Removes the last n variables of the system p.

  -- REQUIRED : n >= Number_of_Unknowns(p(i)), for i in p'range.

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables ( p : in out Poly_Sys; n : in natural32 ) is

  -- DESCRIPTION :
  --   Removes the last n variables of the system p.

  -- REQUIRED : n >= Number_of_Unknowns(p(i)), for i in p'range.

    use Standard_Complex_Polynomials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  function Remove_Embedding
             ( p : Poly_Sys; dim,ns : natural32 ) return Poly_Sys is

  -- DESCRIPTION :
  --   Removes the embedding and extra slack variables from the system p.

  -- ON ENTRY :
  --   p       an embedded polynomial system;
  --   dim     dimension of the solution set used in the embedding;
  --   ns      number of extra slack variables which need to be removed.

  -- REQUIRED :
  --   All slack variables are located as last variables in p.

    res : Poly_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Poly_Sys := Remove_Embedding1(p,dim);
        res1 : Poly_Sys(wrk'range) := wrk;
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  procedure Add_Embedding is

  -- DESCRIPTION :
  --   Embeds a polynomial system, after adding linear equations
  --   and adding extra slack variables to an overconstrained system.

    lq : Link_to_Laur_Sys;
    file : file_type;
    k : natural32 := 0;

  begin
    new_line;
    put_line("Squaring, Slicing, and Embedding a Polynomial System.");
    new_line;
    get(lq);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of slices to add to the system : "); get(k);
    new_line;
    put_line("See the output file for the embedded system...");
    new_line;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Square_Slice_and_Embed(file,lq.all,k);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : Poly_Sys(lq'range)
          := Positive_Laurent_Polynomial_System(lq.all);
      begin
        Square_Slice_and_Embed(file,p,k);
      end;
    end if;
    Close(file);
  end Add_Embedding;

  procedure Remove_Embedding is

  -- DESCRIPTION :
  --   Removes the embedding from a polynomial system.

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    file : file_type;
    k,ns : natural32 := 0;

  begin
    new_line;
    put_line("Removing an Embedding of a Polynomial System.");
    Standard_Read_Embedding(lp,sols,k,ns);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Write_Symbol_Table;
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system...");
    new_line;
    declare
      rp : constant Poly_Sys := Remove_Embedding(lp.all,k,ns);
      nq : constant natural32 := natural32(rp'last);
      nv : constant natural32
         := Standard_Complex_Polynomials.Number_of_Unknowns(rp(rp'first));
      rsols : Solution_List;
    begin
      if k + ns > 0
       then rsols := Remove_Embedding(sols,k+ns);
       else rsols := sols;
      end if;
      put(file,nq,1);
      put(file," ");
      put(file,nv,1);
      new_line(file);
      put(file,rp);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
    end;
    Close(file);
  end Remove_Embedding;

  procedure Parse_Embedding is

  -- DESCRIPTION :
  --   Reads in a polynomial system and parses the symbol table
  --   for the number of embedding variables.

    lp : Link_to_Poly_Sys;
    dim : natural32;

  begin
    new_line;
    put_line("Parsing an embedding from a given polynomial system.");
    get(lp);
    Write_Symbol_Table;
    dim := Count_Embed_Symbols(natural32(lp'last),"zz");
    put("Number of embed symbols : "); put(dim,1); new_line;
    Swap_Symbols_to_End(natural32(lp'last),dim,"zz",lp.all);
    put_line("After swapping embed symbols to the end : ");
    Write_Symbol_Table;
    if dim > 1 then
      Sort_Embed_Symbols(natural32(lp'last),natural32(lp'last)-dim,dim,lp.all);
    end if;
  end Parse_Embedding;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Adding and Removing Embeddings of Polynomial Systems.");
    new_line;
    put_line("Choose one of the following options : ");
    put_line("  1. square, slice, and embed a polynomial system;");
    put_line("  2. remove an embedding from a polynomial system;");
    put_line("  3. parse an embedding from any given system.");
    put("Type 1, 2, or 3 to select your choice : ");
    Ask_Alternative(ans,"123");
    if ans = '1' then
      Add_Embedding;
    elsif ans = '2' then
      Remove_Embedding;
    else
      Parse_Embedding;
    end if;
  end Main;

begin
  Main;
end ts_embed;
