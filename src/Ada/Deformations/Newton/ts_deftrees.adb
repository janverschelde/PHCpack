with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Deflation_Trees;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;
with DoblDobl_Deflation_Trees;
with DoblDobl_Deflation_Trees_io;       use DoblDobl_Deflation_Trees_io;
with QuadDobl_Deflation_Trees;
with QuadDobl_Deflation_Trees_io;       use QuadDobl_Deflation_Trees_io;

procedure ts_deftrees is

-- DESCRIPTION :
--   Tests the creation of deflation trees.

  procedure Create ( nd : in out Standard_Deflation_Trees.Node ) is

    ans : character;
    m : integer32 := 0;

    use Standard_Deflation_Trees;

  begin
    put("Do you wish to create a child at level ");
    put(nd.d,1); put(" ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the number of multipliers : "); get(m);
      Create_Child(nd,m);
      Create(nd.children(m).all);
      Create(nd);
    end if;
  end Create;

  procedure Create ( nd : in out DoblDobl_Deflation_Trees.Node ) is

    ans : character;
    m : integer32 := 0;

    use DoblDobl_Deflation_Trees;

  begin
    put("Do you wish to create a child at level ");
    put(nd.d,1); put(" ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the number of multipliers : "); get(m);
      Create_Child(nd,m);
      Create(nd.children(m).all);
      Create(nd);
    end if;
  end Create;

  procedure Create ( nd : in out QuadDobl_Deflation_Trees.Node ) is

    ans : character;
    m : integer32 := 0;

    use QuadDobl_Deflation_Trees;

  begin
    put("Do you wish to create a child at level ");
    put(nd.d,1); put(" ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the number of multipliers : "); get(m);
      Create_Child(nd,m);
      Create(nd.children(m).all);
      Create(nd);
    end if;
  end Create;

  procedure Interactive_Create
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Deflation_Trees;

    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Node(ne,nv) := Create_Root(p);

  begin
    Create(t);
    put_line("The tree : ");
    Write(standard_output,"output",t);
  end Interactive_Create;

  procedure Interactive_Create
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Deflation_Trees;

    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Node(ne,nv) := Create_Root(p);

  begin
    Create(t);
    put_line("The tree : ");
    Write(standard_output,"output",t);
  end Interactive_Create;

  procedure Interactive_Create
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Deflation_Trees;

    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Node(ne,nv) := Create_Root(p);

  begin
    Create(t);
    put_line("The tree : ");
    Write(standard_output,"output",t);
  end Interactive_Create;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select a level of precision,
  --   reads a polynomial systems and launches the test.

    ans : character;
    stlp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    ddlp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qdlp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the operations in deflation trees ...");
    new_line;
    put_line("MENU for the precision of the coefficients :");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    put_line("Reading the system for the root of the tree.");
    case ans is
      when '0' => get(stlp); Interactive_Create(stlp.all);
      when '1' => get(ddlp); Interactive_Create(ddlp.all);
      when '2' => get(qdlp); Interactive_Create(qdlp.all);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_deftrees;
