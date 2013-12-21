with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
--with Standard_Deflate_Symbols;          use Standard_Deflate_Symbols;
with Standard_Deflation_Trees;          use Standard_Deflation_Trees;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;

procedure ts_deftrees is

  procedure Create ( nd : in out Node ) is

    ans : character;
    m : integer32 := 0;

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

  procedure Interactive_Create ( p : in Poly_Sys ) is

    ne : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    t : Node(ne,nv) := Create_Root(p);

  begin
    Create(t);
    put_line("The tree : ");
    Write(standard_output,"output",t);
  end Interactive_Create;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the operations in deflation trees ...");
    new_line;
    put_line("Reading the system for the root of the tree.");
    get(lp);
    Interactive_Create(lp.all);
  end Main;

begin
  Main;
end ts_deftrees;
