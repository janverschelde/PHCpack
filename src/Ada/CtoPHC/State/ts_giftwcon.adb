with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
with Multprec_Lattice_3d_Facets;
with Multprec_Lattice_3d_Facets_io;     use Multprec_Lattice_3d_Facets_io;
with Multprec_Lattice_4d_Facets;
with Multprec_Lattice_4d_Facets_io;     use Multprec_Lattice_4d_Facets_io;
with Convex_Hull_Methods;
with Multprec_Giftwrap_Container;
with Facets_and_Strings;

procedure ts_giftwcon is

-- DESCRIPTION :
--   Interactive test on the container to store the results of
--   the giftwrapping method in 3 and 4 dimensions.
--   Also the string representation of the facets is shown.

  procedure Test_3d ( n : integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of n points in 3-space.

    A : Multprec_Integer_Matrices.Matrix(1..3,1..n)
      := Convex_Hull_Methods.Random_Data(3,n);
    nbf : natural32;
    fcn : integer32 := 0;
    lft : Multprec_Lattice_3d_Facets.Link_to_3d_Facet;

    use Multprec_Lattice_3d_Facets;

  begin
    put_line("The point configuration :"); put(A);
    put_line("Constructing the convex hull ...");
    Multprec_Giftwrap_Container.Create(A);
    nbf := Multprec_Giftwrap_Container.Number_of_3d_Facets;
    put("Number of facets computed : "); put(nbf,1); new_line;
    loop
      put("Give a facet number (-1 to exit) : "); get(fcn);
      exit when (fcn < 0);
      lft := Multprec_Giftwrap_Container.Facet_3d_Data(natural32(fcn));
      if lft = null then
        put("null pointer for facet number "); put(fcn,1); new_line;
      else
        Write_Facet(A,lft.all);
        put_line("The string representation :");
        put_line(Facets_and_Strings.write(A,lft.all));
      end if;
    end loop; 
  end Test_3d; 

  procedure Test_4d ( n : integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of n points in 4-space.

    A : Multprec_Integer_Matrices.Matrix(1..4,1..n)
      := Convex_Hull_Methods.Random_Data(4,n);
    nbf : natural32;
    fcn : integer32 := 0;
    lft : Multprec_Lattice_4d_Facets.Link_to_4d_Facet;

    use Multprec_Lattice_4d_Facets;

  begin
    put_line("The point configuration :"); put(A);
    put_line("Constructing the convex hull ...");
    Multprec_Giftwrap_Container.Create(A);
    nbf := Multprec_Giftwrap_Container.Number_of_4d_Facets;
    put("Number of facets computed : "); put(nbf,1); new_line;
    loop
      put("Give a facet number (-1 to exit) : "); get(fcn);
      exit when (fcn < 0);
      lft := Multprec_Giftwrap_Container.Facet_4d_Data(natural32(fcn));
      if lft = null then
        put("null pointer for facet number "); put(fcn,1); new_line;
      else
        Write_4D_Facet(A,lft);
        put_line("The string representation :");
        put_line(Facets_and_Strings.write(A,lft.all));
      end if;
    end loop; 
  end Test_4d; 

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select a dimension, 3 or 4,
  --   and then launches the corresponding test.

    ans : character;
    n : integer32 := 0;

  begin
    new_line;
    put_line("Testing the giftwrapping container ...");
    put("What is the dimension, 3 or 4 ? ");
    Ask_Alternative(ans,"34");
    new_line;
    put("Give the number of points : "); get(n);
    if ans = '3'
     then Test_3d(n);
     else Test_4d(n);
    end if;
  end Main;

begin
  Main;
end ts_giftwcon;
