with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Symbol_Table;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Homotopy_Series_Readers;
with Singular_Values_of_Hessians;

procedure ts_hesscrit is

-- DESCRIPTION :
--   Interactive development testing of the Hessian criterion
--   for solution curves defined by polynomial homotopies.

  procedure Standard_Distances
              ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in Standard_Complex_Hessians.Array_of_Hessians;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    eta : double_float;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      eta := Singular_Values_of_Hessians.Standard_Distance(jm,hess,ls.all);
      put("eta : "); put(eta,2); new_line;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Standard_Distances;

  procedure DoblDobl_Distances
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in DoblDobl_Complex_Hessians.Array_of_Hessians;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    eta : double_double;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      eta := Singular_Values_of_Hessians.DoblDobl_Distance(jm,hess,ls.all);
      put("eta : "); put(eta,2); new_line;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end DoblDobl_Distances;

  procedure QuadDobl_Distances
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in QuadDobl_Complex_Hessians.Array_of_Hessians;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    eta : quad_double;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      eta := Singular_Values_of_Hessians.QuadDobl_Distance(jm,hess,ls.all);
      put("eta : "); put(eta,2); new_line;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Distances;

  procedure Standard_Test
              ( nbeq : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in Standard_Homotopy,
  --   computes the array of Hessians.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := Standard_Homotopy.Homotopy_System;
    s : Symbol_Table.Symbol;
    h : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;

    use Singular_Values_of_Hessians;

  begin
    Standard_Jacobian_Hessians_of_Homotopy(jm,h);
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : "); put(p);
    Standard_Distances(jm.all,h.all,sols);
  end Standard_Test;

  procedure DoblDobl_Test
              ( nbeq : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in DoblDobl_Homotopy,
  --   computes the array of Hessians.

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : Symbol_Table.Symbol;
    h : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;

    use Singular_Values_of_Hessians;

  begin
    DoblDobl_Jacobian_Hessians_of_Homotopy(jm,h);
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : "); put(p);
    DoblDobl_Distances(jm.all,h.all,sols);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( nbeq : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in QuadDobl_Homotopy,
  --   computes the array of Hessians.

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : Symbol_Table.Symbol;
    h : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;

    use Singular_Values_of_Hessians;

  begin
    QuadDobl_Jacobian_Hessians_of_Homotopy(jm,h);
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : "); put(p);
    QuadDobl_Distances(jm.all,h.all,sols);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in double precision.

    nbeq : integer32 := 0;
    sols : Standard_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := Standard_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    Standard_Test(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in double double precision.

    nbeq : integer32 := 0;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := DoblDobl_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    DoblDobl_Test(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in quad double precision.

    nbeq : integer32 := 0;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := QuadDobl_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    QuadDobl_Test(nbeq,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches
  --   the corresponding test.

    ans : constant character := Prompt_for_Precision;

  begin
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hesscrit;
