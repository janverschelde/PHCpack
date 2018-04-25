with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_numbers;             use Double_Double_numbers;
with Quad_Double_numbers;               use Quad_Double_numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Query_Matrices;           use Standard_Query_Matrices;
with Standard_Numerical_Rank;           use Standard_Numerical_Rank;
with Standard_Nullity_Matrices;         use Standard_Nullity_Matrices;
with Standard_Multiplicity_Structure;   use Standard_Multiplicity_Structure;
with DoblDobl_Nullity_Matrices;         use DoblDobl_Nullity_Matrices;
with DoblDobl_Query_Matrices;           use DoblDobl_Query_Matrices;
with DoblDobl_Numerical_Rank;           use DoblDobl_Numerical_Rank;
with DoblDobl_Multiplicity_Structure;   use DoblDobl_Multiplicity_Structure;
with QuadDobl_Nullity_Matrices;         use QuadDobl_Nullity_Matrices;
with QuadDobl_Query_Matrices;           use QuadDobl_Query_Matrices;
with QuadDobl_Numerical_Rank;           use QuadDobl_Numerical_Rank;
with QuadDobl_Multiplicity_Structure;   use QuadDobl_Multiplicity_Structure;
 
procedure ts_multip is

-- DESCRIPTION :
--   Testing of the computation of the multiplicity structure.

  procedure Compute_Dimensions_of_Nullity_Matrix is

  -- DESCRIPTION :
  --   Interactive routine to compute dimensions of the nullity matrix.

    nq,nv,nr,nc,k : natural32 := 0;

  begin
    put("Give the number of equations : "); get(nq);
    put("Give the number of variables : "); get(nv);
    loop
      put("Give multiplicity stage (0 to exit) : "); get(k);
      exit when (k=0);
      Standard_Nullity_Matrices.Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      put("  number of rows : "); put(nr,1); new_line;
      put("  number of columns : "); put(nc,1); new_line;
    end loop;
  end Compute_Dimensions_of_Nullity_Matrix;

  procedure Write_Integers ( a : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix as a matrix of integers.

    use Standard_Complex_Numbers;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        put(" "); put(integer32(REAL_PART(a(i,j))),1);
      end loop;
      new_line;
    end loop;
  end Write_Integers;

  procedure Write_Integers ( a : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix as a matrix of integers.

    use DoblDobl_Complex_Numbers;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        put(" "); put(integer32(to_double(REAL_PART(a(i,j)))),1);
      end loop;
      new_line;
    end loop;
  end Write_Integers;

  procedure Write_Integers ( a : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix as a matrix of integers.

    use QuadDobl_Complex_Numbers;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        put(" "); put(integer32(to_double(REAL_PART(a(i,j)))),1);
      end loop;
      new_line;
    end loop;
  end Write_Integers;

  procedure Standard_Evaluation_of_Nullity_Matrix is

  -- DESCRIPTION :
  --   Evaluates the nullity matrix of a polynomial system at zero,
  --   in standard double precision.

    use Standard_Complex_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    f : Link_to_Poly_Sys;
    nq,nv,nr,nc,k : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(f);
    nq := natural32(f'last);
    nv := Number_of_Unknowns(f(f'first));
    put("Give multiplicity stage : "); get(k);
    if k > 0 then
      Standard_Nullity_Matrices.Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva : Matrix(1..integer32(nr),1..integer32(nc));
        z : constant Standard_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => Standard_Complex_Numbers.Create(0.0));
        rnk : natural32;
      begin
        eva := Evaluate_Nullity_Matrix(standard_output,nq,nv,nr,nc,k,f.all,z);
        Show_Matrix(eva);
        Write_Integers(eva);
        rnk := Numerical_Rank(eva,1.0E-8);
        put("The numerical rank : "); put(rnk,1); new_line;
      end;
    end if;
  end Standard_Evaluation_of_Nullity_Matrix;

  procedure DoblDobl_Evaluation_of_Nullity_Matrix is

  -- DESCRIPTION :
  --   Evaluates the nullity matrix of a polynomial system at zero,
  --   in double double precision.

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    f : Link_to_Poly_Sys;
    nq,nv,nr,nc,k : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(f);
    nq := natural32(f'last);
    nv := Number_of_Unknowns(f(f'first));
    put("Give multiplicity stage : "); get(k);
    if k > 0 then
      DoblDobl_Nullity_Matrices.Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva : Matrix(1..integer32(nr),1..integer32(nc));
        zero : constant double_double := create(0.0);
        z : constant DoblDobl_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => DoblDobl_Complex_Numbers.Create(zero));
        rnk : natural32;
      begin
        eva := Evaluate_Nullity_Matrix(standard_output,nq,nv,nr,nc,k,f.all,z);
        Show_Matrix(eva);
        Write_Integers(eva);
        rnk := Numerical_Rank(eva,1.0E-8);
        put("The numerical rank : "); put(rnk,1); new_line;
      end;
    end if;
  end DoblDobl_Evaluation_of_Nullity_Matrix;

  procedure QuadDobl_Evaluation_of_Nullity_Matrix is

  -- DESCRIPTION :
  --   Evaluates the nullity matrix of a polynomial system at zero,
  --   in quad double precision.

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    f : Link_to_Poly_Sys;
    nq,nv,nr,nc,k : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(f);
    nq := natural32(f'last);
    nv := Number_of_Unknowns(f(f'first));
    put("Give multiplicity stage : "); get(k);
    if k > 0 then
      QuadDobl_Nullity_Matrices.Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva : Matrix(1..integer32(nr),1..integer32(nc));
        zero : constant quad_double := create(0.0);
        z : constant QuadDobl_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => QuadDobl_Complex_Numbers.Create(zero));
        rnk : natural32;
      begin
        eva := Evaluate_Nullity_Matrix(standard_output,nq,nv,nr,nc,k,f.all,z);
        Show_Matrix(eva);
        Write_Integers(eva);
        rnk := Numerical_Rank(eva,1.0E-8);
        put("The numerical rank : "); put(rnk,1); new_line;
      end;
    end if;
  end QuadDobl_Evaluation_of_Nullity_Matrix;

  procedure Standard_Computation_of_Multiplicity_Structure is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and solutions and
  --   computes the multiplicity structure in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    f : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions...");
    Standard_System_and_Solutions_io.get(f,sols);
    Driver_to_Multiplicity_Structure(standard_output,f.all,sols);
  end Standard_Computation_of_Multiplicity_Structure;

  procedure DoblDobl_Computation_of_Multiplicity_Structure is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and solutions and
  --   computes the multiplicity structure in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    f : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions...");
    DoblDobl_System_and_Solutions_io.get(f,sols);
    Driver_to_Multiplicity_Structure(standard_output,f.all,sols);
  end DoblDobl_Computation_of_Multiplicity_Structure;

  procedure QuadDobl_Computation_of_Multiplicity_Structure is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and solutions and
  --   computes the multiplicity structure in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    f : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions...");
    QuadDobl_System_and_Solutions_io.get(f,sols);
    Driver_to_Multiplicity_Structure(standard_output,f.all,sols);
  end QuadDobl_Computation_of_Multiplicity_Structure;

  procedure Compute_Basis_of_Dual_Space is

    use Standard_Complex_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    f : Link_to_Poly_Sys;
    nq,nv,k,nr,nc : natural32 := 0;
    tol : double_float := 1.0E-8;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(f);
    nq := natural32(f'last);
    nv := Number_of_Unknowns(f(f'first));
    put("Give the differential order : "); get(k);
    loop
      put("Current tolerance for numerical rank is "); put(tol,3); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans = 'n';
      put("Give new value for tolerance : "); get(tol);
    end loop;
    if k > 0 then
      Standard_Nullity_Matrices.Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva,b : Matrix(1..integer32(nr),1..integer32(nc));
        z : Standard_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => Standard_Complex_Numbers.Create(1.0E-6));
      begin
        eva := Evaluate_Nullity_Matrix(nq,nv,nr,nc,k,f.all,z);
        Orthogonal_Basis(nr,nc,eva,b);
      end;
    end if;
  end Compute_Basis_of_Dual_Space;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Multiplicity structure of an isolated singular zero.");
    new_line;
    put_line("Choose on of the following test routines : ");
    put_line("  1. compute the dimensions of the nullity matrix;");
    put_line("  2. evaluate nullity matrix at zero with doubles;");
    put_line("  3. evaluate nullity matrix at zero with double doubles;");
    put_line("  4. evaluate nullity matrix at zero with quad doubles;");
    put_line("  5. compute the multiplicity structure with doubles;");
    put_line("  6. compute the multiplicity structure with double doubles;");
    put_line("  7. compute the multiplicity structure with quad doubles;");
    put_line("  8. compute a basis for the dual space.");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to choose : ");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Compute_Dimensions_of_Nullity_Matrix;
      when '2' => Standard_Evaluation_of_Nullity_Matrix;
      when '3' => DoblDobl_Evaluation_of_Nullity_Matrix;
      when '4' => QuadDobl_Evaluation_of_Nullity_Matrix;
      when '5' => Standard_Computation_of_Multiplicity_Structure;
      when '6' => DoblDobl_Computation_of_Multiplicity_Structure;
      when '7' => QuadDobl_Computation_of_Multiplicity_Structure;
      when '8' => Compute_Basis_of_Dual_Space;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_multip;
