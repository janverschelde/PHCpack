with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Standard_Query_Matrices;           use Standard_Query_Matrices;
with Standard_Numerical_Rank;           use Standard_Numerical_Rank;
with Standard_Nullity_Matrices;         use Standard_Nullity_Matrices;
with Standard_Multiplicity_Structure;   use Standard_Multiplicity_Structure;
 
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
      Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      put("  number of rows : "); put(nr,1); new_line;
      put("  number of columns : "); put(nc,1); new_line;
    end loop;
  end Compute_Dimensions_of_Nullity_Matrix;

  procedure Write_Integers ( a : in Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix as a matrix of integers.

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        put(" "); put(integer32(REAL_PART(a(i,j))),1);
      end loop;
      new_line;
    end loop;
  end Write_Integers;

  procedure Evaluation_of_Nullity_Matrix is

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
      Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva : Matrix(1..integer32(nr),1..integer32(nc));
        z : constant Standard_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => Create(0.0));
        rnk : natural32;
      begin
        eva := Evaluate_Nullity_Matrix(standard_output,nq,nv,nr,nc,k,f.all,z);
        Show_Matrix(eva);
        Write_Integers(eva);
        rnk := Numerical_Rank(eva,1.0E-8);
        put("The numerical rank : "); put(rnk,1); new_line;
      end;
    end if;
  end Evaluation_of_Nullity_Matrix;

  procedure Computation_of_Multiplicity_Structure is

    f : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions...");
    Standard_System_and_Solutions_io.get(f,sols);
    Driver_to_Multiplicity_Structure(standard_output,f.all,sols);
  end Computation_of_Multiplicity_Structure;

  procedure Compute_Basis_of_Dual_Space is

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
      Dimensions_of_Nullity_Matrix(nq,nv,k,nr,nc);
      declare
        eva,b : Matrix(1..integer32(nr),1..integer32(nc));
        z : Standard_Complex_Vectors.Vector(1..integer32(nv))
          := (1..integer32(nv) => Create(1.0E-6));
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
    put_line("  2. evaluate the nullity matrix for system at zero;");
    put_line("  3. compute the multiplicity structure at a zero;");
    put_line("  4. compute a basis for the dual space.");
    put("Type 1, 2, or 3 to choose : "); Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Compute_Dimensions_of_Nullity_Matrix;
      when '2' => Evaluation_of_Nullity_Matrix;
      when '3' => Computation_of_Multiplicity_Structure;
      when '4' => Compute_Basis_of_Dual_Space;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_multip;
