with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with Symbol_Table;
with DoblDobl_Random_Polynomials;       use DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Randomizers;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;

package body Test_DoblDobl_Coeff_Homotopy is

  procedure Write_Elements ( A,B : in DoblDobl_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put(") = ");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put(") = ");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Write_Elements;

  procedure DoblDobl_Compared_Encapsulated_Eval ( n : in natural32 ) is

    use DoblDobl_Complex_Numbers;

    t : constant double_double := abs(DoblDobl_Random_Numbers.Random);
    ct : constant Complex_Number := Create(t);
    x : constant DoblDobl_Complex_Vectors.Vector(1..integer32(n))
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Homotopy.Eval(x,ct);
    z : constant DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Coefficient_Homotopy.Eval(x,ct);
    A : constant DoblDobl_Complex_Matrices.Matrix(x'range,x'range)
      := DoblDobl_Homotopy.Diff(x,ct);
    B : constant DoblDobl_Complex_Matrices.Matrix(x'range,x'range)
      := DoblDobl_Coefficient_Homotopy.Diff(x,ct);

  begin
    put("A random t : "); put(ct); new_line;
    put_line("A random point : "); put_line(x);
    put_line("-> y = "); put_line(y);
    put_line("-> z = "); put_line(z);
    Write_Elements(A,B);
  end DoblDobl_Compared_Encapsulated_Eval;

  procedure DoblDobl_Random_Systems
              ( n : in integer32;
                p,q : out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    m,d : natural32 := 0;

  begin
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    for i in p'range loop
      p(i) := Random_Sparse_Poly(natural32(n),d,m,0);
      q(i) := Random_Sparse_Poly(natural32(n),d,m,0);
    end loop;
    put_line("-> p = "); put(p);
    put_line("-> q = "); put(q);
  end DoblDobl_Random_Systems;

  procedure DoblDobl_Random_Coefficient_Systems
              ( n : in integer32;
                p,q : out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    m,d : natural32 := 0;

  begin
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    for i in p'range loop
      p(i) := Random_Sparse_Poly(natural32(n),d,m,0);
      q(i) := DoblDobl_Complex_Poly_Randomizers.Complex_Randomize1(p(i));
    end loop;
    put_line("-> p = "); put(p);
    put_line("-> q = "); put(q);
  end DoblDobl_Random_Coefficient_Systems;

  procedure DoblDobl_Compared_Encapsulation_Test is

    use DoblDobl_Complex_Numbers;

    n : natural32 := 0;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
      gamma : constant Complex_Number := DoblDobl_Random_Numbers.Random1;
    begin
      DoblDobl_Random_Systems(integer32(n),p,q);
      DoblDobl_Coefficient_Homotopy.Create(p,q,2,gamma);
      DoblDobl_Homotopy.Create(q,p,2,gamma);
      DoblDobl_Compared_Encapsulated_Eval(n);
    end;
  end DoblDobl_Compared_Encapsulation_Test;

  procedure DoblDobl_Homotopy_Performance ( n,m : natural32 ) is

    use DoblDobl_Complex_Numbers;

    t : double_double;
    ct : Complex_Number;
    x,y : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    A : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(DoblDobl_Random_Numbers.Random);
      ct := Create(t);
      x := DoblDobl_Random_Vectors.Random_Vector(1,integer32(n));
      y := DoblDobl_Homotopy.Eval(x,ct);
      A := DoblDobl_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of dobldobl homotopy");
  end DoblDobl_Homotopy_Performance;

  procedure DoblDobl_Coefficient_Homotopy_Performance ( n,m : natural32 ) is

    use DoblDobl_Complex_Numbers;

    t : double_double;
    ct : Complex_Number;
    x,y : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    A : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(DoblDobl_Random_Numbers.Random);
      ct := Create(t);
      x := DoblDobl_Random_Vectors.Random_Vector(1,integer32(n));
      y := DoblDobl_Coefficient_Homotopy.Eval(x,ct);
      A := DoblDobl_Coefficient_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of coefficient homotopy");
  end DoblDobl_Coefficient_Homotopy_Performance;

  procedure DoblDobl_Performance_Test is

    use DoblDobl_Complex_Numbers;

    n,m : natural32 := 0;

  begin
    put("Give the number of evaluations : "); get(m);
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
      gamma : constant Complex_Number := DoblDobl_Random_Numbers.Random1;
    begin
      DoblDobl_Random_Coefficient_Systems(integer32(n),p,q);
      DoblDobl_Coefficient_Homotopy.Create(p,q,2,gamma);
      DoblDobl_Homotopy.Create(q,p,2,gamma);
    end;
    DoblDobl_Homotopy_Performance(n,m);
    DoblDobl_Coefficient_Homotopy_Performance(n,m);
  end DoblDobl_Performance_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Evaluating a homotopy pair of polynomials ...");
    put_line("  1. compare for double double random systems;");
    put_line("  2. performance of double double homotopies.");
    put("Type 1 or 2 to choose a test : ");
    Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => DoblDobl_Compared_Encapsulation_Test;
      when '2' => DoblDobl_Performance_Test;
      when others => null;
    end case;
  end Main;

end Test_DoblDobl_Coeff_Homotopy;
