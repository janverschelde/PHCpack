with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;       use OctoDobl_Complex_Numbers_io;
with OctoDobl_Random_Numbers;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors_io;       use OctoDobl_Complex_Vectors_io;
with OctoDobl_Random_Vectors;
with Symbol_Table;
with OctoDobl_Random_Polynomials;       use OctoDobl_Random_Polynomials;
with OctoDobl_Complex_Poly_Randomizers;
with OctoDobl_Complex_Poly_Systems_io;  use OctoDobl_Complex_Poly_Systems_io;
with OctoDobl_Homotopy;
with OctoDobl_Coefficient_Homotopy;

package body Test_OctoDobl_Coeff_Homotopy is

  procedure Write_Elements ( A,B : in OctoDobl_Complex_Matrices.Matrix ) is
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

  procedure OctoDobl_Compared_Encapsulated_Eval ( n : in natural32 ) is

    use OctoDobl_Complex_Numbers;

    t : constant octo_double := abs(OctoDobl_Random_Numbers.Random);
    ct : constant Complex_Number := Create(t);
    x : constant OctoDobl_Complex_Vectors.Vector(1..integer32(n))
      := OctoDobl_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant OctoDobl_Complex_Vectors.Vector
      := OctoDobl_Homotopy.Eval(x,ct);
    z : constant OctoDobl_Complex_Vectors.Vector
      := OctoDobl_Coefficient_Homotopy.Eval(x,ct);
    A : constant OctoDobl_Complex_Matrices.Matrix(x'range,x'range)
      := OctoDobl_Homotopy.Diff(x,ct);
    B : constant OctoDobl_Complex_Matrices.Matrix(x'range,x'range)
      := OctoDobl_Coefficient_Homotopy.Diff(x,ct);

  begin
    put("A random t : "); put(ct); new_line;
    put_line("A random point : "); put_line(x);
    put_line("-> y = "); put_line(y);
    put_line("-> z = "); put_line(z);
    Write_Elements(A,B);
  end OctoDobl_Compared_Encapsulated_Eval;

  procedure OctoDobl_Random_Systems
              ( n : in integer32;
                p,q : out OctoDobl_Complex_Poly_Systems.Poly_Sys ) is

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
  end OctoDobl_Random_Systems;

  procedure OctoDobl_Random_Coefficient_Systems
              ( n : in integer32;
                p,q : out OctoDobl_Complex_Poly_Systems.Poly_Sys ) is

    m,d : natural32 := 0;

  begin
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    for i in p'range loop
      p(i) := Random_Sparse_Poly(natural32(n),d,m,0);
      q(i) := OctoDobl_Complex_Poly_Randomizers.Complex_Randomize1(p(i));
    end loop;
    put_line("-> p = "); put(p);
    put_line("-> q = "); put(q);
  end OctoDobl_Random_Coefficient_Systems;

  procedure OctoDobl_Compared_Encapsulation_Test is

    use OctoDobl_Complex_Numbers;

    n : natural32 := 0;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
      gamma : constant Complex_Number := OctoDobl_Random_Numbers.Random1;
    begin
      OctoDobl_Random_Systems(integer32(n),p,q);
      OctoDobl_Coefficient_Homotopy.Create(p,q,2,gamma);
      OctoDobl_Homotopy.Create(q,p,2,gamma);
      OctoDobl_Compared_Encapsulated_Eval(n);
    end;
  end OctoDobl_Compared_Encapsulation_Test;

  procedure OctoDobl_Homotopy_Performance ( n,m : natural32 ) is

    use OctoDobl_Complex_Numbers;

    t : octo_double;
    ct : Complex_Number;
    x,y : OctoDobl_Complex_Vectors.Vector(1..integer32(n));
    A : OctoDobl_Complex_Matrices.Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(OctoDobl_Random_Numbers.Random);
      ct := Create(t);
      x := OctoDobl_Random_Vectors.Random_Vector(1,integer32(n));
      y := OctoDobl_Homotopy.Eval(x,ct);
      A := OctoDobl_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of OctoDobl homotopy");
  end OctoDobl_Homotopy_Performance;

  procedure OctoDobl_Coefficient_Homotopy_Performance ( n,m : natural32 ) is

    use OctoDobl_Complex_Numbers;

    t : octo_double;
    ct : Complex_Number;
    x,y : OctoDobl_Complex_Vectors.Vector(1..integer32(n));
    A : OctoDobl_Complex_Matrices.Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(OctoDobl_Random_Numbers.Random);
      ct := Create(t);
      x := OctoDobl_Random_Vectors.Random_Vector(1,integer32(n));
      y := OctoDobl_Coefficient_Homotopy.Eval(x,ct);
      A := OctoDobl_Coefficient_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of coefficient homotopy");
  end OctoDobl_Coefficient_Homotopy_Performance;

  procedure OctoDobl_Performance_Test is

    use OctoDobl_Complex_Numbers;

    n,m : natural32 := 0;

  begin
    put("Give the number of evaluations : "); get(m);
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
      gamma : constant Complex_Number := OctoDobl_Random_Numbers.Random1;
    begin
      OctoDobl_Random_Coefficient_Systems(integer32(n),p,q);
      OctoDobl_Coefficient_Homotopy.Create(p,q,2,gamma);
      OctoDobl_Homotopy.Create(q,p,2,gamma);
    end;
    OctoDobl_Homotopy_Performance(n,m);
    OctoDobl_Coefficient_Homotopy_Performance(n,m);
  end OctoDobl_Performance_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Evaluating a homotopy pair of polynomials ...");
    put_line("  1. compare for octo double random systems;");
    put_line("  2. performance of octo double homotopies.");
    put("Type 1 or 2 to choose a test : ");
    Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => OctoDobl_Compared_Encapsulation_Test;
      when '2' => OctoDobl_Performance_Test;
      when others => null;
    end case;
  end Main;

end Test_OctoDobl_Coeff_Homotopy;
