with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Blackbox_Continuations;    use Standard_Blackbox_Continuations;
with Planes_and_Polynomials;
with Witness_Sets_io;                    use Witness_Sets_io;

procedure ts_tropawit is

-- DESCRIPTION :
--   Computation of the tropicalization of a witness set,
--   looking for it asymptotics.

  function Special_Position ( n : integer32; p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Given a polynomial representation for a hyperplane in n-space,
  --   the polynomial on return has all coefficients of the hyperplane
   --  set to zero, except the first two.

    res : Poly;
    h : Standard_Complex_Vectors.Vector(0..n)
      := Planes_and_Polynomials.Polynomial(p);

  begin
   -- put_line("Coefficients of slice :");
   -- put_line(h);
    for i in 2..h'last loop
      h(i) := Create(0.0);
    end loop;
   -- put_line("Changed coefficients :");
   -- put_line(h);
    res := Planes_and_Polynomials.Hyperplane(h,1.0E-08);
   -- put("The polynomial hyperplane :");
   -- put_line(res);
    return res;
  end Special_Position;

  function Special_Target ( p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a system to fix the first coordinates for points
  --   on the curve defined by the embedded system p.
  --   Except for the last equation, all polynomials in the
  --   system on return are shared with the given system p.

    res : Poly_Sys(p'range);

  begin
    for i in p'first..p'last-1 loop
      res(i) := p(i);
    end loop;
    res(res'last) := Special_Position(p'last,p(p'last));
    return res;
  end Special_Target;

  procedure Move ( p : in Poly_Sys; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Moves the last hyperplane of p so that the first coordinate
  --   of the witness points are fixed.

    r : constant Poly_Sys(p'range) := Special_Target(p);
    rsols : Solution_List;
   -- one : constant Complex_Number := Create(1.0);
    pd : duration;

  begin
    put("The polynomial hyperplane on input :");
    put_line(p(p'last));
    put("The polynomial hyperplane after the cut:");
    put_line(r(r'last));
    Copy(sols,rsols);
    Set_Continuation_Parameter(rsols,Create(0.0));
    put_line("moving the hyperplane ...");
   -- Black_Box_Polynomial_Continuation(standard_output,r,p,one,rsols,pd);
    Black_Box_Polynomial_Continuation(standard_output,true,r,p,rsols,pd);
  end Move;

  procedure Main is

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : integer32 := 0;

  begin
    new_line;
    put_line("Tropicalization of a witness set...");
    Standard_Read_Embedding(p,sols,natural32(d));
    n := p'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Move(p.all,sols);
  end Main;

begin
  Main;
end ts_tropawit;
