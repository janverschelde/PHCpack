with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Polynomial_Convertors;    use DoblDobl_Polynomial_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Polynomial_Convertors;    use QuadDobl_Polynomial_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Jaco_Matrices;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;       use Lists_of_Integer_Vectors_io;
with Lexicographical_Supports;
with Standard_Polynomial_Flatteners;
with DoblDobl_Polynomial_Flatteners;
with QuadDobl_Polynomial_Flatteners;
with Standard_Jacobian_Evaluations;
with DoblDobl_Jacobian_Evaluations;
with QuadDobl_Jacobian_Evaluations;

procedure ts_speelsys is

-- DESCRIPTION :
--   Development of the efficient evaluation of a system and its matrix
--   of all partial derivatives using the Speelpenning example.
--   The system is assumed to be sparse but all polynomials may share
--   the same support, or share some monomials.

  function Difference_between_Rows
              ( A,B : Standard_Complex_Matrices.Matrix )
              return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of the same number of rows as A and B.
  --   The i-th entry of the vector on returns is the sum of the
  --   differences between the elements on the i-th rows of A and B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2);

    res : Standard_Complex_Vectors.Vector(A'range(1));
    use Standard_Complex_Numbers;

  begin
    for i in A'range(1) loop
      res(i) := Create(0.0);
      for j in A'range(2) loop
        res(i) := res(i) + A(i,j) - B(i,j);
      end loop;
    end loop;
    return res;
  end Difference_between_Rows;

  function Difference_between_Rows
              ( A,B : DoblDobl_Complex_Matrices.Matrix )
              return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of the same number of rows as A and B.
  --   The i-th entry of the vector on returns is the sum of the
  --   differences between the elements on the i-th rows of A and B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2);

    res : DoblDobl_Complex_Vectors.Vector(A'range(1));
    use DoblDobl_Complex_Numbers;

  begin
    for i in A'range(1) loop
      res(i) := Create(integer(0));
      for j in A'range(2) loop
        res(i) := res(i) + A(i,j) - B(i,j);
      end loop;
    end loop;
    return res;
  end Difference_between_Rows;

  function Difference_between_Rows
              ( A,B : QuadDobl_Complex_Matrices.Matrix )
              return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of the same number of rows as A and B.
  --   The i-th entry of the vector on returns is the sum of the
  --   differences between the elements on the i-th rows of A and B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2);

    res : QuadDobl_Complex_Vectors.Vector(A'range(1));
    use QuadDobl_Complex_Numbers;

  begin
    for i in A'range(1) loop
      res(i) := Create(integer(0));
      for j in A'range(2) loop
        res(i) := res(i) + A(i,j) - B(i,j);
      end loop;
    end loop;
    return res;
  end Difference_between_Rows;

  procedure Test_Standard_Jacobian_Evaluation
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Evaluates p and its Jacobian matrix where e contains the
  --   lexicographically sorted list of different exponent vectors of p.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    y : constant Standard_Complex_Vectors.Vector(p'range)
      := Standard_Complex_Poly_SysFun.Eval(p,x);
    z : Standard_Complex_Vectors.Vector(p'range);
    A,B : Standard_Complex_Matrices.Matrix(p'range,1..n);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(A'range(1),A'range(2))
       := Standard_Complex_Jaco_Matrices.Create(p);
    d : Standard_Complex_Vectors.Vector(A'range(1));

    use Standard_Jacobian_Evaluations;

  begin
    A := Standard_Complex_Jaco_Matrices.Eval(jm,x);
    Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    Standard_Jacobian_Evaluation(v,c,k,x,z,B);
    put_line("the first evaluation : "); put_line(y);
    put_line("the second evaluation : "); put_line(z);
    d := Difference_between_Rows(A,B);
    put_line("difference between Jacobian evaluations : ");
    put_line(d);
    Standard_Complex_Jaco_Matrices.Clear(jm);
  end Test_Standard_Jacobian_Evaluation;

  procedure Test_DoblDobl_Jacobian_Evaluation
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Evaluates p and its Jacobian matrix where e contains the
  --   lexicographically sorted list of different exponent vectors of p,
  --   using double double arithmetic.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    c : DoblDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : constant DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    y : constant DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
    z : DoblDobl_Complex_Vectors.Vector(p'range);
    A,B : DoblDobl_Complex_Matrices.Matrix(p'range,1..n);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(A'range(1),A'range(2))
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    d : DoblDobl_Complex_Vectors.Vector(A'range(1));

    use DoblDobl_Jacobian_Evaluations;

  begin
    A := DoblDobl_Complex_Jaco_Matrices.Eval(jm,x);
    DoblDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    DoblDobl_Jacobian_Evaluation(v,c,k,x,z,B);
    put_line("the first evaluation : "); put_line(y);
    put_line("the second evaluation : "); put_line(z);
    d := Difference_between_Rows(A,B);
    put_line("difference between Jacobian evaluations : ");
    put_line(d);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
  end Test_DoblDobl_Jacobian_Evaluation;

  procedure Test_QuadDobl_Jacobian_Evaluation
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Evaluates p and its Jacobian matrix where e contains the
  --   lexicographically sorted list of different exponent vectors of p,
  --   using quad double arithmetic.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : constant QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    y : constant QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    z : QuadDobl_Complex_Vectors.Vector(p'range);
    A,B : QuadDobl_Complex_Matrices.Matrix(p'range,1..n);
    jm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(A'range(1),A'range(2))
       := QuadDobl_Complex_Jaco_Matrices.Create(p);
    d : QuadDobl_Complex_Vectors.Vector(A'range(1));

    use QuadDobl_Jacobian_Evaluations;

  begin
    A := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    QuadDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    QuadDobl_Jacobian_Evaluation(v,c,k,x,z,B);
    put_line("the first evaluation : "); put_line(y);
    put_line("the second evaluation : "); put_line(z);
    d := Difference_between_Rows(A,B);
    put_line("difference between Jacobian evaluations : ");
    put_line(d);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
  end Test_QuadDobl_Jacobian_Evaluation;

  procedure Main is

    use Standard_Complex_Poly_Systems;
    p : Link_to_Poly_Sys;
    s,e : Lists_of_Integer_Vectors.List;
    ans : character;

  begin
    put_line("Reading a polynomial system ..."); get(p);
    put_line("-> your system : "); put(p.all);
    s := Standard_Polynomial_Flatteners.Distinct_Supports(p.all);
    put("number of distinct monomials : ");
    put(Lists_of_Integer_Vectors.Length_Of(s),1); new_line;
    put_line("the distinct exponent vectors : "); put(s);
    e := Lexicographical_Supports.Sort(s);
    put_line("after a lexicographic sort : "); put(e);
    put("standard, double double, or quad double arithmetic ? (s/d/q) ");
    Ask_Alternative(ans,"sdq");
    case ans is
      when 's' => Test_Standard_Jacobian_Evaluation(p.all,e);
      when 'd' => 
        declare
          q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
            := Standard_Poly_Sys_to_DoblDobl_Complex(p.all);
        begin
          Test_DoblDobl_Jacobian_Evaluation(q,e);
          DoblDobl_Complex_Poly_Systems.Clear(q);
        end;
      when 'q' => 
        declare
          q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
            := Standard_Poly_Sys_to_QuadDobl_Complex(p.all);
        begin
          Test_QuadDobl_Jacobian_Evaluation(q,e);
          QuadDobl_Complex_Poly_Systems.Clear(q);
        end;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_speelsys;
