with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Multprec_Integer_Numbers;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Simplex_Systems;           use DoblDobl_Simplex_Systems;
with DoblDobl_Simplex_Solvers;           use DoblDobl_Simplex_Solvers;

package body DoblDobl_Simpomial_Solvers is

-- AUXILIARY ROUTINES :

  function Is_Simplex_System
            ( p : Poly_Sys; nq,nv : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Parses the system p and returns true if it is a simplex system.

  -- ON ENTRY :
  --   p      polynomial system;
  --   nq     number of equations in p;
  --   nv     number of variables in p.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    fail : boolean;

  begin
    Parse(p,nv,A,C,b,fail);
    return not fail;
  end Is_Simplex_System;

  function Is_Simplex_System
            ( p : Laur_Sys; nq,nv : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Parses the system p and returns true if it is a simplex system.

  -- ON ENTRY :
  --   p      Laurent polynomial system;
  --   nq     number of equations in p;
  --   nv     number of variables in p.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    fail : boolean;

  begin
    Parse(p,nv,A,C,b,fail);
    return not fail;
  end Is_Simplex_System;

  function to_Multprec ( A : Standard_Integer64_Matrices.Matrix )
                       return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the matrix A converted to multiprecision format.

    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Multprec_Integer_Numbers.create(integer32(A(i,j)));
      end loop;
    end loop;
    return res;
  end to_Multprec;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in integer32; tol : in double_double;
                sols : out Solution_List; fail,zero_y : out boolean;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance for deciding when a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   sols     solution of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if no solutions with all components unequal zero.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;
    info : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail then
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,info,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,info,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  exception
    when others => put_line("exception raised in Parse_and_Solve...");
                   raise;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in integer32; tol : in double_double;
                rcond : out double_double; sols : out Solution_List;
                fail,zero_y : out boolean;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance used to decide whether a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   rcond    estimate for inverse of the condition number of
  --            the coefficients matrix of p, only if simplex;
  --   sols     solutions of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if p has no solutions with all components unequal zero.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail then
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,rcond,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,rcond,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  exception
    when others => put_line("exception raised in Parse_and_Solve...");
                   put("rcond = "); put(rcond,3); new_line;
                   raise;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in integer32; tol : in double_double;
                rcond : out double_double; sols : out Solution_List;
                fail,zero_y : out boolean; rsum : out double_double;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   sols     solution of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if p has no solutions with all components unequal zero;
  --   rsum     sum of the residuals.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if fail then
      rsum := Double_Double_Numbers.create(1.0);
    else
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,rcond,sols);
        rsum := Sum_Residuals(A,C,b,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,rcond,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in integer32; tol : in double_double;
                sols : out Solution_List; fail,zero_y : out boolean;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   sols     solution of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if p has no solution with all components unequal zero.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;
    info : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail then
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,info,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,info,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in integer32; tol : in double_double;
                rcond : out double_double; sols : out Solution_List;
                fail,zero_y : out boolean;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   rcond    estimate for the inverse of the condition number of
  --            the coefficient matrix of p, only if it is simplex;
  --   sols     solution of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if p has no solutions with all components unequal zero.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail then
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,rcond,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,rcond,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in integer32; tol : in double_double;
                rcond : out double_double; sols : out Solution_List;
                fail,zero_y : out boolean; rsum : out double_double;
                multprec_hermite : in boolean := false ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is simplex,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero;
  --   multprec_hermite indicates whether the Hermite normal form must
  --            be computed with multiprecision arithmetic.

  -- ON RETURN :
  --   sols     solution of p if it is a simplex system;
  --   fail     true if p is not a simplex system;
  --   zero_y   true if p has no solutions with all components unequal zero;
  --   rsum     sum of the residuals.

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;

  begin
    Parse(p,nv,A,C,b,fail);
    if fail then
      rsum := Double_Double_Numbers.create(1.0);
    else
      if not multprec_hermite then
        Solve(A,C,b,tol,zero_y,r,rcond,sols);
        rsum := Sum_Residuals(A,C,b,sols);
      else
        declare
          AA : Multprec_Integer_Matrices.Matrix(1..nv,1..nq) := to_Multprec(A);
        begin
          Solve(AA,C,b,tol,zero_y,r,rcond,sols);
          Multprec_Integer_Matrices.Clear(AA);
        end;
      end if;
    end if;
  end Parse_and_Solve;

-- TARGET ROUTINES :

  function Is_Simplex_System ( p : Poly_Sys ) return boolean is

    use DoblDobl_Complex_Polynomials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    return Is_Simplex_System(p,nq,nv);
  end Is_Simplex_System;

  function Is_Simplex_System ( p : Laur_Sys ) return boolean is

    use DoblDobl_Complex_Laurentials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    return Is_Simplex_System(p,nq,nv);
  end Is_Simplex_System;

  procedure Solve ( p : in Poly_Sys; tol : in double_double;
                    sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Polynomials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,sols,fail,zero_y,multprec_hermite);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_double;
                    sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Laurentials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,sols,fail,zero_y,multprec_hermite);
  end Solve;

  procedure Solve ( p : in Poly_Sys; tol : in double_double;
                    rcond : out double_double; sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Polynomials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,multprec_hermite);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_double;
                    rcond : out double_double;
                    sols : out Solution_List;
                    fail,zero_y : out boolean;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Laurentials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,multprec_hermite);
  end Solve;

  procedure Solve ( p : in Poly_Sys; tol : in double_double;
                    rcond : out double_double; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_double;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Polynomials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,rsum,multprec_hermite);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_double;
                    rcond : out double_double; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_double;
                    multprec_hermite : in boolean := false ) is

    use DoblDobl_Complex_Laurentials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,rsum,multprec_hermite);
  end Solve;

end DoblDobl_Simpomial_Solvers;
