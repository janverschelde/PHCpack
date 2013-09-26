with text_io;                            use text_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Fewnomial_Systems;         use Standard_Fewnomial_Systems;
with Standard_Fewnomial_Solvers;         use Standard_Fewnomial_Solvers;

package body Standard_Sparse_Solvers is

-- AUXILIARY ROUTINES :

  function Is_Fewnomial_System
            ( p : Poly_Sys; nq,nv : natural ) return boolean is

  -- DESCRIPTION :
  --   Parses the system p and returns true if it is fewnomial.

  -- ON ENTRY :
  --   p      polynomial system;
  --   nq     number of equations in p;
  --   nv     number of variables in p.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    fail : boolean;

  begin
    Parse(p,nv,A,C,b,fail);
    return not fail;
  end Is_Fewnomial_System;

  function Is_Fewnomial_System
            ( p : Laur_Sys; nq,nv : natural ) return boolean is

  -- DESCRIPTION :
  --   Parses the system p and returns true if it is fewnomial.

  -- ON ENTRY :
  --   p      Laurent polynomial system;
  --   nq     number of equations in p;
  --   nv     number of variables in p.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    fail : boolean;

  begin
    Parse(p,nv,A,C,b,fail);
    return not fail;
  end Is_Fewnomial_System;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in natural; tol : in double_float;
                sols : out Solution_List; fail,zero_y : out boolean ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance for deciding when a number is zero.

  -- ON RETURN :
  --   sols     solution of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if no solutions with all components unequal zero.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;
    info : integer;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,info,sols);
    end if;
  exception
    when others => put_line("exception raised in Parse_and_Solve...");
                   raise;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in natural; tol : in double_float;
                rcond : out double_float; sols : out Solution_List;
                fail,zero_y : out boolean ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance used to decide whether a number is zero.

  -- ON RETURN :
  --   rcond    estimate for inverse of the condition number of
  --            the coefficients matrix of p, only if fewnomial;
  --   sols     solutions of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if p has no solutions with all components unequal zero.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,rcond,sols);
    end if;
  exception
    when others => put_line("exception raised in Parse_and_Solve...");
                   put("rcond = "); put(rcond,3); new_line;
                   raise;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Poly_Sys; nq,nv : in natural; tol : in double_float;
                rcond : out double_float; sols : out Solution_List;
                fail,zero_y : out boolean; rsum : out double_float ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero.

  -- ON RETURN :
  --   sols     solution of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if p has no solutions with all components unequal zero;
  --   rsum     sum of the residuals.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,rcond,sols);
          rsum := Sum_Residuals(A,C,b,sols);
     else rsum := 1.0;
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in natural; tol : in double_float;
                sols : out Solution_List; fail,zero_y : out boolean ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero.

  -- ON RETURN :
  --   sols     solution of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if p has no solution with all components unequal zero.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;
    info : integer;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,info,sols);
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in natural; tol : in double_float;
                rcond : out double_float; sols : out Solution_List;
                fail,zero_y : out boolean ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero.

  -- ON RETURN :
  --   rcond    estimate for the inverse of the condition number of
  --            the coefficient matrix of p, only if it is fewnomial;
  --   sols     solution of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if p has no solutions with all components unequal zero.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,rcond,sols);
    end if;
  end Parse_and_Solve;

  procedure Parse_and_Solve
              ( p : in Laur_Sys; nq,nv : in natural; tol : in double_float;
                rcond : out double_float; sols : out Solution_List;
                fail,zero_y : out boolean; rsum : out double_float ) is

  -- DESCRIPTION :
  --   Parses the system and returns its solutions if it is fewnomial,
  --   otherwise fail equals true on return.

  -- ON ENTRY :
  --   p        Laurent polynomial system;
  --   nq       number of equations in p;
  --   nv       number of variables in p;
  --   tol      tolerance to decide whether a number is zero.

  -- ON RETURN :
  --   sols     solution of p if it is fewnomial;
  --   fail     true if p is not fewnomial;
  --   zero_y   true if p has no solutions with all components unequal zero;
  --   rsum     sum of the residuals.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : natural;

  begin
    Parse(p,nv,A,C,b,fail);
    if not fail
     then Solve(A,C,b,tol,zero_y,r,rcond,sols);
          rsum := Sum_Residuals(A,C,b,sols);
     else rsum := 1.0;
    end if;
  end Parse_and_Solve;

-- TARGET ROUTINES :

  function Is_Fewnomial_System ( p : Poly_Sys ) return boolean is

    use Standard_Complex_Polynomials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    return Is_Fewnomial_System(p,nq,nv);
  end Is_Fewnomial_System;

  function Is_Fewnomial_System ( p : Laur_Sys ) return boolean is

    use Standard_Complex_Laurentials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    return Is_Fewnomial_System(p,nq,nv);
  end Is_Fewnomial_System;

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    sols : out Solution_List;
                    fail,zero_y : out boolean ) is

    use Standard_Complex_Polynomials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,sols,fail,zero_y);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    sols : out Solution_List;
                    fail,zero_y : out boolean ) is

    use Standard_Complex_Laurentials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,sols,fail,zero_y);
  end Solve;

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean ) is

    use Standard_Complex_Polynomials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    rcond : out double_float;
                    sols : out Solution_List;
                    fail,zero_y : out boolean ) is

    use Standard_Complex_Laurentials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y);
  end Solve;

  procedure Solve ( p : in Poly_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_float ) is

    use Standard_Complex_Polynomials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,rsum);
  end Solve;

  procedure Solve ( p : in Laur_Sys; tol : in double_float;
                    rcond : out double_float; sols : out Solution_List;
                    fail,zero_y : out boolean; rsum : out double_float ) is

    use Standard_Complex_Laurentials;

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));

  begin
    Parse_and_Solve(p,nq,nv,tol,rcond,sols,fail,zero_y,rsum);
  end Solve;

end Standard_Sparse_Solvers;
