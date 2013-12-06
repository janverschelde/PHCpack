with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;       use Standard_Integer_VecVecs_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;      use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;       use Lists_of_Integer_Vectors_io;
with Symbol_Table;
with Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Random_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Lexicographical_Supports;
with Standard_Polynomial_Flatteners;
with DoblDobl_Polynomial_Flatteners;
with QuadDobl_Polynomial_Flatteners;

procedure ts_polyflat is

-- DESCRIPTION :
--   Interactive development of flattening a polynomial system.
--   A flattened representation of a polynomial system consists
--   of a coefficient matrix and a monomial vector.
--   The monomial vector is defined by a vector of exponents.
--   Columns of the coefficient matrix are indexed by the exponents
--   and the rows contain the coefficients of the polynomials in the system.

  procedure Prompt_for_Random_Parameters ( n,d,m : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension n, largest degree d,
  --   and number of monomials m (0 for dense) as the parameters
  --   to generate a random polynomial system.
  --   The symbol table is initialized for n variables.

  begin
    new_line;
    put_line("Generating a random polynomial system ...");
    n := 0; put("Give the dimension : "); get(n);
    Symbol_Table.Init(n);
    d := 0; put("Give the maximum degree : "); get(d);
    m := 0; put("Give the number of monomials (0 for dense) : "); get(m);
  end Prompt_for_Random_Parameters;

  procedure Standard_Random_System
              ( p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                ntp : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.
  --   Returns in ntp the number of terms per polynomial.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..n loop
        if m = 0
         then q(integer32(i))
                := Standard_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(integer32(i))
                := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new Standard_Complex_Poly_Systems.Poly_Sys'(q);
    end;
    ntp := m;
  end Standard_Random_System;

  procedure DoblDobl_Random_System
              ( p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                ntp : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.
  --   Returns in ntp the number of terms per polynomial.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..integer32(n) loop
        if m = 0
         then q(i) := DoblDobl_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(i) := DoblDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(q);
    end;
    ntp := m;
  end DoblDobl_Random_System;

  procedure QuadDobl_Random_System
              ( p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                ntp : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.
  --   Returns in ntp the number of terms per polynomial.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..integer32(n) loop
        if m = 0
         then q(i) := QuadDobl_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(i) := QuadDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(q);
    end;
    ntp := m;
  end QuadDobl_Random_System;

  procedure Dense_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant Standard_Complex_Matrices.Matrix(p'range,v'range)
      := Standard_Polynomial_Flatteners.Coefficient_Matrix(p,v);
    nz : natural32;
    f : double_float;
    

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    Standard_Polynomial_Flatteners.Spy(A,v,nz);
    put("Number of nonzero elements : "); put(nz,1); 
    f := double_float(nz)/double_float(A'last(1)*A'last(2));
    put(" sparse factor : "); put(f,3); new_line;
    Standard_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Sparse_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the flattened sparse representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(s);
    cfv : constant Standard_Integer_VecVecs.VecVec
        := Lexicographical_Supports.Compress(fv);
    nd : constant Lists_of_Integer_Vectors.List 
       := Lexicographical_Supports.Nodes(v);
    nd_cnt : constant natural32 := Lists_of_Integer_Vectors.Length_Of(nd);

  begin
    put("First positive : ");
    put(Lexicographical_Supports.First_Positive(s),1); new_line;
    put_line("The factored supports : "); put(fv);
    put("The list of "); put(nd_cnt,1);
    put_line(" extra nodes : "); put(nd);
    Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    Standard_Polynomial_Flatteners.Test_Eval(p,c,v,fv,cfv,k);
  end Sparse_Coefficients;

  procedure Dense_Coefficients
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant DoblDobl_Complex_Matrices.Matrix(p'range,v'range)
      := DoblDobl_Polynomial_Flatteners.Coefficient_Matrix(p,v);

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    DoblDobl_Polynomial_Flatteners.Spy(A,v);
    DoblDobl_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Sparse_Coefficients
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the flattened sparse representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    c : DoblDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(s);
    cfv : constant Standard_Integer_VecVecs.VecVec
        := Lexicographical_Supports.Compress(fv);

  begin
    DoblDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    DoblDobl_Polynomial_Flatteners.Test_Eval(p,c,v,fv,cfv,k);
  end Sparse_Coefficients;

  procedure Dense_Coefficients
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant QuadDobl_Complex_Matrices.Matrix(p'range,v'range)
      := QuadDobl_Polynomial_Flatteners.Coefficient_Matrix(p,v);

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    QuadDobl_Polynomial_Flatteners.Spy(A,v);
    QuadDobl_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Sparse_Coefficients
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the flattened sparse representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(s);
    cfv : constant Standard_Integer_VecVecs.VecVec
        := Lexicographical_Supports.Compress(fv);

  begin
    QuadDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    QuadDobl_Polynomial_Flatteners.Test_Eval(p,c,v,fv,cfv,k);
  end Sparse_Coefficients;

  procedure Dense_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant Standard_Complex_Matrices.Matrix(p'range,v'range)
      := Standard_Polynomial_Flatteners.Coefficient_Matrix(p,v);
    nz : natural32;
    f : double_float;

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    Standard_Polynomial_Flatteners.Spy(A,v,nz);
    put("Number of nonzero elements : "); put(nz,1);
    f := double_float(nz)/double_float(A'last(1)*A'last(2));
    put(" sparse factor : "); put(f,3); new_line;
    Standard_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Dense_Coefficients
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant DoblDobl_Complex_Matrices.Matrix(p'range,v'range)
      := DoblDobl_Polynomial_Flatteners.Coefficient_Matrix(p,v);

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    DoblDobl_Polynomial_Flatteners.Spy(A,v);
    DoblDobl_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Dense_Coefficients
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(s);
    A : constant QuadDobl_Complex_Matrices.Matrix(p'range,v'range)
      := QuadDobl_Polynomial_Flatteners.Coefficient_Matrix(p,v);

  begin
    put_line("The coefficient matrix : "); put(A,3);
    put_line("The structure of the coefficient matrix : ");
    QuadDobl_Polynomial_Flatteners.Spy(A,v);
    QuadDobl_Polynomial_Flatteners.Test_Eval(p,A,v);
  end Dense_Coefficients;

  procedure Standard_Poly_Flatten_Test is

    use Standard_Polynomial_Flatteners;

    ans : character;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nt,m : natural32;
    sp,ssp : List;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Standard_Random_System(lp,m);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
    end if;
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    new_line;
    put("Is the system sparse ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'n'
     then m := 0;
     else m := nt;
    end if;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    if m = 0
     then Dense_Coefficients(lp.all,ssp);
     else Sparse_Coefficients(lp.all,ssp);
    end if;
  end Standard_Poly_Flatten_Test;

  procedure DoblDobl_Poly_Flatten_Test is

    use DoblDobl_Polynomial_Flatteners;

    ans : character;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nt,m : natural32;
    sp,ssp : List;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      DoblDobl_Random_System(lp,m);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
    end if;
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    new_line;
    put("Is the system sparse ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'n'
     then m := 0;
     else m := nt;
    end if;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    if m = 0
     then Dense_Coefficients(lp.all,ssp);
     else Sparse_Coefficients(lp.all,ssp);
    end if;
  end DoblDobl_Poly_Flatten_Test;

  procedure QuadDobl_Poly_Flatten_Test is

    use QuadDobl_Polynomial_Flatteners;

    ans : character;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nt,m : natural32;
    sp,ssp : List;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      QuadDobl_Random_System(lp,m);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
    end if;
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    new_line;
    put("Is the system sparse ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'n'
     then m := 0;
     else m := nt;
    end if;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    if m = 0
     then Dense_Coefficients(lp.all,ssp);
     else Sparse_Coefficients(lp.all,ssp);
    end if;
  end QuadDobl_Poly_Flatten_Test;

  procedure Standard_Laur_Flatten_Test is

    use Standard_Polynomial_Flatteners;

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    nt : natural32;
    sp,ssp : List;

  begin
    new_line;
    put_line("reading a Laurent polynomial system ...");
    get(lp);
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    Dense_Coefficients(lp.all,ssp);
  end Standard_Laur_Flatten_Test;

  procedure DoblDobl_Laur_Flatten_Test is

    use DoblDobl_Polynomial_Flatteners;

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nt : natural32;
    sp,ssp : List;

  begin
    new_line;
    put_line("reading a Laurent polynomial system ...");
    get(lp);
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    Dense_Coefficients(lp.all,ssp);
  end DoblDobl_Laur_Flatten_Test;

  procedure QuadDobl_Laur_Flatten_Test is

    use QuadDobl_Polynomial_Flatteners;

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nt : natural32;
    sp,ssp : List;

  begin
    new_line;
    put_line("reading a Laurent polynomial system ...");
    get(lp);
    nt := Number_of_Terms(lp.all);
    put("Sum of #terms : "); put(nt,1); new_line;
    sp := Distinct_Supports(lp.all);
    put("#distinct monomials : "); put(Length_Of(sp),1); new_line;
    put_line("The list of distinct monomials : "); put(sp);
    ssp := Lexicographical_Supports.Sort(sp);
    put_line("After lexicographical sort : "); put(ssp);
    Dense_Coefficients(lp.all,ssp);
  end QuadDobl_Laur_Flatten_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test flatteners of polynomial systems :");
    put_line("  1. flatten standard complex polynomial system;");
    put_line("  2. flatten standard complex Laurent polynomial system;");
    put_line("  3. flatten double double complex polynomial system;");
    put_line("  4. flatten double double complex Laurent polynomial system;");
    put_line("  5. flatten quad double complex polynomial system;");
    put_line("  6. flatten quad double complex Laurent polynomial system.");
    put("Type 1, 2, 3, 4, 5, or 6 to make your choice : ");
    Ask_Alternative(ans,"123456");
    case ans is 
      when '1' => Standard_Poly_Flatten_Test;
      when '2' => Standard_Laur_Flatten_Test;
      when '3' => DoblDobl_Poly_Flatten_Test;
      when '4' => DoblDobl_Laur_Flatten_Test;
      when '5' => QuadDobl_Poly_Flatten_Test;
      when '6' => QuadDobl_Laur_Flatten_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_polyflat;
