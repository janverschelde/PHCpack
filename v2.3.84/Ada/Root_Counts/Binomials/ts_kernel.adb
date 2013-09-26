with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Integer_Kernel;
with Standard_Integer64_Kernel;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
with Multprec_Random_Matrices;          use Multprec_Random_Matrices;
with Multprec_Integer_Kernel;

procedure ts_kernel is

-- DESCRIPTION :
--   Interactive development and testing of the calculation of a
--   basis for the kernel of an integer matrix.

  function Is_Zero ( A : Standard_Integer_Matrices.Matrix ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the matrix A equals the zero matrix.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if A(i,j) /= 0
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( A : Standard_Integer64_Matrices.Matrix )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if the matrix A equals the zero matrix.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if A(i,j) /= 0
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( A : Multprec_Integer_Matrices.Matrix )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if the matrix A equals the zero matrix.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if not Equal(A(i,j),0)
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Zero;

  procedure Standard_Compute_and_Test_Kernel
              ( A : in Standard_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Computes the kernel with A and shows the results if output is on.

    r : integer32;
    V : Standard_Integer_Matrices.Link_to_Matrix;

  begin
    Standard_Integer_Kernel.Kernel(A,r,V);
    if output then
      if r = A'last(2) then
        put_line("A matrix of full rank has no kernel.");
      else
        put("The rank of the matrix is "); put(r,1); put_line(".");
        put_line("The kernel : "); put(V.all);
      end if;
    end if;
    bug := false;
    if r < A'last(2) then
      declare
        use Standard_Integer_Matrices;
        Z : constant Matrix := A*V.all;
      begin
        if output
         then put_line("The residuals : "); put(Z);
        end if;
        put("Rank = "); put(r,1); 
        if Is_Zero(Z)
         then put_line(", all residuals are zero, okay.");
         else put_line(", not all residuals are zero, bug!"); bug := true;
        end if;
      end;
    end if;
  end Standard_Compute_and_Test_Kernel;

  procedure Standard_Compute64_and_Test_Kernel
              ( A : in Standard_Integer64_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Computes the kernel with A and shows the results if output is on.

    r : integer32;
    V : Standard_Integer64_Matrices.Link_to_Matrix;

  begin
    Standard_Integer64_Kernel.Kernel(A,r,V);
    if output then
      if r = A'last(2) then
        put_line("A matrix of full rank has no kernel.");
      else
        put("The rank of the matrix is "); put(r,1); put_line(".");
        put_line("The kernel : "); put(V.all);
      end if;
    end if;
    bug := false;
    if r < A'last(2) then
      declare
        use Standard_Integer64_Matrices;
        Z : constant Matrix := A*V.all;
      begin
        if output
         then put_line("The residuals : "); put(Z);
        end if;
        put("Rank = "); put(r,1); 
        if Is_Zero(Z)
         then put_line(", all residuals are zero, okay.");
         else put_line(", not all residuals are zero, bug!"); bug := true;
        end if;
      end;
    end if;
  end Standard_Compute64_and_Test_Kernel;

  procedure Multprec_Compute_and_Test_Kernel
              ( A : in Multprec_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   Computes the kernel with A and shows the results if output is on.

    r : integer32;
    V : Multprec_Integer_Matrices.Link_to_Matrix;

  begin
    Multprec_Integer_Kernel.Kernel(A,r,V);
    if output then
      if r = A'last(2) then
        put_line("A matrix of full rank has no kernel.");
      else
        put("The rank of the matrix is "); put(r,1); put_line(".");
        put_line("The kernel : "); put(V.all);
      end if;
    end if;
    bug := false;
    if r < A'last(2) then
      declare
        use Multprec_Integer_Matrices;
        Z : Matrix(A'range(1),V'range(2)) := A*V.all;
      begin
        if output
         then put_line("The residuals : "); put(Z);
        end if;
        put("Rank = "); put(r,1); 
        if Is_Zero(Z)
         then put_line(", all residuals are zero, okay.");
         else put_line(", not all residuals are zero, bug!"); bug := true;
        end if;
        Multprec_Integer_Matrices.Clear(Z);
      end;
    end if;
    Multprec_Integer_Matrices.Clear(V);
  end Multprec_Compute_and_Test_Kernel;

  procedure Standard_Compute_Kernel ( n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for or generates a k-by-n matrix of exponents
  --   and computes its kernel.

  -- ON ENTRY :
  --   n          the number of variables or the ambient dimension;
  --   k          the number of equations or the expected codimension.

    A : Standard_Integer_Matrices.Matrix(1..k,1..n); 
    A64 : Standard_Integer64_Matrices.Matrix(1..k,1..n); 
    m : integer32 := 0;
    lower,upper : integer32 := 0;
    lower64,upper64 : integer64 := 0;
    ans : character;
    ari64,output,bug : boolean;

  begin
    put("How many random tests ? (0 for interactive) "); get(m);
    put("Use 64-bit arithmetic ? (y/n) "); Ask_Yes_or_No(ans);
    ari64 := (ans = 'y');
    if m > 0 then
      if ari64 then
        put("  give lower bound for exponents : "); get(lower64);
        put("  give upper bound for exponents : "); get(upper64);
      else
        put("  give lower bound for exponents : "); get(lower);
        put("  give upper bound for exponents : "); get(upper);
      end if;
      put("Do you want intermediate output ? (y/n) ");
      Ask_Yes_or_No(ans); output := (ans = 'y');
      for i in 1..m loop
        if ari64
         then A64 := Random_Matrix(natural32(k),natural32(n),lower64,upper64);
         else A := Random_Matrix(natural32(k),natural32(n),lower,upper);
        end if;
        put("-> a random "); put(k,1); put("-by-"); put(n,1);
        put_line(" matrix :");
        if ari64
         then put(A64); Standard_Compute64_and_Test_Kernel(A64,output,bug);
         else put(A);   Standard_Compute_and_Test_Kernel(A,output,bug);
        end if;
        exit when bug;
      end loop;
      if not bug then
        put("Tested "); put(m,1);
        put_line(" random cases successfully.");
      end if;
    else
      put("-> give a "); put(k,1); put("-by-"); put(n,1);
      put_line(" matrix "); 
      if ari64
       then get(A64); Standard_Compute64_and_Test_Kernel(A64,true,bug);
       else get(A);   Standard_Compute_and_Test_Kernel(A,true,bug);
      end if;
    end if;
  end Standard_Compute_Kernel;

  procedure Multprec_Compute_Kernel ( n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for or generates a k-by-n matrix of exponents
  --   and computes its kernel.

  -- ON ENTRY :
  --   n          the number of variables or the ambient dimension;
  --   k          the number of equations or the expected codimension.

    A : Multprec_Integer_Matrices.Matrix(1..k,1..n); 
    m : integer32 := 0;
    lower,upper : integer32 := 0;
    ans : character;
    output,bug : boolean;

  begin
    put("How many random tests ? (0 for interactive) "); get(m);
    if m > 0 then
      put("  give lower bound for exponents : "); get(lower);
      put("  give upper bound for exponents : "); get(upper);
      put("Do you want intermediate output ? (y/n) ");
      Ask_Yes_or_No(ans); output := (ans = 'y');
      for i in 1..m loop
        A := Random_Matrix(natural32(k),natural32(n),lower,upper);
        put("-> a random "); put(k,1); put("-by-"); put(n,1);
        put_line(" matrix :"); put(A);
        Multprec_Compute_and_Test_Kernel(A,output,bug);
        exit when bug;
      end loop;
      if not bug then
        put("Tested "); put(m,1);
        put_line(" random cases successfully.");
      end if;
    else
      put("-> give a "); put(k,1); put("-by-"); put(n,1);
      put_line(" matrix "); get(A);
      Multprec_Compute_and_Test_Kernel(A,true,bug);
    end if;
  end Multprec_Compute_Kernel;

  procedure Main is

    nq,nv : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Computing the kernel of an exponent matrix...");
    put("  give the number of equations : "); get(nq);
    put("  give the number of variables : "); get(nv);
    put("Use multiprecision arithmetic ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Compute_Kernel(nv,nq);
     else Standard_Compute_Kernel(nv,nq);
    end if;
  end Main;

begin
  Main;
end ts_kernel;
