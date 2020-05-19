with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;        use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;     use Standard_Complex_Poly_Functions;
with Evaluation_Differentiation_Errors;

procedure ts_perfhess is

-- DESCRIPTION :
--   Development of better algorithms to compute Hessians.

  function Symbolic
             ( x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Computes the Hessian matrix of the product of the variables in x
  --   symbolically, for testing purposes.

    res : Standard_Complex_Matrices.Matrix(x'range,x'range);
    p,dp,dp2 : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(x'range => 1);
    p := Create(t);
    for i in res'range(1) loop
      dp := Diff(p,i);
      for j in res'range(2) loop
        dp2 := Diff(dp,j);
        res(i,j) := Eval(dp2,x);
        Clear(dp2);
      end loop;
      Clear(dp);
    end loop;
    return res;
  end Symbolic;

  function Algorithmic
             ( x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim >= 2.

    dim : constant integer32 := x'last;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    fwd : Standard_Complex_Vectors.Vector(1..dim-1);
    bck : Standard_Complex_Vectors.Vector(1..dim-2);
    acc : Complex_Number;

  begin
    for k in 1..dim loop
      res(k,k) := Create(0.0);
    end loop;
    if dim = 2 then
      res(1,2) := Create(1.0); res(2,1) := res(1,2);
    elsif dim = 3 then
      res(1,2) := x(3); res(2,1) := res(1,2);
      res(1,3) := x(2); res(3,1) := res(1,3);
      res(2,3) := x(1); res(3,2) := res(2,3);
    else -- dim > 3
      fwd(1) := x(1)*x(2);
      for k in 2..dim-1 loop
        fwd(k) := fwd(k-1)*x(k+1);
      end loop;
      bck(1) := x(dim)*x(dim-1);
      for k in 2..dim-2 loop
        bck(k) := bck(k-1)*x(dim-k);
      end loop;
     -- last element is copy of fwd(dim-3)
      res(dim-1,dim) := fwd(dim-3); res(dim,dim-1) := res(dim-1,dim);
     -- first element is copy of bck(dim-3)
      res(1,2) := bck(dim-3); res(2,1) := res(1,2);
      if dim = 4 then -- special case for all rows
        res(1,3) := x(2)*x(dim);   res(3,1) := res(1,3);
        res(1,4) := x(2)*x(dim-1); res(4,1) := res(1,4);
        res(2,3) := x(1)*x(dim);   res(3,2) := res(2,3);
        res(2,4) := x(1)*x(dim-1); res(4,2) := res(2,4);
      else -- dim > 4
       -- first row is special, starts with x(2) after diagonal
        res(1,3) := x(2)*bck(dim-4); res(3,1) := res(1,3);
        acc := x(2);
        for k in 4..dim-2 loop
          acc := acc*x(k-1);
          res(1,k) := acc*bck(dim-k-1); res(k,1) := res(1,k);
        end loop;
        acc := acc*x(dim-2);
        res(1,dim-1) := acc*x(dim); res(dim-1,1) := res(1,dim-1);
        res(1,dim) := acc*x(dim-1); res(dim,1) := res(1,dim);
       -- second row is special, starts with x(1) after diagonal
        res(2,3) := x(1)*bck(dim-4); res(3,2) := res(2,3);
        acc := x(1);
        for k in 4..dim-2 loop
          acc := acc*x(k-1);
          res(2,k) := acc*bck(dim-k-1); res(k,2) := res(2,k);
        end loop;
        acc := acc*x(dim-2);
        res(2,dim-1) := acc*x(dim); res(dim-1,2) := res(2,dim-1);
        res(2,dim) := acc*x(dim-1); res(dim,2) := res(2,dim);
       -- the row with index dim-2 has a general formula
        res(dim-2,dim-1) := fwd(dim-4)*x(dim);
        res(dim-1,dim-2) := res(dim-2,dim-1);
        res(dim-2,dim) := fwd(dim-4)*x(dim-1);
        res(dim,dim-2) := res(dim-2,dim);
        for row in 3..dim-3 loop  -- row starts with fwd(row-2)
          res(row,row+1) := fwd(row-2)*bck(dim-row-2);
          res(row+1,row) := res(row,row+1);
          acc := fwd(row-2);
          for k in row+2..dim-2 loop
            acc := acc*x(k-1);
            res(row,k) := acc*bck(dim-k-1); res(k,row) := res(row,k);
          end loop;
          acc := acc*x(dim-2);
          res(row,dim-1) := acc*x(dim); res(dim-1,row) := res(row,dim-1);
          res(row,dim) := acc*x(dim-1); res(dim,row) := res(row,dim);
        end loop;
      end if;
    end if;
    return res;
  end Algorithmic;

  procedure Write_Matrix ( A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix A component-wise, with explicit indexing.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put("] : ");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write_Matrix;
 
  procedure Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the computation of the Hessian of a product of dim variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim > 2.

    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    h0 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
       := Symbolic(x);
    h1 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
       := Algorithmic(x);
    err : double_float;

  begin
    put_line("The Hessian computed symbolically :");
    Write_Matrix(h0);
    put_line("The Hessian computed with algorithmic differentiation :");
    Write_Matrix(h1);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(h0,h1);
    put("Sum of errors :"); put(err,3); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension and then generates
  --   as many complex random numbers as the dimension.

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    Test(dim);
  end Main;

begin
  Main;
end ts_perfhess;
