with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;        use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;     use Standard_Complex_Polynomials_io;
-- with Exponent_Indices;
with Evaluation_Differentiation_Errors;
with Standard_Complex_Circuits;
with Standard_Circuit_Makers;

procedure ts_perfhess is

-- DESCRIPTION :
--   Development of better algorithms to compute Hessians.

  function Symbolic
             ( c : Complex_Number; x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Computes the Hessian matrix of the product of the variables in x,
  --   multiplied with the coefficient c, symbolically, for testing purposes.
    res : Standard_Complex_Matrices.Matrix(x'range,x'range);
    p : Poly;
    t : Term;

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(x'range => 1);
    p := Create(t);
    res := Standard_Circuit_Makers.Hessian(p,x);
    Standard_Complex_Polynomials.Clear(t);
    Standard_Complex_Polynomials.Clear(p);
    return res;
  end Symbolic;

  function Symbolic
             ( c : Complex_Number;
               idx : Standard_Integer_Vectors.Vector;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Computes the Hessian matrix of the product of the variables in x,
  --   with indices of the participating variables in the product in idx,
  --   multiplied with the coefficient c, symbolically, for testing purposes.

    res : Standard_Complex_Matrices.Matrix(x'range,x'range);
    p : Poly;
    t : Term;

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(x'range => 0);
    for k in idx'range loop
      t.dg(idx(k)) := 1;
    end loop;
    p := Create(t);
    res := Standard_Circuit_Makers.Hessian(p,x);
    Standard_Complex_Polynomials.Clear(t);
    Standard_Complex_Polynomials.Clear(p);
    return res;
  end Symbolic;

  function Algorithmic
             ( c : Complex_Number; x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
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
      res(k,k) := Create(0.0); -- all diagonal elements are zero
    end loop;
    if dim = 2 then
      res(1,2) := c; res(2,1) := res(1,2);
    elsif dim = 3 then
      res(1,2) := c*x(3); res(2,1) := res(1,2);
      res(1,3) := c*x(2); res(3,1) := res(1,3);
      res(2,3) := c*x(1); res(3,2) := res(2,3);
    else -- dim > 3
      fwd(1) := x(1)*x(2);
      for k in 2..dim-1 loop
        fwd(k) := fwd(k-1)*x(k+1);
      end loop;
      bck(1) := x(dim)*x(dim-1);
      for k in 2..dim-2 loop
        bck(k) := bck(k-1)*x(dim-k);
      end loop;
     -- last element is copy of fwd(dim-3), multiplied with c
      res(dim-1,dim) := c*fwd(dim-3); res(dim,dim-1) := res(dim-1,dim);
     -- first element is copy of bck(dim-3), multiplied with c
      res(1,2) := c*bck(dim-3); res(2,1) := res(1,2);
      if dim = 4 then -- special case for all rows
	acc := c*x(2);
        res(1,3) := acc*x(dim);   res(3,1) := res(1,3);
        res(1,4) := acc*x(dim-1); res(4,1) := res(1,4);
        acc := c*x(1);
        res(2,3) := acc*x(dim);   res(3,2) := res(2,3);
        res(2,4) := acc*x(dim-1); res(4,2) := res(2,4);
      else -- dim > 4
       -- first row is special, starts with x(2) after diagonal
        acc := c*x(2);
        res(1,3) := acc*bck(dim-4); res(3,1) := res(1,3);
        for k in 4..dim-2 loop
          acc := acc*x(k-1);
          res(1,k) := acc*bck(dim-k-1); res(k,1) := res(1,k);
        end loop;
        acc := acc*x(dim-2);
        res(1,dim-1) := acc*x(dim); res(dim-1,1) := res(1,dim-1);
        res(1,dim) := acc*x(dim-1); res(dim,1) := res(1,dim);
       -- second row is special, starts with x(1) after diagonal
        acc := c*x(1);
        res(2,3) := acc*bck(dim-4); res(3,2) := res(2,3);
        for k in 4..dim-2 loop
          acc := acc*x(k-1);
          res(2,k) := acc*bck(dim-k-1); res(k,2) := res(2,k);
        end loop;
        acc := acc*x(dim-2);
        res(2,dim-1) := acc*x(dim); res(dim-1,2) := res(2,dim-1);
        res(2,dim) := acc*x(dim-1); res(dim,2) := res(2,dim);
       -- the row with index dim-2 has a general formula
        acc := c*fwd(dim-4);
        res(dim-2,dim-1) := acc*x(dim);
        res(dim-1,dim-2) := res(dim-2,dim-1);
        res(dim-2,dim) := acc*x(dim-1);
        res(dim,dim-2) := res(dim-2,dim);
        for row in 3..dim-3 loop  -- row starts with fwd(row-2)
          acc := c*fwd(row-2);
          res(row,row+1) := acc*bck(dim-row-2);
          res(row+1,row) := res(row,row+1);
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

  function Algorithmic
             ( c : Complex_Number;
               idx : Standard_Integer_Vectors.Vector;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   and with idx the indices of participating variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim >= 2.

    dim : constant integer32 := x'last;
    size : constant integer32 := idx'last;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    fwd : Standard_Complex_Vectors.Vector(1..size-1);
    bck : Standard_Complex_Vectors.Vector(1..size-2);
    acc : Complex_Number;

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    if size = 2 then
      res(idx(1),idx(2)) := c;
      res(idx(2),idx(1)) := res(idx(1),idx(2));
    elsif size = 3 then
      res(idx(1),idx(2)) := c*x(idx(3));
      res(idx(2),idx(1)) := res(idx(1),idx(2));
      res(idx(1),idx(3)) := c*x(idx(2));
      res(idx(3),idx(1)) := res(idx(1),idx(3));
      res(idx(2),idx(3)) := c*x(idx(1));
      res(idx(3),idx(2)) := res(idx(2),idx(3));
    else -- size > 3
      fwd(1) := x(idx(1))*x(idx(2));
      for k in 2..size-1 loop
        fwd(k) := fwd(k-1)*x(idx(k+1));
      end loop;
      bck(1) := x(idx(size))*x(idx(size-1));
      for k in 2..size-2 loop
        bck(k) := bck(k-1)*x(idx(size-k));
      end loop;
     -- last element is copy of fwd(size-3), multiplied with c
      res(idx(size-1),idx(size)) := c*fwd(size-3);
      res(idx(size),idx(size-1)) := res(idx(size-1),idx(size));
     -- first element is copy of bck(size-3), multiplied with c
      res(idx(1),idx(2)) := c*bck(size-3);
      res(idx(2),idx(1)) := res(idx(1),idx(2));
      if size = 4 then -- special case for all rows
        acc := c*x(idx(2));
        res(idx(1),idx(3)) := acc*x(idx(size));
        res(idx(3),idx(1)) := res(idx(1),idx(3));
        res(idx(1),idx(4)) := acc*x(idx(size-1));
        res(idx(4),idx(1)) := res(idx(1),idx(4));
        acc := c*x(idx(1));
        res(idx(2),idx(3)) := acc*x(idx(size));
        res(idx(3),idx(2)) := res(idx(2),idx(3));
        res(idx(2),idx(4)) := acc*x(idx(size-1));
        res(idx(4),idx(2)) := res(idx(2),idx(4));
      else -- size > 4
       -- first row is special, starts with x(idx(2)) after diagonal
        acc := c*x(idx(2));
        res(idx(1),idx(3)) := acc*bck(size-4);
        res(idx(3),idx(1)) := res(idx(1),idx(3));
        for k in 4..size-2 loop
          acc := acc*x(idx(k-1));
          res(idx(1),idx(k)) := acc*bck(size-k-1);
          res(idx(k),idx(1)) := res(idx(1),idx(k));
        end loop;
        acc := acc*x(idx(size-2));
        res(idx(1),idx(size-1)) := acc*x(idx(size));
        res(idx(size-1),idx(1)) := res(idx(1),idx(size-1));
        res(idx(1),idx(size)) := acc*x(idx(size-1));
        res(idx(size),idx(1)) := res(idx(1),idx(size));
       -- second row is special, starts with x(idx(1)) after diagonal
        acc := c*x(idx(1));
        res(idx(2),idx(3)) := acc*bck(size-4);
        res(idx(3),idx(2)) := res(idx(2),idx(3));
        for k in 4..size-2 loop
          acc := acc*x(idx(k-1));
          res(idx(2),idx(k)) := acc*bck(size-k-1);
          res(idx(k),idx(2)) := res(idx(2),idx(k));
        end loop;
        acc := acc*x(idx(size-2));
        res(idx(2),idx(size-1)) := acc*x(idx(size));
        res(idx(size-1),idx(2)) := res(idx(2),idx(size-1));
        res(idx(2),idx(size)) := acc*x(idx(size-1));
        res(idx(size),idx(2)) := res(idx(2),idx(size));
       -- the row with index size-2 has a general formula
        acc := c*fwd(size-4);
        res(idx(size-2),idx(size-1)) := acc*x(idx(size));
        res(idx(size-1),idx(size-2)) := res(idx(size-2),idx(size-1));
        res(idx(size-2),idx(size)) := acc*x(idx(size-1));
        res(idx(size),idx(size-2)) := res(idx(size-2),idx(size));
        for row in 3..size-3 loop  -- row starts with fwd(row-2)
          acc := c*fwd(row-2);
          res(idx(row),idx(row+1)) := acc*bck(size-row-2);
          res(idx(row+1),idx(row)) := res(idx(row),idx(row+1));
          for k in row+2..size-2 loop
            acc := acc*x(idx(k-1));
            res(idx(row),idx(k)) := acc*bck(size-k-1);
            res(idx(k),idx(row)) := res(idx(row),idx(k));
          end loop;
          acc := acc*x(idx(size-2));
          res(idx(row),idx(size-1)) := acc*x(idx(size));
          res(idx(size-1),idx(row)) := res(idx(row),idx(size-1));
          res(idx(row),idx(size)) := acc*x(idx(size-1));
          res(idx(size),idx(row)) := res(idx(row),idx(size));
        end loop;
      end if;
    end if;
    return res;
  end Algorithmic;

  procedure Algorithmic
             ( H : in out Standard_Complex_Matrices.Matrix;
               c : in Complex_Number;
               idx : in Standard_Integer_Vectors.Vector;
               x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   and with idx the indices of participating variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- ON ENTRY :
  --   H       the current Hessian matrix,
  --           initialized with zero if called for the first time.
  --   c       coefficient of the term in the circuit;
  --   idx     index of the participating variables;
  --   x       values for all variables.

  -- ON RETURN :
  --   H       updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last >= 2.

    sz : constant integer32 := idx'last;
    fwd : Standard_Complex_Vectors.Vector(1..sz-1);
    bck : Standard_Complex_Vectors.Vector(1..sz-2);
    acc : Complex_Number;

  begin
    if sz = 2 then
      H(idx(1),idx(2)) := H(idx(1),idx(2)) + c;
    elsif sz = 3 then
      H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*x(idx(3));
      H(idx(1),idx(3)) := H(idx(1),idx(3)) + c*x(idx(2));
      H(idx(2),idx(3)) := H(idx(2),idx(3)) + c*x(idx(1));
    else -- sz > 3
      fwd(1) := x(idx(1))*x(idx(2));
      for k in 2..sz-1 loop
        fwd(k) := fwd(k-1)*x(idx(k+1));
      end loop;
      bck(1) := x(idx(sz))*x(idx(sz-1));
      for k in 2..sz-2 loop
        bck(k) := bck(k-1)*x(idx(sz-k));
      end loop;
     -- last element is copy of fwd(sz-3), multiplied with c
      H(idx(sz-1),idx(sz)) := H(idx(sz-1),idx(sz)) + c*fwd(sz-3);
     -- first element is copy of bck(sz-3), multiplied with c
      H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*bck(sz-3);
      if sz = 4 then -- special case for all rows
        acc := c*x(idx(2));
        H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*x(idx(sz));
        H(idx(1),idx(4)) := H(idx(1),idx(4)) + acc*x(idx(sz-1));
        acc := c*x(idx(1));
        H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*x(idx(sz));
        H(idx(2),idx(4)) := H(idx(2),idx(4)) + acc*x(idx(sz-1));
      else -- sz > 4
       -- first row is special, starts with x(idx(2)) after diagonal
        acc := c*x(idx(2));
        H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*bck(sz-4);
        for k in 4..sz-2 loop
          acc := acc*x(idx(k-1));
          H(idx(1),idx(k)) := H(idx(1),idx(k)) + acc*bck(sz-k-1);
        end loop;
        acc := acc*x(idx(sz-2));
        H(idx(1),idx(sz-1)) := H(idx(1),idx(sz-1)) + acc*x(idx(sz));
        H(idx(1),idx(sz)) := H(idx(1),idx(sz)) + acc*x(idx(sz-1));
       -- second row is special, starts with x(idx(1)) after diagonal
        acc := c*x(idx(1));
        H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*bck(sz-4);
        for k in 4..sz-2 loop
          acc := acc*x(idx(k-1));
          H(idx(2),idx(k)) := H(idx(2),idx(k)) + acc*bck(sz-k-1);
        end loop;
        acc := acc*x(idx(sz-2));
        H(idx(2),idx(sz-1)) := H(idx(2),idx(sz-1)) + acc*x(idx(sz));
        H(idx(2),idx(sz)) := H(idx(2),idx(sz)) + acc*x(idx(sz-1));
       -- the row with index sz-2 has a general formula
        acc := c*fwd(sz-4);
        H(idx(sz-2),idx(sz-1)) := H(idx(sz-2),idx(sz-1)) + acc*x(idx(sz));
        H(idx(sz-2),idx(sz)) := H(idx(sz-2),idx(sz)) + acc*x(idx(sz-1));
        for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
          acc := c*fwd(rw-2);
          H(idx(rw),idx(rw+1)) := H(idx(rw),idx(rw+1)) + acc*bck(sz-rw-2);
          for k in rw+2..sz-2 loop
            acc := acc*x(idx(k-1));
            H(idx(rw),idx(k)) := H(idx(rw),idx(k)) + acc*bck(sz-k-1);
          end loop;
          acc := acc*x(idx(sz-2));
          H(idx(rw),idx(sz-1)) := H(idx(rw),idx(sz-1)) + acc*x(idx(sz));
          H(idx(rw),idx(sz)) := H(idx(rw),idx(sz)) + acc*x(idx(sz-1));
        end loop;
      end if;
    end if;
  end Algorithmic;

  function Algorithmic
             ( c : Standard_Complex_Circuits.Circuit;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the Hessian matrix of the circuit c at x.

    dim : constant integer32 := x'last;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for k in 1..c.nbr loop
      declare
        idx : constant Standard_Integer_Vectors.Vector := c.xps(k).all;
           -- := Exponent_Indices.Exponent_Index(c.xps(k).all);
      begin
        put("update for "); put(idx); put_line(" ...");
        Algorithmic(res,c.cff(k),idx,x);
      end;
    end loop;
    for i in 2..dim loop
      for j in 1..(i-1) loop
        res(i,j) := res(j,i);
      end loop;
    end loop;
    return res;
  end Algorithmic;

  procedure Test_Product ( dim : in integer32; size : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Tests the computation of the Hessian of a product of dim variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim > 2.

    c : constant Complex_Number := Standard_Random_Numbers.Random1;
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    h0,h1 : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    err : double_float;

  begin
    if size = 0 then
      h0 := Symbolic(c,x);
      h1 := Algorithmic(c,x);
    else
      declare
        idx : constant Standard_Integer_Vectors.Vector(1..size)
            := Standard_Circuit_Makers.Random_Indices(dim,size);
      begin
        put("The indices : "); put(idx); new_line;
        h0 := Symbolic(c,idx,x);
        h1 := Algorithmic(c,idx,x);
      end;
    end if;
    put_line("The Hessian computed symbolically :");
    Standard_Circuit_Makers.Write_Matrix(h0);
    put_line("The Hessian computed with algorithmic differentiation :");
    Standard_Circuit_Makers.Write_Matrix(h1);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(h0,h1);
    put("Sum of errors :"); put(err,3); new_line;
  end Test_Product;

  procedure Test_Circuit ( dim,nbr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables and nbr of terms
  --   to test the computation of the Hessian.

    c : constant Standard_Complex_Circuits.Circuit
      := Standard_Circuit_Makers.Random_Complex_Circuit(nbr,dim);
    p : constant Standard_Complex_Polynomials.Poly
      := Standard_Circuit_Makers.Make_Polynomial(c,true);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    h0 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
       := Standard_Circuit_Makers.Hessian(p,x);
    h1 : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
       := Algorithmic(c,x);
    err : double_float;

  begin
    new_line;
    put_line("The polynomial : "); put(p); new_line;
    put_line("The Hessian computed symbolically :");
    Standard_Circuit_Makers.Write_Matrix(h0);
    put_line("The Hessian computed algorithmically :");
    Standard_Circuit_Makers.Write_Matrix(h1);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(h0,h1);
    put("Sum of errors :"); put(err,3); new_line;
  end Test_Circuit;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension and then generates
  --   as many complex random numbers as the dimension.

    dim,size,nbr : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("MENU to test the algorithmic Hessian computations :");
    put_line("  0. on a product of variables");
    put_line("  1. on a random circuit");
    put("Type 0 or 1 to select the test : "); Ask_Alternative(ans,"01");
    case ans is
      when '0' =>
        new_line;
        put("Indexed product of variables ? (y/n) "); Ask_Yes_or_No(ans);
        if ans /= 'y' then
          new_line;
          put("Give the dimension : "); get(dim);
          Test_Product(dim);
        else
          put("Give the size of the product : "); get(size);
          put("Give the dimension (> "); put(size,1); put(") : "); get(dim);
          Test_Product(dim,size);
        end if;
      when '1' =>
        new_line;
        put("Give the dimension : "); get(dim);
        put("Give the number of terms : "); get(nbr);
        Test_Circuit(dim,nbr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfhess;
