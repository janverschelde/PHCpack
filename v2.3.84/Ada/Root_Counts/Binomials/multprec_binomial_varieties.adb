with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
with Multprec_Integer_Linear_Solvers;
with Multprec_Integer_Matrix_Inverse;
with Multprec_Smith_Normal_Form;
with Multprec_Integer_Kernel;

package body Multprec_Binomial_Varieties is

-- STAGE I

  procedure Cone_of_Tropisms 
               ( A : in Multprec_Integer_Matrices.Matrix;
                 rank : out integer32;
                 V : out Multprec_Integer_Matrices.Link_to_Matrix ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    B : Multprec_Integer_Matrices.Matrix(1..k,1..n) 
      := Multprec_Integer_Matrices.Transpose(A);

  begin
    Multprec_Integer_Kernel.Kernel(B,rank,V);
    Multprec_Integer_Matrices.Clear(B);
  end Cone_of_Tropisms;

  procedure Check_Inner_Products
               ( A : in Multprec_Integer_Matrices.Matrix;
                 V : in Multprec_Integer_Matrices.Matrix;
                 output : in boolean; bug : out boolean ) is

    use Multprec_Integer_Numbers;
    sum,acc : Integer_Number;

  begin
    bug := false;
    for j in V'range(2) loop    -- all inner products of j-th tropism
      for k in A'range(2) loop  -- with all columns of A
        sum := Create(integer(0));
        for i in A'range(1) loop
          acc := A(i,k)*V(i,j);
          Add(sum,acc); Clear(acc);
        end loop;
        if output
         then put(" "); Multprec_Integer_Numbers_io.put(sum);
        end if;
        bug := not Equal(sum,0);
        Clear(sum);
        exit when bug;
      end loop;
      if output
       then new_line;
      end if;
    end loop;
  end Check_Inner_Products;

  procedure Check_Rank
               ( V : in Multprec_Integer_Matrices.Matrix;
                 rank : in integer32; output : in boolean;
                 bug : out boolean ) is

    n : constant integer32 := V'last(1);
    d : constant integer32 := V'last(2);
    r : integer32;

  begin
    if rank /= n-d then
      bug := true;
    else
      r := Multprec_Integer_Linear_Solvers.Rank(V);
      if output
       then put("Rank of tropisms cone : "); put(r,1); new_line;
      end if;
      bug := not (r = d);
    end if;
  end Check_Rank;

  procedure Check_Cone 
               ( A : in Multprec_Integer_Matrices.Matrix;
                 V : in Multprec_Integer_Matrices.Matrix;
                 rank : in integer32; output : in boolean;
                 bug : out boolean ) is
  begin
    Check_Inner_Products(A,V,output,bug);
    if not bug
     then Check_Rank(V,rank,output,bug);
    end if;
  end Check_Cone;

  procedure Expected_Dimension
               ( A,V : in Multprec_Integer_Matrices.Matrix;
                 rnk : in integer32; output : in boolean;
                 dim : out integer32 ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    d,m : integer32;

  begin
    d := n - rnk;
    if output then
      if rnk = k then
        put_line("The rank equals the codimension.");
        put("=> For nonzero coefficients, we have a "); put(d,1);
        put_line("-dimensional solution set.");
      else
        put_line("The rank does not equal the codimension.");
      end if;
      if d > 0 then
        put_line("The cone of tropisms is spanned by"); put(V); 
        m := V'last(2);
        put("The expected dimension is "); put(d,1);
        put(", found "); put(m,1); put_line(" tropisms.");
        dim := d;
      else
        dim := 0;
      end if;
    end if;
  end Expected_Dimension;

-- STAGE II

  procedure UV_Coordinate_Transformation
              ( U,V : in Multprec_Integer_Matrices.Matrix;
                M : out Multprec_Integer_Matrices.Matrix ) is

    Uinv : Multprec_Integer_Matrices.Matrix(U'range(1),U'range(2));
    Vinv : Multprec_Integer_Matrices.Matrix(V'range(1),V'range(2));
    E : Multprec_Integer_Matrices.Matrix(V'range(1),V'range(2));
    use Multprec_Integer_Numbers,Multprec_Integer_Matrices;
    d : Integer_Number;

  begin
    Multprec_Integer_Matrix_Inverse.Inverse(U,d,Uinv);
    Clear(d);
    for i in Uinv'range(1) loop
      for j in Uinv'range(2) loop
        Copy(Uinv(i,j),E(i,j));
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        E(i,j) := Create(integer(0));
      end loop;
    end loop;
    for i in Uinv'last(1)+1..E'last(1) loop
      for j in Uinv'range(2) loop
        E(i,j) := Create(integer(0));
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        if i = j
         then E(i,j) := Create(integer(1));
         else E(i,j) := Create(integer(0));
        end if;
      end loop;
    end loop;
    Clear(Uinv);
    Multprec_Integer_Matrix_Inverse.Inverse(V,d,Vinv);
    Clear(d);
    M := E*Vinv;
    Clear(Vinv); Clear(E);
  end UV_Coordinate_Transformation;

  function Diagonal_Product 
               ( A : Multprec_Integer_Matrices.Matrix )
               return Multprec_Integer_Numbers.Integer_Number is

    use Multprec_Integer_Numbers;
    res : Integer_Number := Create(integer(1));
    acc : Integer_Number;

  begin
    for i in A'range(1) loop
      exit when (i > A'last(2));
      acc := res*A(i,i);
      Mul(res,acc); Clear(acc);
    end loop;
    return res;
  end Diagonal_Product;

  procedure Unimodular_Coordinate_Transformation
               ( C : in Multprec_Integer_Matrices.Matrix;
                 U,T,V,M : out Multprec_Integer_Matrices.Matrix;
                 p : out Multprec_Integer_Numbers.Integer_Number;
                 uv,fail : out boolean ) is

    use Multprec_Integer_Numbers;
    d : Integer_Number;

  begin
    T := Multprec_Integer_Matrices.Transpose(C);
    Multprec_Smith_Normal_Form.Diagonalize(U,T,V);
    p := Diagonal_Product(T);
    if Multprec_Integer_Matrix_Inverse.Is_Identity(U) then
      Multprec_Integer_Matrix_Inverse.Inverse(V,d,M);
      Clear(d);
      fail := false; uv := false;
    elsif Equal(p,1) then
      UV_Coordinate_Transformation(U,V,M);
      fail := false; uv := true;
    else
      fail := true; uv := false;
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
               ( C : in Multprec_Integer_Matrices.Matrix;
                 output : in boolean; fail : out boolean;
                 M : out Multprec_Integer_Matrices.Matrix ) is

    T : Multprec_Integer_Matrices.Matrix(C'range(2),C'range(1))
      := Multprec_Integer_Matrices.Transpose(C); 
    U : Multprec_Integer_Matrices.Matrix(T'range(1),T'range(1))
      := Multprec_Smith_Normal_Form.Identity(natural32(T'last(1)));
    V : Multprec_Integer_Matrices.Matrix(T'range(2),T'range(2))
      := Multprec_Smith_Normal_Form.Identity(natural32(T'last(2)));
    use Multprec_Integer_Numbers,Multprec_Integer_Matrices;
    p : Integer_Number;
    uv : boolean;

  begin
    Unimodular_Coordinate_Transformation(C,U,T,V,M,p,uv,fail);
    if output then
      put_line("The Smith normal form of the tropisms :"); put(T);
      put("The product of the diagonal of the Smith form : ");
      Multprec_Integer_Numbers_io.put(p); put_line(".");
      if fail then
        put_line("No rational parameterization possible.");
      else
        if not uv then
          put_line("U is the identity matrix => M = V^-1.");
        else
          put_line("U is not the identity matrix.");
          put("Smith form diagonal product = 1 => M = E V^-1,");
          put_line(" E contains U^-1.");
        end if;
        put_line("The unimodular transformation M :"); put(M);
      end if;
    end if;
    Clear(T); Clear(U); Clear(V); Clear(p);
  end Unimodular_Coordinate_Transformation;

  function Is_Zero_Row
               ( A : Multprec_Integer_Matrices.Matrix; i : integer32 )
               return boolean is

    use Multprec_Integer_Numbers;

  begin
    for j in A'range(2) loop
      if not Equal(A(i,j),0)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero_Row;

  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Multprec_Integer_Matrices.Matrix;
                 dim : in integer32; output : in boolean;
                 nzr : out natural32; fail : out boolean ) is

    use Multprec_Integer_Matrices;
    p : Matrix(M'range(1),A'range(2)) := M*A;

  begin
    if output then
      put("Unimodular transformation M multiplied by exponent matrix A, ");
      put_line("M*A :"); put(p);
    end if;
    fail := false;
    nzr := 0;
    for i in 1..dim loop
      if not Is_Zero_Row(p,i) then
        if output then
          put("Row "); put(i,1); put_line(" of M*A is nonzero, bug!");
          fail := true;
        end if;
      else
        nzr := nzr + 1;
      end if;
    end loop;
    for i in dim+1..p'last(1) loop
      if Is_Zero_Row(p,i)
       then nzr := nzr + 1;
      end if;
    end loop;
    if output
     then put("The number of nonzero rows in M*A : "); put(nzr,1); new_line;
    end if;
    Clear(p);
  end Test_Unimodular_Coordinate_Transformation;

-- STAGE III

  procedure Upper_Transformed_Exponents
               ( A,M : in Multprec_Integer_Matrices.Matrix;
                 dim : in integer32; output : in boolean;
                 U : out Multprec_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector ) is

    use Multprec_Integer_Matrices;
    p : Matrix(M'range(1),A'range(2)) := M*A;
    W : Matrix(1..M'last(1)-dim,A'range(2));

  begin
    for i in W'range(1) loop
      for j in W'range(2) loop
        Multprec_Integer_Numbers.Copy(p(i+dim,j),W(i,j));
      end loop;
    end loop;
    Clear(p);
    Multprec_Integer_Linear_Solvers.Upper_Triangulate(W);
    if output
     then put_line("The upper triangulated M*A : "); put(W);
    end if;
    Multprec_Integer_Kernel.Pivots_in_Upper(W,rank,pivots);
    if output then
      put("The rank of the matrix M*A : "); put(rank,1); new_line;
      put("The pivots : "); put(pivots); new_line; 
    end if;
    U := W;
  end Upper_Transformed_Exponents;

  function Product_of_Pivots
              ( U : Multprec_Integer_Matrices.Matrix;
                pivots : Standard_Integer_Vectors.Vector ) 
              return Multprec_Integer_Numbers.Integer_Number is

    use Multprec_Integer_Numbers;
    res : Integer_Number := Create(integer(1));
    acc : Integer_Number;

  begin
    for i in pivots'range loop
      exit when i > U'last(1);
      acc := res*U(i,pivots(i));
      Mul(res,acc); Clear(acc);
    end loop;
    return res;
  end Product_of_Pivots;

end Multprec_Binomial_Varieties;
