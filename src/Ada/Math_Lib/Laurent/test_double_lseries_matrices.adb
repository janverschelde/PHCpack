with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Double_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;
with Double_Linear_Laurent_Solvers;     use Double_Linear_Laurent_Solvers;

package body Test_Double_Lseries_Matrices is

  procedure Matrix_Matrix_Product
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Clead : out Standard_Integer_Matrices.Matrix;
                Ccffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Arow,Brow,Crow : Standard_Complex_VecVecs.Link_to_VecVec;
    Acff,Bcff,Ccff : Standard_Complex_Vectors.Link_to_Vector;
    ze,ewrk,eprd : integer32;
    zc,cwrk,cprd : Standard_Complex_Vectors.Vector(0..deg);

  begin
    for i in 1..nrows loop           -- initialize C to zero
      Crow := Ccffs(i); 
      for j in 1..ncols loop
        Clead(i,j) := 0;
        Ccff := Crow(j);             -- coefficients of C(i,j)
        for k in 0..deg loop
          Ccff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
    end loop;
    for i in 1..nrows loop
      Arow := Acffs(i);
      Crow := Ccffs(i);
      for j in 1..ncols loop
        ewrk := 0;              -- accumulates leading exponent of C(i,j)
        for k in 0..deg loop
          cwrk(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        for k in 1..ncols loop
          Acff := Arow(k);                   -- Acff is A(i,k)
          Brow := Bcffs(k); Bcff := Brow(j); -- Bcff is B(k,j)
          Double_Laurent_Series.Multiply
            (deg,Alead(i,k),Blead(k,j),Acff.all,Bcff.all,eprd,cprd);
         -- eprd is the leading exponent of A(i,k)*B(k,j)
         -- cprd has the coefficients of A(i,k)*B(k,j)
          Double_Laurent_Series.Add(deg,ewrk,eprd,cwrk,cprd,ze,zc);
          ewrk := ze;          -- leading exponent of accumulator 
          for L in 0..deg loop -- copy coefficients of accumulator
            cwrk(L) := zc(L);
          end loop;
        end loop;
        Clead(i,j) := ewrk;
        Ccff := Crow(j);
        for k in 0..deg loop
          Ccff(k) := cwrk(k);
        end loop;
      end loop;
    end loop;
  end Matrix_Matrix_Product;

  procedure Write_Difference
              ( deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Arow,Brow : Standard_Complex_VecVecs.Link_to_VecVec;
    Acff,Bcff : Standard_Complex_Vectors.Link_to_Vector;
    ze : integer32;
    zc : Standard_Complex_Vectors.Vector(0..deg);
    nrm : double_float := 0.0;

  begin
    for i in Alead'range(1) loop
      Arow := Acffs(i);
      Brow := Bcffs(i);
      for j in Alead'range(2) loop
        Acff := Arow(j);
        Bcff := Brow(j);
        Double_Laurent_Series.Subtract
          (deg,Alead(i,j),Blead(i,j),Acff.all,Bcff.all,ze,zc);
        put("D("); put(i,1); put(","); put(j,1); put_line(") :");
        Double_Laurent_Series.Write(ze,zc);
        for k in 0..deg loop
          nrm := nrm + Standard_Complex_Numbers.AbsVal(zc(k));
        end loop;
      end loop;
    end loop;
    put("The difference sum : "); put(nrm); new_line;
  end Write_Difference;

  procedure Plain_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Llead : out Standard_Integer_Matrices.Matrix;
                Lcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Ulead : out Standard_Integer_Matrices.Matrix;
                Ucffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Lrow,Urow,Arow : Standard_Complex_VecVecs.Link_to_VecVec;
    Lcff,Ucff,Acff : Standard_Complex_Vectors.Link_to_Vector;
    ze,ewrk,eprd : integer32;
    zc,cwrk,cprd : Standard_Complex_Vectors.Vector(0..deg);

  begin
    for i in 1..nrows loop -- initialization of L and U
      Urow := Ucffs(i);
      for j in 1..(i-1) loop -- set the lower triangular part of U to zero
        Ulead(i,j) := 0;
        Ucff := Urow(j);     -- coefficients of U(i,j)
        for k in 0..deg loop
          Ucff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
      Llead(i,i) := 0;
      Lrow := Lcffs(i);      -- set the diagonal of L to one
      Lcff := Lrow(i);       -- coefficients of L(i,i)
      Lcff(0) := Standard_Complex_Numbers.Create(1.0);
      for k in 1..deg loop
        Lcff(k) := Standard_Complex_Numbers.Create(0.0);
      end loop;
      for j in (i+1)..ncols loop -- make upper triangular part of L zero
        Llead(i,j) := 0;
        Lcff := Lrow(j);     -- coefficients of L(i,j)
        for k in 0..deg loop
          Lcff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
    end loop;
    for j in 1..ncols loop
      for i in 1..j loop              -- row reduction assigns to U
        Arow := Acffs(i); Acff := Arow(j);    -- U(i,j) := A(i,j)
        ewrk := Alead(i,j);
        for k in 0..deg loop
          cwrk(k) := Acff(k);                 -- accumulates U(i,j)
        end loop;
        for k in 1..(i-1) loop        -- U(i,j) := U(i,j) - L(i,k)*U(k,j)
          Lrow := Lcffs(i); Lcff := Lrow(k); -- Lcff is L(i,k)
          Urow := Ucffs(k); Ucff := Urow(j); -- Ucff is U(k,j)
          Double_Laurent_Series.Multiply
            (deg,Llead(i,k),Ulead(k,j),Lcff.all,Ucff.all,eprd,cprd);
         -- eprd is the leading exponent of L(i,k)*U(k,j)
         -- cprd has the coefficients of L(i,k)*U(k,j)
          Double_Laurent_Series.Subtract(deg,ewrk,eprd,cwrk,cprd,ze,zc);
          ewrk := ze;          -- leading exponent of U(i,j) - L(i,k)*U(k,j)
          for L in 0..deg loop -- copy coefficients of U(i,j) - L(i,k)*U(k,j)
            cwrk(L) := zc(L);
          end loop;
        end loop;
        Ulead(i,j) := ewrk;                -- copy accumulator to U(i,j)
        Urow := Ucffs(i); Ucff := Urow(j); -- Ucff is U(i,j)
        for k in 0..deg loop
          Ucff(k) := cwrk(k);
        end loop;
      end loop;
      for i in j+1..nrows loop   -- row reduction assigns to L
        Arow := Acffs(i); Acff := Arow(j);      -- L(i,j) := A(i,j)
        ewrk := Alead(i,j);
        for k in 0..deg loop
          cwrk(k) := Acff(k);                   -- accumulates L(i,j)
        end loop;
        for k in 1..(j-1) loop -- L(i,j) := L(i,j) - L(i,k)*U(k,j)
          Lrow := Lcffs(i); Lcff := Lrow(k); -- Lcff is L(i,k)
          Urow := Ucffs(k); Ucff := Urow(j); -- Ucff is U(k,j)
          Double_Laurent_Series.Multiply
            (deg,Llead(i,k),Ulead(k,j),Lcff.all,Ucff.all,eprd,cprd);
         -- eprd is the leading exponent of L(i,k)*U(k,j)
         -- cprd has the coefficients of L(i,k)*U(k,j)
          Double_Laurent_Series.Subtract(deg,ewrk,eprd,cwrk,cprd,ze,zc);
          ewrk := ze;          -- leading exponent of L(i,j) - L(i,k)*U(k,j)
          for L in 0..deg loop -- copy coefficients of L(i,j) - L(i,k)*U(k,j)
            cwrk(L) := zc(L);
          end loop;
        end loop;
        Llead(i,j) := ewrk;                -- copy accumulator to L(i,j)
        Lrow := Lcffs(i); Lcff := Lrow(j); -- Lcff is L(i,j)
        for k in 0..deg loop
          Lcff(k) := cwrk(k);
        end loop;
        Urow := Ucffs(j); Ucff := Urow(j); -- Ucff is U(j,j) 
        Double_Laurent_Series.Divide
          (deg,Llead(i,j),Ulead(j,j),Lcff.all,Ucff.all,ze,zc,cwrk);
        Llead(i,j) := ze;    -- leading exponent of L(i,j)/U(j,j)
        for k in 0..deg loop -- copy coefficients of L(i,j)/U(j,j)
          Lcff(k) := zc(k);
        end loop;
      end loop;
    end loop;
  end Plain_LU_Factorization;

  procedure Test_Plain_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Llead,Ulead,Plead : Standard_Integer_Matrices.Matrix(1..nrows,1..ncols);
    Lcffs,Ucffs,Pcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Complex_VecVecVecs.Allocate(Lcffs,1,nrows,1,ncols,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Ucffs,1,nrows,1,ncols,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Pcffs,1,nrows,1,ncols,0,deg);
    put_line("Computing a plain LU factorization ...");
    Plain_LU_Factorization
      (nrows,ncols,deg,Alead,Acffs,Llead,Lcffs,Ulead,Ucffs);
    put_line("Computing the product of L with U ...");
    Matrix_Matrix_Product
      (nrows,ncols,deg,Llead,Lcffs,Ulead,Ucffs,Plead,Pcffs);
    Write(Llead,Lcffs,"L");
    Write(Ulead,Ucffs,"U");
    Write(Plead,Pcffs,"P");
    Write(Alead,Acffs,"A");
    Write_Difference(deg,Alead,Acffs,Plead,Pcffs);
  end Test_Plain_LU_Factorization;

  procedure Lower_Triangular_Part
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Llead : out Standard_Integer_Matrices.Matrix;
                Lcffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Lrow,Arow : Standard_Complex_VecVecs.Link_to_VecVec;
    Lcff,Acff : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in 1..nrows loop
      Lrow := Lcffs(i);
      Arow := Acffs(i);
      for j in 1..(i-1) loop
        Llead(i,j) := Alead(i,j);
        Lcff := Lrow(j);
        Acff := Arow(j);
        for k in 0..deg loop
          Lcff(k) := Acff(k);
        end loop;
      end loop;
      Llead(i,i) := 0;
      Lcff := Lrow(i);
      Lcff(0) := Standard_Complex_Numbers.Create(1.0);
      for k in 1..deg loop
        Lcff(k) := Standard_Complex_Numbers.Create(0.0);
      end loop;
      for j in (i+1)..ncols loop
        Llead(i,j) := 0;
        Lcff := Lrow(j);
        for k in 0..deg loop
          Lcff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
    end loop;
  end Lower_Triangular_Part;

  procedure Upper_Triangular_Part
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Ulead : out Standard_Integer_Matrices.Matrix;
                Ucffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Urow,Arow : Standard_Complex_VecVecs.Link_to_VecVec;
    Ucff,Acff : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in 1..nrows loop
      Urow := Ucffs(i);
      for j in 1..(i-1) loop
        Ulead(i,j) := 0;
        Ucff := Urow(j);
        for k in 0..deg loop
          Ucff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
      Arow := Acffs(i);
      for j in i..ncols loop
        Ulead(i,j) := Alead(i,j);
        Ucff := Urow(j);
        Acff := Arow(j);
        for k in 0..deg loop
          Ucff(k) := Acff(k);
        end loop;
      end loop;
    end loop;
  end Upper_Triangular_Part;

  procedure Permute
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : in Standard_Integer_Vectors.Vector;
                Blead : out Standard_Integer_Matrices.Matrix;
                Bcffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Arow,Brow : Standard_Complex_VecVecs.Link_to_VecVec;
    Acff,Bcff : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in 1..nrows loop
      Arow := Acffs(pivots(i));
      Brow := Bcffs(i);
      for j in 1..ncols loop
        Acff := Arow(j);
        Bcff := Brow(j);
        for k in 0..deg loop
          Bcff(k) := Acff(k);
        end loop;
        Blead(i,j) := Alead(pivots(i),j);
      end loop;
    end loop;
  end Permute;

  procedure Copy
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : out Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Arow,Brow : Standard_Complex_VecVecs.Link_to_VecVec;
    Acff,Bcff : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in 1..nrows loop
      Arow := Acffs(i);
      Brow := Bcffs(i);
      for j in 1..ncols loop
        Blead(i,j) := Alead(i,j);
        Acff := Arow(j);
        Bcff := Brow(j);
        for k in 0..deg loop
          Bcff(k) := Acff(k);
        end loop;
      end loop;
    end loop;
  end Copy;

  procedure Test_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

    Llead,Ulead,Plead : Standard_Integer_Matrices.Matrix(1..nrows,1..ncols);
    Lcffs,Ucffs,Pcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    pivots : Standard_Integer_Vectors.Vector(1..ncols);
    A2lead,Blead : Standard_Integer_Matrices.Matrix(1..nrows,1..ncols);
    A2cffs,Bcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Complex_VecVecVecs.Allocate(A2cffs,1,nrows,1,ncols,0,deg);
    Copy(nrows,ncols,deg,Alead,Acffs,A2lead,A2cffs);
    Standard_Complex_VecVecVecs.Allocate(Lcffs,1,nrows,1,ncols,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Ucffs,1,nrows,1,ncols,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Pcffs,1,nrows,1,ncols,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Bcffs,1,nrows,1,ncols,0,deg);
    put_line("Computing a LU factorization ...");
    LU_Factorization(nrows,ncols,deg,Alead,Acffs,pivots);
    Lower_Triangular_Part(nrows,ncols,deg,Alead,Acffs,Llead,Lcffs);
    Upper_Triangular_Part(nrows,ncols,deg,Alead,Acffs,Ulead,Ucffs);
    put_line("Computing the product of L with U ...");
    Matrix_Matrix_Product
      (nrows,ncols,deg,Llead,Lcffs,Ulead,Ucffs,Plead,Pcffs);
    Write(Llead,Lcffs,"L");
    Write(Ulead,Ucffs,"U");
    Write(Plead,Pcffs,"P");
    put("The pivots : "); put(pivots); new_line;
    Permute(nrows,ncols,deg,A2lead,A2cffs,pivots,Blead,Bcffs);
    Write(Blead,Bcffs,"permutedA");
    Write_Difference(deg,Blead,Bcffs,Plead,Pcffs);
  end Test_LU_Factorization;

  procedure Test ( nrows,ncols,deg,low,upp : in integer32;
                   lower,upper : in boolean ) is

    nbrows : constant natural32 := natural32(nrows);
    nbcols : constant natural32 := natural32(ncols);
    Alead : Standard_Integer_Matrices.Matrix(1..nrows,1..ncols);
    xlead : Standard_Integer_Vectors.Vector(1..ncols);
    Acffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    xcffs,ycffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ylead : Standard_Integer_Vectors.Vector(1..ncols);
    blead : Standard_Integer_Vectors.Vector(1..nrows);
    bcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ans : character;

  begin
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,nrows,1,ncols,0,deg);
    if not lower and not upper then
      Random_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    elsif lower then
      Random_Lower_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    else -- upper must be true
      Random_Upper_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    end if;
    put_line("The matrix of leading exponents :"); put(Alead,1);
    put("A "); put(nrows,1); put("-by-"); put(ncols,1);
    put_line(" matrix of Laurent series : "); Write(Alead,Acffs);
    Random_Vector(ncols,deg,low,upp,xlead,xcffs);
    put("The vector of leading exponents :"); put(xlead,1); new_line;
    put("A "); put(ncols,1); put_line("-vector of Laurent series :");
    Write(xlead,xcffs,"x");
    Allocate_Series_Coefficients(nrows,deg,bcffs);
    Matrix_Vector_Product(deg,Alead,Acffs,xlead,xcffs,blead,bcffs);
    put_line("The product of the matrix with the vector :");
    Write(blead,bcffs,"b");
    Allocate_Series_Coefficients(ncols,deg,ycffs);
    if lower then
      Forward_Substitution(deg,Alead,Acffs,blead,bcffs,ylead,ycffs);
      put_line("The computed solution :");
      Write(ylead,ycffs,"y");
    elsif upper then
      Backward_Substitution(deg,Alead,Acffs,blead,bcffs,ylead,ycffs);
      put_line("The computed solution :");
      Write(ylead,ycffs,"y");
    else
      new_line;
      put("Apply pivoting in the LU factorization ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Test_LU_Factorization(nrows,ncols,deg,Alead,Acffs);
       else Test_Plain_LU_Factorization(nrows,ncols,deg,Alead,Acffs);
      end if;
    end if;
  end Test;

  procedure Specific_Test is

    deg : constant integer32 := 2;
    Alead,Llead,Ulead : Standard_Integer_Matrices.Matrix(1..2,1..2);
    Acffs,Lcffs,Ucffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    Arow : Standard_Complex_VecVecs.Link_to_VecVec;
    Acff : Standard_Complex_Vectors.Link_to_Vector;
    pivots : Standard_Integer_Vectors.Vector(1..2);

  begin
    Alead(1,1) := 1; Alead(1,2) := 0;
    Alead(2,1) := 1; Alead(2,2) := 0;
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,2,1,2,0,deg);
    Arow := Acffs(1);
    Acff := Arow(1); -- Acff is A(1,1) = t
    Acff(0) := Standard_Complex_Numbers.Create(1.0);
    for k in 1..deg loop
      Acff(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    Acff := Arow(2); -- Acff is A(1,2) = 1
    Acff(0) := Standard_Complex_Numbers.Create(1.0);
    for k in 1..deg loop
      Acff(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    Arow := Acffs(2);
    Acff := Arow(1); -- Acff is A(2,1) = t
    Acff(0) := Standard_Complex_Numbers.Create(1.0);
    for k in 1..deg loop
      Acff(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    Acff := Arow(2); -- Acff is A(2,2) = 0
    Acff(0) := Standard_Complex_Numbers.Create(0.0);
    for k in 1..deg loop
      Acff(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    Write(Alead,Acffs,"A");
    Standard_Complex_VecVecVecs.Allocate(Lcffs,1,2,1,2,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Ucffs,1,2,1,2,0,deg);
    LU_Factorization(2,2,deg,Alead,Acffs,pivots);
    put("The pivots :"); put(pivots); new_line;
    Lower_Triangular_Part(2,2,deg,Alead,Acffs,Llead,Lcffs);
    Upper_Triangular_Part(2,2,deg,Alead,Acffs,Ulead,Ucffs);
    Write(Llead,Lcffs,"L");
    Write(Ulead,Ucffs,"U");
  end Specific_Test;

  function Seed_Prompt return integer32 is

    ans : character;
    seed : integer32 := 0;

  begin
    new_line;
    put("Fixed seed ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      seed := Standard_Random_Numbers.Get_Seed;
    else
      put("Give the seed : "); get(seed);
      Standard_Random_Numbers.Set_Seed(natural32(seed));
    end if;
    return seed;
  end Seed_Prompt;

  procedure Determinant_Test is

    deg,low,upp,seed : integer32 := 0;
    Alead,Llead,Ulead : Standard_Integer_Matrices.Matrix(1..2,1..2);
    Acffs,Lcffs,Ucffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    Arow1,Arow2 : Standard_Complex_VecVecs.Link_to_VecVec;
    pivots : Standard_Integer_Vectors.Vector(1..2);
    xdetlead,left,right : integer32 := 0;
    xdetcff0 : Complex_Number;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the lower bound on the leading exponent : "); get(low);
    put("Give the upper bound on the leading exponent : "); get(upp);
    seed := Seed_Prompt;
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,2,1,2,0,deg);
    Random_Matrix(2,2,low,upp,Alead,Acffs);
    Write(Alead,Acffs,"A");
    put_line("The leading exponents : "); put(Alead);
   -- xdetlead := min(Alead(1,1) + Alead(2,2),Alead(2,1) + Alead(1,2));
    left := Alead(1,1) + Alead(2,2);
    right := Alead(2,1) + Alead(1,2);
    put("Leading exponent of determinant : "); put(xdetlead,1); new_line;
    Arow1 := Acffs(1);
    Arow2 := Acffs(2);
    if left = right then
      xdetcff0 := Arow1(1)(0)*Arow2(2)(0) - Arow2(1)(0)*Arow1(2)(0);
      xdetlead := left;
    elsif left < right then 
      xdetcff0 := Arow1(1)(0)*Arow2(2)(0);
      xdetlead := left;
    else
      xdetcff0 := -Arow2(1)(0)*Arow1(2)(0);
      xdetlead := right;
    end if;
    Standard_Complex_VecVecVecs.Allocate(Lcffs,1,2,1,2,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Ucffs,1,2,1,2,0,deg);
    LU_Factorization(2,2,deg,Alead,Acffs,pivots);
    put_line("The matrix after the LU factorization :");
    Write(Alead,Acffs,"A");
    put("The pivots :"); put(pivots); new_line;
    Lower_Triangular_Part(2,2,deg,Alead,Acffs,Llead,Lcffs);
    Upper_Triangular_Part(2,2,deg,Alead,Acffs,Ulead,Ucffs);
    Write(Llead,Lcffs,"L");
    Write(Ulead,Ucffs,"U");
    new_line;
    put_line("Now we do the determinant test ...");
    declare
      ydetlead : integer32 := 0;
      ydetcff : Standard_Complex_Vectors.Vector(0..deg);
    begin
      Double_Laurent_Series.Multiply
        (deg,Ulead(1,1),Ulead(2,2),Ucffs(1)(1).all,Ucffs(2)(2).all,
         ydetlead,ydetcff);
      put("xdetlead : "); put(xdetlead,1); new_line;
      put("ydetlead : "); put(ydetlead,1); 
      if xdetlead /= ydetlead then
        put_line("  mismatched leading exponents of determinant!");
      else
        put_line("  matching leading exponents of determinant.");
        put("xdet(0) : "); put(xdetcff0); new_line;
        put("ydet(0) : "); put(ydetcff(0)); new_line;
      end if;
    end;
    new_line;
    put("The seed used : "); put(seed,1); new_line;
  end Determinant_Test;

  procedure Random_Test is

    deg,nrows,ncols,low,upp,seed : integer32 := 0;
    ans : character;
    lower,upper : boolean;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the lower bound on the leading exponent : "); get(low);
    put("Give the upper bound on the leading exponent : "); get(upp);
    put("Give the number of rows : "); get(nrows);
    put("Give the number of columns : "); get(ncols);
    new_line;
    if nrows /= ncols then
      lower := false;
    else
      put("Lower triangular matrix ? (y/n) "); Ask_Yes_or_No(ans);
      lower := (ans = 'y');
      if lower then
        upper := false;
      else
        put("Upper triangular matrix ? (y/n) "); Ask_Yes_or_No(ans);
        upper := (ans = 'y');
      end if;
    end if;
    seed := Seed_Prompt;
    Test(nrows,ncols,deg,low,upp,lower,upper);
    new_line;
    put("The seed used : "); put(seed,1); new_line;
  end Random_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for testing the Laurent matrix operations :");
    put_line("  1. LU factorization of a specific 2-by-2 matrix");
    put_line("  2. determinant of a random 2-by-2 matrix");
    put_line("  3. LU factorization of a random general matrix");
    put("Type 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Specific_Test;
      when '2' => Determinant_Test;
      when '3' => Random_Test;
      when others => null;
    end case;
  end Main;

end Test_Double_Lseries_Matrices;
