with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural64_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vecvecs;
with Standard_Complex_Matrices;        use Standard_Complex_Matrices;
with Standard_Complex_VecMats;         use Standard_Complex_VecMats;
with Standard_Jacobian_Trees;          use Standard_Jacobian_Trees;

package Standard_Deflation_Matrices is

-- DESCRIPTION :
--   The operations provided in this package assign values to blocks
--   in the Jacobian matrices which arise in the deflation algorithm.

-- STORAGE THRESHOLD :
--   Deflation has the tendency to generate rather large matrices
--   and the GNAT GPL 2009 compiler has the unfortunate behaviour
--   to run in an infinite loop when matrices get too large.
--   As a crude solution, we check for a storage threshold:

  storage_threshold : constant natural32 := 2**16;

  function Number_of_Columns
              ( d,nv,R1 : Standard_Natural_Vectors.Vector;
                m : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of columns in the result of the
  --   application of the derivative operator d to the matrix A(m).

  -- ON ENTRY :
  --   d        derivative operator;
  --   nv       nv(i) equals the number of variables in A(i);
  --   R1       R1(i) is number of multiplier variables in A(i);
  --   m        current deflation matrix where d is applied to.

  -- REQUIRED : m <= nv'last = R1'last.

  function Zero_Matrix ( rows,cols : natural32 ) return Matrix;

  -- DESCRIPTION :
  --   Raises STORAGE_ERROR if rows*cols exceeds the storage threshold,
  --   else returns a zero matrix with #rows in rows and #columns in cols.

  generic
    with procedure Loop_Body ( index : in Standard_Natural_Vectors.Vector );
  procedure Multi_Loop
              ( s : in natural32; d,n : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Executes an s-loop, consisting of s loops nested into each other.
  --   The s is the sum of the elements in the derivative operator d.
  --   The number of iterations in the i-th loop is determined by n(i).
  --   The index in the Loop_Body is an s-vector, with the values of the
  --   the running variables.

  -- REQUIRED : Sum(d) = s, d'range = n'range, index'range = 1..s.

-- MULTIPLIERS :

  function Multiply ( B : VecMat; x : Standard_Complex_VecVecs.VecVec )
                    return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the result of the multiplication of the matrices in B
  --   with all the multipliers in x.

  procedure Multiply ( A : in out Link_to_Matrix; r,c : in integer32;
                       JM,B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Does an inline multiplication of the matrix JM with B,
  --   starting assignments at row r and column c in the matrix A.

  procedure Multiply ( A : in out Link_to_Matrix; r,c : in integer32;
                       JM : in Link_to_Matrix;
                       Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Assigns the result of JM*Bl into the matrix A,
  --   starting at row r, and into column c.

  procedure Multiply ( A : in out Link_to_Matrix; r,c,m : in integer32;
                       JM : in Link_to_Matrix;
                       Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Assigns the result of JM*Bl into the matrix A,
  --   starting at row r, and into m columns starting at c.

  procedure Alternating_Permute
              ( A : in out Link_to_Matrix;
                row,col,nbrows,nbcols,nbvars : in integer32 );

  procedure Alternating_Permute
              ( file : in file_type; A : in out Link_to_Matrix;
                row,col,nbrows,nbcols,nbvars : in integer32 );

  -- DESCRIPTION :
  --   This routine must be executed after multiplication of Jacobian
  --   matrices with a matrix B, in case it represents the derivatives
  --   with respect to lambda.

  -- ON ENTRY :
  --   A        deflation matrix;
  --   row      starting row in A for the data to permute;
  --   col      starting column in A for the data to permute;
  --   nbrows   number of rows in A which will be permuted;
  --   nbcols   number of columns in A which will be permuted;
  --   nbvars   number of multiplier variables.

  -- ON RETURN :
  --   A        deflation matrix with columns alternatingly permuted.

  procedure One_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  procedure One_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  procedure Multi_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k,s : in integer32;
                d,R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  procedure Multi_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k,s : in integer32;
                d,R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs the right multiplication on the deflation matrix JM
  --   with the vector Bl and stores the result in A, starting at (r,c).
  --   A deflation matrix is a derivative with respect to the multipliers.
  --   This multiplication respects the 2-column block format of JM.
  --   One_* assumes is that JM is a first-order derivative.
  --   Multi_* assumes is that JM is a second-order derivative.

  -- ON ENTRY :
  --   file     to write diagnostics on;
  --   A        Jacobian matrix A(m) at stage m in the deflation;
  --   m        number of column blocks in A(m);
  --   r        start row to position result of the multiplication;
  --   c        start column for result of multiplication;
  --   k        JM is k-th child in the tree;
  --   s        order of the derivative, is omitted, then s = 1;
  --   d        derivative operator of JM;
  --   R0       R0(i) is number of multipliers at stage i in deflation,
  --            R0(0) is number of original variables;
  --   JM       Jacobian matrix for right-multiplication with Bl;
  --   Bl       multiplication of B(m) with lambda(m).

  -- ON RETURN :
  --   A        matrix with updated entries.

  procedure One_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix );

  procedure One_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Similar as the One_Right_Multiply_Deflation from above,
  --   except that instead of a vector, we multiply with a matrix B.

  procedure One_Right_Multiply
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  procedure One_Right_Multiply
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs the right multiplication of Bl with a Jacobian matrix
  --   and stores the results in A, starting at position (r,c).
  --   One_* assumes JM is a first-order derivative matrix.
  --   Two_* assumes JM is a second-order derivative matrix.

  -- ON ENTRY :
  --   file     to write diagnostics or intermediate results to;
  --   A        Jacobian matrix A(m) at stage m in the deflation;
  --   m        number of column blocks in A(m);
  --   r        start row to position result of the multiplication;
  --   c        start column for result of multiplication;
  --   k        JM is k-th child in the tree;
  --   R0       R0(i) is number of multipliers at stage i in deflation,
  --            R0(0) is number of original variables;
  --   JM       Jacobian matrix for right-multiplication with Bl;
  --   Bl       multiplication of B(m) with lambda(m).

  -- ON RETURN :
  --   A        matrix with updated entries.

  procedure One_Right_Multiply
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix );

  procedure One_Right_Multiply
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Similar to the other "Right_Multiply", except we use a matrix B,
  --   instead of a vector to multiply the Jacobian matrices with.

-- ASSIGNMENTS :

  procedure Assign_Scaling_Coefficients
              ( A : in out Link_to_Matrix;
                h : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Assigns the values of h into the lower right corner of A.

  procedure Assign_Lower_Jacobian_Matrices
              ( A : in out Link_to_Matrix; row,col : in integer32;
                jms : in Link_to_VecMat;
                Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Assigns values of the Jacobian matrices in jms to A, after
  --   multiplication with Bl, starting at position (row,col) in A.

  -- REQUIRED : restricted to first-order (or lower) derivatives in jms.

  procedure Assign_Lower_Jacobian_Matrices
              ( A : in out Link_to_Matrix; row,col : in integer32;
                jms : in Link_to_VecMat; B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Assigns values of the Jacobian matrices in jms, after
  --   multiplication with B, starting at position (row,col) of A.

  -- REQUIRED : restricted to first-order (or lower) derivatives in jms.

  procedure Assign_from_Jacobian_Matrices
              ( A : in out Link_to_Matrix; col : in out integer32;
                jms : in Link_to_VecMat; nbcols : in natural32 );

  -- DESCRIPTION :
  --   Assigns values of the Jacobian matrices in jms to A,
  --   starting at the top row, and column number equal to col.
  --   The number of columns in a nonvoid Jacobian matrix is nbcols.
  --   On return is the next free column in A.

  -- REQUIRED : jms contains either the original Jacobian matrix,
  --            or its first derivatives.

  procedure Assign_from_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Assigns the values of the Jacobian matrices in jms to A,
  --   starting at the given row, and column number col.
  --   Before assigning, the matrices are multiplied with B.
  --   On return is the next free column in A.

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 k,n : in integer32 );

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 k,n : in integer32 );

  -- DESCRIPTION :
  --   Assign the values of all Jacobian matrices stored in jms into A,
  --   the matrices are stored next to each other.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   A         deflation matrix;
  --   row       row where to start placing the values in A;
  --   col       column to start placing the values in A;
  --   jms       derivatives of Jacobian matrix of order k;
  --   monkeys   hash table of monomial codes;
  --   k         order of the derivative of the Jacobian;
  --   n         number of original variables in the system.

  -- ON RETURN :
  --   A         matrix with all Jacobian matrices of order k
  --             assigned to it, starting at (row,col);
  --   col       updated column index.

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector;
                 k,n : in integer32 );

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector;
                 k,n : in integer32 );

  -- DESCRIPTION :
  --   Multiplies all Jacobian matrices with the vector Bl before
  --   assigning them to the proper places in the matrix A.

  -- ON ENTRY :
  --   file      to write diagnostics and intermediate output;
  --   A         deflation matrix;
  --   row       row where to start placing the values in A;
  --   col       column to start placing the values in A;
  --   jms       derivatives of Jacobian matrix of order k;
  --   monkeys   hash table of monomial codes;
  --   Bl        productive of random matrix with multiplier vector;
  --   k         order of the derivative of the Jacobian;
  --   n         number of original variables in the system.

  -- ON RETURN :
  --   A         matrix with all Jacobian matrices of order k
  --             multiplied with Bl assigned, starting at (row,col);
  --   col       updated column index.

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 B : in Link_to_Matrix; k,n : in integer32 );

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 B : in Link_to_Matrix; k,n : in integer32 );

  -- DESCRIPTION :
  --   Multiplies all Jacobian matrices with the matrix B before
  --   assigning them to the proper places in the matrix A.

  -- ON ENTRY :
  --   file      to write diagnostics and intermediate output;
  --   A         deflation matrix;
  --   row       row where to start placing the values in A;
  --   col       column to start placing the values in A;
  --   jms       derivatives of Jacobian matrix of order k;
  --   monkeys   hash table of monomial codes;
  --   B         random matrix to multiply Jacobians with;
  --   k         order of the derivative of the Jacobian;
  --   n         number of original variables in the system.

  -- ON RETURN :
  --   A         matrix with all Jacobian matrices of order k
  --             multiplied with B assigned, starting at (row,col);
  --   col       updated column index.

  procedure Assign_Root_Child
               ( A : in out Link_to_Matrix; d0,m : in natural32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix );

  procedure Assign_Root_Child
               ( file : in file_type;
                 A : in out Link_to_Matrix; d0,m : in natural32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Assigns the first child v0 to the proper position in A(m),
  --   where A is the matrix at the m-th stage in the deflation.
  --   This is a subroutine of the second Assign_Children,
  --   see the explanation of the parameters below.

  -- REQUIRED : v0 /= null.

  procedure Assign_Jacobian_Child
               ( A : in out Link_to_Matrix; m : in integer32; 
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector;
                 ind : in integer32; col : in out integer32 );

  procedure Assign_Jacobian_Child
               ( file : in file_type;
                 A : in out Link_to_Matrix; m : in integer32; 
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector;
                 ind : in integer32; col : in out integer32 );

  -- DESCRIPTION :
  --   Assigns into A a Jacobian matrix of the original system.

  procedure Assign_Children
               ( A : in out Link_to_Matrix; m : in integer32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v : in VecMat; B : in Link_to_Matrix );

  procedure Assign_Children
               ( file : in file_type;
                 A : in out Link_to_Matrix; m : in integer32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v : in VecMat; B : in Link_to_Matrix );

  -- DESCRIPTION :
  --   Multiplies the matrices in v with the matrix B
  --   and assigns the results of this multiplication into the matrix A.

  procedure Assign_Children
               ( A : in out Link_to_Matrix; nq,d0,m,order : in natural32;
                 chtp : in Standard_Natural_VecVecs.VecVec;
                 R0 : in Standard_Natural_Vectors.Vector; v : in VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector );

  procedure Assign_Children
               ( file : in file_type;
                 A : in out Link_to_Matrix; nq,d0,m,order : in natural32;
                 chtp : in Standard_Natural_VecVecs.VecVec;
                 R0 : in Standard_Natural_Vectors.Vector; v : in VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   The vector v contains the values of the deflation matrices
  --   of the children 0 to m (v'range = 0..m).  We assign these
  --   values to the proper position in the matrix A.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   A         Jacobian matrix of deflation at stage m;
  --   nq        number of equations in original system;
  --   d0        order of the derivative in the original variables;
  --   m         number of the current deflation stage;
  --   order     order of the derivatives of the children;
  --   chtp      derivative operator for each child;
  --   R0        number of variables added in every deflation stage,
  --             R0(0) = number of original variables,
  --             R0(i) = number of multipliers added at stage i;
  --   v         values of children of derivative with the m groups;
  --   monkeys   monomial hashing keys to derivatives of Jacobian matrices;
  --   jrt       remember table with all derivates of Jacobian matrices;
  --   B         random matrix B(m) used in m-th deflation stage;
  --   Bl        product of B(m) with lambda(m), the m-th multipliers.

  -- ON RETURN :
  --   A         matrix updated with values v into A.

end Standard_Deflation_Matrices;
