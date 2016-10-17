with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vecvecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Poly_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Brackets;                          use Brackets;

package Moving_Flag_Homotopies is

-- DESCRIPTION :
--   A moving flag homotopy introduces new coefficients in a general flag
--   as defined by a generalizing sequence in a checker game.
--   Different versions of the procedures with the same name are for
--   different precision: standard double, double double, or quad double;
--   and either stay silent or write extra output to file.

  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy with one flag,
  --   in standard double, double double, or quad double precision,
  --   without intermediate output to file.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection condition for the fixed flag;
  --   vf      coordinates of a general flag to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy with one flag,
  --   in standard double, double double, or quad double precision,
  --   with intermediate output written to file.

  -- ON ENTRY :
  --   file    for intermediate output and diagnostics;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection condition for the fixed flag;
  --   vf      coordinates of a general flag to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in DoblDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in QuadDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space x,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   x       represents a moving k-plane in n-space;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in DoblDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in QuadDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space x,
  --   using a more efficient representation for the Schubert problem,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   x       represents a moving k-plane in n-space;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a degenerate moving flag, equal to the identity matrix.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a degenerate moving flag, equal to the identity matrix,
  --   using a more efficient formulation for the Schubert problem,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in Standard_Complex_Matrices.Matrix;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in DoblDobl_Complex_Matrices.Matrix;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in QuadDobl_Complex_Matrices.Matrix;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a particular moving flag, computed
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   mf      coordinates of the moving flag;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf, where the solution
  --           plane x is multiplied by mf, as mf*x.

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in Standard_Complex_Matrices.Matrix;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in DoblDobl_Complex_Matrices.Matrix;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in QuadDobl_Complex_Matrices.Matrix;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a particular moving flag,
  --   using a more efficient formulation for the Schubert problem,
  --   in standard double, double double, or quad double precision.
  --   These versions are silent, no output is written.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   mf      coordinates of the moving flag;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf, where the solution
  --           plane x is multiplied by mf, as mf*x.

  procedure Flag_Conditions
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a flag intermediate between the identity and opposite flag,
  --   in standard double, double double, or quad double precision,
  --   with intermediate output written to file.

  -- ON ENTRY :
  --   file    for diagnostic output;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed;
  --   mf      coordinates of the opposite flag;
  --   nf      current accumulated flag, moving to the opposite flag.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a particular moving flag,
  --   using a more efficient formulation for the Schubert problem,
  --   in standard double, double double, or quad double precision.
  --   These versions write diagnostic output to file.

  -- ON ENTRY :
  --   file    for diagnostic output;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   mf      coordinates of the moving flag;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf, where the solution
  --           plane x is multiplied by mf, as mf*x.

  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy for a moving flag,
  --   capable to keep track of multiple intersection conditions,
  --   in standard double, double double, or quad double precision.
  --   These versions are silent, do not write diagnostic output.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy for a moving flag,
  --   capable to keep track of multiple intersection conditions.

  -- ON ENTRY :
  --   file    for intermediate output and diagnostics;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns start*(1-t) + target*t for continuation parameter t
  --   as a matrix of nv+1 unknowns, where t is at the last position,
  --   in standard double, double double, or quad double precision.

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in Standard_Complex_Matrices.Matrix;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the flag conditions for a k-plane in n-space with
  --   localization pattern determined by (p,rows,cols) that has
  --   to meet the Cheater homotopy flag under the conditions cond,
  --   in standard double, double double, or quad double precision.

  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f all polynomial equations that must be satisfied
  --   to meet all conditions by the brackets in cond,
  --   in standard double, double double, or quad double precision,
  --   for the matrices in start and target.

  -- REQUIRED : conds'range = start'range = target'range.

end Moving_Flag_Homotopies;
