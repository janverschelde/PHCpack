with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Brackets;                          use Brackets;
with Localization_Posets;               use Localization_Posets;

package Deformation_Posets is

-- DESCRIPTION :
--   This package is the continuous analogue of "Localization_Posets".

-- DATA STRUCTURES :

  type Array_of_VecMats is array ( integer32 range <> ) of Link_to_VecMat;
  type Link_to_Array_of_VecMats is access Array_of_VecMats;
  type Array_of_Array_of_VecMats is
    array ( integer32 range <> ) of Link_to_Array_of_VecMats;

  type Duration_Array is array ( integer32 range <> ) of duration;

-- CREATORS :

  function Create ( index_poset : Array_of_Array_of_Nodes )
                  return Array_of_Array_of_VecMats;

  -- DESCRIPTION :
  --   The array on return mirrors the index poset in the following sense:
  --     res(i)(j)'length = poset(i)(j).roco, if roco > 0,
  --     res(i)(j) = null, if roco = 0.

-- SELECTOR :

  function Empty ( poset : Array_of_Array_of_VecMats;
                   level,label : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if poset(level)(label) is a vector of null pointers.

-- ANALOGUES TO ROOT COUNTERS :

  procedure Solve ( n : in natural32;
                    poset : in out Array_of_Array_of_VecMats;
                    nd : in Node; planes : in VecMat;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );
  procedure Solve ( file : in file_type; n : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Computes the p-planes that intersect with the given m-planes.
  --   This is the solver for hypersurface intersection conditions.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input planes in general position;
  --   report       indicates if reporting path tracker is needed;
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains all p-planes that
  --                intersect nontrivially with planes(1..nd.level);
  --   npaths       number of paths traced at each level;
  --   timings      CPU user time at each level.

  procedure One_Solve
                  ( file : in file_type; n : in natural32; cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Computes the p-planes that intersect with the given planes of
  --   dimension (m+1-cod(i)), for i in cod'range.
  --   There is only one moving equation in the homotopies.

  -- REQUIRED :
  --   The localization poset has been built either by consistently
  --   incrementing the top or by decrementing the bottom pivots.
  --   The co-dimension conditions are not interlaced.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   cod          co-dimensions of input planes, dim(L_i) = m+1-cod(i)
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input planes L_i in general position;
  --   report       indicates if reporting path tracker is needed;
  --   outlog       flag to write homotopies on file if set to true. 

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains all p-planes that
  --                intersect nontrivially with planes(1..nd.level);
  --   npaths       number of paths traced at each level;
  --   timings      CPU user time at each level.

  procedure Solve ( file : in file_type; n : in natural32; cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Computes the p-planes that intersect with the given planes of
  --   dimension (m+1-cod(i)), for i in cod'range.
  --   The co-dimension conditions are interlaced here.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   cod          co-dimensions of input planes, dim(L_i) = m+1-cod(i)
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input planes L_i in general position;
  --   report       indicates if reporting path tracker is needed;
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains all p-planes that
  --                intersect nontrivially with planes(1..nd.level);
  --   npaths       number of paths traced at each level;
  --   timings      CPU user time at each level.

  procedure Solve ( n,q,nb : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );
  procedure Solve ( file : in file_type; n,q,nb : in natural32;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Computes the q-curve that produces p-planes in n-dimensional space
  --   that intersect with the given m-planes at specific s-values.
  --   This is the quantum analogue to the hypersurface solver.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   q            degree of the map;
  --   nb           number of solution maps wanted;
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input m-planes in general position;
  --   s            interpolation points to meet the input m-planes;
  --   report       indicates if reporting path tracker is needed.
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains coefficients
  --                of nb curves that intersect nontrivially with
  --                planes(1..nd.level) at the specified s-values.
  --   npaths       number of paths traced at each level;
  --   timings      CPU user times at each level.

  procedure One_Solve
                  ( file : in file_type;
                    n,q,nb : in natural32; cnt : in out natural32;
                    cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   This is the quantum analogue to the other "One_Solve".
  --   Computes the p-planes that intersect with the given planes of
  --   dimension (m+1-cod(i)), for i in cod'range.
  --   There is only one moving equation in the homotopies.

  -- REQUIRED :
  --   The localization poset has been built either by consistently
  --   incrementing the top or by decrementing the bottom pivots.
  --   The co-dimension conditions are not interlaced.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   q            degree of the map;
  --   nb           number of maps wanted;
  --   cnt          counts the number of maps already computed;
  --   cod          co-dimensions of input planes, dim(L_i) = m+1-cod(i)
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input planes L_i in general position;
  --   s            interpolation points to meet the input planes;
  --   report       indicates if reporting path tracker is needed;
  --   outlog       flag to write homotopies on file if set to true. 

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains all p-planes that
  --                intersect nontrivially with planes(1..nd.level);
  --   npaths       number of paths traced at each level;
  --   timings      CPU user time at each level.

  procedure Solve ( file : in file_type;
                    n,q,nb : in natural32; cnt : in out natural32;
                    cod : in Bracket;
                    poset : in out Array_of_Array_of_VecMats; nd : in Node;
                    planes : in VecMat; s : in Standard_Complex_Vectors.Vector;
                    report,outlog : in boolean;
                    npaths : in out Standard_Natural_Vectors.Vector;
                    timings : in out Duration_Array );

  -- DESCRIPTION :
  --   Computes the q-curve that produces p-planes in n-dimensional space
  --   that intersect with the given m-planes at specific s-values.
  --   This is the quantum analogue to the general solver.

  -- ON ENTRY :
  --   file         output file for intermediate results and logs;
  --   n            dimension of the working space, equals m+p;
  --   q            degree of the map;
  --   nb           number of maps wanted;
  --   cnt          counts the number of maps already computed;
  --   cod          co-dimensions of input planes, dim(L_i) = m+1-cod(i);
  --   poset        deformation poset, properly created from index poset;
  --   nd           target root node from the localization poset;
  --   planes       input m-planes in general position;
  --   s            interpolation points to meet the input planes;
  --   report       indicates if reporting path tracker is needed.
  --   outlog       flag to write homotopies on file if set to true.

  -- ON RETURN :
  --   poset        poset(nd.level)(nd.label) contains coefficients
  --                of all curves that intersect nontrivially with
  --                planes(1..nd.level) at the specified s-values.
  --   npaths       number of paths traced at each level;
  --   timings      CPU user times at each level.

-- DESTRUCTORS :

  procedure Clear ( avm : in out Array_of_VecMats );
  procedure Clear ( avm : in out Link_to_Array_of_VecMats );
  procedure Clear ( avm : in out Array_of_Array_of_VecMats );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory.

end Deformation_Posets;
