with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions;

package Mixed_Volume_Computation is

-- DESCRIPTION :
--   This package offers a number of routines for the computation
--   of the mixed volume of a system of polynomial equations.

-- UTILITIES :

  procedure Compute_Mixture ( supports : in out Array_of_Lists;
                              mix,perms : out Link_to_Vector ); 
  -- DESCRIPTION : 
  --   Computes the type of mixture of the supports of a system.

  -- ON ENTRY :
  --   supports  the supports of a polynomial system.

  -- ON RETURN :
  --   supports  a permuted array of supports, so that the same
  --             supports stand all toghether;
  --   mix       mix(k) indicates number of occurrencies of the kth support;
  --   perms     perms(k) gives the place of the kth support,
  --             after permutation to make supports correspond with mix.

  function Compute_Index ( k : integer32; mix : Vector ) return integer32;

  -- DESCRIPTION :
  --   Given k, an entry in the supports, the number this function returns
  --   indicates the number of different support, w.r.t. the type of  mixture.

  function Compute_Permutation
               ( n : integer32; mix : Vector; supports : Array_of_Lists )
               return Link_to_Vector;

  -- DESCRIPTION :
  --   Given the type of mixture and the support, the permutation vector
  --   will be computed.

  -- ON RETURN :
  --   perms     perms(k) gives the place of the kth support,
  --             after permutation to make supports correspond with mix.

  function Typed_Lists ( mix : Vector; points : Array_of_Lists )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns a tuple of lists where each list occurs only once,
  --   according to the given type of mixture.

  function Permute ( p : Poly_Sys; perm : Link_to_Vector ) return Poly_Sys;
  function Permute ( p : Laur_Sys; perm : Link_to_Vector ) return Laur_Sys;
  function Permute ( supports : Array_of_Lists ; perm : Link_to_Vector )
                   return Array_of_Lists;

  -- DESCRIPTION :
  --   Permutes the polynomials in the system or the supports,
  --   according to the vector perm.

-- COMPUTE MIXED SUBDIVISION :

  procedure Mixed_Coherent_Subdivision
               ( n : in integer32; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean; low,upp : in Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision );

  -- DESCRIPTION :
  --   Given a set of points and a lifting function,
  --   a subdivision of the polytope will be computed.

  -- ON ENTRY :
  --   n         the dimension of the vector space;
  --   mix       mix(k) indicates how many times the kth point set occurs;
  --   points    an array of all different point sets;
  --   linear    indicates wether a linear lifting should be used;
  --   lift      an array of lifting polynomials or an m dimensional
  --             array of vectors, where the length of the kth vector
  --             must equal the length of the kth support,
  --             when nonlinear, otherwise the length equals n.
  --   low,upp   lower and upper bounds for random lifting.

  -- ON RETURN :
  --   lifted    the lifted points which can later be used for lifting
  --             the polynomial system;
  --   nbsucc    the number of successful face-face combinations that
  --             have been computed;
  --   nbfail    the number of unsuccessful face-face combinations;
  --   mixsub    the mixed subdivision of the polytope, defined as
  --             lower hull of the lifted points.

-- MIXED VOLUME COMPUTATION, GIVEN A SUBDIVISION :

  function Mixed_Volume
               ( n : integer32; mix : Vector;
                 mic : Integer_Mixed_Subdivisions.Mixed_Cell;
                 multprec_hermite : boolean := false )
               return natural32;
  function Mixed_Volume 
               ( n : integer32; mix : Vector;
                 mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 multprec_hermite : boolean := false )
               return natural32;

  -- DESCRIPTION :
  --   Computes the mixed volume based on a mixed cell and subdivision.
  --   When the cells are not fine enough, they will be refined but will
  --   be lost after returning the result.

  procedure Mixed_Volume 
               ( n : in integer32; mix : in Vector; 
                 mic : in out Integer_Mixed_Subdivisions.Mixed_Cell;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );
  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector;
                 mic : in out Floating_Mixed_Subdivisions.Mixed_Cell;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );
  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector; 
                 mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );
  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector; 
                 mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Computes the mixed volume based on a mixed cell and subdivision.
  --   When the cells are not fine enough, they will be refined by lifting.
  --   The refinement is stored in the subdivision field of the cells.

-- MIXED VOLUME COMPUTATIONS, GIVEN THE SUPPORTS :

  function Mixed_Volume ( n : integer32; supports : Array_of_Lists )
			return natural32;

  function Mixed_Volume ( file : file_type; n : integer32;
                          supports : Array_of_Lists ) return natural32;

  function Mixed_Volume ( n : integer32; mix : Vector;
                          supports : Array_of_Lists ) return natural32;

  function Mixed_Volume ( file : file_type; n : integer32; mix : Vector;
                          supports : Array_of_Lists ) return natural32;

  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector;
                 supports : in Array_of_Lists; lifted : out Array_of_Lists;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );

  procedure Mixed_Volume
               ( file : in file_type; n : in integer32; mix : in Vector;
                 supports : in Array_of_Lists;
                 lifted : out Array_of_Lists;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   All these routines compute the mixed volume of support lists.

  -- ON ENTRY :
  --   file      if specified, then the mixed subdivision will be 
  --             written on file;
  --   n         the dimension of the system;
  --   mix       mix(k) is the number of times the kth support occurs;
  --   supports  the supports of a system of n polynomials in n unknowns;
  --   multprec_hermite indicates whether multiprecision arithmetic must
  --             be use to compute the Hermite normal form.

  -- ON RETURN :
  --   lifted    array of listed points;
  --   mixsub    mixed subdivision used;
  --   mv        the mixed volume.

end Mixed_Volume_Computation;
