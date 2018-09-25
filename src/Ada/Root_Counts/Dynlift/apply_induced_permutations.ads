with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Apply_Induced_Permutations is

  procedure Compute_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the induced permutation for p.

  -- ON ENTRY :
  --   sup      supports of a polynomial system;
  --   r        number of distinct supports;
  --   perm     permutation of the supports;
  --   mtype    type of mixture, computed by MixedVol;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   p        permuted polynomials.

  procedure Compute_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the induced permutation for p.

  -- ON ENTRY :
  --   sup      supports of a polynomial system;
  --   stlb     stable lifting bound;
  --   r        number of distinct supports;
  --   perm     permutation of the supports;
  --   mtype    type of mixture, computed by MixedVol;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   p        permuted polynomials.

  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Computes the induced permutation and applies it to p.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   r        number of distinct supports;
  --   perm     permutation of the supports;
  --   mtype    type of mixture, computed by MixedVol;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   p        permuted polynomials.

-- WITH STABLE LIFTING BOUND :

  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );
  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Computes the induced permutation and applies it to p.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   stlb     stable lifting bound;
  --   r        number of distinct supports;
  --   perm     permutation of the supports;
  --   mtype    type of mixture, computed by MixedVol;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   p        permuted polynomials.

end Apply_Induced_Permutations;
