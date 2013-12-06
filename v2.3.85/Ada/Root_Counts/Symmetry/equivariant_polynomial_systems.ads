with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Symmetry_Group;                     use Symmetry_Group;

package Equivariant_Polynomial_Systems is

-- DESCRIPTION :
--   This package contains procedures for, given a group representation V,
--   to compute the associated representation W, of a (G,V,W)-symmetric
--   polynomial system.

  procedure Act ( v : in List_of_Permutations; s : in Poly_Sys;
                  w : in out List_of_Permutations;
                  fail,inva,equi : out boolean );

  -- DESCRIPTION :
  --   Each permutation of the list v will be applied on the system s;
  --   the list w contains the results of each permutation.

  -- ON ENTRY :
  --   v          a list of permutations;
  --   s          a polynomial system.

  -- ON RETURN :
  --   w          a list of Natural_Vectors x
  --               x(i) = j, where j indicates the index of the
  --               resulting polynomial of s, after permutation,
  --              if j = n+1, then the permuted polynomial did not belong to s
  --              and fail will be true on return;
  --   fail       true if the system is not (G,V,W)-symmetric, false otherwise;
  --   inva       true, if every polynomial in the system remains invariant,
  --              i.e.: Permute(s(i),p) = s(i), with i in s'range,
  --              for every permutation p, false otherwise;
  --   equi       true, if v = w, false otherwise.

  function Symmetric ( s : Poly_Sys; v,w : List_of_Permutations )
                     return boolean;

  -- DESCRIPTION :
  --   This routine returns true if s is (G,V,W)-symmetric,
  --   it returns false when s is (G,V,W)-symmetric.

  -- ON ENTRY :
  --   s          a polynomial system;
  --   v,w        representations of the group.

end Equivariant_Polynomial_Systems;
