with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;

package Supported_Subsystems is

-- DESCRIPTION :
--   This package collects operations to select a subsystem from a
--   polynomial system, given a subset of its supports.
--   The supports may be given as integer or floating-point lists.

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Natural_Vectors.Vector )
                    return boolean;
  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Integer_Vectors.Vector )
                    return boolean;

  -- DESCRIPTION :
  --   The vectors v1 and v2 are equal if and only if all elements of
  --   the vector v2 converted into double_float equal v1.
  --   If v1 is longer than v2, then these elements will be ignored.

  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Natural_Vectors.Vector ) return boolean;
  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector v belongs to the list s.

  function Select_Terms ( p : Standard_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns those terms of the polynomial p whose exponent vectors
  --   occur in the list s.

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns those terms of the polynomials in p whose exponent vectors
  --   occur in the lists s.  If the type of mixture m is omitted,
  --   then s'range must equal p'range, otherwise the type of mixture
  --   is taken into account when selecting terms from p.

end Supported_Subsystems;
