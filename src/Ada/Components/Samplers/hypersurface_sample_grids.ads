with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Standard_Stacked_Sample_Grids;     use Standard_Stacked_Sample_Grids;

package Hypersurface_Sample_Grids is

-- DESCRIPTION :
--   This package provide tools to construct grids of points
--   sampled from hypersurfaces.

-- FORMAT CONVERSIONS :

  function Create ( b,v : Vector ) return VecVec;

  -- DESCRIPTION :
  --   Returns the implicit representation of the line b+t*v
  --   as a vector of hyperplanes defining the same line.

  function Create ( b,v : Vector; t : Complex_Number ) return Solution;
  function Create ( b,v : Vector; t : Complex_Number ) return Standard_Sample;

  -- DESCRIPTION :
  --   Given a point on the line b+t*v, the representation of the
  --   sample as a solution or standard sample is returned.

  function Create ( b,v,t : Vector ) return Standard_Sample_List;

  -- DESCRIPTION :
  --   Returns the list of witness points on the line b+t*v,
  --   represented as a list of standard samples.

-- THE SAMPLING MACHINE :

  procedure Initialize ( p : in Poly );

  -- DESCRIPTION :
  --   Initializes the sampling machine with the polynomial p.

  function Sample ( b,v,t : Vector ) return Standard_Sample_List;
  function Sample ( file : file_type; output : boolean;
                    b,v,t : Vector ) return Standard_Sample_List;

  -- DESCRIPTION :
  --   Generates a new random line b1 + t*v1 and computes a set of witness
  --   points on it, starting with the witness points on b + t*v.
  --   Eventually writes intermediate results to file, with diagnostics
  --   of the predictor-corrector stage if output is set to true.

  function Parallel_Sample ( b,v,t : Vector ) return Standard_Sample_List;
  function Parallel_Sample ( file : file_type; output : boolean;
                             b,v,t : Vector ) return Standard_Sample_List;

  -- DESCRIPTION :
  --   Same as Sample, except that now v remains the same.

  function Parallel_Sample1 ( b,v,t : Vector ) return Standard_Sample_List;
  function Parallel_Sample1 ( file : file_type; output : boolean;
                              b,v,t : Vector ) return Standard_Sample_List;

  -- DESCRIPTION :
  --   Same as Parallel_Sample, except that only the constant of the
  --   first hyperplane of the n-1 hyperplanes defined by b+t*v changes.

  function Sample ( b,v,t : Vector; k : integer32 )
                  return Array_of_Standard_Sample_Lists;
  function Sample ( file : file_type; output : boolean;
                    b,v,t : Vector; k : integer32 )
                  return Array_of_Standard_Sample_Lists;

  -- DESCRIPTION :
  --   Returns an array of sample lists of range 0..k, at 0 we find the
  --   list created from b,v,t, and the entries 1..k contain new samples.
  --   Eventually writes intermediate results to file, with diagnostics
  --   of the predictor-corrector stage if output is set to true.

  function Parallel_Sample ( b,v,t : Vector; k : integer32 )
                           return Array_of_Standard_Sample_Lists;
  function Parallel_Sample ( file : file_type; output : boolean;
                             b,v,t : Vector; k : integer32 )
                           return Array_of_Standard_Sample_Lists;

  -- DESCRIPTION :
  --   Same as Sample, except that now v remains the same.

  function Parallel_Sample1 ( b,v,t : Vector; k : integer32 )
                            return Array_of_Standard_Sample_Lists;
  function Parallel_Sample1 ( file : file_type; output : boolean;
                              b,v,t : Vector; k : integer32 )
                            return Array_of_Standard_Sample_Lists;

  -- DESCRIPTION :
  --   Same as Parallel_Sample, except that only the constant of the
  --   first hyperplane of the n-1 hyperplanes defined by b+t*v changes.

  function Sample ( b,v,t : Vector ) return Stacked_Sample_Grid;
  function Sample ( file : file_type;
                    b,v,t : Vector ) return Stacked_Sample_Grid;

  -- DESCRIPTION :
  --   Returns a grid of sample to interpolate a polynomial in n variables
  --   whose witness points on b + t*v are encoded in the vector t.
  --   This stacked sample grid works for the Newton-Taylor form
  --   of the interpolating polynomial.

  function Full_Sample ( b,v,t : Vector ) return Stacked_Sample_Grid;
  function Full_Sample ( file : file_type;
                         b,v,t : Vector ) return Stacked_Sample_Grid;

  -- DESCRIPTION :
  --   This procedure create the full grid, as required by the basic
  --   version of the trace form of the interpolating polynomial.

  procedure Clear;

  -- DESCRIPTION :
  --   Destroys the internal data structures.

end Hypersurface_Sample_Grids;
