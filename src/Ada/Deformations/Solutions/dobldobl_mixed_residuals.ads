with text_io;                            use text_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;

package DoblDobl_Mixed_Residuals is

-- DESCRIPTION :
--   A mixed residual of a polynomial system mixes relative and absolute
--   criteria to decide whether a point is a solution of the system.
--   The computations happen in double double precision.

  function AbsVal ( z : Vector ) return Vector;

  -- DESCRIPTION :
  --   The i-th component of the vector on return contains the
  --   absolute value (or radius) of z(i).

  function AbsVal ( t : DoblDobl_Complex_Polynomials.Term )
                  return DoblDobl_Complex_Polynomials.Term;
  function AbsVal ( t : DoblDobl_Complex_Laurentials.Term )
                  return DoblDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Returns the same term as t, but with a coefficient
  --   equal to the absolute value of t.cf.

  function AbsVal ( p : DoblDobl_Complex_Polynomials.Poly )
                  return DoblDobl_Complex_Polynomials.Poly;
  function AbsVal ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function AbsVal ( p : DoblDobl_Complex_Laurentials.Poly )
                  return DoblDobl_Complex_Laurentials.Poly;
  function AbsVal ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the polynomials with the same terms as p, but with 
  --   coefficients equal to their absolute values.

  function Residual ( val,avl : Complex_Number ) return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual, given in val the value of a polynomial
  --   and in avl the value of the same polynomial with absolute coefficients
  --   at the root with absolute coordinates.

  function Residual ( val,avl : Vector ) return double_double;

  -- DESCRIPTION :
  --   Given in val the value of a polynomial system at a point and
  --   given in avl the value of the same system with absolute value
  --   coefficient at the point with absolute coordinates,
  --   returns the average of the mixed residuals of the components.

  procedure Residual ( val,avl : in Complex_Number;
                       vpz,vap,res : out double_double );

  -- DESCRIPTION :
  --   Evaluates the formula for the mixed residual.

  -- ON ENTRY :
  --   val      value of a polynomial at a point;   
  --   avl      value of the same polynmial with absolute coefficients,
  --            at the point with absolute coordinates.

  -- ON RETURN :
  --   vpz      the radius of the number val;
  --   vap      the radius of the number avl;
  --   res      the value of vpz/(vap+1.0), the mixed residual.

  procedure Write_Residuals
              ( file : in file_type;
                vaz,vpz,vap,res : in double_double );

  -- DESCRIPTION :
  --   Writes the mixed residual res and the auxiliary vaz, vpz, vap
  --   to the output file.

  procedure Residual ( file : in file_type;
                       pol,abp : in DoblDobl_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( pol,abp : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( pol,abp : in DoblDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in DoblDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( pol,abp : in DoblDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in DoblDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_double );

  -- DESCRIPTION :
  --   Computes the mixed residual, with the return of all auxiliary values.
  --   If a file is provided, then one line is written to file.

  -- ON ENTRY :
  --   pol       a polynomial in several variables;
  --   abp       the same as pol, but with absolute values of coefficients,
  --             abp = AbsVal(pol);
  --   z         values for the variables, z'first = 1, and z'last equals
  --             the number of variables in pol.

  -- ON RETURN :
  --   abz       AbsVal(z), radii of the coordinates of z;
  --   vaz       the max norm of the vector abz;
  --   vpz       radius of the value of pol at z;
  --   vap       radius of the value of abp at abz;
  --   res       the mixed residual is vpz/(vap+1.0).

  function Residual ( pol,abp : DoblDobl_Complex_Polynomials.Poly;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Polynomials.Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laurentials.Poly;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Laurentials.Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).  If a file is provided, then one line
  --   for every polynomial is written to the output file.

  function Residual ( pol,abp : DoblDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_double;
  function Residual ( file : file_type;
                      pol,abp : DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).  If a file is provided, then one line
  --   for every polynomial is written to the output file.

  function Mixed_Residual
              ( valres,absres : in DoblDobl_Complex_Vectors.Vector )
              return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual for the evaluation in valres of the
  --   convolutions at a point, and for absres the evaluated convolutions
  --   with radii as coefficients at the absolute values of the point.

end DoblDobl_Mixed_Residuals;
