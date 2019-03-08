with text_io;                            use text_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;

package QuadDobl_Mixed_Residuals is

-- DESCRIPTION :
--   A mixed residual of a polynomial system mixes relative and absolute
--   criteria to decide whether a point is a solution of the system.
--   The computations happen in quad double precision.

  function AbsVal ( z : Vector ) return Vector;

  -- DESCRIPTION :
  --   The i-th component of the vector on return contains the
  --   absolute value (or radius) of z(i).

  function AbsVal ( t : QuadDobl_Complex_Polynomials.Term )
                  return QuadDobl_Complex_Polynomials.Term;
  function AbsVal ( t : QuadDobl_Complex_Laurentials.Term )
                  return QuadDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Returns the same term as t, but with a coefficient
  --   equal to the absolute value of t.cf.

  function AbsVal ( p : QuadDobl_Complex_Polynomials.Poly )
                  return QuadDobl_Complex_Polynomials.Poly;
  function AbsVal ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function AbsVal ( p : QuadDobl_Complex_Laurentials.Poly )
                  return QuadDobl_Complex_Laurentials.Poly;
  function AbsVal ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the polynomials with the same terms as p, but with 
  --   coefficients equal to their absolute values.

  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( pol,abp : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( pol,abp : in QuadDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( pol,abp : in QuadDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );
  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double );

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

  function Residual ( pol,abp : QuadDobl_Complex_Polynomials.Poly;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type; 
                      pol,abp : QuadDobl_Complex_Polynomials.Poly;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Laurentials.Poly;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laurentials.Poly;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).  If a file is provided, then one line
  --   for every polynomial is written to the output file.

  function Residual ( pol,abp : QuadDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return quad_double;
  function Residual ( pol,abp : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return quad_double;
  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).  If a file is provided, then one line
  --   for every polynomial is written to the output file.

end QuadDobl_Mixed_Residuals;
