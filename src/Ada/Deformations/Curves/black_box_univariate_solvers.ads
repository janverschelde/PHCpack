with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with Multprec_Floating_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Multprec_Complex_Laurentials;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Black_Box_Univariate_Solvers is

-- DESCRIPTION :
--   Forms an interface to the black box univariate solver in PHCpack,
--   to the Durand-Kerner method, also known as the method of Weierstrass.

  function Create_Solution_List
             ( z : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution_List;
  function Create_Solution_List
             ( z : DoblDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Create_Solution_List
             ( z : QuadDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Solutions.Solution_List;
  function Create_Solution_List
             ( z : Multprec_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the root of a degree one equation in the format of a list.

  function Create_Solution_List
               ( z : Standard_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return Standard_Complex_Solutions.Solution_List;
  function Create_Solution_List
               ( z : DoblDobl_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return DoblDobl_Complex_Solutions.Solution_List;
  function Create_Solution_List
               ( z : QuadDobl_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return QuadDobl_Complex_Solutions.Solution_List;
  function Create_Solution_List
               ( z : Multprec_Complex_Vectors.Vector;
                 err,rco,res : Multprec_Floating_Vectors.Vector )
               return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the roots and residuals in the format of a solution list.

  function Coefficient_Vector
              ( d : natural32; p : Standard_Complex_Polynomials.Poly )
              return Standard_Complex_Vectors.Vector;
  function Coefficient_Vector
              ( d : natural32; p : DoblDobl_Complex_Polynomials.Poly )
              return DoblDobl_Complex_Vectors.Vector;
  function Coefficient_Vector
              ( d : natural32; p : QuadDobl_Complex_Polynomials.Poly )
              return QuadDobl_Complex_Vectors.Vector;
  function Coefficient_Vector
              ( d : natural32; p : Multprec_Complex_Polynomials.Poly )
              return Multprec_Complex_Vectors.Vector;

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : d is the degree of p.

  function Coefficient_Vector
              ( mind,maxd : integer32;
                p : Standard_Complex_Laurentials.Poly )
              return Standard_Complex_Vectors.Vector;
  function Coefficient_Vector
              ( mind,maxd : integer32;
                p : Multprec_Complex_Laurentials.Poly )
              return Multprec_Complex_Vectors.Vector;

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : mind is the minimal and maxd the maximal degree of p.

-- MAIN TARGET ROUTINES :

  procedure Black_Box_Durand_Kerner
              ( p : in Standard_Complex_Polynomials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( p : in Multprec_Complex_Polynomials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( p : in Standard_Complex_Laurentials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( p : in Multprec_Complex_Laurentials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Standard_Complex_Laurentials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Multprec_Complex_Laurentials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies the method of Durand-Kerner to the polynomial p,
  --   writing all results to file and returning the solutions
  --   in the list sols on return.
  --   For the multiprecision version, the working precision is set
  --   by the value of size.

  -- REQUIRED : Number_of_Unknowns(p) = 1.

end Black_Box_Univariate_Solvers;
