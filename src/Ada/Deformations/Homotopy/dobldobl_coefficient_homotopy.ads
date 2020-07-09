with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;         use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;

package DoblDobl_Coefficient_Homotopy is

-- DESCRIPTION :
--   A coefficient homotopy for (1-t)*f + t*g is beneficial if f and g
--   share many common monomials.  In that case, the evaluation of the
--   shared monomials happens only once.

-- PART I : for a pair of polynomials

  function Coefficients ( p : Poly ) return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector of the polynomial p.

  function Labeled_Coefficients ( p : Poly; real : boolean ) return Poly;

  -- DESCRIPTION :
  --   Replaces the k-th coefficient in p by k if real or k*I otherwise.

  function Index_of_Labels
             ( c : DoblDobl_Complex_Vectors.Vector; real : boolean )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The k-th entry in the vector on return contains the index to the
  --   position of k in c if real, or the position of k*I in c otherwise.

  procedure Evaluated_Coefficients
             ( cff : in out DoblDobl_Complex_Vectors.Vector;
               cp,cq : in DoblDobl_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               t : in double_double );
  procedure Evaluated_Coefficients
             ( cff : in out DoblDobl_Complex_Vectors.Vector;
               cp,cq : in DoblDobl_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               k : in natural32;
               gamma : in DoblDobl_Complex_Vectors.Vector;
               t : in Complex_Number );
  procedure Evaluated_Coefficients
             ( cff : in DoblDobl_Complex_Vectors.Link_to_Vector;
               cp,cq : in DoblDobl_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               t : in double_double );
  procedure Evaluated_Coefficients
             ( cff : in DoblDobl_Complex_Vectors.Link_to_Vector;
               cp,cq : in DoblDobl_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               k : in natural32;
               gamma : in DoblDobl_Complex_Vectors.Vector;
               t : in Complex_Number );

  -- DESCRIPTION :
  --   Evaluates the coefficients of the polynomial (1-t)*p + t*q.

  -- ON ENTRY :
  --   cff     large enough for all coefficients of (1-t)*p + t*q,
  --           must be properly initialized to zero;
  --   cp      coefficient vector of p;
  --   cq      coefficient vector of q;
  --   ip      index of the labels of coefficients of p in h:
  --           ip(k) locates in ch the k-th coefficient of p,
  --           in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq      index of the labels of coefficients of q in h:
  --           ip(k) locates in ch the k-th coefficient of q,
  --           in particular: ch(iq(k)) := cq(k) sets h to q;
  --   k       power of t and (1-t) is one if omitted;
  --   gamma   gamma constants are one if omitted;
  --   t       some double float, typically in [0,1].

  -- ON RETURN :
  --   cff     evaluated coefficients of (1-t)*p + t*q.

  function Evaluated_Coefficients
             ( nbcff : integer32;
               cp,cq : DoblDobl_Complex_Vectors.Vector;
               ip,iq : Standard_Integer_Vectors.Vector;
               t : double_double ) return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..nbcff with the coefficients of
  --   the polynomial h = (1-t)*p + t*q, with coefficient vector ch.

  -- ON ENTRY :
  --   nbcff   total number of coefficients in h;
  --   cp      coefficient vector of p;
  --   cq      coefficient vector of q;
  --   ip      index of the labels of coefficients of p in h:
  --           ip(k) locates in ch the k-th coefficient of p,
  --           in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq      index of the labels of coefficients of q in h:
  --           ip(k) locates in ch the k-th coefficient of q,
  --           in particular: ch(iq(k)) := cq(k) sets h to q;
  --   t       some double float, typically in [0,1].

-- PART II : for a pair of polynomial systems

  function Coefficients
             ( p : Poly_Sys ) return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the polynomials in the system p.

  function Labeled_Coefficients
             ( p : Poly_Sys; real : boolean ) return Poly_Sys;

  -- DESCRIPTION :
  --   Replaces the k-th coefficient in p by k if real or k*I otherwise.

  function Index_of_Labels
             ( c : DoblDobl_Complex_VecVecs.VecVec; real : boolean )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The k-th entry in the vector e(i) on return contains the index to the
  --   position of k in c(i) if real, or the position of k*I in c(i) otherwise.

  procedure Evaluated_Coefficients
             ( cff : in out DoblDobl_Complex_VecVecs.VecVec;
               cp,cq : in DoblDobl_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               t : in double_double );
  procedure Evaluated_Coefficients
             ( cff : in out DoblDobl_Complex_VecVecs.VecVec;
               cp,cq : in DoblDobl_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               k : in natural32;
               gamma : in DoblDobl_Complex_Vectors.Vector;
               t : in Complex_Number );

  -- DESCRIPTION :
  --   Evaluates the coefficients of the system (1-t)*p + t*q.

  -- ON ENTRY :
  --   cff     large enough for all coefficients of (1-t)*p + t*q,
  --           must be properly initialized to zero;
  --   cp      coefficient vector of p;
  --   cq      coefficient vector of q;
  --   ip      index of the labels of coefficients of p in h:
  --           ip(k) locates in ch the k-th coefficient of p,
  --           in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq      index of the labels of coefficients of q in h:
  --           ip(k) locates in ch the k-th coefficient of q,
  --           in particular: ch(iq(k)) := cq(k) sets h to q;
  --   k       power of t and (1-t) is one if omitted;
  --   gamma   gamma constants are one if omitted;
  --   t       some double float, typically in [0,1].

  -- ON RETURN :
  --   cff     evaluated coefficients of (1-t)*p + t*q.

-- PART III : encapsulation

  procedure Create ( p,q : in Poly_Sys; k : in natural32;
                     a : in Complex_Number );

  -- DESCRIPTION :
  --   The following artificial-parameter homotopy is constructed:
  --     H(x,t) = a * ((1 - t)^k) * p + (t^k) * q.
  --   Note that p is the start system and q the target!

  function Number_of_Equations return integer32;

  -- DESCRIPTION :
  --   Returns the number of equations in the homotopy,
  --   returns -1 if there is no homotopy created.

  function All_Start_Coefficients return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the start equations
  --   in the stored coefficient homotopy.

  -- REQUIRED : Number_of_Equations /= -1.

  function Start_Coefficients
             ( idx : integer32 )
             return DoblDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the start equations of index idx
  --   in the stored coefficient homotopy.

  -- REQUIRED : Number_of_Equations /= -1.

  function All_Target_Coefficients return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the target equations
  --   in the stored coefficient homotopy.

  -- REQUIRED : Number_of_Equations /= -1.

  function Target_Coefficients
             ( idx : integer32 )
             return DoblDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the target equations of index idx
  --   in the stored coefficient homotopy.

  -- REQUIRED : Number_of_Equations /= -1.

  function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   The homotopy is evaluated in x and t and a vector is returned.

  function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                  t : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   The homotopy is differentiated to x and the Jacobian matrix
  --   of H(x,t) is returned.

  procedure Clear;

  -- DESCRIPTION :
  --   The homotopy is cleared.

end DoblDobl_Coefficient_Homotopy;
