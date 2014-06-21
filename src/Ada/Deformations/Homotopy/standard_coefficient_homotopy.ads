with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Standard_Coefficient_Homotopy is

-- DESCRIPTION :
--   A coefficient homotopy for (1-t)*f + t*g is beneficial if f and g
--   share many common monomials.  In that case, the evaluation of the
--   shared monomials happens only once.

-- PART I : for a pair of polynomials

  function Coefficients ( p : Poly ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector of the polynomial p.

  function Labeled_Coefficients ( p : Poly; real : boolean ) return Poly;

  -- DESCRIPTION :
  --   Replaces the k-th coefficient in p by k if real or k*I otherwise.

  function Index_of_Labels
             ( c : Standard_Complex_Vectors.Vector; real : boolean )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The k-th entry in the vector on return contains the index to the
  --   position of k in c if real, or the position of k*I in c otherwise.

  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_Vectors.Vector;
               cp,cq : in Standard_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               t : in double_float );
  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_Vectors.Link_to_Vector;
               cp,cq : in Standard_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               t : in double_float );

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
  --   t       some double float, typically in [0,1].

  -- ON RETURN :
  --   cff     evaluated coefficients of (1-t)*p + t*q.

  function Evaluated_Coefficients
             ( nbcff : integer32;
               cp,cq : Standard_Complex_Vectors.Vector;
               ip,iq : Standard_Integer_Vectors.Vector;
               t : double_float ) return Standard_Complex_Vectors.Vector;

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
             ( p : Poly_Sys ) return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the polynomials in the system p.

  function Labeled_Coefficients
             ( p : Poly_Sys; real : boolean ) return Poly_Sys;

  -- DESCRIPTION :
  --   Replaces the k-th coefficient in p by k if real or k*I otherwise.

  function Index_of_Labels
             ( c : Standard_Complex_VecVecs.VecVec; real : boolean )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The k-th entry in the vector e(i) on return contains the index to the
  --   position of k in c(i) if real, or the position of k*I in c(i) otherwise.

  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_VecVecs.VecVec;
               cp,cq : in Standard_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               t : in double_float );

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
  --   t       some double float, typically in [0,1].

  -- ON RETURN :
  --   cff     evaluated coefficients of (1-t)*p + t*q.
  --
end Standard_Coefficient_Homotopy;
