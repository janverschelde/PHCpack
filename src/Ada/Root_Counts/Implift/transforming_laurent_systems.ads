with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with Standard_Integer32_Transformations; 
 use Standard_Integer32_Transformations;

package Transforming_Laurent_Systems is

-- DESCRIPTION :
--   This package offers some routines for transforming Laurent polynomials.

  procedure Shift ( p : in out Standard_Complex_Laurentials.Poly );
  function  Shift ( p : Standard_Complex_Laurentials.Poly )
                  return Standard_Complex_Laurentials.Poly;
  procedure Shift ( p : in out DoblDobl_Complex_Laurentials.Poly );
  function  Shift ( p : DoblDobl_Complex_Laurentials.Poly )
                  return DoblDobl_Complex_Laurentials.Poly;
  procedure Shift ( p : in out QuadDobl_Complex_Laurentials.Poly );
  function  Shift ( p : QuadDobl_Complex_Laurentials.Poly )
                  return QuadDobl_Complex_Laurentials.Poly;

  procedure Shift ( L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Shift ( L : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  procedure Shift ( L : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  function  Shift ( L : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  procedure Shift ( L : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );
  function  Shift ( L : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Shifts the support of the polynomial so that the constant term
  --   belongs to p.
  --   This Shift does not change the term order in p!

  procedure Transform ( t : in Transfo;
                        p : in out Standard_Complex_Laurentials.Poly );
  function  Transform ( t : Transfo;
                        p : Standard_Complex_Laurentials.Poly )
                      return Standard_Complex_Laurentials.Poly;

  procedure Transform ( t : in Transfo;
                        L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Transform ( t : Transfo;
                        L : Standard_Complex_Laur_Systems.Laur_Sys )
                      return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION : Application of the transformation t.

  function Maximal_Support
             ( p : Standard_Complex_Laurentials.Poly;
               v : Vector ) return integer32;
  function Maximal_Support
             ( p : Standard_Complex_Laurentials.Poly; 
               v : Link_to_Vector ) return integer32;

  -- DESCRIPTION :
  --   Computes the value of the supporting function of the polynomial p,
  --   for the direction v.

  procedure Face ( i,m : in integer32;
                   p : in out Standard_Complex_Laurentials.Poly );
  function  Face ( i,m : integer32;
                   p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Poly;

  procedure Face ( i,m : in integer32;
                   L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Face ( i,m : integer32;
                   L : Standard_Complex_Laur_Systems.Laur_Sys )
                 return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   returns only the terms t for which deg(t,xi) = m.

  procedure Face ( v : in Vector; m : in integer32;
                   p : in out Standard_Complex_Laurentials.Poly );
  function  Face ( v : Vector; m : integer32;
                   p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Poly;

  procedure Face ( v,m : in Vector;
                   L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Face ( v,m : Vector;
                   L : Standard_Complex_Laur_Systems.Laur_Sys )
                 return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Only the terms for which for the degrees d the following holds
  --    < d , v > = m, are left.

  procedure Reduce ( i : in integer32;
                     p : in out Standard_Complex_Laurentials.Poly );
  function  Reduce ( i : integer32;
                     p : Standard_Complex_Laurentials.Poly )
                   return Standard_Complex_Laurentials.Poly;

  procedure Reduce ( i : in integer32;
                     L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Reduce ( i : integer32;
                     L : Standard_Complex_Laur_Systems.Laur_Sys )
                   return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   The i-th component out of every monomial will be removed,
  --   so that the polynomials will have an unknown less.

  procedure Insert ( i,d : in integer32;
                     p : in out Standard_Complex_Laurentials.Poly );
  function  Insert ( i,d : integer32;
                     p : Standard_Complex_Laurentials.Poly )
                   return Standard_Complex_Laurentials.Poly;

  procedure Insert ( i,d : in integer32;
                     L : in out Standard_Complex_Laur_Systems.Laur_Sys );
  function  Insert ( i,d : integer32;
                     L : Standard_Complex_Laur_Systems.Laur_Sys )
                   return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   The i-th component of each monomial will be inserted,
  --   using the value d.

end Transforming_Laurent_Systems;
