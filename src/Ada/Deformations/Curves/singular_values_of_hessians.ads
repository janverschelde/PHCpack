with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Hessians;

package Singular_Values_of_Hessians is

-- DESCRIPTION :
--   Given symbolic definitions of Hessian matrices,
--   evaluates the Hessians at numerical vectors and
--   returns the singular values.

   function Standard_Singular_Values
              ( h : Standard_Complex_Hessians.Link_to_Hessian;
                x : Standard_Complex_Vectors.Vector )
              return  Standard_Complex_Vectors.Vector;
   function DoblDobl_Singular_Values
              ( h : DoblDobl_Complex_Hessians.Link_to_Hessian;
                x : DoblDobl_Complex_Vectors.Vector )
              return  DoblDobl_Complex_Vectors.Vector;
   function QuadDobl_Singular_Values
              ( h : QuadDobl_Complex_Hessians.Link_to_Hessian;
                x : QuadDobl_Complex_Vectors.Vector )
              return  QuadDobl_Complex_Vectors.Vector;

   -- DESCRIPTION :
   --   Given a Hessian h and a numerical vector v,
   --   returns the singular values computed in double precision,
   --   or double double precision, or quad double precision.

end Singular_Values_of_Hessians;
