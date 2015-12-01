with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Exponent_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Single_Polyhedral_Trackers is

-- DESCRIPTION :
--   This package offers routines to track one single path defined
--   by a polyhedral homotopy, for use in multitasking,
--   in double, double double, or quad double precision.

  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in Standard_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in Standard_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in DoblDobl_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure Track_Path
               ( mxt : in Standard_Integer_Vectors.Vector;
                 lft : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 nor : in Standard_Floating_Vectors.Link_to_Vector;
                 cff : in QuadDobl_Complex_VecVecs.VecVec;
                 dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                 cft : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                 epv : in Exponent_Vectors.Exponent_Vectors_Array;
                 hom : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                 ejf : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                 jmf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                 sol : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Calls the path tracker on the solution sol.

  -- ON ENTRY :
  --   mxt      type of mixture of the supports;
  --   lft      lifted supports;
  --   cff      coefficients in the polyhedral homotopy;
  --   dpw      workspace for scaled powers of the polyhedral homotopy,
  --            computed and scaled here for the particular cell;
  --   cft      coefficients of the polyhedral homotopy that allow
  --            to evaluate at the continuation parameter t;
  --   epv      exponents vectors array;
  --   hom      Laurent homotopy to evaluate a coefficient system;
  --   ejf      Jacobian matrix that corresponds to the homotopy hom;
  --   jmf      multiplication factors in the Jacobian matrix ejf;
  --   sol      pointer to the start solution.

  -- ON RETURN :
  --   sol      the solution at the end of the solution path.

end Single_Polyhedral_Trackers;
