with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;

package DoblDobl_Parameter_Solutions is

-- DESCRIPTION :
--   The two major operations specific to solutions of systems with
--   natural parameters are the separation of the variables from the
--   parameters and the joining of variables with parameters.

  function Select_Variables
              ( s : Solution; nv : integer32;
                v : Standard_Integer_Vectors.Vector ) return Solution;

  -- DESCRIPTION :
  --   Returns the solution obtained by selecting those nv variables
  --   from s which are indexed by the vector v.

  function Select_Variables
              ( s : Solution_List; nv : integer32;
                v : Standard_Integer_Vectors.Vector ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the solutions obtained by selecting those nv variables
  --   from s which are indexed by the vector v.

  function Join_Variables
              ( s : Solution; n : integer32;
                vars,pars : Standard_Integer_Vectors.Vector;
                val_pars : DoblDobl_Complex_Vectors.Vector ) return Solution;

  -- DESCRIPTION :
  --   Joins the solution and values for the parameters into one
  --   solution of vectors of n entries long.

  -- ON ENTRY :
  --   s        a solution;
  --   n        dimension of the solution on return;
  --   vars     index to the variables of s in the solution on return;
  --   pars     index to the parameters in the solution on return;
  --   val_pars contains the values for the parameters.

  function Join_Variables
              ( s : Solution_List; n : integer32;
                vars,pars : Standard_Integer_Vectors.Vector;
                val_pars : DoblDobl_Complex_Vectors.Vector )
              return Solution_List;

  -- DESCRIPTION :
  --   Joins the solutions and values for the parameters into a
  --   list of solutions.

end DoblDobl_Parameter_Solutions;
