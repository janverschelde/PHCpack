with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package m_Homogeneous_Start_Systems is

-- DESCRIPTION :
--   the purpose of this package is to provide a routine
--   for the automatic construction of a m-homomogeneous
--   start system, given a partition of the set of unknowns.

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in Partition );

  -- ON ENTRY :
  --   p           polynomial system;
  --   z           partition of the set of unknowns.

  -- ON RETURN :
  --   The data managed by Standard_Linear_Product_System
  --   now contains a linear-product start system.

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in Partition;
                   q : out Poly_Sys; qsols : in out Solution_List );
  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in Partition;
                   q : out Poly_Sys; rq : out Prod_Sys;
                   qsols : in out Solution_List );

  -- ON ENTRY :
  --   p           polynomial system;
  --   z           partition of the set of unknowns.

  -- ON RETURN :
  --   q           an m-homogeneous start system;
  --   rq          start system in linear-product format;
  --   qsols       the solutions of the start system q.

end m_Homogeneous_Start_Systems;
