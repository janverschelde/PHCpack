with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_VecVecs;
with TripDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;

package Test_Series_Solutions is

-- DESCRIPTION :
--   Test on the development of series as solutions of polynomial systems.

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies Newton's method in standard double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies Newton's method in double double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies Newton's method in triple double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies Newton's method in quad double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in TripDobl_Complex_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in triple double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.


  procedure Standard_Test_at_Zero_Order ( echelon : in boolean );
  procedure DoblDobl_Test_at_Zero_Order ( echelon : in boolean );
  procedure TripDobl_Test_at_Zero_Order ( echelon : in boolean );
  procedure QuadDobl_Test_at_Zero_Order ( echelon : in boolean );

  -- DESCRIPTION :
  --   Tests in double, double double, trip double,
  --   or quad double precision,
  --   prompting for a system and its solutions,
  --   starting Newton's method at zero-th order series.

  procedure Standard_Test_at_Series ( echelon : in boolean );
  procedure DoblDobl_Test_at_Series ( echelon : in boolean );
  procedure TripDobl_Test_at_Series ( echelon : in boolean );
  procedure QuadDobl_Test_at_Series ( echelon : in boolean );

  -- DESCRIPTION :
  --   Tests in double, double double, trip double,
  --   or quad double precision,
  --   prompting for a system and its solutions,
  --   starting Newton's method at a series.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the setup, the working precision,
  --   and then launches the corresponding test.

end Test_Series_Solutions;
