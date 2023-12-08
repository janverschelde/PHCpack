with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_VecMats;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_VecMats;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with PentDobl_Speelpenning_Convolutions;
with OctoDobl_Speelpenning_Convolutions;
with DecaDobl_Speelpenning_Convolutions;
with HexaDobl_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;
with DoblDobl_Coefficient_Convolutions;
with QuadDobl_Coefficient_Convolutions;

package Multitasked_AlgoDiff_Convolutions is

-- DESCRIPTION :
--   Offers algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return Standard_Floating_VecVecVecs.VecVecVec;

  -- DESCRIPTION :
  --   Returns work space for nbt tasks to evaluate coefficient circuits
  --   of dimension dim at power series of degree deg in double precision.

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return Standard_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return DoblDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return TripDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return QuadDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return PentDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return OctoDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return DecaDobl_Speelpenning_Convolutions.VecVecVec;
  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return HexaDobl_Speelpenning_Convolutions.VecVecVec;

  -- DESCRIPTION :
  --   Returns work space for nbt tasks to evaluate circuits of
  --   dimension dim at power series of degree deg,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, deca double, or hexa double precision.

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Computes the power table at x.
  --   Evaluates and differentiates the coefficient convolution circuits 
  --   in c at x with multitasking in double precision.

  -- REQUIRED :
  --   The power table rpwt, ipwt is allocated for the degrees in mxe.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   rx       real parts of the coefficients of power series;
  --   ix       imaginary parts of the coefficients of power series;
  --   mxe      exponent maxima over all dimensions;
  --   rpwt     real parts of the power table;
  --   ipwt     imaginary parts of the power table;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x;
  --   static   flag to apply static load balancing,
  --            which by default is turned off because less efficient,
  --            by default, dynamic load balancing is applied;
  --   output   true for the writing of intermediate output.

  procedure DoblDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Computes the power table at x.
  --   Evaluates and differentiates the coefficient convolution circuits 
  --   in c at x with multitasking in double double precision.

  -- REQUIRED :
  --   The power table *pwt is allocated for the degrees in mxe.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   rhx      real high parts of the coefficients of power series;
  --   ihx      imaginary high parts of the coefficients of power series;
  --   rlx      real low parts of the coefficients of power series;
  --   ilx      imaginary low parts of the coefficients of power series;
  --   mxe      exponent maxima over all dimensions;
  --   rhpwt    real high parts of the power table;
  --   ihpwt    imaginary low parts of the power table;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x;
  --   static   flag to apply static load balancing,
  --            which by default is turned off because less efficient,
  --            by default, dynamic load balancing is applied;
  --   output   true for the writing of intermediate output.

  procedure QuadDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Computes the power table at x.
  --   Evaluates and differentiates the coefficient convolution circuits 
  --   in c at x with multitasking in quad double precision.

  -- REQUIRED :
  --   The power table *pwt is allocated for the degrees in mxe.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   xr       real parts of the coefficients of power series;
  --   xr       imaginary parts of the coefficients of power series;
  --   rlx      real low parts of the coefficients of power series;
  --   ilx      imaginary low parts of the coefficients of power series;
  --   mxe      exponent maxima over all dimensions;
  --   rpwt     real parts of the power table;
  --   ipwt     imaginary parts of the power table;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x;
  --   static   flag to apply static load balancing,
  --            which by default is turned off because less efficient,
  --            by default, dynamic load balancing is applied;
  --   output   true for the writing of intermediate output.

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure DoblDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure TripDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in TripDobl_Speelpenning_Convolutions.Circuits;
                x : in TripDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in TripDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in TripDobl_Complex_VecVecs.VecVec;
                vm : in TripDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure QuadDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure PentDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in PentDobl_Speelpenning_Convolutions.Circuits;
                x : in PentDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in PentDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in PentDobl_Complex_VecVecs.VecVec;
                vm : in PentDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure OctoDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in OctoDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in OctoDobl_Complex_VecVecs.VecVec;
                vm : in OctoDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure DecaDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in DecaDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in DecaDobl_Complex_VecVecs.VecVec;
                vm : in DecaDobl_Complex_VecMats.VecMat;
                output : in boolean := false );
  procedure HexaDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in HexaDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in HexaDobl_Complex_VecVecs.VecVec;
                vm : in HexaDobl_Complex_VecMats.VecMat;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Computes the power table at x.
  --   Evaluates and differentiates the convolution circuits in c at x
  --   with multitasking, in double, double double, triple double,
  --   quad double, octo double, deca double, or hexa double precision.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   x        coefficients of power series;
  --   mxe      exponent maxima over all dimensions;
  --   pwt      power table, allocated for the degrees in mxe;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x;
  --   output   true for the writing of intermediate output.

end Multitasked_AlgoDiff_Convolutions;
