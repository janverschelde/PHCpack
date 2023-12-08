with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with TripDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with PentDobl_Complex_VecVecs;
with OctoDobl_Complex_VecVecs;
with DecaDobl_Complex_VecVecs;
with HexaDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with PentDobl_Speelpenning_Convolutions;
with OctoDobl_Speelpenning_Convolutions;
with DecaDobl_Speelpenning_Convolutions;
with HexaDobl_Speelpenning_Convolutions;

package Multitasked_Power_Newton is

-- DESCRIPTION :
--   Applies Newton's method on power series to apply Fabry's theorem
--   in double, double double, triple double, quad double, penta double,
--   octo double, deca double, and hexa double precision,
--   with multitasking for shared memory parallel computers.

  procedure Standard_Run
              ( nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_float; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure Standard_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_float; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure DoblDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure DoblDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure TripDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out triple_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure TripDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out triple_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure QuadDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out quad_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure QuadDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out quad_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure PentDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in PentDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out penta_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure PentDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in PentDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out penta_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure OctoDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in OctoDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out octo_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure OctoDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in OctoDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out octo_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure DecaDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DecaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out deca_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure DecaDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DecaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out deca_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure HexaDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in HexaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in HexaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out hexa_double; 
                output : in boolean := false;
                verbose : in boolean := true );
  procedure HexaDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in HexaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in HexaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out hexa_double; 
                output : in boolean := false;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs Newton's method with nbt tasks 
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- ON ENTRY :
  --   file     optional output file;
  --   nbt      the number of tasks;
  --   dim      number of equations and variables;
  --   maxit    maximum number of iterations;
  --   s        system of convolution circuits;
  --   scf      its leading coefficients are solution vector;
  --   tol      tolerance on the magnitude of the update vector;
  --   estco    if true, then the condition number is estimated;
  --   output   if true and verbose, then Newton's method is verbose;
  --   verbose  if true, then each step has some output.

  -- ON RETURN :
  --   fail     true if the tolerance was not reached;
  --   info     if not estco, if nonzero, then LU encountered a zero pivot;
  --   nbrit    number of iterations;
  --   rcond    if estco, estimate for the inverse condition;
  --   absdx    magnitude of the update vector.

end Multitasked_Power_Newton;
