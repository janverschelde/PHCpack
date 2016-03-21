with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package Root_Refining_Parameters is

-- DESCRIPTION :
--   This packages offers default choices for the root refining parameters,
--   as well as interactive procedures for the user to set the values.

  procedure Standard_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean );
  procedure DoblDobl_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean );
  procedure QuadDobl_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural32; deflate,wout : out boolean );
  procedure Multprec_Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out Floating_Number;
                maxit,deci : out natural32; deflate,wout : out boolean );

  -- DESCRIPTION :
  --  Defines the default values for the root refining parameters,
  --  for double, double double, quad double, and multiprecision arithmetic.

  -- ON RETURN :
  --  epsxa     precision for correction on x;
  --  epsfa     precision for residual;
  --  tolsing   tolerance on inverse condition numbers;
  --  maxit     maximal number of Newton iterations;
  --  deflate   if deflation is wanted;
  --  wout      if intermediate output is wanted.

  procedure Standard_Put_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in double_float;
                maxit : in natural32; deflate,wout : in boolean );
  procedure Multprec_Put_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in Floating_Number;
                maxit,deci : in natural32; deflate,wout : in boolean );

  -- DESCRIPTION :
  --   Writes the parameters for the root refiner on file.

  procedure Standard_Menu_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in out double_float;
                maxit : in out natural32; deflate,wout : in out boolean );
  procedure Multprec_Menu_Root_Refining_Parameters
              ( file : in file_type;
                epsxa,epsfa,tolsing : in out Floating_Number;
                maxit,deci : in out natural32; deflate,wout : in out boolean );

  -- DESCRIPTION :
  --   The user can set the parameters of the root refiner via menus.

end Root_Refining_Parameters;
