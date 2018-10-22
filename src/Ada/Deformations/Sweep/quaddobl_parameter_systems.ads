with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;        use QuadDobl_Complex_Solutions;

package QuadDobl_Parameter_Systems is

-- DESCRIPTION :
--   Any system with parameters is a coefficient-parameter homotopy.
--   Operations specific to such systems concern the definition of
--   the parameters, as interactively determined by the user.
--   Once the system has been solved for particular values of the
--   parameters, we need to substitute the parameters by these values
--   for root refinement.

  procedure Sort ( v : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   The labels of the parameters must be sorted in increasing order,
  --   otherwise, the substitute will not work.

  procedure Read_Solution_Parameters
              ( infile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 );

  -- DESCRIPTION :
  --   Scans the infile for the solutions for the system p.
  --   This procedure is called by Read_Parameter_Homotopy.

  -- REQUIRED : number of variables > number of equations,
  --   the symbol table is already initiliazed by the reading
  --   of the polynomial which define the parameter homotopy.

  -- ON ENTRY :
  --   infile   file with system already read from;
  --   p        a polynomial system with parameters.

  -- ON RETURN :
  --   sols     solutions for particular values of the parameters;
  --   nb_equ   number of equations in the system;
  --   nb_unk   number of unknowns in the system;
  --   nb_par   number of parameters: nb_unk - nb_equ.

  procedure Read_Solution_Parameters
              ( infile : in file_type; outfile : out file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 );

  -- DESCRIPTION :
  --   Scans the infile for the solutions for the system p.
  --   This procedure is called by Read_Parameter_Homotopy.
  --   Prompts the user for the name of an output file.

  -- REQUIRED : number of variables > number of equations

  -- ON ENTRY :
  --   infile   file with system already read from;
  --   p        a polynomial system with parameters.

  -- ON RETURN :
  --   sols     solutions for particular values of the parameters;
  --   nb_equ   number of equations in the system;
  --   nb_unk   number of unknowns in the system;
  --   nb_par   number of parameters: nb_unk - nb_equ.

  procedure Read_Parameter_Homotopy
              ( lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 );

  -- DESCRIPTION :
  --   Reads a coefficient-paremeter homotopy with solutions for 
  --   some values of the parameters.

  -- ON RETURN :
  --   lp       a polynomial system;
  --   sols     solutions for particular values of the parameters;
  --   nb_equ   number of equations in the system;
  --   nb_unk   number of unknowns in the system;
  --   nb_par   number of parameters: nb_unk - nb_equ.

  procedure Read_Parameter_Homotopy
              ( outfile : out file_type;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 );

  -- DESCRIPTION :
  --   Reads a coefficient-paremeter homotopy with solutions for 
  --   some values of the parameters.

  -- ON RETURN :
  --   outfile  file created for output;
  --   lp       a polynomial system;
  --   sols     solutions for particular values of the parameters;
  --   nb_equ   number of equations in the system;
  --   nb_unk   number of unknowns in the system;
  --   nb_par   number of parameters: nb_unk - nb_equ.

  function Define_Parameters ( nb_equ,nb_unk,nb_par : integer32 )
                             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector with indices to the parameters,
  --   as defined interactively by the user.
  --   The vector on return is sorted in increasing order.

  function Complement ( n : integer32; v : Standard_Integer_Vectors.Vector )
                      return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the complement of the elements in v with respect to
  --   the set of natural numbers 1,2,..,n.

  -- NOTE: the complement of the parameters are those variables which
  --   depend on the natural parameters via the equations of the system.

  -- REQUIRED : the entries in v are sorted in increasing order.

  function Substitute
              ( t : QuadDobl_Complex_Polynomials.Term;
                pars : Standard_Integer_Vectors.Vector;
                vals : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Term;

  -- DESCRIPTION :
  --   Replaces the parameters in t by the given values in vals.
  --   The parameters are defined by the indices in pars.

  function Substitute
              ( p : QuadDobl_Complex_Polynomials.Poly;
                pars : Standard_Integer_Vectors.Vector;
                vals : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Replaces the parameters defined by the indices in pars
  --   by the values in vals in the polynomial p.

  function Substitute
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pars : Standard_Integer_Vectors.Vector;
                vals : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Replaces the parameters defined by the indices in pars
  --   by the values in vals in the polynomial p.

end QuadDobl_Parameter_Systems;
