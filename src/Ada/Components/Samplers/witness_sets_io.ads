with text_io;                            use text_io;
with Symbol_Table;                       use Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Witness_Sets_io is

-- DESCRIPTION :
--   A witness set contains an overdetermined polynomial system
--   which is made square by the addition of slack variables.
--   The number of slack variables -- usually starting with "zz" --
--   equals the dimension of the algebraic set.
--   If the system on input has more equations than unknowns,
--   then extra "ss" variables are added to make the system square.
--   This package offers utilities to recognize these "zz" variables.

  procedure Write_Symbol_Table ( file : in file_type );
 
  -- DESCRPITION :
  --   Writes the content of the symbol table in the current order.

  procedure Add_Slack_Symbols ( k : in natural32 );

  -- DESCRIPTION :
  --   Adds k new symbols of the form "ssi" to the symbol table,
  --   where i is a number ranging from 1 to k.

  procedure Add_Embed_Symbols ( k : in natural32 );

  -- DESCRIPTION :
  --   Adds k new symbols of the form "zzi" to the symbol table,
  --   where i is a number ranging from 1 to k.

  procedure Add_Extra_Symbols ( k : in natural32 );

  -- DESCRIPTION :
  --   When variables are missing from the input, simply because
  --   they are free, the user must provide these symbols.
  --   This routine writes the symbol table and prompts the user
  --   for k extra symbols.

  function Count_Embed_Symbols
              ( sbt : Array_of_Symbols; s : string ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of symbols in sbt which start with the same
  --   characters as in the string s.

  function Count_Embed_Symbols ( n : natural32; s : string ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of symbols starting with the characters in s.

  -- ON ENTRY :
  --   n        number of variables, n <= Symbol_Table.Number;
  --   s        common start of variables, typical choices are
  --            "zz" for slack variables in the embedding, and
  --            "ss" for slack variables for overdetermined systems.

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Swaps the k symbols whose initial characters match s completely,
  --   to the end of the symbol table, marked by n.
  --   Also swaps the variables in the polynomial system p.

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );
  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Swaps the k symbols whose initial characters match s completely,
  --   to the end of the symbol table, marked by n.
  --   The swaps also apply to the system p and its solutions sols.

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Sorts the embedded symbols in ascending order.

  -- ON ENTRY :
  --   nv       total number of variables;
  --   n        dimension before the embedding;
  --   d        number of embed symbols, typically: n + d = nv;
  --   f        embedded polynomial system.

  -- ON RETURN :
  --   f        system with embedded variables sorted in ascending order.

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );
  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Sorts the embedded symbols in ascending order.

  -- ON ENTRY :
  --   nv       total number of variables;
  --   n        dimension before the embedding;
  --   d        number of embed symbols, typically: n + d = nv;
  --   f        embedded polynomial system;
  --   sols     solutions of f.

  -- ON RETURN :
  --   f        system with embedded variables sorted in ascending order;
  --   sols     solutions with sorted embedded variables.

  procedure Standard_Read_Embedding
               ( file : in file_type;
                 lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure Standard_Read_Embedding
               ( file : in file_type;
                 lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );

  -- DESCRIPTION :
  --   Reads the embedded system and its solutions from file.
  --   Ensures the added variables (marked with zz symbols) occur last.

  -- ON RETURN :
  --   file      the writing starts from file;
  --   lp        pointer to the embedded polynomial system;
  --   sols      solutions of lp;
  --   dim       number of added variables in the embedding.

  procedure Standard_Read_Embedding
               ( file : in file_type;
                 lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure Standard_Read_Embedding
               ( file : in file_type;
                 lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( file : in file_type;
                 lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );

  -- DESCRIPTION :
  --   Reads the embedded system and its solutions from file.
  --   Ensures the added variables (marked with zz symbols) occur last.

  -- ON RETURN :
  --   file      the writing starts from file;
  --   lp        pointer to the embedded polynomial system;
  --   sols      solutions of lp;
  --   dim       number of added variables in the embedding;
  --   nsl       number of slack variables used to square the system.

  procedure Standard_Read_Embedding
               ( lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure Standard_Read_Embedding
               ( lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );

  -- DESCRIPTION :
  --   Reads the embedded system and its solutions from file.
  --   Ensures the added variables (marked with zz symbols) occur last.

  -- ON RETURN :
  --   lp        pointer to the embedded polynomial system;
  --   sols      solutions of lp;
  --   dim       number of added variables in the embedding.

  procedure Standard_Read_Embedding
               ( lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure Standard_Read_Embedding
               ( lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 dim,nsl : out natural32 );

  -- DESCRIPTION :
  --   Reads the embedded system and its solutions from file.
  --   Ensures the added variables (marked with zz symbols) occur last.

  -- ON RETURN :
  --   lp        pointer to the embedded polynomial system;
  --   sols      solutions of lp;
  --   dim       number of added variables in the embedding;
  --   nsl       number of slack variables used to square the system.

  procedure Standard_Read_Embedding
               ( name : in string;
                 lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( name : in string;
                 lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( name : in string;
                 lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );

  -- DESCRIPTION :
  --   Attempts to open the file with name in the string filename to read
  --   a witness set defined by an ordinary polynomial system.
  --   If an exception occurs, the user is prompted for a file.

  -- ON ENTRY :
  --   filename  name of a file with a witness set.

  -- ON RETURN :
  --   lp        sliced and embedded polynomial system;
  --   sols      list of generic points on the slices;
  --   dim       dimension of the solution component, or in case of
  --             one polynomial, this is the number of variables.

  procedure Standard_Read_Embedding
               ( name : in string;
                 lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure DoblDobl_Read_Embedding
               ( name : in string;
                 lp : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );
  procedure QuadDobl_Read_Embedding
               ( name : in string;
                 lp : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                 sols : out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 );

  -- DESCRIPTION :
  --   Attempts to open the file with name in the string filename to read
  --   a witness set defined by a Laurent polynomial system.
  --   If an exception occurs, the user is prompted for a file.

  -- ON ENTRY :
  --   filename  name of a file with a witness set.

  -- ON RETURN :
  --   lp        sliced and embedded polynomial system;
  --   sols      list of generic points on the slices;
  --   dim       dimension of the solution component, or in case of
  --             one polynomial, this is the number of variables.

  procedure Get_Multprec_System 
               ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mpsys : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                 size,dim : in natural32 );

  -- DESCRIPTION :
  --   Depending on the answer of the user, the mpsys on return is read in
  --   from file or is the conversion of the given stsys.
  --   Note that stsys is an embedded system!

  procedure Determine_Order
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Determine_Order
               ( p : in out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Determine_Order
               ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Determine_Order
               ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Determine_Order
               ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Determine_Order
               ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   The user can interactively determine the order of the variables.
  --   The effect of this routine is that the symbol table and the
  --   internal presentation of the polynomial system is changed.
  --   The external representation of the system does not change.

  procedure Determine_Order
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The user can interactively determine the order of the variables.
  --   The effect of this routine is that the symbol table and the
  --   internal presentation of the polynomial system is changed.
  --   The external representation of the system does not change.
  --   The analogue operations are performed on the solution list.

  function Square_and_Embed
               ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : natural32 )
               return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Performs squaring, slicing and embedding with k additional slices
  --   and slack variables on the system p.  In case the original system
  --   is underdetermined, then padding with zero equations happens.
  --   For overdetermined systems, the symbol table is extended with
  --   slack variables "ss", and with k embeddings variables "zz".

end Witness_Sets_io;
