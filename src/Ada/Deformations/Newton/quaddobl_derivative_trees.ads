with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;     use QuadDobl_Complex_Poly_Systems;

package QuadDobl_Derivative_Trees is

-- DESCRIPTION :
--   This package allows to store all partial derivatives
--   of the polynomials in a system, in double double precision.

-- DATA STRUCTURES :

  type Node;
  type Link_to_Node is access Node;
  type Array_of_Nodes is array ( integer32 range <> ) of Link_to_Node;

  type Node ( n : integer32 ) is record   -- n is #variables in p
    p : Poly;                             -- polynomial in n variables
    d : Array_of_Nodes(1..n);             -- d(i) is derivative of p wrt i
  end record;

-- CONSTRUCTORS :

  function Initialize ( p : Poly ) return Node;

  -- DESCRIPTION :
  --   Returns a node with a copy of p in its p field
  --   and all partial derivatives in d initialized to null.

  function Initialize ( p : Poly_Sys ) return Array_of_Nodes;

  -- DESCRIPTION :
  --   Applies the other Initialize from above to all polynomials
  --   in the system p.  The array on return has range p'range.

  procedure Create ( nd : in out Node );

  -- DESCRIPTION :
  --   Creates all partial derivatives of the polynomial in nd.p
  --   and stores them in the appropriate entries of nd.d.

  function Create ( p : Poly ) return Node;

  -- DESCRIPTION :
  --   Returns a node with a copy of p in its p field
  --   and all partial derivatives in its d field,
  --   completely recursively created.

  function Create ( p : Poly_Sys ) return Array_of_Nodes;

  -- DESCRIPTION :
  --   Applies the Create from above to all polynomials in p.
  --   The array on return has range p'range.

-- SELECTOR (and updator) :

  procedure Derivative ( nd : in out Node; v : in Vector; dp : out Poly );

  -- DESCRIPTION :
  --   Returns in dp the derivative of nd.p, derived v(i) times
  --   w.r.t. the i-th variable.  Whenever the derivatives have not yet
  --   been computed the node nd is adjusted with new derivatives.

-- ENUMERATORS :

  generic
    with procedure Process ( nd : in node; continue : out boolean );
  procedure Enumerate ( nd : in Node );

  -- DESCRIPTION :
  --   Applies Process to the node and continues with calling
  --   the Enumerate_Nodes below if continue is not set to false
  --   in the procedure Process.

  generic
    with procedure Process ( nd : in node; continue : out boolean );
  procedure Enumerate_Nodes ( nd : in Array_of_Nodes );

  -- DESCRIPTION :
  --   Applies the procedure Process to every nonempty node in nd.
  --   The enumeration stops when the continue in the Process
  --   has been set to false.

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node );
  procedure Clear ( nd : in out Link_to_Node );
  procedure Clear ( nd : in out Array_of_Nodes );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the nodes.

end QuadDobl_Derivative_Trees;
