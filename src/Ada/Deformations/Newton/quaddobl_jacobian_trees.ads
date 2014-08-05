with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecMats;          use QuadDobl_Complex_VecMats;
with QuadDobl_Complex_Jaco_Matrices;    use QuadDobl_Complex_Jaco_Matrices;

package QuadDobl_Jacobian_Trees is

-- DESCRIPTION :
--   This package allows to store all Jacobian matrices of a system
--   of m polynomial equations in n unknowns with double double
--   complex coefficients.

-- DATA STRUCTURES :

  type Node;
  type Link_to_Node is access Node;
  type Array_of_Nodes is array ( integer32 range <> ) of Link_to_Node;

  type Node ( n : integer32 ) is record   -- n is number of variables
    a : Link_to_Jaco_Mat;                 -- Jacobian matrix a
    d : Array_of_Nodes(1..n);             -- d(i) is derivative of a wrt i
  end record;

  type Eval_Node;
  type Link_to_Eval_Node is access Eval_Node;
  type Array_of_Eval_Nodes is
    array ( integer32 range <> ) of Link_to_Eval_Node;

  type Eval_Node ( n : integer32 ) is record -- n is number of variables
    a : Link_to_Jaco_Mat;                    -- Jacobian matrix
    f : Link_to_Eval_Jaco_Mat;               -- Horner form of Jacobi matrix
    d : Array_of_Eval_Nodes(1..n);           -- a(i) is derivative wrt i
  end record;

  type Jacobian_Remember_Table is
    array ( integer32 range <> ) of Link_to_VecMat;

-- CONSTRUCTORS :

  function Initialize ( a : Jaco_Mat ) return Node;
  function Initialize ( a : Jaco_Mat ) return Eval_Node;

  -- DESCRIPTION :
  --   Returns a node with in its a field a pointer to the given a
  --   and all partial derivatives in d initialized to null.

  procedure Create ( nd : in out Node );
  procedure Create ( nd : in out Eval_Node );

  -- DESCRIPTION :
  --   Creates all partial derivatives of the Jacobian matrix in nd.a
  --   and stores them in the appropriate entries of nd.d.

  function Create ( a : Jaco_Mat ) return Node;
  function Create ( a : Jaco_Mat ) return Eval_Node;

  -- DESCRIPTION :
  --   Returns a node with a link to a in its a field
  --   and all partial derivatives in its d field,
  --   completely recursively created.

  function Diff ( a : Jaco_Mat; k : integer32 ) return Jaco_Mat;
  function Diff ( a : Link_to_Jaco_Mat; k : integer32 )
                return Link_to_Jaco_Mat;

  -- DESCRIPTION :
  --   Computes the derivative of every polynomial in the matrix a
  --   with respect to the k-th variable.

-- SELECTORS :

  function Degree ( a : Jaco_Mat; k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the maximal degree of all polynomials in the matrix a
  --   in the k-th variable.

  procedure Derivative ( nd : in out Node;
                         v : in Standard_Natural_Vectors.Vector;
                         da : out Link_to_Jaco_Mat );
  procedure Derivative ( nd : in out Eval_Node;
                         v : in Standard_Natural_Vectors.Vector;
                         da : out Link_to_Jaco_Mat;
                         df : out Link_to_Eval_Jaco_Mat );

  -- DESCRIPTION :
  --   Returns in da the derivative of nd.a, derived v(i) times
  --   w.r.t. the i-th variable.  Whenever the derivatives have not yet
  --   been computed the node nd is adjusted with new derivatives.
  --   Note that the matrices in da and df are just links, it is shared with
  --   the data in the tree under the node nd.

  procedure Derivative ( nd : in out Node; k : in integer32 );
  procedure Derivative ( nd : in out Eval_Node; k : in integer32 );

  -- DESCRIPTION :
  --   Updates the node with all k-th derivatives of the matrix in nd.a.

  procedure Dimensions ( v : in Link_to_VecMat; rows,cols : out integer32 );

  -- DESCRIPTION :
  --   Returns number of rows and columns in the first nonvoid matrix in v.

-- ENUMERATORS :

  generic
    with procedure Process ( nd : in Node; continue : out boolean );
  procedure Enumerate_Nodes ( nd : in Node );
  generic
    with procedure Process ( nd : in Eval_Node; continue : out boolean );
  procedure Enumerate_Eval_Nodes ( nd : in Eval_Node );

  -- DESCRIPTION :
  --   Applies Process to the node and continues with calling
  --   the Enumerate_Nodes below if continue is not set to false
  --   in the procedure Process.

-- EVALUATORS :

  procedure Create_Remember_Derivatives
               ( a : in Jaco_Mat; k : in integer32;
                 nd : out Link_to_Eval_Node );

  -- DESCRIPTION :
  --   Creates a remember table for all k-th order derivatives of the
  --   Jacobian matrix in a.

  -- ON ENTRY :
  --   a         Jacobian matrix of a polynomial system;
  --   k         order of the derivative.

  -- ON RETURN :
  --   nd        remember table for Horner forms of all derivatives
  --             of the Jacobian matrix a, up to order k.

  function Evaluate_Jacobian_Remember_Table
               ( nd : Link_to_Eval_Node; k,n : integer32;
                 x : QuadDobl_Complex_Vectors.Vector )
               return Jacobian_Remember_Table;

  -- DESCRIPTION :
  --   Returns a remember table for the derivatives of the Jacobian
  --   matrix (up to order k) evaluated at x.

  -- ON ENTRY :
  --   nd        remember table for all Jacobian matrices up to order k;
  --   k         order of the deflation matrices;
  --   n         number of variables in the Jacobian matrices;
  --   x         values for the variables in the Jacobian matrices.

  -- ON RETURN :
  --   An array of range 0..k with vectors of matrices,
  --   the i-th contains all derivatives of order i.

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node );
  procedure Clear ( nd : in out Eval_Node );
  procedure Clear ( nd : in out Link_to_Node );
  procedure Clear ( nd : in out Link_to_Eval_Node );
  procedure Clear ( nd : in out Array_of_Nodes );
  procedure Clear ( nd : in out Array_of_Eval_Nodes );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the nodes.

  procedure Clear ( jrt : in out Jacobian_Remember_Table );

  -- DESCRIPTION :
  --   Clears all values in the Jacobian remember table.

end QuadDobl_Jacobian_Trees;
