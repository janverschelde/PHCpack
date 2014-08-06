with text_io;                          use text_io;
with Generic_Lists;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural64_VecVecs;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;        use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;         use QuadDobl_Complex_VecMats;
with QuadDobl_Complex_Poly_SysFun;     use QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;   use QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Jacobian_Trees;          use QuadDobl_Jacobian_Trees;

package QuadDobl_Evaluate_Deflation is

-- DESCRIPTION :
--   This package facilitates the efficient evaluation of the Jacobian
--   matrices which arise in Newton's method with deflation.

-- DATA STRUCTURES :

  type Eval_Tree;
  type Link_to_Eval_Tree is access Eval_Tree;
  type Array_of_Eval_Trees is
    array ( integer32 range <> ) of Link_to_Eval_Tree;

  type Eval_Tree ( k,m : integer32 ) is record
    key : natural32;                              -- order of node in tree
    d : Standard_Natural_Vectors.Vector(0..k);    -- derivative applied to A(m)
    c : Array_of_Eval_Trees(0..m);                -- children needed for A(m)
    e : Standard_Integer_Vectors.Vector(0..m);    -- labels for children
    v : QuadDobl_Complex_Matrices.Link_to_Matrix; -- value of A(m)
  end record;

  package List_of_Nodes is new Generic_Lists(Link_to_Eval_Tree);
  type Node_List is new List_of_Nodes.List;

-- OPERATIONS :

  function Create ( k : natural32 ) return Link_to_Eval_Tree;

  -- DESCRIPTION :
  --   Unwinds the application of multipliers in A(k) in the tree.
  --   Only generates nodes which do not yet occur.

  function Is_In ( t : Eval_Tree; d : Standard_Natural_Vectors.Vector;
                   m : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the derivative operator d applied to A(m)
  --   already belongs to the tree t.

  function Key_In ( t : Eval_Tree; d : Standard_Natural_Vectors.Vector;
                    m,max_key : integer32 ) return integer32;

  -- DESCRIPTION :
  --   If d applied to A(m) belongs to a node of t with key <= max_key,
  --   then this key is returned, otherwise -1 is returned.

  function Key_In ( l : Node_List; d : Standard_Natural_Vectors.Vector;
                    m,max_key : integer32 ) return integer32;

  -- DESCRIPTION :
  --   If d applied to A(m) belongs to a node of t with key <= max_key,
  --   then this key is returned, otherwise -1 is returned.

  function Look_Up ( t : Link_to_Eval_Tree; key : integer32 )
                   return Link_to_Eval_Tree;

  -- DESCRIPTION :
  --   Searches in the tree for a node with the given key.

  function Is_Leaf ( t : Eval_Tree ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all children of t are null.

  function Node_Count ( evt : Eval_Tree ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nodes in the tree.

  function Different_Node_Count ( evt : Eval_Tree ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of different nodes in the tree.

-- ENUMERATOR to unwind the multipliers in the deflation matrix :

  generic

    with procedure Multiplier_Derivative
          ( d : in Standard_Natural_Vectors.Vector; i,depth : in natural32;
            zero : in boolean; continue : out boolean );

    -- DESCRIPTION :
    --   Called in enumeration with new multiplier operator.

    -- ON ENTRY :
    --   d        derivative operator d;
    --   i        operator d is applied to A(i);
    --   depth    depth in the enumeration tree;
    --   zero     when true, instead of the derivative operator,
    --            zero should be printed.

    -- ON RETURN :
    --  continue  when false, the enumeration stops,
    --            otherwise it continues.

  procedure Enumerate_Multiplier_Derivatives ( k : in natural32 );

  -- DESCRIPTION :
  --   Recursively unwinds all derivatives in A(k), the Jacobian matrix
  --   of the k-th stage in the deflation.  Each time a new operator is
  --   created, the procedure Multiplier_Derivative is called.

  procedure Update_Stack
               ( sd : in out Standard_Natural_VecVecs.VecVec;
                 sk : in out Standard_Natural_Vectors.Vector;
                 sz : in out natural32;
                 d : in Standard_Natural_Vectors.Vector;
                 k : in natural32; found : out boolean );

  -- DESCRIPTION :
  --   Updates the stack of symbols with d and k.
  --   This procedure maintains the remember table for an efficient
  --   execution of the enumerations of the multiplier derivatives.

  -- REQUIRED : sz < sd'last.

  -- ON ENTRY :
  --   sd        stack of derivative operators;
  --   sk        corresponding k's: sd(i) applies to A(sk(i));
  --   sz        current size of the stack, maximal index,
  --             first-time call is with sz = sd'first-1;
  --   d         new derivative operator;
  --   k         d applies to A(k).

  -- ON RETURN :
  --  sd         unchanged if found, otherwise d is added to sd;
  --  sk         unchanged if found, otherwise k is added to sk;
  --  sz         unchanged if found, otherwise incremented by one;
  --  found      true if (d,k) belonged already to (sd,sk),
  --             false otherwise.

-- EVALUATORS :

  function Eval ( f : Eval_Poly_Sys; jm : Eval_Jaco_Mat;
                  t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec; k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : QuadDobl_Complex_VecVecs.VecVec )
                 return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of the k-th deflation of the system f at x.

  -- ON ENTRY :
  --   f         original system of nq(0) equations in nv(0) variables;
  --   jm        Jacobian matrix of the original system;
  --   t         directed acyclic graph with derivative operators;
  --   nd        remember table for Jacobian matrices up to order k,
  --             created with the procedure Create_Remember_Derivatives;
  --   monkeys   monomial key codes for every Jacobian matrix;
  --   k         order of the deflation matrix;
  --   nv        nv(i) is the number of columns in the i-th deflation matrix,
  --             nv is short for the number of variables;
  --   nq        nq(i) is the number of row in i-th deflation matrix,
  --             nq is short for the number of equations;
  --   R1        R1(i) is the number of multipliers in i-th stage,
  --             R1 is short for the rank plus one;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers.

  -- ON RETURN :
  --   vector of range 1..nq(k) with the values of the k-th deflated system.  

  function Eval ( t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec; k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : QuadDobl_Complex_VecVecs.VecVec )
                return Matrix;

  function Eval ( file : file_type;
                  t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec; k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : QuadDobl_Complex_VecVecs.VecVec )
                return Matrix;

  -- DESCRIPTION :
  --   Returns the value of the k-th deflation matrix at x.

  -- ON ENTRY :
  --   file      for writing diagnostics to file;
  --   t         directed acyclic graph with derivative operators;
  --   nd        remember table for Jacobian matrices up to order k,
  --             created with the procedure Create_Remember_Derivatives;
  --   monkeys   monomial key codes for every Jacobian matrix;
  --   k         order of the deflation matrix;
  --   nv        nv(i) is the number of columns in the i-th deflation matrix,
  --             nv is short for the number of variables;
  --   nq        nq(i) is the number of row in i-th deflation matrix,
  --             nq is short for the number of equations;
  --   R1        R1(i) is the number of multipliers in i-th stage,
  --             R1 is short for the rank plus one;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers.

  -- ON RETURN :
  --   matrix with nq(k) rows and nv(k) columns with the values of the
  --   Jacobian matrix of the k-th deflated system.

-- DESTRUCTORS :

  procedure Clear ( nl : in out Node_List );
  procedure Clear ( evt : in out Eval_Tree );
  procedure Clear ( evt : in out Link_to_Eval_Tree );

  -- DESCRIPTION :
  --   Deallocation of memory.

end QuadDobl_Evaluate_Deflation;
