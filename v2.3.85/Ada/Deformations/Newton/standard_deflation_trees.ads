with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Deflation_Trees is

-- DESCRIPTION :
--   A deflation tree stores several deflations of one polynomial system.
--   The root of the tree is the original system.
--   The level in the tree corresponds to the number of deflations.
--   A node in the tree can have as many children as the number
--   of variables, the position of a child in the array equals the
--   number of multipliers in the deflation.

  type Node;
  type Link_to_Node is access Node;
  type Array_of_Nodes is array ( integer32 range <> ) of Link_to_Node;

  type Node ( ne,nv : integer32 ) is record
    d : natural32;      -- number of deflations
    s : Poly_Sys(1..ne);
    f : Eval_Poly_Sys(1..ne);
    jm : Jaco_Mat(1..ne,1..nv);
    jf : Eval_Jaco_Mat(1..ne,1..nv);
    children : Array_of_Nodes(1..nv);
    sols,last : Solution_List;
  end record;

  function Create_Root ( p : Poly_Sys ) return Node;

  -- DESCRIPTION :
  --   Returns all the fields of the node, created with p.

  procedure Create_Child ( nd : in out Node; m : in integer32 );

  -- DESCRIPTION :
  --   Creates a child of the current node, adding m multipliers.

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node );
  procedure Clear ( nd : in out Link_to_Node );
  procedure Clear ( nd : in out Array_of_Nodes );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the nodes.

end Standard_Deflation_Trees;
