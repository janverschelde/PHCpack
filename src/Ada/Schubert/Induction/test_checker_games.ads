with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Brackets;                          use Brackets;
with Bracket_Monomials;                 use Bracket_Monomials;
with Checker_Boards;                    use Checker_Boards;
with Checker_Posets;                    use Checker_Posets;

package Test_Checker_Games is

-- DESCRIPTION :
--   Tests the operations on checker posets to count all solutions
--   to Schubert problems with checker games.

  procedure Interactive_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user to enter a permutation of n elements
  --   and calculates the descending and rising checker.

  procedure Specialization_Order ( n : in integer32);

  -- DESCRIPTION :
  --   Shows all the moves in the specialization order,
  --   starting from the identity permutation on n elements.

  procedure Show_Specializing_Moves ( n : in integer32 );

  -- DESCRIPTION :
  --   Shows all moves of n black checkers in the specializing order.

  procedure Show_Generalizing_Moves ( n : in integer32 );

  -- DESCRIPTION :
  --   Shows all moves of n black checkers in the generalizing order.

  procedure Show_All_Moves ( n : in natural32 );

  -- DESCRIPTION :
  --   On a board with n black checkers we can either play the
  --   specialization or the generalization of the moves.

  procedure Initialize_Checkerboard
               ( p,rows,cols : out Vector; b : out Board );

  -- DESCRIPTION :
  --   Initialization of the checkerboard for k-planes in n-space.
  --   Checks the happiness of the white checkers and forces the
  --   user to re-enter the rows and columns of the white checkers
  --   in case the white checkers are not happy.

  -- ON RETURN :
  --   p         permutation defines the n black checkers;
  --   rows      rows of the k white checkers;
  --   cols      columns of the k white checkers;
  --   b         a checkerboard with n black and k white checkers.

  procedure Test_White_Moves ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Tests the calculation of the moves of the white checkers.

  procedure Group_Black_Checkers ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Groups the black checkes into four zones.

  procedure Write_Children ( p : in Vector; nd : in Node );

  -- DESCRIPTION :
  --   Writes the children of the node lnd, also displaying the
  --   configuration of black checkers encoded in the permutation p.

  procedure Pick_a_Node ( n,k : in integer32; ps : in Poset );

  -- DESCRIPTION :
  --   Prompts the user for coordinates (level and position) of a node
  --   in the poset.  All relevant information of that node is displayed.

  procedure Create_Poset ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Creates a poset for k-planes in n-space.

  procedure Create_Intersection_Poset
              ( n,m : in integer32; p,rows,cols : in Vector );

  -- DESCRIPTION :
  --   Creates an intersection poset for at most m intersection
  --   conditions on k-planes in n-space.  The initial configuration
  --   of black checkers is in the permutation p and the positions of
  --   the white checkers is defined by rows and cols.

  procedure Count_Planes ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Counts the number of k-planes in n-space which satisfy
  --   certain user given Schubert conditions.

  procedure Create_Intersection_Poset
              ( n : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Makes the intersection poset to resolve a general intersection 
  --   condition in n-space, as defined by the conditions in the monomial
  --   of brackets.

  procedure Resolve_Intersection_Condition ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a bracket monomial on a k-plane in n-space
  --   and resolves the entered Schubert intersection condition.

  procedure Write_Paths_in_Poset ( ps : in Poset );

  -- DESCRIPTION :
  --   Enuerations all paths defined by the poset.

  procedure Walk_to_First_Parent ( ps : in Poset );

  -- DESCRIPTION :
  --   Walks from leaves to the root of the poset,
  --   using only the first parent.

  procedure Enumerate_Paths_in_Poset ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a configuration of white checkers,
  --   creates a poset and writes all paths through the poset.

  procedure Flip ( b : in out Bracket );

  -- DESCRIPTION :
  --   Flips the order of the numbers in the bracket.

  function Partition_to_Bracket
             ( k,n : integer32; p : Bracket ) return Bracket;

  -- DESCRIPTION :
  --   Turns a partition p into a bracket, for a k-plane in n-space.

  function Random_Partition ( k,n : integer32 ) return Bracket;

  -- DESCRIPTION :
  --   Returns a decreasing sequence of k numbers with entries of at
  --   most n-k-1, generated uniformly with the first number at least 1.

  function Random_Complement ( k,n : integer32; p : Bracket ) return Bracket;

  -- DESCRIPTION :
  --   Returns an increasing list q of k numbers with the restriction that
  --   p[i] + q[i] < n-k for all i, the last element is at least 1.

  function Difference ( k,n : integer32; p,q : Bracket ) return Bracket;

  -- DECRIPTION :
  --   Returns the bracket with i-th entry n - k - p(i) - q(i).

  procedure Sort_into_Partition ( b : in out Bracket );

  -- DESCRIPTION :
  --   In a partition, the conditions appear in decreasing order,
  --   so the numbers in b will be sorted accordingly.

  function Random_Redistribute
             ( k,n : integer32; b : Bracket ) return Bracket;

  -- DESCRIPTION :
  --   Takes the sum of the conditions imposed by the bracket
  --   and redistributes this sum randomly into a partition.

  procedure Random_Triplet ( k,n : in integer32; roco : out Natural_Number );

  -- DECRIPTION :
  --   Returns 3 conditions on a k-plane in n-space with formal root count
  --   returned in roco.

  procedure Generate_Intersection_Problem
              ( k,n,m : in integer32; roco : out Natural_Number );

  -- DESCRIPTION :
  --   Generates m conditions on a k-plane in n-space
  --   and returns the final sum as root count in roco.

  procedure Random_Intersection_Problem ( k,n : integer32 );

  -- DESCRIPTION :
  --   Generates a random intersection problem on k-planes in n-space.

  procedure NotAbove ( k,n : in integer32 );

  -- DESCRIPTION :
  --   Checks for brackets not above a given bracket.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the ambient dimension n, for a test, and depending 
  --   on the test selection also for the dimension k of the planes.

end Test_Checker_Games;
