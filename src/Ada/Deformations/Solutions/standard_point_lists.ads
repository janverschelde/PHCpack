with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Generic_Lists;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Point_Lists is

-- DESCRIPTION :
--   This package defines lists of points in the plane, obtained by
--   projecting standard complex solution lists of polynomial systems.

-- DATA STRUCTURES :

  type Point is record
    label : integer32;
    x,y : double_float;
  end record;

  type Link_to_Point is access Point;

  package List_of_Points is new Generic_Lists(Link_to_Point);
  type Point_List is new List_of_Points.List;

-- HASH FUNCTION :

  function Complex_Inner_Product
             ( x,y : Standard_Complex_Vectors.Vector )
	     return Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex inner product, multiplying the components of x
  --   with the conjugated components of y and returning the sum.

-- CREATORS :

  function Create ( sv,h1,h2 : Standard_Complex_Vectors.Vector;
                    label : integer32 ) return Point;
  function Create ( ls : Link_to_Solution;
                    h1,h2 : Standard_Complex_Vectors.Vector;
                    label : integer32 ) return Point;

  -- DESCRIPTION :
  --   Applies the hash function to create a point.

  -- ON ENTRY :
  --   sv        solution vector of a polynomial system;
  --   ls        pointer to a solution of a polynomial system;
  --   h1        hash function to define the x-coordinate of the point;
  --   h2        hash function to define the y-coordinate of the point;
  --   label     label used to indentify the point.

  procedure Append ( first,last : in out Point_List;
                     h1,h2 : in Standard_Complex_Vectors.Vector;
                     label : in integer32;
                     sv : in Standard_Complex_Vectors.Vector );
  procedure Append ( first,last : in out Point_List;
                     h1,h2 : in Standard_Complex_Vectors.Vector;
                     label : in integer32; ls : in Link_to_Solution );

  -- DESCRIPTION :
  --   Applies the hash function defined by h1 and h2 to the solution 
  --   vector ls.v and appends the point to the list first.

  -- ON ENTRY :
  --   first     pointer to the first element in the point list;
  --   last      pointer to the last element in the point list;
  --   h1,h2     defines a hash function to create a point;
  --   label     identification label for the solution;
  --   sv        solution vector of a polynomial system;
  --   ls        pointer to a solution.
 
  -- ON RETURN :
  --   first     updated pointer to the first element in the point list;
  --   last      updated pointer to the last element in the point list.

-- SELECTOR :

  procedure Center ( pl : in Point_List; cx,cy : out double_float );

  -- DESCRIPTION :
  --   Returns in (cx,cy) the coordinates of the center of the point list,
  --   computed as the average of all coordinates of the points in the list.

-- SORTING :

  function "<" ( lp1,lp2 : Link_to_Point ) return boolean;

  -- DESCRIPTION :
  --   Returns true if lp1 < lp2.

  procedure Swap ( lp1,lp2 : in out Link_to_Point );

  -- DESCRIPTION :
  --   The points lp1 and lp2 are swapped.

  procedure Sort ( pl : in out Point_List );

  -- DESCRIPTION :
  --   Sorts the list in increasing order, first on x, then on y.

  procedure Insert ( pl : in out Point_List; pt : in Link_to_Point );

  -- DESCRIPTION :
  --   For a given list pl, sorted in increasing order,
  --   the point pt is inserted into pl so that pl remains sorted.

  function Equal ( lp1,lp2 : Link_to_Point;
                   tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the x and y coordinate are within tol distance.

  procedure Insert_no_Duplicates
              ( pl : in out Point_List; pt : in Link_to_Point;
                tol : in double_float; lbl : out integer32 );

  -- DESCRIPTION :
  --   For a given list pl, sorted in increasing order,
  --   the point pt is inserted into pl so that pl remains sorted.
  --   If the point pt is equal to any other point q in pl
  --   (with respect to the tolerance), then the point is not
  --   inserted and on return the value of the label lbl is q.label.
  --   If the point pt is inserted, the lbl equals pt.label on return.

  procedure Insert_with_Duplicates
              ( pl : in out Point_List; pt : in Link_to_Point;
                tol : in double_float; cnt : out integer32;
                ptpl  : out Point_List );

  -- DESCRIPTION :
  --   For a given list pl, sorted in increasing order,
  --   the point pt is inserted into pl so that pl remains sorted.
  --   If the point pt is equal to any other point q in pl
  --   (with respect to the tolerance), then the point is
  --   inserted and on return the value of the cnt counts
  --   the number of occurrences of the point pt in the list.

  -- ON ENTRY :
  --   pl       list of hashed coordinates as points;
  --   pt       hash point;
  --   tol      tolerance to decide whether two points are equal.

  -- ON RETURN :
  --   cnt      number of occurrences of the point pt in the list
  --            is an upper bound on the multiplicity;
  --   ptpl     points to the place in pl where pt was inserted.

  generic
    with procedure Report ( lp1,lp2 : in Link_to_Point );
  procedure Clusters ( pl : in Point_List; tol : in double_float );

  -- DECRIPTION :
  --   Given a sorted point list and a tolerance, this procedure calls
  --   the procedure Report for every pair of points equal up to tol.

-- DESTRUCTORS :

  procedure Clear ( lp : in out Link_to_Point );

  -- DESCRIPTION :
  --   Deallocation of the memory occupied by lp.

  procedure Shallow_Clear ( pl : in out Point_List );
  procedure Deep_Clear ( pl : in out Point_List );

  -- DESCRIPTION :
  --   Deallocation of the memory occupied by the variables.

end Standard_Point_Lists;
