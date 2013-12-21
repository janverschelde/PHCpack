with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;

package Standard_Integer32_Simplices is

-- DESCRIPTION :
--   The simplices are facets on the lower hull of a lifted polytope
--   spanned by standard 32-bit integer vectors.

-- DATA STRUCTURE :

  type Simplex is private;
  Null_Simplex : constant Simplex;  -- empty simplex

-- CREATORS :

  function Create ( x : VecVec ) return Simplex;

  -- DESCRIPTION :
  --   Given the vertices, the simplex will be created.

  procedure Update ( s : in out Simplex; x : in Link_to_Vector;
                     k : in integer32 );

  -- DESCRIPTION :
  --   Creates a neighboring simplex by pivoting the kth point in the
  --   simplex by the given point x.

  procedure Update ( s : in out Simplex; x : in Link_to_Vector;
                     pos : in Vector );

  -- DESCRIPTION :
  --   Creates neighboring simplices by pivoting the points in the
  --   simplex by the given point x.

  -- ON ENTRY :
  --   s             simplex, volume(s) > 0;
  --   x             a certain point to be considered for update;
  --   pos           the position vector of x.
 
  -- ON RETURN :
  --   s             updated, so that the neighbors of s contain the new
  --                 simplices.

  procedure Update_One ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector );
  procedure Update_One ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector; news : out Simplex );


  -- DESCRIPTION :
  --   Creates at most one new simplex that contains x.
  --   It will be assumed that there are no internal points, so no more
  --   than one new simplex will be constructed.  When x already belongs
  --   to  s or to one of its neighbor, then nothing will be done.
  --   The specifications of s,x and pos are the same as the other
  --   Update operation, except for the new parameter news wich returns
  --   the new simplex that contains x.

  generic
    with procedure Process ( news : in Simplex; cont : out boolean );
  procedure Update_All ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector; ancestor : in Simplex );

  -- DESCRIPTION :
  --   Computes all new simplices that can be derived from the initial
  --   simplex s and the new point x.
  --   The ancestor simplex prevents the walk from turning back.

  procedure Connect ( s1,s2 : in out Simplex );
 
  -- DESCRIPTION :
  --   Computes the intersection of the two simplices.
  --   If they are neighbors to each other, then the appropriate
  --   connections will be made.

  procedure Flatten ( s : in out Simplex );

  -- DESCRIPTION :
  --   All lifting values will become zero and the inner normal changes
  --   into (0,0,..,0,1).

-- SELECTORS :

  function Dimension ( s : Simplex ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of points in the simplex.

  function Normal ( s : Simplex ) return Vector;

  -- DESCRIPTION :
  --   Returns the normal to the given simplex.

  function Is_Flat ( s : Simplex ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all components of the normal equal zero,
  --   except the last one, which must be equal to one.
  --   Returns false otherwise.

  function Vertices ( s : Simplex ) return VecVec;

  -- DESCRIPTION :
  --   Returns the vertices that span the simplex.

  function Vertex ( s : Simplex; k : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the kth vector that spans the simplex.

  function Is_Vertex ( s : Simplex; x : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x is one of the vertices that spans the simplex.

  function Equal ( s1,s2 : Simplex ) return boolean;

  -- DESCRIPTION :
  --   returns true if the representations of both simplices are the same.

  function Index ( s : Simplex; x : Vector ) return integer32;

  -- DESCRIPTION :
  --   k := Index(s,x);
  --   if k = 0 then Is_Vertex(s,x) = false,
  --   else Neighbor(s,k) is the neighboring simplex of s by pivoting x.

  function Neighbor ( s : Simplex; k : integer32 ) return Simplex;
  function Neighbor ( s : Simplex; k : integer32; pos : Vector ) return Simplex;

  -- DESCRIPTION :
  --   Returns the neighbor of the kth point of the given simplex.
  --   If the position vector is supplied, then it will only return
  --   a non-empty simplex if that is allowed by the position.

  function Position ( s : Simplex; x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Computes the position of the given vector w.r.t. that simplex,
  --   i.e., if l = Position(pt,s), then sum of l(i)*s(i) + l(l'last)*pt = 0,
  --   where s(i) = ith point that spans the simplex, and the sum
  --   of all entries in l equals 1.

  function Is_In ( s : Simplex; x : Vector ) return boolean;
  function Is_In ( pos : Vector ) return boolean;
  function Is_In_All ( s : Simplex; x : Vector ) return boolean;
  function Is_In_All ( s : Simplex; x : Vector ) return Simplex;
  function Is_In_All ( s : Simplex; x,pos : Vector ) return boolean;
  function Is_In_All ( s : Simplex; x,pos : Vector ) return Simplex;

  -- DESCRIPTION :
  --   Returns true if the point x is contained in the convex hull of s.
  --   The function runs more efficiently if the position vector is
  --   already available, otherwise it will be computed.
  --   With Is_In_All, also all neighbors of the given simplex are checked.
  --   Is_In_All(s,x) returns true if x is contained in s or in one of its
  --   neigbors, or in one of the neighbors of the neigbors...
  --   Either the simplex which contains x or the empty simplex is returned.

  generic

    with procedure Process_Neighbor
                      ( nei : in out Simplex; k : in integer32; 
                        continue : out boolean );

  procedure Neighbors ( s : in out Simplex; x : in Vector );

  -- DESCRIPTION :
  --   Computes all neighbors (and neighbors of ...) of the simplex s,
  --   that lie closest to x, i.e. with no other simplices in between.
  --   This procedure implements a walk from the simplex s to
  --   a given point x.
  --   As the parameters are of type 'in out', simplices are
  --   allowed to be updated.

  -- To the user defined procedure Process_Neighbor the following
  -- information is passed:

  -- ON ENTRY :
  --   nei          a neighboring simplex close to x;
  --   k            index for point in nei that can be replaced by s.

  -- ON RETURN :
  --   nei          updated simplex;
  --   continue     if true then more neighboring simplices may be delivered;
  --                if false then the iteration will stop.

  function Volume ( s : Simplex ) return natural32;

   -- DESCRIPTION :
   --   Returns n! times the volume of the simplex.

-- DESTRUCTORS :

  procedure Destroy_Neighbor  ( s : in out Simplex; k : in integer32 );
  procedure Destroy_Neighbors ( s : in out Simplex );

  -- DESCRIPTION :
  --   Sets the kth neighbor of the simplex to the empty simplex.
  --   If k is not specified, then all neighbors will be set to empty.

  procedure Clear_Neighbor  ( s : in out Simplex; k : in integer32 );
  procedure Clear_Neighbors ( s : in out Simplex );

  -- DESCRIPTION :
  --   Clears the kth neighbor of the simplex s.
  --   If k is not specified, then all neighbors will be cleared.

  procedure Clear ( s : in out Simplex );
   
  -- DESCRIPTION :
  --   Makes the allocated memory space available.

private

  type Simplex_Rep ( n : integer32 );
  type Simplex is access Simplex_Rep;
  Null_Simplex : constant Simplex := null;

end Standard_Integer32_Simplices;
