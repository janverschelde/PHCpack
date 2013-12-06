with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Integer64_VecMats;         use Standard_Integer64_VecMats;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Lists_of_Integer64_Vectors;
with Standard_Lattice_3d_Facets;         use Standard_Lattice_3d_Facets;

package Standard_Pretropical_Facets is

-- DESCRIPTION :
--   A pretropical facet is determined by a facet pretropism,
--   a vector perpendicular to a tuple of facets.
--   A genuine pretropical edge is determined by an edge pretropism,
--   perpendicular to a tuple of edges not adjacent to a pretropical facet.

-- converting supports into matrix formats :

  function List2Matrix ( A : Lists_of_Integer_Vectors.List ) return Matrix;

  -- DESCRIPTION :
  --   Returns the representation of the list A as a matrix,
  --   the columns of the matrix on return are the vectors in A.

  -- REQUIRED : not Is_Null(A).

  function Lists2VecMat ( A : Array_of_Lists ) return VecMat;

  -- DESCRIPTION :
  --   Returns the matrix representations of the lists of A.
  --   The columns of the matrices on return contain the coordinates
  --   of the points in the lists.

-- computing pretropical facets :

  function Facet_Pretropisms
              ( f,g : Facet_3d_List ) return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Returns the set of inner facet normals common to both f and g.

  function Facet_Pretropisms
              ( f : Array_of_Facet_3d_Lists )
              return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   A facet pretropism is a vector perpendicular to a tuple of facets,
  --   each facet in the tuple from a different polytope.

  function Pretropical_Facets
               ( f : Array_of_Facet_3d_Lists;
                 v : Lists_of_Integer64_Vectors.List ) 
               return Array_of_Facet_3d_Lists;

  -- DESCRIPTION :
  --   Given the tuple of complete facet lists in f and the list of facet
  --   pretropisms in v, this function returns the pretropical facets.

  function Is_Subset ( p,s : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if every element of s also occurs in p.

  function On_Pretropical_Edge
                ( A,B : Matrix; e,v : Standard_Integer64_Vectors.Vector )
               return boolean;

  -- DESCRIPTION :
  --   The vector v is on a pretropical edge e if its supports on A and B
  --   are entirely contained in the points supported by f.

  function On_Pretropical_Edge
                ( A,B : Matrix; e : Lists_of_Integer64_Vectors.List;
                  v : Standard_Integer64_Vectors.Vector )
               return boolean;

  -- DESCRIPTION :
  --   Returns true if v is on a pretropical edge defined by one of
  --   the inner normals in the list e.

  function Subset_Count
               ( p,s : Standard_Integer_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of indices of s that belong to p.

  function On_Pretropical_Facet
                ( A,B : Matrix; f,v : Standard_Integer64_Vectors.Vector )
               return boolean;

  -- DESCRIPTION :
  --   The vector v is on the pretropical facet defined by normal f if the
  --   sets supported by f and v share at least two points for both A and B.
  --   Returns true if v is on the pretropical facet, false otherwise.

  function On_Pretropical_Facet
                ( A,B : Matrix; f : Lists_of_Integer64_Vectors.List;
                  v : Standard_Integer64_Vectors.Vector )
               return boolean;

  -- DESCRIPTION :
  --   Returns true if v is on a pretropical facet defined by one of
  --   the inner normals in the list f.

  function On_Pretropical_Facet
                ( A : VecMat; f : Array_of_3d_Facets;
                  v : Standard_Integer64_Vectors.Vector )
                return boolean;

  -- DESCRIPTION :
  --   Returns true if the support of tuple A in the v direction
  --   contains at least two or more points of every facet in f. 

  function On_Pretropical_Facet
                ( A : VecMat; f : Array_of_Facet_3d_Lists;
                  v : Standard_Integer64_Vectors.Vector )
                return boolean;

  -- DESCRIPTION :
  --   A vector v is on a pretropical facet if it supports more than one
  --   point of each facet in the tuple of pretropical facets.

  -- REQUIRED : each list in f has the same length
  --   and the i-th facet belongs to the i-th tuple.

-- computing pretropical edges via Minkowski sum :

  function Minkowski_Sum ( A,B : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix with column range 1..A'last(2)*B'last(2),
  --   containing all possible combinations of sums of columns of A
  --   with columns of B.

  function Is_Mixed_Facet_Normal
               ( A,B : Matrix; v : Standard_Integer64_Vectors.Vector )
               return boolean;

  -- DESCRIPTION :
  --   Returns true if v supports at least two points for A and B.

  function Mixed_Facet_Normals
              ( A,B : Matrix; f : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of inner normals to mixed facets,
  --   spanned by at least one edge of A and one edge of B.

  function Edge_Tropisms_by_Sum
              ( A,B : Matrix ) return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Computes all tropisms to (A,B) via their Minkowski sum.
  --   If (A,B) has facet tropisms, then they will be returned as well.

-- pretropical edges of two supports via wrapping :

  procedure Check_Edge
              ( a,b,u,v : in Standard_Integer64_Vectors.Vector;
                uk,vk : out integer64 );

  -- DESCRIPTION :
  --   Computes coefficients uk and vk with u and v, so w = uk*u + vk*v 
  --   is perpendicular to the edge spanned by a and b.

  function Is_Tropism 
              ( A,B : Matrix;
                v : Standard_Integer64_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Computes the points v supports of A and B and returns true
  --   if the number of supported points equals at least 2 for both A and B.

  procedure Wrap_Edge_for_Tropism
              ( A,B : in Matrix; ip,iq : in integer32;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List;
                fail : out boolean );
  procedure Wrap_Edge_for_Tropism
              ( A,B : in Matrix; ip,iq : in integer32;
                fpt : in Lists_of_Integer64_Vectors.List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List;
                fail : out boolean );

  -- DESCRIPTION :
  --   An edge spanned by (ip,iq) of B defines a tropism if it has an inner
  --   normal spanned by the vectors u and v.  If fpt is given on entry,
  --   the tropism is checked for being on a pretropical facet.

  -- ON ENTRY :
  --   A        support of the first polytope;
  --   B        support of the second polytope;
  --   ip       index to a vertex of B with inner normal spanned by u and v;
  --   iq       index to a vertex of B connected to ip via an edge;
  --   fpt      list of inner normals to pretropical facets;
  --   u        inner normal to a facet of A;
  --   v        inner normal to neighbor facet of A.
  --   trp      list of current tropisms;
  --   trp_last pointer to last element in trp.

  -- ON RETURN :
  --   trp      updated list of tropisms;
  --   trp_last pointer to last element in trp;
  --   fail     true if (ip,iq) did not define a tropism,
  --            false if (ip,iq) defines a tropism.

  procedure Wrap_Edges
              ( A,B : in Matrix; f : in Facet_3d_List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List );
  procedure Wrap_Edges
              ( A,B : in Matrix; f : in Facet_3d_List;
                fpt : in Lists_of_Integer64_Vectors.List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List );

  -- DESCRIPTION :
  --   Applies the giftwrapping idea to find a tropism of (A,B).
  --   The version without fpt does not take facet pretropisms into account.

  -- ON ENTRY :
  --   A        support of the first polytope;
  --   B        support of the second polytope;
  --   f        list of facets for B;
  --   fpt      list of normals to facet pretropisms;
  --   u        inner normal to a facet of A;
  --   v        inner normal to neighbor facet of A.
  --   trp      list of current tropisms;
  --   trp_last pointer to last element in trp.

  -- ON RETURN :
  --   trp      updated list of tropisms;
  --   trp_last pointer to last element in trp.

  procedure Visit_Edges
              ( A,B : in Matrix; f : in Link_to_3d_Facet; g : in Facet_3d_List;
                w : in Lists_of_Integer64_Vectors.List;
                t,t_last : in out Lists_of_integer64_Vectors.List;
                cnt : in out natural32 );
  procedure Visit_Edges
              ( A,B : in Matrix; f : in Link_to_3d_Facet; g : in Facet_3d_List;
                fpt,w : in Lists_of_Integer64_Vectors.List;
                t,t_last : in out Lists_of_integer64_Vectors.List;
                cnt : in out natural32 );

  -- DESCRIPTION :
  --   Visits all edges of the facet f of the support B,
  --   avoiding edges neighbor to already visited facets.
  --   If fpt is given, then facet pretropisms are not reported.

  -- ON ENTRY :
  --   A        support of the first polytope;
  --   B        support of the second polytope;
  --   f        current facet visited;
  --   g        list of facets of B;
  --   fpt      iist of facet pretropisms;
  --   w        list of already visited facets;
  --   t        list of current pretropisms;
  --   t_last   pointer to last element of t;
  --   cnt      current count of visited edges.

  -- ON RETURN :
  --   t        updated list of pretropisms;
  --   t_last   pointer to last element of t;
  --   cnt      updated count of visited edges.

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix; f,g : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Computes all tropisms to edges spanned by A and B.
  --   The facets of A and B are respectively stored in f and g.
  --   This function does not consider pretropical facets and will return
  --   redundant normals in case A and B are not in general position.

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix; fpt : Lists_of_Integer64_Vectors.List;
                f,g : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Computes all tropisms to edges spanned by A and B, taking into
  --   account the normals to pretropical facets given in fpt.
  --   The facets of A and B are respectively stored in f and g.

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix ) return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Computes all tropisms perpendicular to one edge of A
  --   and one edge of B via giftwrapping.

  procedure Edge_Tropisms_by_Wrap
              ( A,B : in Matrix;
                fpt,ept : out Lists_of_Integer64_Vectors.List );

  -- DESCRIPTION :
  --   Applies wrapping to compute face and edge pretropisms.

  -- ON ENTRY :
  --   A        support set of the first polytope;
  --   B        support set of the second polytope.

  -- ON RETURN :
  --   fpt      facet pretropisms;
  --   ept      pretropisms to edges.

-- computing pretropical edges :

  procedure Wrap ( A : in VecMat; g : in Array_of_Facet_3d_Lists;
                   u,v : in Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   Attempts to compute a positive combination of the vectors u and v
  --   so the support of this vector is an edge of A.

  procedure Report_Edges
              ( A : in VecMat; f : in Link_to_3d_Facet;
                g : in Array_of_Facet_3d_Lists;
                v,w : in Lists_of_Integer64_Vectors.List;
                cnt : in out natural32 );

  -- DESCRIPTION :
  --   Reports all edges of the facet f, avoiding edges that are neighbors
  --   to facets with inner normals in v or in w.
  --   Updates the counter with every good edge.

  -- ON ENTRY :
  --   A        tuple of supports;
  --   f        facet to enumerate edges from;
  --   g        tuples of pretropical facets;
  --   v        facet pretropisms;
  --   w        inner normals to already visited facets;
  --   cnt      current count of the number of good edges.

  -- ON RETURN :
  --   cnt      updated counter of good edges.

  procedure Edge_Pretropisms
              ( A : in VecMat; f,g : in Array_of_Facet_3d_Lists;
                fpt : in Lists_of_Integer64_Vectors.List );

  -- DESCRIPTION :
  --   Computes pretropisms perpendicular to edges of each polytope.

  -- ON ENTRY :
  --   f        lists of facets to each polytope;
  --   g        tuples of pretropical facets;
  --   fpt      inner normals to facet pretropisms.

end Standard_Pretropical_Facets;
