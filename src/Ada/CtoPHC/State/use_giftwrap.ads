with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_giftwrap ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer;
                        vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Gateway to the giftwrapping method in PHCpack.

-- ON ENTRY :
--   job   =   1 : convex hull in the plane, on input are
--                 in a: the number of characters in the string representation
--                       of the point configuration,
--                 in b: as many integers as there are characters in the
--                       string representation of the point configuration,
--                 on return are 
--                 in a: the number of characters of the representation of
--                       the string representation of the convex hull
--                 in b: is a string representation of a tuple:
--                 the vertices and the inner normals, both represented
--                 as string representations of point configurations.
--   job   =   2 : convex hull in 3-space or 4-space, on input are
--                 in a: the number of characters in the string representation
--                       of the point configuration,
--                 in b: as many integers as there are characters in the
--                       string representation of the point configuration,
--                 this operation initializes the giftwrap container with
--                 the list of facets of the convex hull.
--   job   =   3 : if given in a either 3 or 4, then the number of facets
--                 of a convex hull in 3-space or 4-space respectively is
--                 returned in b.
--   job   =   4 : if given in a either 3 or 4 and in b the label of a facet,
--                 then the string representation of the facet is returned
--                 in b, with in a the number of characters in the string.
--   job   =   5 : clear the list of facets for the convex hull in 3-space.
--   job   =   6 : clear the list of facets for the convex hull in 4-space.
--   job   =   7 : returns in a the number of characters in the string
--                 representation of a support set of the first
--                 polynomial in the Laurent systems container;
--   job   =   8 : given in a the number computed in job 7,
--                 and with sufficient space allocated in b,
--                 returns the string representation of the support in b;
--   job   =   9 : deallocates the string computed in job 7;
--   job   =  10 : replaces the system in the Laurent systems container
--                 with its initial form as defined by the inner normal
--                 with coordinates in b and number of variables in a.
--                 The normal is given as a string of a Python tuple.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
