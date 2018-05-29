with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with C_Integer_Arrays;                   use C_Integer_Arrays;
with C_Double_Arrays;                    use C_Double_Arrays;

package DEMiCs_Algorithm is

-- DESCRIPTION :
--   This package provides an interface to DEMiCs,
--   to compute all mixed cells by dynamic enumeration.
--   DEMiCs was developed by Tomohiko Mizutani, Akiko Takeda, and
--   Masakazu Kojima and licensed under GNU GPL Version 2 or higher.

  procedure Extract_Supports 
              ( p : in Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );
  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the supports and computes the type of mixture.

  procedure Extract_Supports 
              ( p : in Poly_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );
  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the supports and computes the type of mixture.
  --   In addition, returns the permutation to put the same supports
  --   in consecutive order.

  procedure Add_Artificial_Origin
              ( dim : in integer32;
                sup : in out Lists_of_Integer_Vectors.List;
                added : out boolean );

  -- DESCRIPTION :
  --   If the origin does not belong to the list of points sup,
  --   then the origin is appended to sup and added is true on return.
  --   The dimension dim on entry equals the dimension of the points.

  procedure Add_Artificial_Origins
              ( dim : in integer32;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                nbadd : out integer32;
                added : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Add artificial origins to the supports for the stable mixed volume.

  -- ON ENTRY :
  --   dim      dimension of the points in the supports;
  --   sup      supports of a polynomial system.

  -- ON RETURN :
  --   nbadd    number of artificial origins added;
  --   added    vector of sup'range to indicate where the artificial origins
  --            are in sup: sup(k) = 1 then added to k-th support, 0 if not.

  function Mixture_Type
              ( mix : Standard_Integer_Vectors.Vector )
              return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the type of mixture as an array suitable to pass to C.

  function Cardinalities
              ( mix : Standard_Integer_Vectors.Vector;
                sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
              return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the length of each support list in sup
  --   as an array suitable to pass to C.
  --   The type of mixture is given in mix.

  function Number_of_Points
              ( mix : Standard_Integer_Vectors.Vector;
                sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
              return integer32;

  -- DESCRIPTION :
  --   Returns the total number of points in the supports sup.

  function Coordinates
              ( nbr : integer32;
                mix : Standard_Integer_Vectors.Vector;
                sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
              return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the length of each support list in sup
  --   as an array suitable to pass to C.

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns random lifting values for each of the support sets,
  --   as arranged in the order of their mixture, defined by mix.

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               stlb : double_float;
               added : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns random lifting values for each of the support sets,
  --   assigning the stable lifting bound to the artifial origins,
  --   added to compute the stable mixed volume.

  -- ON ENTRY :
  --   mix     type of mixture;
  --   sup     supports of a polynomial system;
  --   stlb    stable lifting bound, for the stable mixed volume;
  --   added   if added(k), then the artificial origin has been
  --           added to the k-th support set, as returned by the
  --           procedure Add_Artificial_Origins.

  function Random_Lifting ( nbr : integer32 ) return C_Double_Array;

  -- DESCRIPTION :
  --   Given in nbr is the total number of points in the supports.
  --   On return is a vector of range 0..nbr-1 with random doubles.
  --   The random doubles are in the interval [0, 1].

  function Size ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                return integer32;

  -- DESCRIPTION :
  --   Returns the total number of elements in v.
 
  function Flatten ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                   return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector with all values in v.

  function Copy_Lifting
              ( lif : Standard_Floating_Vectors.Vector )
              return C_Double_Array;

  -- DESCRIPTION :
  --   Returns a C double array of range 0..lif'last-1
  --   which contains the copied lifting values of lif.

  function run_demics
              ( v,d,r : integer32; mtp,crd,pts : C_Integer_Array )
              return integer32;
  pragma import(C, run_demics, "demicsrun");

  -- DESCRIPTION :
  --   Interface to the C++ function demicsrun.
  --   In this version, demics generates the random lifting values.

  -- ON ENTRY :
  --   v        0 or 1 whether silent or verbose;
  --   d        dimension;
  --   r        number of distinct support sets;
  --   mtp      mixture type;
  --   crd      cardinalities of each support set;
  --   pts      coordinates of the points in the each support.

  function fly_demics
              ( v,d,r : integer32; mtp,crd,pts : C_Integer_Array;
                lif : C_Double_Array ) return integer32;
  pragma import(C, fly_demics, "demicsfly");

  -- DESCRIPTION :
  --   Interface to the C++ function demicsrun,
  --   with given lifting values, which better be random numbers.

  -- ON ENTRY :
  --   v        0 or 1 whether silent or verbose;
  --   d        dimension;
  --   r        number of distinct support sets;
  --   mtp      mixture type;
  --   crd      cardinalities of each support set;
  --   pts      coordinates of the points in the each support;
  --   lif      random lifting values for each point,
  --            in an array of size equal to the total number of points.

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs on the given supports, of type of mixture mix.

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs on the given supports, of type of mixture mix.
  --   If stable, then the stable mixed volume will be computed
  --   and stlb is the stable lifting bound.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   stable   if the stable mixed volume is wanted;
  --   stlb     stable lifting bound if stable, otherwise 0.0;
  --   verbose  for extra output.

  -- ON RETURN :
  --   sup      supports with artificial origins added if stable.

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifting : in Standard_Floating_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs on the given supports, of type of mixture mix,
  --   with given random lifting values.

  -- ON ENTRY :
  --   mix      type of mixture of the supports.
  --   supports supports of a polynomial system;
  --   lifting  a vector of range 1..Number_of_Points(mix.all,supports)
  --            with sufficiently random values for the lifting;
  --   verbose  if additional output is desired.

  procedure Show_Output;

  -- DESCRIPITON :
  --   Shows the output stored in DEMiCs_Output_Data.

  function Extract_Indices ( s : string ) return string;

  -- DESCRIPTION :
  --   Scans s to the first colon and then returns the rest of s,
  --   all characters of s after the first colon.

  procedure Process_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : out Mixed_Subdivision;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Converts the data computed by DEMiCs to lifted supports
  --   and a mixed cell configuration.

  -- ON ENTRY :
  --   dim      dimension of the points before lifting;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag to indicate if extra output is wanted.

  -- ON RETURN :
  --   lif      lifted support sets;
  --   mcc      a regular mixed-cell configuration.

end DEMiCs_Algorithm;
