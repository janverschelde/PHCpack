with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_VecVecs;            use Standard_Natural_VecVecs;

package Monodromy_Partitions is

-- DESCRIPTION :
--   This package helps in the management of partitions created by
--   the monodromy breakup algorithms.

  function Init_Factors ( d : natural32 ) return Link_to_VecVec;

  -- DESCRIPTION :
  --   Initializes the data structure to hold the factors of a polynomial
  --   of degree d as d singletons.

  function Map ( t1,t2 : Standard_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector;
  function Map ( t1,t2 : DoblDobl_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector;
  function Map ( t1,t2 : QuadDobl_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   If mp = Map(t1,t2,tol), then t1(i) = t2(mp(i)).  The tolerance tol is
  --   used to decide whether two complex numbers are equal to each other.

  procedure Write_Map ( file : in file_type;
                        map : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the map on file.

  function Is_In ( v : Standard_Natural_Vectors.Vector;
                   i : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if v(k) = i for some k, returns false otherwise.

  function Is_Connected ( deco : Link_to_VecVec; i,j : in natural32 ) 
                        return boolean;

  -- DESCRIPTION :
  --   Two points indexed by i and j are connected if they belong to
  --   the same vector in deco.

  procedure Connect ( deco : in Link_to_VecVec;
                      i,j : in natural32 );

  -- DESCRIPTION :
  --   Updates the connectivity information about i and j in deco.

  procedure Add_Map ( deco : in Link_to_VecVec; nb : in out natural32;
                      map : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Adds the information of the map to the decomposition.
  --   The number of factors "nb" is updated accordingly.

  function Number_of_Factors ( deco : in VecVec ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nonempty entries of deco.

  procedure Write_Factors ( file : in file_type; deco : in VecVec );
  procedure Write_Factors ( file : in file_type; deco : in VecVec;
                            m : in Standard_Natural_Vectors.vector );

  -- DESCRIPTION :
  --   Writes the irreducible factors in the decomposition to the file,
  --   eventually followed by the multiplity of each factor if m contains
  --   the multiplicity of each witness point.

  procedure Assign_Legend ( deco : in Link_to_VecVec;
                            legend : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Everywhere in deco, i will be replaced by legend(i).

  procedure Merge ( deco : in Link_to_VecVec;
                    subdeco : in Link_to_VecVec );

  -- DESCRIPTION :
  --   Copies the connectivity information in subdeco over to deco.

  procedure Remove_Empty_Entries ( deco : in out Link_to_VecVec );

  -- DESCRIPTION :
  --   Removes empty entries from the decomposition structure.

end Monodromy_Partitions;
