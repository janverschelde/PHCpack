with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_VecVecs;

package Write_Witness_Solutions is

  procedure Standard_Write ( topdim,lowdim : in natural32 );

  -- DESCRIPTION :
  --   Writes the number of solutions at each dimension,
  --   computed in standard double precision,
  --   starting at the top dimension topdim.

  procedure DoblDobl_Write ( topdim,lowdim : in natural32 );

  -- DESCRIPTION :
  --   Writes the number of solutions at each dimension,
  --   computed in double double precision,
  --   starting at the top dimension topdim.

  procedure QuadDobl_Write ( topdim,lowdim : in natural32 );

  -- DESCRIPTION :
  --   Writes the number of solutions at each dimension,
  --   computed in quad double precision,
  --   starting at the top dimension topdim.

  procedure Write_Counts
              ( filter,factor : in boolean;
                pc,fc : in Standard_Natural_VecVecs.Link_to_VecVec;
                idx : in Standard_Natural_VecVecs.Link_to_Array_of_VecVecs );

  -- DESCRIPTION :
  --   Writes the path counts pc, and if filter, also the filter counts fc,
  --   and if factor, also the indices idx to the irreducible factors.

end Write_Witness_Solutions;
