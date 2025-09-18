with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;

package demics_ftest is

-- DESCRIPTION :
--   Defines classes to store floating-point data of linear programs.
--   Translation initiated by g++ -c -fdump-ada-spec fTest.h.

  package class_theData is

    type theData;
    type Link_to_theData is access theData;
    type Array_of_theData is array ( integer32 range <> ) of Link_to_theData;
    type Link_to_Array_of_theData is access Array_of_theData;

    type theData is record
      row : integer32;
      col : integer32;
      termS : integer32;
      next : Link_to_theData;
      flag : integer32;
      polyDim : integer32;
      nbN : integer32;
      nfN : integer32;
      artV : integer32;
      pivOutNum : integer32;
      fIdx : integer32;
      sw : integer32;
      invB : Standard_Floating_Vectors.Link_to_Vector;
      transMat : Standard_Floating_Vectors.Link_to_Vector;
      transRed : Standard_Floating_Vectors.Link_to_Vector;
      p_sol : Standard_Floating_Vectors.Link_to_Vector;
      d_sol : Standard_Floating_Vectors.Link_to_Vector;
      redVec : Standard_Floating_Vectors.Link_to_Vector;
      basisIdx : Standard_Integer_Vectors.Link_to_Vector;
      nbIdx : Standard_Integer_Vectors.Link_to_Vector;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      rIdx : Standard_Integer_Vectors.Link_to_Vector;
      pivOutList : Standard_Integer_Vectors.Link_to_Vector;
      pivOutCheck : Standard_Integer_Vectors.Link_to_Vector;
      invB_ptr : Standard_Floating_Vectors.Link_to_Vector;
      transMat_ptr : Standard_Floating_Vectors.Link_to_Vector;
      transRed_ptr : Standard_Floating_Vectors.Link_to_Vector;
      p_sol_ptr : Standard_Floating_Vectors.Link_to_Vector;
      d_sol_ptr : Standard_Floating_Vectors.Link_to_Vector;
      redVec_ptr : Standard_Floating_Vectors.Link_to_Vector;
      basisIdx_ptr : Standard_Integer_Vectors.Link_to_Vector;
      nbIdx_ptr : Standard_Integer_Vectors.Link_to_Vector;
      nf_pos_ptr : Standard_Integer_Vectors.Link_to_Vector;
      nodeLabel : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    function new_theData return theData;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_theData ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Deallocates those fields in this that are not ending in _ptr.

    procedure create ( this : in Link_to_theData;
                       ori_row : in integer32; ori_col : in integer32;
                       ori_termS : in integer32; ori_polyDim : in integer32;
                       vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates space for the fields in the this record.

    procedure joint ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Assign all this.*_ptr values to the corresponding this.* values.

    procedure iJoint ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Assigns _ptr values for transMat, transRed, redVec, and nbIdx.

    procedure mJoint ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Assigns this.nf_pos to this.nf_pos_ptr.

    procedure clear ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Sets values in this to zero.

    procedure clear_transMat ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Sets all values of the this.transMat to zero.

    procedure put_info ( this : in Link_to_theData;
                         repIdx : in integer32; idx2 : out integer32;
                         lNbN : out integer32; lNfN : out integer32 );

    -- DESCRIPTION :
    --   Sets idx2 to this.rIdx(repIdx), sets lNbN to this.nbN, and
    --   sets lNfN to this.nfN.

    function invB_out ( this : Link_to_theData;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float;

    -- DESCRIPTION :
    --   Returns the value of this.invB as defined by the rowIdx
    --   and the colIdx.

    function transMat_out ( this : Link_to_theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float;

    -- DESCRIPTION :
    --   Returns the value of this.transMat as defined by the rowIdx
    --   and the colIdx.

    function invB_ptr_out ( this : Link_to_theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float;

    -- DESCRIPTION :
    --   Returns the value of this.invB_ptr as defined by the rowIdx
    --   and the colIdx.

    function transMat_ptr_out ( this : Link_to_theData;
                                rowIdx : integer32; colIdx : integer32 )
                              return double_float;

    -- DESCRIPTION :
    --   Returns the value of this.transMat_ptr as defined by the rowIdx
    --   and the colIdx.

    procedure info_p_sol ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.p_sol.

    procedure info_d_sol ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.d_sol.

    procedure info_invB ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.invB.

    procedure info_transMat ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.transMat.

    procedure info_transRed ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.transRed.

    procedure info_basisIdx ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.basisIdx.

    procedure info_nf_pos ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.nf_pos.

    procedure info_nbIdx ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.nbIdx.

    procedure info_redVec ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.redVec.

    procedure info_rIdx ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Write the numbers stored in this.rIdx.

    procedure info_pivOutIdx ( this : in Link_to_theData );

    -- DESCRPIPTION :
    --   Writes the numbers in this.pivOutCheck and this.pivOutList.

    procedure info_p_sol_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.p_sol_ptr.

    procedure info_d_sol_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.d_sol_ptr.

    procedure info_invB_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.invB_ptr.

    procedure info_transMat_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.transMat_ptr.

    procedure info_transRed_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers in this.transRed_ptr.

    procedure info_basisIdx_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.info_basisIdx_ptr.

    procedure info_nf_pos_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.info_nf_pos_ptr.

    procedure info_nbIdx_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in nbIdx_ptr.

    procedure info_redVec_ptr ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.redVec_ptr.

    procedure info_fIdx ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the value of this.fIdx + 1, followed by a new line.

    procedure info_node ( this : in Link_to_theData );

    -- DESCRIPTION :
    --   Writes the values in this.nodeLabel, each augmented by one.

  end class_theData;

  use class_theData;

  package class_ftData is

    type ftData is record
      dim : integer32;
      elemNum : integer32;
      cur : Link_to_theData;
      parent : Link_to_theData;
      limit : Link_to_theData;
      head : Link_to_theData;
      last : Link_to_theData;
    end record;

    type Link_to_ftData is access ftData;
    type Array_of_ftData is array ( integer32 range <> ) of Link_to_ftData;
    type Link_to_Array_of_ftData is access Array_of_ftData;

    function new_ftData return ftData;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_ftData ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Empty destructor.

    procedure clear ( lftd : in out Link_to_Array_of_ftData );

    -- DESCRIPTION :
    --   Deallocates the pointer to the array.

    procedure create_elem
                ( this : in Link_to_ftData;
                  row : in integer32; col : in integer32;
                  termS : in integer32; polyDim : in integer32;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Makes a new data element in this.cur.

    procedure add_elem ( this : in Link_to_ftData;
                         vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Appends this.cur to the list this.head.

    procedure mark ( this : in Link_to_ftData;
                     vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Sets this.limit to this.last.

    procedure clear ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls the clear procedure to all elements of this.head.

    procedure clear_transMat ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls the clear_transMat procedure to all elements of this.head.

    procedure delete_cur ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls delete_theData on this.cur.

    procedure delete_all ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Deletes all elements of this.head and sets then all
    --   fields to zero or to null.

    procedure delete_addedElem ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Deletes all elements of this.limit.

    procedure init_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Sets the pointers this.parent and this.cur to this.head.

    procedure make_init_data
                ( this : in Link_to_ftData;
                  termSumNum : in integer32; supN : in integer32;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Initializes the indices in this.cur.

    procedure next_data ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Sets the this.parent to this.parent.next if not null.

    procedure copy ( this : in Link_to_ftData;
                     col : in integer32; pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Copies values of p_sol, d_sol, basisIdx, and nbIdx
    --   from corresponding *_ptr fields from pre_data into this.cur.

    procedure get_ptr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Sets the values for p_sol, d_sol, basisIdx, nbIdx in this.cur
    --   from corresponding fields in pre_data, after deallocating
    --   the fields in this.cur.

    procedure create_rIdx
                 ( this : in Link_to_ftData;
                   preNbN : in integer32; repIdx : in integer32;
                   candIdx : in Standard_Integer_Vectors.Link_to_Vector );

    -- DESCRIPTION :
    --   Sets the values of this.cur.rIdx.

    procedure init_info ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Sets this.pivOutcheck and this.transRef to zero.

    procedure get_nbIdx_rIdx
                ( this : in Link_to_ftData;
                  preNbN : in integer32; repIdx : in integer32;
                  candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  reTermS : in integer32; pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Gets the values of this.cur.nbIdx and this.cur.rIdx
    --   from the input parameters.

    procedure iCopy ( this : in Link_to_ftData;
                      preNbN : in integer32; nfN : in integer32;
                      repIdx : in integer32; -- termS : in integer32;
                      reTermS : in integer32;
                      candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                      pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Copies values into nbIdx, rIdx, invB of this.cur
    --   from the input arguments.
    --   The original input parameter termS is not referenced.

    procedure iGetPtr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Copies values of _ptr fields from pre_data into this.cur.

    procedure output ( this : in Link_to_ftData;
                       repIdx : in integer32; idx2 : out integer32;
                       nbN : out integer32; nfN : out integer32 );

    -- DESCRIPTION :
    --   Assigns to the output variables values from this.parent.

    procedure decrease_nfN ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Decrements the value of nfN of this.cur by one.

    procedure copy_rIdx
                ( this : in Link_to_ftData;
                  pre_data : in Link_to_theData; termS : in integer32 );

    -- DESCRIPTION :
    --   Copies termS values for pre_data.rIdx into this.cur.

    procedure copy_pivOutIdx ( this : in Link_to_ftData;
                               pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Copies indices from pre_data into this.cur.pivOutCheck
    --   and this.cur.pivoutList.

    procedure get_nf_pos
                ( this : in Link_to_ftData; pre_data : in Link_to_theData;
                  nfN : in integer32; idx2 : in integer32 );

    -- DESCRIPTION :
    --   Gets values from pre_data.nf_pos_ptr into this.cur.nf_pos.

    procedure mCopy ( this : in Link_to_ftData;
                      nbN : in integer32; nfN : in integer32;
                      idx2 : in integer32; termS : in integer32;
                      pre_data : in Link_to_theData );

    procedure mGetPtr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData );

    -- DESCRIPTION :
    --   Copies into this.cur the values p_sol, d_sol, basisIdx, nbIdx
    --   of pre_data after clearing the fields in this.cur.

    procedure put_sup ( this : in Link_to_ftData; sup : out integer32 );

    -- DESCRIPTION :
    --   Returns in sup the first nodeLabel of this.parent.

    procedure info_parent_nbN_nfN ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Writes the values of nbN ad nfN of this.parent.

    procedure info_parent_p_sol ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_p_sol on this.parent.

    procedure info_parent_d_sol ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_d_sol on this.parent.

    procedure info_parent_invB ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_invB on this.parent.

    procedure info_parent_transMat ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transMat on this.parent.

    procedure info_parent_transRed ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transRed on this.parent.

    procedure info_parent_basisIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_basisIdx on this.parent.

    procedure info_parent_nf_pos ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nf_pos on this.parent.

    procedure info_parent_nbIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nbIdx on this.parent.

    procedure info_parent_redVec ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_redVec on this.parent.

    procedure info_parent_rIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_rIdx on this.parent.

    procedure info_parent_pivOutIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_pivOutIdx on this.parent.

    procedure info_parent_p_sol_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_p_sol_ptr on this.parent.

    procedure info_parent_d_sol_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_d_sol_ptr on this.parent.

    procedure info_parent_invB_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_invB_ptr on this.parent.

    procedure info_parent_transMat_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transMat_ptr on this.parent.

    procedure info_parent_transRed_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transRed_ptr on this.parent.

    procedure info_parent_basisIdx_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_basisIdx_ptr on this.parent.

    procedure info_parent_nf_pos_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nf_pos_ptr on this.parent.

    procedure info_parent_nbIdx_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nbIdx_ptr on this.parent.

    procedure info_parent_redVec_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_redVec_ptr on this.parent.

    procedure info_parent ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls all info_* procedures on the this.parent fields.

    procedure info_parent_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls all info_*_ptr procedures on the this.parent fields.

    procedure info_parent_node ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_node on this.parent.

    procedure info_cur_nbN_nfN ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Writes the values of nbN ad nfN of this.cur.

    procedure info_cur_p_sol ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_p_sol on this.cur.

    procedure info_cur_d_sol ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_d_sol on this.cur.

    procedure info_cur_invB ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_invB on this.cur.

    procedure info_cur_transMat ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transMat on this.cur.

    procedure info_cur_transRed ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transRed on this.cur.

    procedure info_cur_basisIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_basisIdx on this.cur.

    procedure info_cur_nf_pos ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nf_pos on this.cur.

    procedure info_cur_nbIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nbIdx on this.cur.

    procedure info_cur_redVec ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_redVec on this.cur.

    procedure info_cur_rIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_rIdx on this.cur.

    procedure info_cur_pivOutIdx ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_pivOutIdx on this.cur.

    procedure info_cur_p_sol_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_p_sol_ptr on this.cur.

    procedure info_cur_d_sol_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_d_sol_ptr on this.cur.

    procedure info_cur_invB_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_invB_ptr on this.cur.

    procedure info_cur_transMat_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transMat_ptr on this.cur.

    procedure info_cur_transRed_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_transRed_ptr on this.cur.

    procedure info_cur_basisIdx_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_basisIdx_ptr on this.cur.

    procedure info_cur_nf_pos_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nf_pos_ptr on this.cur.

    procedure info_cur_nbIdx_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_nbIdx_ptr on this.cur.

    procedure info_cur_redVec_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_redVec_ptr on this.cur.

    procedure info_cur ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls all info_* procedures on the this.cur fields.

    procedure info_cur_ptr ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls all info_*_ptr procedures on the this.cur fields.

    procedure info_cur_node ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_node on this.cur.

    procedure info_all_node ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls info_node on all elements in the list this.head.

    procedure info_all_cur ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Calls all info_ procedures on the all elements in this.head.

    procedure info_all_nodeNum ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Writes the index counter of all elemens in the list this.head.

    procedure info_numElem ( this : in Link_to_ftData );

    -- DESCRIPTION :
    --   Counts and writes the number of elements in the list this.head.

  end class_ftData;

  use class_ftData;

  package class_lvData is

    type lvData is record
      dim : integer32;
      length : integer32;
      termMax : integer32;
      mRepN : Standard_Integer_Vectors.Link_to_Vector;
      mFeaIdx : Standard_Integer_VecVecs.Link_to_VecVec;
      mFea : Standard_Integer_Vectors.Link_to_Vector;
      fTest : Link_to_Array_of_ftData;
      Node : Link_to_ftData;
    end record;

    type Link_to_lvData is access lvData;
    type Array_of_lvData is array ( integer32 range <> ) of Link_to_lvData;
    type Link_to_Array_of_lvdata is access Array_of_lvData;

    function new_lvData return lvData;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_lvData ( this : in Link_to_lvData );

    -- DESCRIPTION :
    --   Dealocates the space occupied by lvData.

    procedure clear ( lvd : in out Link_to_Array_of_lvData );

    -- DESCRIPTION :
    --   Deallocates the pointer to the array.

    procedure create ( this : in Link_to_lvData; depth : in integer32;
                      -- supN : in integer32; dim : in integer32;
                       ori_length : in integer32; ori_termMax : in integer32;
                       vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates data for the lvData object.
    --   The original input arguments supN and dim are not references.

    procedure get_info
                ( this : in Link_to_lvData;
                  g_mRepN : out Standard_Integer_Vectors.Link_to_Vector;
                  g_mFeaIdx : out Standard_Integer_VecVecs.Link_to_VecVec;
                  g_mFea : out Standard_Integer_Vectors.Link_to_Vector );

    -- DESCRIPTION :
    --   Assigns to the three output parameters this.mRepN, this.mFeaIdx,
    --   and this.mFea respectively.

    procedure init_ptr ( this : in Link_to_lvData );

    -- DESCRIPTION :
    --   Calls init_ptr on each element on this.fTest.

    procedure info_mFea ( this : in Link_to_lvData );

    -- DESCRIPTION :
    --   Writes the values stored in this.mFea and this.mRepN.

  end class_lvData;

end demics_ftest;
