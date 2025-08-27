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

    type theData is record
      row : integer32;
      col : integer32;
      termS : integer32;
      next : access theData;
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

    type Link_to_theData is access theData;
    type Array_of_theData is array ( integer range <> ) of Link_to_theData;
    type Link_to_Array_of_theData is access Array_of_theData;

    function new_theData return theData;

    procedure Delete_theData ( this : access theData );

    procedure create ( this : access theData;
                       ori_row : in integer32; ori_col : in integer32;
                       ori_termS : in integer32; ori_polyDim : in integer32 );

    procedure joint ( this : access theData );

    procedure iJoint ( this : access theData );

    procedure mJoint ( this : access theData );

    procedure clear ( this : access theData );

    procedure clear_transMat ( this : access theData );

    procedure put_info ( this : access theData; repIdx : in integer32;
                         idx2 : out integer32; lNbN : out integer32;
                         lNfN : out integer32 );

    function invB_out ( this : access theData;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float;

    function transMat_out ( this : access theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float;

    function invB_ptr_out ( this : access theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float;

    function transMat_ptr_out ( this : access theData;
                                rowIdx : integer32; colIdx : integer32 )
                              return double_float;

    procedure info_p_sol ( this : access theData );

    procedure info_d_sol ( this : access theData );

    procedure info_invB ( this : access theData );

    procedure info_transMat ( this : access theData );

    procedure info_transRed ( this : access theData );

    procedure info_basisIdx ( this : access theData );

    procedure info_nf_pos ( this : access theData );

    procedure info_nbIdx ( this : access theData );

    procedure info_redVec ( this : access theData );

    procedure info_rIdx ( this : access theData );

    procedure info_pivOutIdx ( this : access theData );

    procedure info_p_sol_ptr ( this : access theData );

    procedure info_d_sol_ptr ( this : access theData );

    procedure info_invB_ptr ( this : access theData );

    procedure info_transMat_ptr ( this : access theData );

    procedure info_transRed_ptr ( this : access theData );

    procedure info_basisIdx_ptr ( this : access theData );

    procedure info_nf_pos_ptr ( this : access theData );

    procedure info_nbIdx_ptr ( this : access theData );

    procedure info_redVec_ptr ( this : access theData );

    procedure info_fIdx ( this : access theData );

    procedure info_node ( this : access theData );

  end class_theData;

  use class_theData;

  package Class_ftData is

    type ftData is record
      dim : integer32;
      elemNum : integer32;
      cur : access theData;
      parent : access theData;
      limit : access theData;
      head : access theData;
      last : access theData;
    end record;

    type Link_to_ftData is access ftData;
    type Array_of_ftData is array ( integer32 range <> ) of Link_to_ftData;
    type Link_to_Array_of_ftdata is access Array_of_ftData;

    function new_ftData return ftData;

    procedure delete_ftData ( this : access ftData );

    procedure create_elem
      ( this : access ftData;
        row : in integer32; col : in integer32;
        termS : in integer32; polyDim : in integer32 );

    procedure add_elem ( this : access ftData );

    procedure mark ( this : access ftData );

    procedure clear ( this : access ftData );

    procedure clear_transMat ( this : access ftData );

    procedure delete_cur ( this : access ftData );

    procedure delete_all ( this : access ftData );

    procedure delete_addedElem ( this : access ftData );

    procedure init_ptr ( this : access ftData );

    procedure make_init_data
      ( this : access ftData;
        termSumNum : in integer32; supN : in integer32;
        termS : in integer32; reTermS : in integer32 );

    procedure next_data ( this : access ftData );

    procedure copy ( this : access ftData;
                     col : in integer32; pre_data : access theData );

    procedure get_ptr ( this : access ftData;
                        pre_data : access theData );

    procedure create_rIdx
      ( this : access ftData; nbN : in integer32; repIdx : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector );

    procedure init_info ( this : access ftData );

    procedure get_nbIdx_rIdx
      ( this : access ftData;
        preNbN : in integer32; repIdx : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector;
        reTermS : in integer32; pre_data : access theData );

    procedure iCopy
      ( this : access ftData;
        nbN : in integer32; nfN : in integer32;
        repIdx : in integer32; termS : in integer32;
        reTermS : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector;
        pre_data : access theData );

    procedure iGetPtr ( this : access ftData;
                        pre_data : access theData );

    procedure output
      ( this : access ftData;
        repIdx : in integer32;
        idx2 : out integer32; nbN : out integer32; nfN : out integer32 );

    procedure decrease_nfN ( this : access ftData );

    procedure copy_rIdx
      ( this : access ftData;
        pre_data : access theData; termS : in integer32 );

    procedure copy_pivOutIdx ( this : access ftData;
                               pre_data : access theData );

    procedure get_nf_pos
      ( this : access ftData; pre_data : access theData;
        nfN : in integer32; idx2 : in integer32 );

    procedure mCopy ( this : access ftData;
                      nbN : in integer32; nfN : in integer32;
                      idx2 : in integer32; termS : in integer32;
                      pre_data : access theData );

    procedure mGetPtr ( this : access ftData;
                        pre_data : access theData );

    procedure put_sup ( this : access ftData; sup : out integer32 );

    procedure info_parent_nbN_nfN ( this : access ftData );

    procedure info_parent_p_sol ( this : access ftData );

    procedure info_parent_d_sol ( this : access ftData );

    procedure info_parent_invB ( this : access ftData );

    procedure info_parent_transMat ( this : access ftData );

    procedure info_parent_transRed ( this : access ftData );

    procedure info_parent_basisIdx ( this : access ftData );

    procedure info_parent_nf_pos ( this : access ftData );

    procedure info_parent_nbIdx ( this : access ftData );

    procedure info_parent_redVec ( this : access ftData );

    procedure info_parent_rIdx ( this : access ftData );

    procedure info_parent_pivOutIdx ( this : access ftData );

    procedure info_parent_p_sol_ptr ( this : access ftData );

    procedure info_parent_d_sol_ptr ( this : access ftData );

    procedure info_parent_invB_ptr ( this : access ftData );

    procedure info_parent_transMat_ptr ( this : access ftData );

    procedure info_parent_transRed_ptr ( this : access ftData );

    procedure info_parent_basisIdx_ptr ( this : access ftData );

    procedure info_parent_nf_pos_ptr ( this : access ftData );

    procedure info_parent_nbIdx_ptr ( this : access ftData );

    procedure info_parent_redVec_ptr ( this : access ftData );

    procedure info_parent_pivOutIdx_ptr ( this : access ftData );

    procedure info_parent ( this : access ftData );

    procedure info_parent_ptr ( this : access ftData );

    procedure info_parent_node ( this : access ftData );

    procedure info_cur_nbN_nfN ( this : access ftData );

    procedure info_cur_p_sol ( this : access ftData );

    procedure info_cur_d_sol ( this : access ftData );

    procedure info_cur_invB ( this : access ftData );

    procedure info_cur_transMat ( this : access ftData );

    procedure info_cur_transRed ( this : access ftData );

    procedure info_cur_basisIdx ( this : access ftData );

    procedure info_cur_nf_pos ( this : access ftData );

    procedure info_cur_nbIdx ( this : access ftData );

    procedure info_cur_redVec ( this : access ftData );

    procedure info_cur_rIdx ( this : access ftData );

    procedure info_cur_pivOutIdx ( this : access ftData );

    procedure info_cur_p_sol_ptr ( this : access ftData );

    procedure info_cur_d_sol_ptr ( this : access ftData );

    procedure info_cur_invB_ptr ( this : access ftData );

    procedure info_cur_transMat_ptr ( this : access ftData );

    procedure info_cur_transRed_ptr ( this : access ftData );

    procedure info_cur_basisIdx_ptr ( this : access ftData );

    procedure info_cur_nf_pos_ptr ( this : access ftData );

    procedure info_cur_nbIdx_ptr ( this : access ftData );

    procedure info_cur_redVec_ptr ( this : access ftData );

    procedure info_cur ( this : access ftData );

    procedure info_cur_ptr ( this : access ftData );

    procedure info_cur_node ( this : access ftData );

    procedure info_all_node ( this : access ftData );

    procedure info_all_cur ( this : access ftData );

    procedure info_all_nodeNum ( this : access ftData );

    procedure info_numElem ( this : access ftData );

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
      fTest : access ftData;
      Node : access ftData;
    end record;

    function new_lvData return lvData;

    procedure delete_lvData ( this : access lvData );

    procedure create
      ( this : access lvData;
        depth : in integer32; supN : in integer32; dim : in integer32;
        ori_length : in integer32; ori_termMax : in integer32 );

    procedure get_info
      ( this : access lvData;
        g_mRepN : out Standard_Integer_VecVecs.Link_to_VecVec;
        g_mFeaIdx : out Standard_Integer_VecVecs.Link_to_VecVec;
        g_mFea : out Standard_Integer_VecVecs.Link_to_VecVec );

    procedure init_ptr ( this : access lvData );

    procedure info_mFea ( this : access lvData );

  end class_lvData;

end demics_ftest;
