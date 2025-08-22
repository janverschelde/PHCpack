package body demics_ftest is

  package body class_theData is

    function new_theData return theData is

      res : theData;

    begin
      return res;
    end new_theData;

    procedure delete_theData ( this : access theData ) is
    begin
      null;
    end delete_theData;

    procedure create ( this : access theData;
                       ori_row : in integer32; ori_col : in integer32;
                       ori_termS : in integer32;
                       ori_polyDim : in integer32 ) is
    begin
      null;
    end create;

    procedure joint ( this : access theData ) is
    begin
      null;
    end joint;

    procedure iJoint ( this : access theData ) is
    begin
      null;
    end iJoint;

    procedure mJoint ( this : access theData ) is
    begin
      null;
    end mJoint;

    procedure clear ( this : access theData ) is
    begin
      null;
    end clear;

    procedure clear_transMat ( this : access theData ) is
    begin
      null;
    end clear_transMat;

    procedure put_info ( this : access theData; repIdx : in integer32;
                         idx2 : out integer32; lNbN : out integer32;
                         lNfN : out integer32 ) is
    begin
      null;
    end put_info;

    function invB_out ( this : access theData;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return 0.0;
    end invB_out;

    function transMat_out ( this : access theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float is
    begin
      return 0.0;
    end transMat_out;

    function invB_ptr_out ( this : access theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float is
    begin
      return 0.0;
    end invB_ptr_out;

    function transMat_ptr_out ( this : access theData;
                                rowIdx : integer32; colIdx : integer32 )
                              return double_float is
    begin
      return 0.0;
    end transMat_ptr_out;

    procedure info_p_sol ( this : access theData ) is
    begin
      null;
    end info_p_sol;

    procedure info_d_sol ( this : access theData ) is
    begin
      null;
    end info_d_sol;

    procedure info_invB ( this : access theData ) is
    begin
      null;
    end info_invB;

    procedure info_transMat ( this : access theData ) is
    begin
      null;
    end info_transMat;

    procedure info_transRed ( this : access theData ) is
    begin
      null;
    end info_transRed;

    procedure info_basisIdx ( this : access theData ) is
    begin
      null;
    end info_basisIdx;

    procedure info_nf_pos ( this : access theData ) is
    begin
      null;
    end info_nf_pos;

    procedure info_nbIdx ( this : access theData ) is
    begin
      null;
    end info_nbIdx;

    procedure info_redVec ( this : access theData ) is
    begin
      null;
    end info_redVec;

    procedure info_rIdx ( this : access theData ) is
    begin
      null;
    end info_rIdx;

    procedure info_pivOutIdx ( this : access theData ) is
    begin
      null;
    end info_pivOutIdx;

    procedure info_p_sol_ptr ( this : access theData ) is
    begin
      null;
    end info_p_sol_ptr;

    procedure info_d_sol_ptr ( this : access theData ) is
    begin
      null;
    end info_d_sol_ptr;

    procedure info_invB_ptr ( this : access theData ) is
    begin
      null;
    end info_invB_ptr;

    procedure info_transMat_ptr ( this : access theData ) is
    begin
      null;
    end info_transMat_ptr;

    procedure info_transRed_ptr ( this : access theData ) is
    begin
      null;
    end info_transRed_ptr;

    procedure info_basisIdx_ptr ( this : access theData ) is
    begin
      null;
    end info_basisIdx_ptr;

    procedure info_nf_pos_ptr ( this : access theData ) is
    begin
      null;
    end info_nf_pos_ptr;

    procedure info_nbIdx_ptr ( this : access theData ) is
    begin
      null;
    end info_nbIdx_ptr;

    procedure info_redVec_ptr ( this : access theData ) is
    begin
      null;
    end info_redVec_ptr;

    procedure info_fIdx ( this : access theData ) is
    begin
      null;
    end info_fIdx;

    procedure info_node ( this : access theData ) is
    begin
      null;
    end info_node;

  end class_theData;

  package body Class_ftData is

    function new_ftData return ftData is

      res : ftData;

    begin
      return res;
    end new_ftData;

    procedure delete_ftData ( this : access ftData ) is
    begin
      null;
    end delete_ftData;

    procedure create_elem
      ( this : access ftData;
        row : in integer32; col : in integer32;
        termS : in integer32; polyDim : in integer32 ) is

    begin
      null;
    end create_elem;

    procedure add_elem ( this : access ftData ) is
    begin
      null;
    end add_elem;

    procedure mark ( this : access ftData ) is
    begin
      null;
    end mark;

    procedure clear ( this : access ftData ) is
    begin
      null;
    end clear;

    procedure clear_transMat ( this : access ftData ) is
    begin
      null;
    end clear_transMat;

    procedure delete_cur ( this : access ftData ) is
    begin
      null;
    end delete_cur;

    procedure delete_all ( this : access ftData ) is
    begin
      null;
    end delete_all;

    procedure delete_addedElem ( this : access ftData ) is
    begin
      null;
    end delete_addedElem;

    procedure init_ptr ( this : access ftData ) is
    begin
      null;
    end init_ptr;

    procedure make_init_data
      ( this : access ftData;
        termSumNum : in integer32; supN : in integer32;
        termS : in integer32; reTermS : in integer32 ) is

    begin
      null;
    end make_init_data;

    procedure next_data ( this : access ftData ) is
    begin
      null;
    end next_data;

    procedure copy ( this : access ftData;
                     col : in integer32; pre_data : access theData ) is
    begin
      null;
    end copy;

    procedure get_ptr ( this : access ftData;
                        pre_data : access theData ) is
    begin
      null;
    end get_ptr;

    procedure create_rIdx
      ( this : access ftData; nbN : in integer32; repIdx : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector ) is

    begin
      null;
    end create_rIdx;

    procedure init_info ( this : access ftData ) is
    begin
      null;
    end init_info;

    procedure get_nbIdx_rIdx
      ( this : access ftData;
        preNbN : in integer32; repIdx : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector;
        reTermS : in integer32; pre_data : access theData ) is

    begin
      null;
    end get_nbIdx_rIdx;

    procedure iCopy
      ( this : access ftData;
        nbN : in integer32; nfN : in integer32;
        repIdx : in integer32; termS : in integer32;
        reTermS : in integer32;
        candIdx : in Standard_Integer_Vectors.Link_to_Vector;
        pre_data : access theData ) is

    begin
      null;
    end iCopy;

    procedure iGetPtr ( this : access ftData; pre_data : access theData ) is
    begin
      null;
    end iGetPtr;

    procedure output
      ( this : access ftData;
        repIdx : in integer32;
        idx2 : out integer32; nbN : out integer32; nfN : out integer32 ) is

    begin
      null;
    end output;

    procedure decrease_nfN ( this : access ftData ) is
    begin
      null;
    end decrease_nfN;

    procedure copy_rIdx
      ( this : access ftData; pre_data : access theData;
        termS : in integer32 ) is

    begin
      null;
    end copy_rIdx;

    procedure copy_pivOutIdx ( this : access ftData;
                               pre_data : access theData ) is
    begin
      null;
    end copy_pivOutIdx;

    procedure get_nf_pos
      ( this : access ftData; pre_data : access theData;
        nfN : in integer32; idx2 : in integer32 ) is

    begin
      null;
    end get_nf_pos;

    procedure mCopy ( this : access ftData;
                      nbN : in integer32; nfN : in integer32;
                      idx2 : in integer32; termS : in integer32;
                      pre_data : access theData ) is
    begin
      null;
    end mCopy;

    procedure mGetPtr ( this : access ftData;
                        pre_data : access theData ) is
    begin
      null;
    end mGetPtr;

    procedure put_sup ( this : access ftData; sup : out integer32 ) is
    begin
      null;
    end put_sup;

    procedure info_parent_nbN_nfN ( this : access ftData ) is
    begin
      null;
    end info_parent_nbN_nfN;

    procedure info_parent_p_sol ( this : access ftData ) is
    begin
      null;
    end info_parent_p_sol;

    procedure info_parent_d_sol ( this : access ftData ) is
    begin
      null;
    end info_parent_d_sol;

    procedure info_parent_invB ( this : access ftData ) is
    begin
      null;
    end info_parent_invB;

    procedure info_parent_transMat ( this : access ftData ) is
    begin
      null;
    end info_parent_transMat;

    procedure info_parent_transRed ( this : access ftData ) is
    begin
      null;
    end info_parent_transRed;

    procedure info_parent_basisIdx ( this : access ftData ) is
    begin
      null;
    end info_parent_basisIdx;

    procedure info_parent_nf_pos ( this : access ftData ) is
    begin
      null;
    end info_parent_nf_pos;

    procedure info_parent_nbIdx ( this : access ftData ) is
    begin
      null;
    end info_parent_nbIdx;

    procedure info_parent_redVec ( this : access ftData ) is
    begin
      null;
    end info_parent_redVec;

    procedure info_parent_rIdx ( this : access ftData ) is
    begin
      null;
    end info_parent_rIdx;

    procedure info_parent_pivOutIdx ( this : access ftData ) is
    begin
      null;
    end info_parent_pivOutIdx;

    procedure info_parent_p_sol_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_p_sol_ptr;

    procedure info_parent_d_sol_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_d_sol_ptr;

    procedure info_parent_invB_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_invB_ptr;

    procedure info_parent_transMat_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_transMat_ptr;

    procedure info_parent_transRed_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_transRed_ptr;

    procedure info_parent_basisIdx_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_basisIdx_ptr;

    procedure info_parent_nf_pos_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_nf_pos_ptr;

    procedure info_parent_nbIdx_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_nbIdx_ptr;

    procedure info_parent_redVec_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_redVec_ptr;

    procedure info_parent_pivOutIdx_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_pivOutIdx_ptr;

    procedure info_parent ( this : access ftData ) is
    begin
      null;
    end info_parent;

    procedure info_parent_ptr ( this : access ftData ) is
    begin
      null;
    end info_parent_ptr;

    procedure info_parent_node ( this : access ftData ) is
    begin
      null;
    end info_parent_node;

    procedure info_cur_nbN_nfN ( this : access ftData ) is
    begin
      null;
    end info_cur_nbN_nfN;

    procedure info_cur_p_sol ( this : access ftData ) is
    begin
      null;
    end info_cur_p_sol;

    procedure info_cur_d_sol ( this : access ftData ) is
    begin
      null;
    end info_cur_d_sol;

    procedure info_cur_invB ( this : access ftData ) is
    begin
      null;
    end info_cur_invB;

    procedure info_cur_transMat ( this : access ftData ) is
    begin
      null;
    end info_cur_transMat;

    procedure info_cur_transRed ( this : access ftData ) is
    begin
      null;
    end info_cur_transRed;

    procedure info_cur_basisIdx ( this : access ftData ) is
    begin
      null;
    end info_cur_basisIdx;

    procedure info_cur_nf_pos ( this : access ftData ) is
    begin
      null;
    end info_cur_nf_pos;

    procedure info_cur_nbIdx ( this : access ftData ) is
    begin
      null;
    end info_cur_nbIdx;

    procedure info_cur_redVec ( this : access ftData ) is
    begin
      null;
    end info_cur_redVec;

    procedure info_cur_rIdx ( this : access ftData ) is
    begin
      null;
    end info_cur_rIdx;

    procedure info_cur_pivOutIdx ( this : access ftData ) is
    begin
      null;
    end info_cur_pivOutIdx;

    procedure info_cur_p_sol_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_p_sol_ptr;

    procedure info_cur_d_sol_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_d_sol_ptr;

    procedure info_cur_invB_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_invB_ptr;

    procedure info_cur_transMat_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_transMat_ptr;

    procedure info_cur_transRed_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_transRed_ptr;

    procedure info_cur_basisIdx_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_basisIdx_ptr;

    procedure info_cur_nf_pos_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_nf_pos_ptr;

    procedure info_cur_nbIdx_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_nbIdx_ptr;

    procedure info_cur_redVec_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_redVec_ptr;

    procedure info_cur ( this : access ftData ) is
    begin
      null;
    end info_cur;

    procedure info_cur_ptr ( this : access ftData ) is
    begin
      null;
    end info_cur_ptr;

    procedure info_cur_node ( this : access ftData ) is
    begin
      null;
    end info_cur_node;

    procedure info_all_node ( this : access ftData ) is
    begin
      null;
    end info_all_node;

    procedure info_all_cur ( this : access ftData ) is
    begin
      null;
    end info_all_cur;

    procedure info_all_nodeNum ( this : access ftData ) is
    begin
      null;
    end info_all_nodeNum;

    procedure info_numElem ( this : access ftData ) is
    begin
      null;
    end info_numElem;

  end class_ftData;

  package body class_lvData is

    function new_lvData return lvData is

      res : lvData;

    begin
      return res;
    end new_lvData;

    procedure delete_lvData ( this : access lvData ) is
    begin
      null;
    end delete_lvData;

    procedure create
      ( this : access lvData;
        depth : in integer32; supN : in integer32; dim : in integer32;
        ori_length : in integer32; ori_termMax : in integer32 ) is

    begin
      null;
    end create;

    procedure get_info
      ( this : access lvData;
        g_mRepN : out Standard_Integer_VecVecs.Link_to_VecVec;
        g_mFeaIdx : out Standard_Integer_VecVecs.Link_to_VecVec;
        g_mFea : out Standard_Integer_VecVecs.Link_to_VecVec ) is

    begin
      null;
    end get_info;

    procedure init_ptr ( this : access lvData ) is
    begin
      null;
    end init_ptr;

    procedure info_mFea ( this : access lvData ) is
    begin
      null;
    end info_mFea;

  end class_lvData;

end demics_ftest;
