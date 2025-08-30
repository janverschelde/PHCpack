package body demics_simplex is

  package body class_supportSet is

    function new_supportSet return supportSet is

      res : supportSet;

    begin
      return res;
    end new_supportSet;

    procedure delete_supportSet ( this : in Link_to_supportSet ) is
    begin
      null;
    end delete_supportSet;

    procedure allocSupp
                ( this : in Link_to_supportSet;
                  data : demics_input_data.class_dataSet.dataSet;
                  level : in integer32;
                  num : in integer32;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      null;
    end allocSupp;

    procedure allocAux
                ( this : in Link_to_supportSet;
                  data : in demics_input_data.class_dataSet.dataSet ) is
    begin
      null;
    end allocAux;

    procedure supMat_in ( this : in Link_to_supportSet;
                          rowIdx : in integer32;
                          colIdx : in integer32;
                          elem : in double_float ) is
    begin
      null;
    end supMat_in;

    procedure supMat_neg ( this : in Link_to_supportSet;
                           rowIdx : in integer32;
                           colIdx : in integer32 ) is
    begin
      null;
    end supMat_neg;

    function supMat_out ( this : Link_to_supportSet;
                          rowIdx : integer32;
                          colIdx : integer32 ) return double_float is
    begin
      return 0.0;
    end supMat_out;

    function redVal ( this : Link_to_supportSet;
                      d_sol : Standard_Floating_Vectors.Link_to_Vector;
                      idx : integer32;
                      ii : integer32 ) return double_float is
    begin
      return 0.0;
    end redVal;

    procedure info_sup ( this : in Link_to_supportSet ) is
    begin
      null;
    end info_sup;

    procedure info_costVec ( this : in Link_to_supportSet ) is
    begin
      null;
    end info_costVec;

  end class_supportSet;

  package body class_simplex is

-- relation table

    function checkFrIdx ( this : Link_to_simplex ) return integer32 is
    begin
      return 0;
    end checkFrIdx;

    procedure elimFrIdx ( this : in Link_to_simplex;
                          sub_pivOutIdx : in integer32 ) is
    begin
      null;
    end elimFrIdx;

-- phase 1

    procedure reMakeNonBasisIdx ( this : in Link_to_simplex;
                                  reTermS : in integer32 ) is
    begin
      null;
    end reMakeNonBasisIdx;

    procedure reMakeNonBasisIdx_tab ( this : in Link_to_simplex ) is
    begin
      null;
    end reMakeNonBasisIdx_tab;

    procedure elimArt ( this : in Link_to_simplex;
                        depth : in integer32; preNbN : in integer32;
                        termS : in integer32; reTermS : in integer32;
                        iter : in out integer32 ) is
    begin
      null;
    end elimArt;

    procedure calRedCost ( this : in Link_to_simplex;
                           pivInIdx : in integer32;
                           redCost : out double_float ) is
    begin
      null;
    end calRedCost;

    procedure isZeroDirEle ( this : in Link_to_simplex;
                             termS : in integer32; idx : in integer32;
                             preNbN : in integer32;
                             sub_pivInIdx : out integer32;
                             result : out integer32 ) is
    begin
      null;
    end isZeroDirEle;

    procedure IP_vec_mat ( this : in Link_to_simplex ) is
    begin
      null;
    end IP_vec_mat;

-- reduced cost

    procedure reducedCost_tab_p1
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_tab_p1;

    procedure reducedCost_tab
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_tab;

    procedure reducedCost_p1
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_p1;

    procedure reducedCost
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost;

    procedure reducedCost_Bland
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_Bland;

    procedure reducedCost_mFst
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_mFst;

    procedure reducedCost_iFst
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float;
                  termS : in integer32;
                  reTermS : in integer32;
                  preNbN : in integer32; flag : out integer32 ) is
    begin
      null;
    end reducedCost_iFst;

    procedure extend_nbIdx ( this : in Link_to_simplex;
                             cIdx : in integer32;
                             pre_pivInIdx : in integer32;
                             pre_pivOutIdx : in integer32;
                             pre_length : in integer32;
                             reTermS : in integer32;
                             cnt : in out integer32 ) is
    begin
      null;
    end extend_nbIdx;

    procedure extend_nbIdx_comp
                ( this : in Link_to_simplex;
                  non_basisIdx : out integer32;
                  cIdx : in integer32;
                  pre_pivInIdx : in integer32;
                  pre_pivOutIdx : in integer32;
                  pre_length : in integer32;
                  reTermS : in integer32;
                  cnt : in out integer; flag : out integer32 ) is
    begin
      null;
    end extend_nbIdx_comp;

    procedure getIdx ( this : in Link_to_simplex;
                       level : out integer32;
                       idx : out integer32;
                       idx2 : out integer32;
                       ii : out integer32;
                       d_nbIdx : in out integer32 ) is
    begin
      null;
    end getIdx;

-- ratio test

    procedure ratioTest
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 ) is
    begin
      null;
    end ratioTest;

    procedure ratioTest_artFst
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 ) is
    begin
      null;
    end ratioTest_artFst;

    procedure ratioTest_art
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 ) is
    begin
      null;
    end ratioTest_art;

    procedure ratioTest_art_Bland
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 ) is
    begin
      null;
    end ratioTest_art_Bland;

    function ratioTest_frIdx ( this : Link_to_simplex;
                               pivInIdx : integer32 ) return integer32 is
    begin
      return 0;
    end ratioTest_frIdx;

    procedure IP_mat_vec ( this : in Link_to_simplex;
                           pivInIdx : in integer32 ) is
    begin
      null;
    end IP_mat_vec;

    procedure IP_mat_vec_fst ( this : in Link_to_simplex;
                               pivInIdx : in integer32 ) is
    begin
      null;
    end IP_mat_vec_fst;

    procedure update_p1_d_sol ( this : in Link_to_simplex;
                                pivInIdx : in integer32;
                                sub_pivOutIdx : in integer32 ) is
    begin
      null;
    end update_p1_d_sol;

    procedure modify_p_sol ( this : in Link_to_simplex;
                             pivInIdx : in integer32 ) is
    begin
      null;
    end modify_p_sol;

    procedure calElem ( this : in Link_to_simplex; idx : in integer32 ) is
    begin
      null;
    end calElem;

-- create new basis and nonbasis

    procedure createNewBandN_tab
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float ) is
    begin
      null;
    end createNewBandN_tab;

    procedure createNewBandN_p1
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 ) is
    begin
      null;
    end createNewBandN_p1;

    procedure createNewBandN
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 ) is
    begin
      null;
    end createNewBandN;

    procedure createNewBandN_iFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 ) is
    begin
      null;
    end createNewBandN_iFst;

    procedure createNewBandN_mFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 ) is
    begin
      null;
    end createNewBandN_mFst;

    procedure createNewBandN_art
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 ) is
    begin
      null;
    end createNewBandN_art;

    procedure invB_in ( this : in Link_to_simplex;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ) is
    begin
      null;
    end invB_in;

    function invB_out ( this : Link_to_simplex;
                        rowIdx : integer32;
                        colIdx : integer32 ) return double_float is
    begin
      return 0.0;
    end invB_out;

    function transMat_out ( this : Link_to_simplex;
                            rowIdx : integer32;
                            colIdx : integer32 ) return double_float is
    begin
      return 0.0;
    end transMat_out;

    procedure supp_in ( this : in Link_to_simplex;
                        lvl : in integer32;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ) is
    begin
      null;
    end supp_in;

    function supp_out ( this : Link_to_simplex;
                        lvl : integer32;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return 0.0;
    end supp_out;

    function isZero ( this : Link_to_simplex;
                      val : double_float ) return integer32 is
    begin
      return 0;
    end isZero;

    procedure info_p_sol ( this : in Link_to_simplex ) is
    begin
      null;
    end info_p_sol;

    procedure info_d_sol ( this : in Link_to_simplex ) is
    begin
      null;
    end info_d_sol;

    procedure info_p1_d_sol ( this : in Link_to_simplex ) is
    begin
      null;
    end info_p1_d_sol;

    procedure info_invB ( this : in Link_to_simplex ) is
    begin
      null;
    end info_invB;

    procedure info_transMat ( this : in Link_to_simplex ) is
    begin
      null;
    end info_transMat;

    procedure info_transRed ( this : in Link_to_simplex ) is
    begin
      null;
    end info_transRed;

    procedure info_basisIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_basisIdx;

    procedure info_nf_pos ( this : in Link_to_simplex ) is
    begin
      null;
    end info_nf_pos;

    procedure info_nbIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_nbIdx;

    procedure info_rIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_rIdx;

    procedure info_redVec ( this : in Link_to_simplex ) is
    begin
      null;
    end info_redVec;

    procedure info_dir ( this : in Link_to_simplex ) is
    begin
      null;
    end info_dir;

    procedure info_frIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_frIdx;

    procedure info_candIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_candIdx;

    procedure info_repIdx ( this : in Link_to_simplex ) is
    begin
      null;
    end info_repIdx;

    procedure info_oriSup ( this : in Link_to_simplex ) is
    begin
      null;
    end info_oriSup;

    function new_simplex return simplex is

      res : simplex;

    begin
      return res;
    end new_simplex;

    procedure delete_simplex ( this : in Link_to_simplex ) is
    begin
      null;
    end delete_simplex;

    procedure get_iNbN_nfN
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData;
                  lNbN : in integer32;
                  lNfN : in integer32 ) is
    begin
      null;
    end get_iNbN_nfN;

    procedure get_mNbN_nfN
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                ) is
    begin
      null;
    end get_mNbN_nfN;

    procedure get_repIdx_candIdx
                ( this : in Link_to_simplex;
                  ori_candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_repIdx : in integer32 ) is
    begin
      null;
    end get_repIdx_candIdx;

    procedure get_parent
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end get_parent;

    procedure get_cur
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                ) is
    begin
      null;
    end get_cur;

    procedure get_res ( this : in Link_to_simplex; 
                        iData : in demics_ftest.class_ftData.Link_to_ftData
                      ) is
    begin
      null;
    end get_res;

    procedure get_pivOutNum
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                ) is
    begin
      null;
    end get_pivOutNum;

    procedure get_nbN_nfN ( this : in Link_to_simplex;
                            ori_nbN : in integer32;
                            ori_nfN : in integer32 ) is
    begin
      null;
    end get_nbN_nfN;

    procedure get_p_sol
                ( this : in Link_to_simplex;
                  ori_p_sol : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      null;
    end get_p_sol;

    procedure get_d_sol
                ( this : in Link_to_simplex;
                  ori_d_sol : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      null;
    end get_d_sol;

    procedure get_basisIdx
                ( this : in Link_to_simplex;
                  ori_basisIdx : in Standard_Integer_Vectors.Link_to_Vector
                ) is
    begin
      null;
    end get_basisIdx; 

    procedure get_nf_pos
                ( this : in Link_to_simplex;
                  ori_nf_pos : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      null;
    end get_nf_pos;

    procedure get_nbIdx
                ( this : in Link_to_simplex;
                  ori_nbIdx : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      null;
    end get_nbIdx;

    procedure get_invB
                ( this : in Link_to_simplex;
                  invB : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      null;
    end get_invB;

    procedure get_frIdx ( this : in Link_to_simplex;
                          ori_frIdx : in integer32 ) is
    begin
      null;
    end get_frIdx;

    procedure copy_p1_d_sol
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end copy_p1_d_sol;

    procedure copy_eye
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                ) is
    begin
      null;
    end copy_eye;

    procedure allocateAndIni
                ( this : in Link_to_simplex;
                  data : demics_input_data.class_dataSet.dataSet;
                  ori_firIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  seedNum : in integer32;
                  ori_output : in integer32 ) is
    begin
      null;
    end allocateAndIni;

-- for relation table

    procedure tSolLP ( this : in Link_to_simplex;
                       iter : in out integer32;
                       mode : in integer32; flag : out integer32 ) is
    begin
      null;
    end tSolLP;

-- for phase 1 and 2

    procedure fSolLP ( this : in Link_to_simplex;
                       termS : in integer32;
                       reTermS : in integer32;
                       iter : in out integer32; flag : out integer32 ) is
    begin
      null;
    end fSolLP;

-- iCheck

    procedure fstRed_candIdx
                 ( this : in Link_to_simplex;
                   curInif : in demics_iTest.class_inifData.Link_to_inifData;
                   mCandIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                   pivInIdx : out integer32;
                   sub_pivInIdx : out integer32 ) is
    begin
      null;
    end fstRed_candIdx;

    procedure cal_redVec
                ( this : in Link_to_simplex;
                  termS : in integer32;
                  reTermS : in integer32;
                  fst_pivInIdx : in integer32;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                ) is
    begin
      null;
    end cal_redVec;

    function put_redCost
               ( this : Link_to_simplex;
                 fst_pivInIdx : integer32 ) return double_float is
    begin
      return 0.0;
    end put_redCost;

-- iCheck_art

    procedure solLP_art ( this : in Link_to_simplex;
                          depth : in integer32;
                          idx_one : in integer32;
                          fst_pivIn : in integer32;
                          preNbN : in integer32;
                          termS : in integer32;
                          reTermS : in integer32;
                          iter : in out integer32; flag : out integer32 ) is
    begin
      null;
    end solLP_art;

    procedure solLP_art_Bland ( this : in Link_to_simplex;
                                pivInIdx : in integer32;
                                sub_pivInIdx : in integer32;
                                pivOutIdx : in integer32;
                                sub_pivOutIdx : in integer32;
                                redFlag : in integer32;
                                theta : in double_float;
                                redCost : in double_float;
                                termS : in integer32;
                                reTermS : in integer32;
                                iter : in out integer32;
                                flag : out integer32 ) is
    begin
      null;
    end solLP_art_Bland;

-- for mLP

    procedure solLP ( this : in Link_to_simplex;
                      depth : in integer32;
                      fst_pivInIdx : in integer32;
                      fst_sub_pivInIdx : in integer32;
                      fst_redCost : in double_float;
                      mode : in integer32;
                      termS : in integer32;
                      reTermS : in integer32;
                      preNbN : in integer32;
                      iter : in out integer; flag : out integer32 ) is
    begin
      null;
    end solLP;

    procedure solLP_Bland ( this : in Link_to_simplex;
                            pivInIdx : in integer32;
                            sub_pivInIdx : in integer32;
                            pivOutIdx : in integer32;
                            sub_pivOutIdx : in integer32;
                            redFlag : in integer32;
                            theta : in double_float;
                            redCost : in double_float;
                            termS : in integer32;
                            reTermS : in integer32;
                            iter : in out integer32; flag : out integer32 ) is
    begin
      null;
    end solLP_Bland;

    procedure initIter ( this : in Link_to_simplex;
                         mode : in integer32;
                         fst_pivInIdx : in integer32;
                         fst_sub_pivInIdx : in integer32;
                         fst_redCost : in double_float;
                         redFlag : out integer32;
                         pivInIdx : out integer32;
                         sub_pivInIdx : out integer32;
                         pivOutIdx : out integer32;
                         sub_pivOutIdx : out integer32;
                         theta : out double_float;
                         redCost : out double_float;
                         termS : in integer32;
                         reTermS : in integer32;
                         preNbN : in integer32; flag : out integer32 ) is
    begin
      null;
    end initIter;

    procedure calMixedVol ( this : in Link_to_simplex;
                            lv : in demics_fTest.class_lvData.Link_to_lvData;
                            sp : in Standard_Integer_Vectors.Link_to_Vector;
                            supN : in integer32 ) is
    begin
      null;
    end calMixedVol;

    function lu ( this : Link_to_simplex;
                  n : integer32;
                  a : Standard_Floating_Vectors.Vector ) return double_float is
    begin
      return 0.0;
    end lu;

    function matinv ( this : Link_to_simplex;
                      n : integer32;
                      a : Standard_Floating_Vectors.Link_to_Vector;
                      a_inv : Standard_Floating_Vectors.Link_to_Vector )
                    return double_float is
    begin
      return 0.0;
    end matinv;

    function put_elem_supp ( this : Link_to_simplex;
                             lvl : integer32;
                             idx : integer32;
                             row : integer32;
                             col : integer32 ) return double_float is
    begin
      return 0.0;
    end put_elem_supp;

    procedure mult_elem_supp ( this : in Link_to_simplex;
                               lvl : in integer32;
                               idx : in integer32;
                               row : in integer32;
                               col : in integer32 ) is
    begin
      null;
    end mult_elem_supp;

    procedure check_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  depth : in integer32 ) is
    begin
      null;
    end check_dirRed;

    procedure dbg_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  nextInif : in demics_itest.class_inifData.Link_to_inifData;
                  depth : in integer32 ) is
    begin
      null;
    end dbg_dirRed;

    procedure info_mv ( this : in Link_to_simplex ) is
    begin
      null;
    end info_mv;

    procedure info_allSup ( this : in Link_to_simplex ) is
    begin
      null;
    end info_allSup;

    procedure info_allCostVec ( this : in Link_to_simplex ) is
    begin
      null;
    end info_allCostVec;

    procedure info_lifting ( this : in Link_to_simplex ) is
    begin
      null;
    end info_lifting;

    procedure info_simplexData ( this : in Link_to_simplex ) is
    begin
      null;
    end info_simplexData;

  end class_simplex;

end demics_simplex;
