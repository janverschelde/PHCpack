with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with demics_input_data;
with demics_itest;
with demics_ftest;

package demics_simplex is

-- DESCRIPTION :
--   Defines the simple algorithm as used in DEMiCs.
--   Translation initiated by g++ -c -fdump-ada-spec simplex.h.

  package class_supportSet is

    type supportSet is record
      row : integer32;
      col : integer32;
      supMat : Standard_Floating_Vectors.Link_to_Vector;
      costVec : Standard_Floating_Vectors.Link_to_Vector;
    end record;

    type Link_to_supportSet is access supportSet;
    type Array_of_supportSets is
      array ( integer32 range <> ) of Link_to_supportSet;
    type Link_to_Array_of_supportSets is access Array_of_supportSets;

    function new_supportSet return supportSet;

    procedure delete_supportSet ( this : access supportSet );

    procedure allocSupp
                ( this : access supportSet;
                  data : access demics_input_data.class_dataSet.dataSet;
                  level : in integer32;
                  num : in integer32;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector );

    procedure allocAux
                ( this : access supportSet;
                  data : access demics_input_data.class_dataSet.dataSet );

    procedure supMat_in ( this : access supportSet;
                          rowIdx : in integer32;
                          colIdx : in integer32;
                          elem : in double_float );

    procedure supMat_neg ( this : access supportSet;
                           rowIdx : in integer32;
                           colIdx : in integer32 );

    function supMat_out ( this : access supportSet;
                          rowIdx : integer32;
                          colIdx : integer32 ) return double_float;

    function redVal ( this : access supportSet;
                      d_sol : Standard_Floating_Vectors.Link_to_Vector;
                      idx : integer32;
                      ii : integer32 ) return double_float;

    procedure info_sup ( this : access supportSet );

    procedure info_costVec ( this : access supportSet );

  end class_supportSet;

  use class_supportSet;

  package class_simplex is

    type simplex is record
      dim : integer32;
      supN : integer32;
      termSumNum : integer32;
      repIdx : integer32;
      candIdx : Standard_Integer_Vectors.Link_to_Vector;
      firIdx : Standard_Integer_Vectors.Link_to_Vector;
      termSet : Standard_Integer_Vectors.Link_to_Vector;
      termStart : Standard_Integer_Vectors.Link_to_Vector;
      re_termStart : Standard_Integer_Vectors.Link_to_Vector;
      output : integer32;
      mixedVol : double_float;
      mixedCell : integer32;
      ip : Standard_Integer_Vectors.Link_to_Vector;
      weight : Standard_Floating_Vectors.Link_to_Vector;
      vol : Standard_Floating_Vectors.Link_to_Vector;
      eye : Standard_Floating_Vectors.Link_to_Vector;
      nbN : integer32;
      nfN : integer32;
      artV : integer32;
      pivOutNum : integer32;
      frIdx :  integer32;
      Supp : Link_to_Array_of_supportSets;
      oriSupp : Standard_Floating_VecVecs.Link_to_VecVec;
      invB : Standard_Floating_Vectors.Link_to_Vector;  -- row oriented
      transMat : Standard_Floating_Vectors.Link_to_Vector; -- row oriented
      transRed : Standard_Floating_Vectors.Link_to_Vector;
      p_sol : Standard_Floating_Vectors.Link_to_Vector;
      d_sol : Standard_Floating_Vectors.Link_to_Vector;
      p1_d_sol : Standard_Floating_Vectors.Link_to_Vector;
      fst_d_sol : Standard_Floating_Vectors.Link_to_Vector;
      aux_cvec : Standard_Floating_Vectors.Link_to_Vector;
      dir : Standard_Floating_Vectors.Link_to_Vector;
      fst_redVec : Standard_Floating_Vectors.Link_to_Vector;
      redVec : Standard_Floating_Vectors.Link_to_Vector;
      basisIdx : Standard_Integer_Vectors.Link_to_Vector;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      nbIdx : Standard_Integer_Vectors.Link_to_Vector;
      rIdx : Standard_Integer_Vectors.Link_to_Vector;
      pivOutList : Standard_Integer_Vectors.Link_to_Vector;
      pivOutCheck : Standard_Integer_Vectors.Link_to_Vector;
      tmp_newInvB : Standard_Floating_Vectors.Link_to_Vector;
      tmp_transMat : Standard_Floating_Vectors.Link_to_Vector;
      nIdx : Standard_Integer_Vectors.Link_to_Vector;
      pre_p_sol : Standard_Floating_Vectors.Link_to_Vector;
      pre_d_sol : Standard_Floating_Vectors.Link_to_Vector;
      pre_redVec : Standard_Floating_Vectors.Link_to_Vector;
      pre_basisIdx : Standard_Integer_Vectors.Link_to_Vector;
      pre_nbIdx : Standard_Integer_Vectors.Link_to_Vector;
      pre_nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      pre_invB : Standard_Floating_Vectors.Link_to_Vector;
      pre_transMat : Standard_Floating_Vectors.Link_to_Vector;
      pre_transRed : Standard_Floating_Vectors.Link_to_Vector;
      lifting : Standard_Floating_Vectors.Link_to_Vector;
    end record;

-- relation table

    function checkFrIdx ( this : access simplex ) return integer32;

    procedure elimFrIdx ( this : access simplex;
                          sub_pivOutIdx : in integer32 );

-- phase 1

    procedure reMakeNonBasisIdx ( this : access simplex;
                                  reTermS : in integer32 );

    procedure reMakeNonBasisIdx_tab ( this : access simplex );

    procedure elimArt ( this : access simplex;
                        depth : in integer32; preNbN : in integer32;
                        termS : in integer32; reTermS : in integer32;
                        iter : in out integer32 );

    procedure calRedCost ( this : access simplex;
                           pivInIdx : in integer32;
                           redCost : out double_float );

    procedure isZeroDirEle ( this : access simplex;
                             termS : integer32; idx : integer32;
                             preNbN : integer32;
                             sub_pivInIdx : out integer32;
                             result : out integer32 );

    -- isZeroDirEle was a function returning (TRUE) or (FALSE)
    -- with a side effect: assigning to sub_pivInIdx

    procedure IP_vec_mat ( this : access simplex );

-- reduced cost

    procedure reducedCost_tab_p1
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab_p1 was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_tab
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_p1
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab_p1 was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost was defined as a function, assigning to its arguments
    -- as side effects and returning a flag value

    procedure reducedCost_Bland
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_Bland was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_mFst
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_mFst was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_iFst
                ( this : access simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float;
                  termS : in integer32;
                  reTermS : in integer32;
                  preNbN : in integer32; flag : out integer32 );

    -- reducedCost_iFst was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure extend_nbIdx ( this : access simplex;
                             cIdx : in integer32;
                             pre_pivInIdx : in integer32;
                             pre_pivOutIdx : in integer32;
                             pre_length : in integer32;
                             reTermS : in integer32;
                             cnt : in out integer32 );

    procedure extend_nbIdx_comp
                ( this : access simplex;
                  non_basisIdx : out integer32;
                  cIdx : in integer32;
                  pre_pivInIdx : in integer32;
                  pre_pivOutIdx : in integer32;
                  pre_length : in integer32;
                  reTermS : in integer32;
                  cnt : in out integer; flag : out integer32 );

    -- extend_nbIdx_comp was declared as a function, returning a flag,
    -- assigning to non_basisIdx and cnt as side effects

    procedure getIdx ( this : access simplex;
                       level : out integer32;
                       idx : out integer32;
                       idx2 : out integer32;
                       ii : out integer32;
                       d_nbIdx : in out integer32 );

-- ratio test

    procedure ratioTest
                ( this : access simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_artFst
                ( this : access simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_artFst was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_art
                ( this : access simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_art was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_art_Bland
                ( this : access simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_art_Bland was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    function ratioTest_frIdx ( this : access simplex;
                               pivInIdx : integer32 ) return integer32;

    procedure IP_mat_vec ( this : access simplex;
                           pivInIdx : in integer32 );

    procedure IP_mat_vec_fst ( this : access simplex;
                               pivInIdx : in integer32 );

    procedure update_p1_d_sol ( this : access simplex;
                                pivInIdx : in integer32;
                                sub_pivOutIdx : in integer32 );

    procedure modify_p_sol ( this : access simplex;
                             pivInIdx : in integer32 );

    procedure calElem ( this : access simplex; idx : in integer32 );

-- create new basis and nonbasis

    procedure createNewBandN_tab
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float );

    procedure createNewBandN_p1
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_iFst
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_mFst
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_art
                ( this : access simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure invB_in ( this : access simplex;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float );

    function invB_out ( this : access simplex;
                        rowIdx : integer32;
                        colIdx : integer32 ) return double_float;

    function transMat_out ( this : access simplex;
                            rowIdx : integer32;
                            colIdx : integer32 ) return double_float;

    procedure supp_in ( this : access simplex;
                        lvl : in integer32;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ); 

    function supp_out ( this : access simplex;
                        lvl : integer32;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float;

    function isZero ( this : access simplex;
                      val : double_float ) return integer32;

    procedure info_p_sol ( this : access simplex );

    procedure info_d_sol ( this : access simplex );

    procedure info_p1_d_sol ( this : access simplex );

    procedure info_invB ( this : access simplex );

    procedure info_transMat ( this : access simplex );

    procedure info_transRed ( this : access simplex );

    procedure info_basisIdx ( this : access simplex );

    procedure info_nf_pos ( this : access simplex );

    procedure info_nbIdx ( this : access simplex );

    procedure info_rIdx ( this : access simplex );

    procedure info_redVec ( this : access simplex );

    procedure info_dir ( this : access simplex );

    procedure info_frIdx ( this : access simplex );

    procedure info_candIdx ( this : access simplex );

    procedure info_repIdx ( this : access simplex );

    procedure info_oriSup ( this : access simplex );

    function new_simplex return simplex;

    procedure delete_simplex ( this : access simplex );

    procedure get_iNbN_nfN
                ( this : access simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData;
                  lNbN : in integer32;
                  lNfN : in integer32 );

    procedure get_mNbN_nfN
                ( this : access simplex;
                  parent : access demics_ftest.class_theData.theData;
              cur : in demics_ftest.class_theData.Link_to_Array_of_theData );

    procedure get_repIdx_candIdx
                ( this : access simplex;
                  ori_candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_repIdx : in integer32 );

    procedure get_parent
                ( this : access simplex;
                  parent : access demics_ftest.class_theData.theData );

    procedure get_cur
                ( this : access simplex;
              cur : in demics_ftest.class_theData.Link_to_Array_of_theData );

    procedure get_res ( this : access simplex; 
                        iData : access demics_ftest.class_ftData.ftData );

    procedure get_pivOutNum
                ( this : access simplex;
              cur : in demics_ftest.class_theData.Link_to_Array_of_theData );

    procedure get_nbN_nfN ( this : access simplex;
                            ori_nbN : in integer32;
                            ori_nfN : in integer32 );

    procedure get_p_sol
                ( this : access simplex;
                  ori_p_sol : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_d_sol
                ( this : access simplex;
                  ori_d_sol : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_basisIdx
                ( this : access simplex;
                  ori_basisIdx : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_nf_pos
                ( this : access simplex;
                  ori_nf_pos : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_nbIdx
                ( this : access simplex;
                  ori_nbIdx : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_invB
                ( this : access simplex;
                  invB : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_frIdx ( this : access simplex;
                          ori_frIdx : in integer32 );

    procedure copy_p1_d_sol
                ( this : access simplex;
                  cur : access demics_ftest.class_theData.theData );

    procedure copy_eye
                ( this : access simplex;
              cur : in demics_ftest.class_theData.Link_to_Array_of_theData );

    procedure allocateAndIni
                ( this : access simplex;
                  data : access demics_input_data.class_dataSet.dataSet;
                  ori_firIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  seedNum : in integer32;
                  ori_output : in integer32 );

-- for relation table

    procedure tSolLP ( this : access simplex;
                       iter : in out integer32;
                       mode : in integer32; flag : out integer32 );

    -- tSolLP was declared as a function, assigning to iter
    -- as a side effect, and returning a flag

-- for phase 1 and 2

    procedure fSolLP ( this : access simplex;
                       termS : in integer32;
                       reTermS : in integer32;
                       iter : in out integer32; flag : out integer32 );

    -- fSolLP was declared as a function, assigning to iter
    -- as a side effect, and returning a flag

-- iCheck

    procedure fstRed_candIdx
                 ( this : access simplex;
                   curInif : access demics_iTest.class_inifData.inifData;
                   mCandIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                   pivInIdx : out integer32;
                   sub_pivInIdx : out integer32 );

    procedure cal_redVec
                ( this : access simplex;
                  termS : in integer32;
                  reTermS : in integer32;
                  fst_pivInIdx : in integer32;
              cur : in demics_ftest.class_theData.Link_to_Array_of_theData );

    function put_redCost
               ( this : access simplex;
                 fst_pivInIdx : integer32 ) return double_float;

-- iCheck_art

    procedure solLP_art ( this : access simplex;
                          depth : in integer32;
                          idx_one : in integer32;
                          fst_pivIn : in integer32;
                          preNbN : in integer32;
                          termS : in integer32;
                          reTermS : in integer32;
                          iter : in out integer32; flag : out integer32 );

    -- solLP_art was declared as a function, updating iter as a side effect,
    -- and returning a flag

    procedure solLP_art_Bland ( this : access simplex;
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
                                flag : out integer32 );

    -- solLP_art_bland was declared as a function, updating iter
    -- as a side effect, and returning a flag

-- for mLP

    procedure solLP ( this : access simplex;
                      depth : in integer32;
                      fst_pivInIdx : in integer32;
                      fst_sub_pivInIdx : in integer32;
                      fst_redCost : in double_float;
                      mode : in integer32;
                      termS : in integer32;
                      reTermS : in integer32;
                      preNbN : in integer32;
                      iter : in out integer; flag : out integer32 );

    -- solLP was declared as a function, updating iter as a side effect,
    -- and returning a flag value

    procedure solLP_Bland ( this : access simplex;
                            pivInIdx : in integer32;
                            sub_pivInIdx : in integer32;
                            pivOutIdx : in integer32;
                            sub_pivOutIdx : in integer32;
                            redFlag : in integer32;
                            theta : in double_float;
                            redCost : in double_float;
                            termS : in integer32;
                            reTermS : in integer32;
                            iter : in out integer32; flag : out integer32 );

    -- solLP_Bland was declared as a function, updating iter 
    -- as a side effect, and returning a flag value

    procedure initIter ( this : access simplex;
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
                         preNbN : in integer32; flag : out integer32 );

    -- initIter was declared as a function, returning a flag,
    -- with many assignments to the arguments as side effects

    procedure calMixedVol ( this : access simplex;
                            lv : access demics_fTest.class_lvData.lvData;
                            sp : Standard_Integer_Vectors.Link_to_Vector;
                            supN : in integer32 );

    function lu ( this : access simplex;
                  n : integer32;
                  a : Standard_Floating_Vectors.Vector ) return double_float;

    function matinv ( this : access simplex;
                      n : integer32;
                      a : Standard_Floating_Vectors.Link_to_Vector;
                      a_inv : Standard_Floating_Vectors.Link_to_Vector )
                    return double_float;

    function put_elem_supp ( this : access simplex;
                             lvl : integer32;
                             idx : integer32;
                             row : integer32;
                             col : integer32 ) return double_float;

    procedure mult_elem_supp ( this : access simplex;
                               lvl : in integer32;
                               idx : in integer32;
                               row : in integer32;
                               col : in integer32 );

    procedure check_dirRed
                ( this : access simplex;
                  parent : access demics_ftest.class_theData.theData;
                  depth : in integer32 );

    procedure dbg_dirRed
                ( this : access simplex;
                  parent : access demics_ftest.class_theData.theData;
                  nextInif : access demics_itest.class_inifData.inifData;
                  depth : in integer32 );

    procedure info_mv ( this : access simplex );

    procedure info_allSup ( this : access simplex );

    procedure info_allCostVec ( this : access simplex );

    procedure info_lifting ( this : access simplex );

    procedure info_simplexData ( this : access simplex );

  end class_simplex;

end demics_simplex;
