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

    type VecVec_of_supportSets is
      array ( integer32 range <> ) of Link_to_Array_of_supportSets;
    type Link_to_VecVec_of_supportSets is access VecVec_of_supportSets;

    function new_supportSet return supportSet;

    -- DESCRIPTION :
    --   Returns a supportSet record with zero and null values.

    procedure delete_supportSet ( this : in Link_to_supportSet );

    -- DESCRIPTION :
    --   Deallocates the supMath and costVec fields of the record.

    procedure allocSupp
                ( this : in Link_to_supportSet;
                  data : in demics_input_data.class_dataSet.dataSet;
                  level : in integer32;
                  num : in integer32;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates and initializes the support set.

    procedure allocAux
                ( this : in Link_to_supportSet;
                  data : in demics_input_data.class_dataSet.dataSet;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocation and initialization of the last support set.

    procedure supMat_in ( this : in Link_to_supportSet;
                          rowIdx : in integer32;
                          colIdx : in integer32;
                          elem : in double_float );

    -- DESCRIPTION :
    --   Sets the value in this.supMat defined by rowIdx and colIdx to elem.

    procedure supMat_neg ( this : in Link_to_supportSet;
                           rowIdx : in integer32;
                           colIdx : in integer32 );

    -- DESCRIPTION :
    --   Flips the sign of the value in this.supMat,
    --   as defined by rowIdx and colIdx.

    function supMat_out ( this : Link_to_supportSet;
                          rowIdx : integer32;
                          colIdx : integer32 ) return double_float;

    -- DESCRIPTION :
    --   Returns the value in this.supMat as defined by rowIdx and colIdx.

    function redVal ( this : Link_to_supportSet;
                      d_sol : Standard_Floating_Vectors.Link_to_Vector;
                      idx : integer32;
                      ii : integer32 ) return double_float;

    -- DESCRIPTION :
    --   Returns the value of the solution.

    procedure info_sup ( this : in Link_to_supportSet );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.supMat.

    procedure info_costVec ( this : in Link_to_supportSet );

    -- DESCRIPTION :
    --   Writes the numbers stored in this.costVec.

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
      Supp : Link_to_VecVec_of_supportSets;
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

    type Link_to_simplex is access simplex;

-- relation table

    function checkFrIdx ( this : Link_to_simplex ) return integer32;

    procedure elimFrIdx ( this : in Link_to_simplex;
                          sub_pivOutIdx : in integer32 );

-- phase 1

    procedure reMakeNonBasisIdx ( this : in Link_to_simplex;
                                  reTermS : in integer32 );

    procedure reMakeNonBasisIdx_tab ( this : in Link_to_simplex );

    procedure elimArt ( this : in Link_to_simplex;
                        depth : in integer32; preNbN : in integer32;
                        termS : in integer32; reTermS : in integer32;
                        iter : in out integer32 );

    procedure calRedCost ( this : in Link_to_simplex;
                           pivInIdx : in integer32;
                           redCost : out double_float );

    procedure isZeroDirEle ( this : in Link_to_simplex;
                             termS : in integer32; idx : in integer32;
                             preNbN : in integer32;
                             sub_pivInIdx : out integer32;
                             result : out integer32 );

    -- isZeroDirEle was a function returning (TRUE) or (FALSE)
    -- with a side effect: assigning to sub_pivInIdx

    procedure IP_vec_mat ( this : in Link_to_simplex );

-- reduced cost

    procedure reducedCost_tab_p1
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab_p1 was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_tab
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_p1
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_tab_p1 was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost was defined as a function, assigning to its arguments
    -- as side effects and returning a flag value

    procedure reducedCost_Bland
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_Bland was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_mFst
                ( this : in Link_to_simplex;
                  enterIdx : out integer32;
                  sub_enterIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float; flag : out integer32 );

    -- reducedCost_mFst was defined as a function, assigning to its
    -- arguments as side effects and returning a flag value

    procedure reducedCost_iFst
                ( this : in Link_to_simplex;
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

    procedure extend_nbIdx ( this : in Link_to_simplex;
                             cIdx : in integer32;
                             pre_pivInIdx : in integer32;
                             pre_pivOutIdx : in integer32;
                             pre_length : in integer32;
                             reTermS : in integer32;
                             cnt : in out integer32 );

    procedure extend_nbIdx_comp
                ( this : in Link_to_simplex;
                  non_basisIdx : out integer32;
                  cIdx : in integer32;
                  pre_pivInIdx : in integer32;
                  pre_pivOutIdx : in integer32;
                  pre_length : in integer32;
                  reTermS : in integer32;
                  cnt : in out integer; flag : out integer32 );

    -- extend_nbIdx_comp was declared as a function, returning a flag,
    -- assigning to non_basisIdx and cnt as side effects

    procedure getIdx ( this : in Link_to_simplex;
                       level : out integer32;
                       idx : out integer32;
                       idx2 : out integer32;
                       ii : out integer32;
                       d_nbIdx : in out integer32 );

-- ratio test

    procedure ratioTest
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_artFst
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_artFst was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_art
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_art was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    procedure ratioTest_art_Bland
                ( this : in Link_to_simplex;
                  redFlag : in integer32;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32;
                  sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 );

    -- ratioTest_art_Bland was declared as a function, returning a flag,
    -- assigning to pivOutIdx, sub_pivOutIdx, and theta as side effects

    function ratioTest_frIdx ( this : Link_to_simplex;
                               pivInIdx : integer32 ) return integer32;

    procedure IP_mat_vec ( this : in Link_to_simplex;
                           pivInIdx : in integer32 );

    procedure IP_mat_vec_fst ( this : in Link_to_simplex;
                               pivInIdx : in integer32 );

    procedure update_p1_d_sol ( this : in Link_to_simplex;
                                pivInIdx : in integer32;
                                sub_pivOutIdx : in integer32 );

    procedure modify_p_sol ( this : in Link_to_simplex;
                             pivInIdx : in integer32 );

    procedure calElem ( this : in Link_to_simplex; idx : in integer32 );

-- create new basis and nonbasis

    procedure createNewBandN_tab
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float );

    procedure createNewBandN_p1
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_iFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_mFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  theta : in double_float;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure createNewBandN_art
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32;
                  sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : in double_float;
                  termS : in integer32;
                  reTermS : in integer32 );

    procedure invB_in ( this : in Link_to_simplex;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float );

    function invB_out ( this : Link_to_simplex;
                        rowIdx : integer32;
                        colIdx : integer32 ) return double_float;

    function transMat_out ( this : Link_to_simplex;
                            rowIdx : integer32;
                            colIdx : integer32 ) return double_float;

    procedure supp_in ( this : in Link_to_simplex;
                        lvl : in integer32;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ); 

    function supp_out ( this : Link_to_simplex;
                        lvl : integer32;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float;

    function isZero ( this : Link_to_simplex;
                      val : double_float ) return integer32;

    procedure info_p_sol ( this : in Link_to_simplex );

    procedure info_d_sol ( this : in Link_to_simplex );

    procedure info_p1_d_sol ( this : in Link_to_simplex );

    procedure info_invB ( this : in Link_to_simplex );

    procedure info_transMat ( this : in Link_to_simplex );

    procedure info_transRed ( this : in Link_to_simplex );

    procedure info_basisIdx ( this : in Link_to_simplex );

    procedure info_nf_pos ( this : in Link_to_simplex );

    procedure info_nbIdx ( this : in Link_to_simplex );

    procedure info_rIdx ( this : in Link_to_simplex );

    procedure info_redVec ( this : in Link_to_simplex );

    procedure info_dir ( this : in Link_to_simplex );

    procedure info_frIdx ( this : in Link_to_simplex );

    procedure info_candIdx ( this : in Link_to_simplex );

    procedure info_repIdx ( this : in Link_to_simplex );

    procedure info_oriSup ( this : in Link_to_simplex );

    function new_simplex return simplex;

    -- DESCRIPTION :
    --   Returns a simplex record with zero and null values.

    procedure delete_simplex ( this : in Link_to_simplex );

    procedure get_iNbN_nfN
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData;
                  lNbN : in integer32;
                  lNfN : in integer32 );

    procedure get_mNbN_nfN
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                );

    procedure get_repIdx_candIdx
                ( this : in Link_to_simplex;
                  ori_candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_repIdx : in integer32 );

    procedure get_parent
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData );

    procedure get_cur
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                );

    procedure get_res ( this : in Link_to_simplex; 
                        iData : in demics_ftest.class_ftData.Link_to_ftData );

    procedure get_pivOutNum
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                );

    procedure get_nbN_nfN ( this : in Link_to_simplex;
                            ori_nbN : in integer32;
                            ori_nfN : in integer32 );

    procedure get_p_sol
                ( this : in Link_to_simplex;
                  ori_p_sol : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_d_sol
                ( this : in Link_to_simplex;
                  ori_d_sol : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_basisIdx
                ( this : in Link_to_simplex;
                  ori_basisIdx : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_nf_pos
                ( this : in Link_to_simplex;
                  ori_nf_pos : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_nbIdx
                ( this : in Link_to_simplex;
                  ori_nbIdx : in Standard_Integer_Vectors.Link_to_Vector );

    procedure get_invB
                ( this : in Link_to_simplex;
                  invB : in Standard_Floating_Vectors.Link_to_Vector );

    procedure get_frIdx ( this : in Link_to_simplex;
                          ori_frIdx : in integer32 );

    procedure copy_p1_d_sol
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData );

    procedure copy_eye
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                );

    procedure allocateAndIni
                ( this : in Link_to_simplex;
                  data : in demics_input_data.class_dataSet.dataSet;
                  ori_firIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  seedNum : in integer32; ori_output : in integer32;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates and initializes the simplex record.
    --   Generates random lifting values, setting the seed to seedNum.

-- for relation table

    procedure tSolLP ( this : in Link_to_simplex;
                       iter : in out integer32;
                       mode : in integer32; flag : out integer32 );

    -- tSolLP was declared as a function, assigning to iter
    -- as a side effect, and returning a flag

-- for phase 1 and 2

    procedure fSolLP ( this : in Link_to_simplex;
                       termS : in integer32;
                       reTermS : in integer32;
                       iter : in out integer32; flag : out integer32 );

    -- fSolLP was declared as a function, assigning to iter
    -- as a side effect, and returning a flag

-- iCheck

    procedure fstRed_candIdx
                 ( this : in Link_to_simplex;
                   curInif : in demics_iTest.class_inifData.Link_to_inifData;
                   mCandIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                   pivInIdx : out integer32;
                   sub_pivInIdx : out integer32 );

    procedure cal_redVec
                ( this : in Link_to_simplex;
                  termS : in integer32;
                  reTermS : in integer32;
                  fst_pivInIdx : in integer32;
                  cur : in demics_ftest.class_theData.Link_to_Array_of_theData
                );

    function put_redCost
               ( this : Link_to_simplex;
                 fst_pivInIdx : integer32 ) return double_float;

-- iCheck_art

    procedure solLP_art ( this : in Link_to_simplex;
                          depth : in integer32;
                          idx_one : in integer32;
                          fst_pivIn : in integer32;
                          preNbN : in integer32;
                          termS : in integer32;
                          reTermS : in integer32;
                          iter : in out integer32; flag : out integer32 );

    -- solLP_art was declared as a function, updating iter as a side effect,
    -- and returning a flag

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
                                flag : out integer32 );

    -- solLP_art_bland was declared as a function, updating iter
    -- as a side effect, and returning a flag

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
                      iter : in out integer; flag : out integer32 );

    -- solLP was declared as a function, updating iter as a side effect,
    -- and returning a flag value

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
                            iter : in out integer32; flag : out integer32 );

    -- solLP_Bland was declared as a function, updating iter 
    -- as a side effect, and returning a flag value

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
                         preNbN : in integer32; flag : out integer32 );

    -- initIter was declared as a function, returning a flag,
    -- with many assignments to the arguments as side effects

    procedure calMixedVol ( this : in Link_to_simplex;
                            lv : in demics_fTest.class_lvData.Link_to_lvData;
                            sp : in Standard_Integer_Vectors.Link_to_Vector;
                            supN : in integer32 );

    function lu ( this : Link_to_simplex;
                  n : integer32;
                  a : Standard_Floating_Vectors.Vector ) return double_float;

    function matinv ( this : Link_to_simplex;
                      n : integer32;
                      a : Standard_Floating_Vectors.Link_to_Vector;
                      a_inv : Standard_Floating_Vectors.Link_to_Vector )
                    return double_float;

    function put_elem_supp ( this : Link_to_simplex;
                             lvl : integer32;
                             idx : integer32;
                             row : integer32;
                             col : integer32 ) return double_float;

    procedure mult_elem_supp ( this : in Link_to_simplex;
                               lvl : in integer32;
                               idx : in integer32;
                               row : in integer32;
                               col : in integer32 );

    procedure check_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  depth : in integer32 );

    procedure dbg_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  nextInif : in demics_itest.class_inifData.Link_to_inifData;
                  depth : in integer32 );

    procedure info_mv ( this : in Link_to_simplex );

    procedure info_allSup ( this : in Link_to_simplex );

    procedure info_allCostVec ( this : in Link_to_simplex );

    procedure info_lifting ( this : in Link_to_simplex );

    procedure info_simplexData ( this : Link_to_simplex );

  end class_simplex;

end demics_simplex;
