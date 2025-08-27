with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with demics_input_data;
with demics_itest;
with demics_ftest;
with demics_reltab;
with demics_simplex;

package demics_mvc is

-- DESCRIPTION :
--   Enumerates all mixed cells with DEMiCs.
--   Translation initiated by g++ -c -fdump-ada-spec mvc.h.

  package class_mvc is

    type mvc is record
      dim : integer32;
      supN : integer32;
      row : integer32;
      col : integer32;
      termSumNum : integer32;
      termMax : integer32;
      maxLength : integer32;
      total_iter : double_float;
      total_feasLP : double_float;
      total_LPs : double_float;
      total_1PT : double_float;
      total_2PT : double_float;
      total_triLPs_mLP : double_float;
      total_unbLP_tab : double_float;
      lvl_1PT : Standard_Floating_Vectors.Link_to_Vector;
      lvl_2PT : Standard_Floating_Vectors.Link_to_Vector;
      actNode : Standard_Floating_Vectors.Link_to_Vector;
      mfNum : Standard_Integer_Vectors.Link_to_Vector;
      firIdx : Standard_Integer_Vectors.Link_to_Vector;
      repN : Standard_Integer_Vectors.Link_to_Vector;
      negIdx : Standard_Integer_VecVecs.Link_to_VecVec;
      termSet : Standard_Integer_Vectors.Link_to_Vector;
      termStart : Standard_Integer_Vectors.Link_to_Vector;
      re_termStart : Standard_Integer_Vectors.Link_to_Vector;
      supType : Standard_Integer_Vectors.Link_to_Vector;
      mRepN : Standard_Integer_Vectors.Link_to_Vector;
      mFeaIdx : Standard_Integer_VecVecs.Link_to_VecVec;
      mFea : Standard_Integer_Vectors.Link_to_Vector;
      trNeg : Standard_Integer_VecVecs.Link_to_VecVec;
      sp : Standard_Integer_Vectors.Link_to_Vector;
      candIdx : Standard_Integer_Vectors.Link_to_Vector;
      trMat : Standard_Floating_Vectors.Link_to_Vector;
      table : Standard_Integer_Vectors.Link_to_Vector;
      the_Reltab : demics_reltab.class_reltab.reltab;
      the_Simplex : demics_simplex.class_simplex.simplex;
      lv : access demics_ftest.class_lvData.lvData;
      iLv : access demics_itest.class_iLvData.iLvData;
    end record;

    procedure getMemory
                ( this : access mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32 );

    procedure initMemoryCheck
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32 );

    procedure memoryCheck
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32 );

    procedure get_candIdx
                ( this : access mvc;
                  curInif : access demics_itest.class_inifData.inifData );

    function chooseSup
                ( this : access mvc;
                  depth : integer32;
                  curNode : access demics_ftest.class_theData.theData;
                  curInif : access demics_itest.class_inifData.inifData;
                  nextInif : access demics_itest.class_inifData.inifData )
                return integer32;

    procedure fUpdateDirRed
                ( this : access mvc;
                  curInif : demics_itest.class_inifData.Array_of_inifData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.Class_theData.theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 );

    procedure updateDirRed
                ( this : access mvc;
                  curInif : demics_itest.class_inifData.Array_of_inifData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 );

    function findUnbDir
                ( this : access mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32;

    function findUnbDir_art
                ( this : access mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) return integer32;

    function checkDir
                ( this : access mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : access demics_itest.class_uData.uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32;

    function checkDir_art
                ( this : access mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : access demics_itest.class_uData.uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32;

    procedure skipPtr
                ( this : access mvc;
                  curr : demics_itest.class_uData.Link_to_Array_of_uData;
                  fHead : demics_itest.class_uData.Link_to_Array_of_uData );

    procedure get_tuple_index
                ( this : access mvc;
                  node : access demics_ftest.class_ftData.ftData;
                  data : access demics_ftest.class_ftData.ftData;
                  length : in integer32 );

    procedure dbg_init_transMat
                ( this : access mvc;
                  curNode : access demics_ftest.class_theData.theData );

    procedure dbg_transMat
                ( this : access mvc;
                  preNode : access demics_ftest.class_theData.theData;
                  curNode : access demics_ftest.class_theData.theData );

    procedure check_transMat
                ( this : access mvc;
                  preNode : access demics_ftest.class_theData.theData;
                  curNode : access demics_ftest.class_theData.theData );

    procedure check_init_transRed
                ( this : access mvc;
                  curNode : access demics_ftest.class_theData.theData );

    function checkSign_red
                ( this : access mvc;
                  curRed : double_float;
                  tarRed : double_float ) return integer32;

    function checkNonNeg_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32;

    function checkNonPos_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32;

    function checkZero_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float) return integer32;

    function table_out
                ( this : access mvc;
                  row : integer32; col : integer32 ) return integer32;

    procedure info_neg
                ( this : access mvc;
                  termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec );

    procedure info_sp ( this : access mvc; depth : in integer32 );

    procedure info_parent_node ( this : access mvc; depth : in integer32 );

    procedure info_tuple
                ( this : access mvc;
                  lvl : in integer32; depth : in integer32 );

    procedure info_all_dirRed
                ( this : access mvc;
                  depth : in integer32;
                  node : access demics_ftest.class_ftData.ftData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData );

    procedure info_mFea ( this : access mvc; length : in integer32 );

    procedure info_firIdx ( this : access mvc; length : in integer32 );

    procedure info_fIdx
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData );

    procedure info_candIdx ( this : access mvc );

    procedure info_elemNum
                ( this : access mvc;
                  length : in integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  node : demics_ftest.class_ftData.ftData );

    procedure info_prop_elemNum
                ( this : access mvc;
                  length : in integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  node : demics_ftest.class_ftData.ftData );

    procedure info_table ( this : access mvc );

    function new_mvc return mvc;

    procedure delete_mvc ( this : access mvc );

    procedure allocateAndIni
                ( this : access mvc;
                  data : access demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32 );

    procedure initFeasTest ( this : access mvc; depth : in integer32 );

    procedure initCheck
                ( this : access mvc;
                  depth : in integer32;
                  data : access demics_ftest.class_ftData.ftData );

    procedure initLP
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  depth : in integer32; idx : in integer32;
                  feaNum : in out integer32 );

    function feasTest
                ( this : access mvc;
                  depth : integer32;
                  parent : access demics_ftest.class_theData.theData )
                return integer32;

    procedure upFeasTest
                ( this : access mvc; depth : in out integer32;
                  flag : out integer32 );

    -- upFeasTest was declared as a function, assigning to depth
    -- as a side effect, returning a flag

    procedure findMixedCell
                ( this : access mvc;
                  depth : in integer32;
                  parent : access demics_ftest.class_theData.theData );

    procedure findAllMixedCells ( this : access mvc; depth : in integer32 );

    function iCheck
                ( this : access mvc;
                  depth : integer32;
                  parent : access demics_ftest.class_theData.theData;
                  data : access demics_ftest.class_ftData.ftData;
                  inifData : access demics_itest.class_inifData.inifData )
                return integer32;

    procedure iLP ( this : access mvc;
                    parent : access demics_ftest.class_theData.theData;
                    data : access demics_ftest.class_ftData.ftData;
                    depth : in integer32;
                    idx_one : in integer32;
                    fst_pivInIdx : in integer32;
                    sub_fst_pivInIdx : in integer32;
                    preNbN : in integer32;
                    feaNum : in out integer32 );

    procedure iLP_Art
                ( this : access mvc;
                  parent : access demics_ftest.class_theData.theData;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32;
                  idx_one : in integer32;
                  fst_pivInIdx : in integer32;
                  sub_fst_pivInIdx : in integer32;
                  preNbN : in integer32;
                  feaNum : in out integer32 );

    procedure findNode
                ( this : access mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  flag : out integer32 );

    -- findNode was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag

    procedure findNextNode
                ( this : access mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  Data : access demics_ftest.class_ftData.ftData;
                  flag : out integer32 );

    -- findNode was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag

    procedure findUpNode
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  pre : demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  cur : demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  lvl : in out integer32;
                  polyDim : in integer32;
                  depth : in integer32 );

    procedure mLP ( this : access mvc;
                    Pre : access demics_ftest.class_ftData.ftData;
                    Cur : access demics_ftest.class_ftData.ftData;
                    data : access demics_ftest.class_ftData.ftData;
                    repIdx : in Standard_Integer_Vectors.Link_to_Vector;
                    feaIdx : in Standard_Integer_Vectors.Link_to_Vector;
                    tarIdx : in integer32;
                    mRepN : in Standard_Integer_Vectors.Link_to_Vector;
                    totalN : in integer32;
                    depth : in integer32;
                    feaNum : in out integer32;
                    lvl : in integer32;
                    length : in integer32; flag : out integer32 );

    -- mLP was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag

    function checkBasis
                ( this : access mvc;
                  target : access demics_ftest.class_theData.theData;
                  sub_sIdx : integer32 ) return integer32;

    function checkAnotherBasis
                ( this : access mvc;
                  repIdx : integer32; dist : integer32;
                  target : demics_ftest.class_theData.Link_to_Array_of_theData
                ) return integer32;

    procedure get_firIdx
                ( this : access mvc;
                  data_a : demics_ftest.class_ftData.ftData;
                  data_b : demics_ftest.class_ftData.ftData;
                  sn : in integer32; lvl : in integer32 );

    procedure info_cpuTime
                ( this : access mvc;
                  cpuTime_start : double_float;
                  cpuTime_end : double_float );

    procedure info_final ( this : access mvc );

    procedure enum ( this : access mvc );

  end class_mvc;

end demics_mvc;
