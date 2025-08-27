package body demics_mvc is

  package body class_mvc is

    procedure getMemory
                ( this : access mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32 ) is
    begin
      null;
    end getMemory;

    procedure initMemoryCheck
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32 ) is
    begin
      null;
    end initMemoryCheck;

    procedure memoryCheck
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32 ) is
    begin
      null;
    end memoryCheck;

    procedure get_candIdx
                ( this : access mvc;
                  curInif : access demics_itest.class_inifData.inifData ) is
    begin
      null;
    end get_candIdx;

    function chooseSup
                ( this : access mvc;
                  depth : integer32;
                  curNode : access demics_ftest.class_theData.theData;
                  curInif : access demics_itest.class_inifData.inifData;
                  nextInif : access demics_itest.class_inifData.inifData )
                return integer32 is
    begin
      return 0;
    end chooseSup;

    procedure fUpdateDirRed
                ( this : access mvc;
                  curInif : demics_itest.class_inifData.Array_of_inifData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.Class_theData.theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) is
    begin
      null;
    end fUpdateDirRed;

    procedure updateDirRed
                ( this : access mvc;
                  curInif : demics_itest.class_inifData.Array_of_inifData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) is
    begin
      null;
    end updateDirRed;

    function findUnbDir
                ( this : access mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32 is
    begin
      return 0;
    end findUnbDir;

    function findUnbDir_art
                ( this : access mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : access demics_ftest.class_theData.theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) return integer32 is
    begin
      return 0;
    end findUnbDir_art;

    function checkDir
                ( this : access mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : access demics_itest.class_uData.uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32 is
    begin
      return 0;
    end checkDir;

    function checkDir_art
                ( this : access mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : access demics_itest.class_uData.uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32 is
    begin
      return 0;
    end checkDir_art;

    procedure skipPtr
                ( this : access mvc;
                  curr : demics_itest.class_uData.Link_to_Array_of_uData;
                  fHead : demics_itest.class_uData.Link_to_Array_of_uData ) is
    begin
      null;
    end skipPtr;

    procedure get_tuple_index
                ( this : access mvc;
                  node : access demics_ftest.class_ftData.ftData;
                  data : access demics_ftest.class_ftData.ftData;
                  length : in integer32 ) is
    begin
      null;
    end get_tuple_index;

    procedure dbg_init_transMat
                ( this : access mvc;
                  curNode : access demics_ftest.class_theData.theData ) is
    begin
      null;
    end dbg_init_transMat;

    procedure dbg_transMat
                ( this : access mvc;
                  preNode : access demics_ftest.class_theData.theData;
                  curNode : access demics_ftest.class_theData.theData ) is
    begin
      null;
    end dbg_transMat;

    procedure check_transMat
                ( this : access mvc;
                  preNode : access demics_ftest.class_theData.theData;
                  curNode : access demics_ftest.class_theData.theData ) is
    begin
      null;
    end check_transMat;

    procedure check_init_transRed
                ( this : access mvc;
                  curNode : access demics_ftest.class_theData.theData ) is
    begin
      null;
    end check_init_transRed;

    function checkSign_red
                ( this : access mvc;
                  curRed : double_float;
                  tarRed : double_float ) return integer32 is
    begin
      return 0;
    end checkSign_red;

    function checkNonNeg_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkNonNeg_dir;

    function checkNonPos_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkNonPos_dir;

    function checkZero_dir
                ( this : access mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkZero_dir;

    function table_out
                ( this : access mvc;
                  row : integer32; col : integer32 ) return integer32 is
    begin
      return 0;
    end table_out;

    procedure info_neg
                ( this : access mvc;
                  termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec ) is
    begin
      null;
    end info_neg;

    procedure info_sp ( this : access mvc; depth : in integer32 ) is
    begin
      null;
    end info_sp;

    procedure info_parent_node ( this : access mvc; depth : in integer32 ) is
    begin
      null;
    end info_parent_node;

    procedure info_tuple
                ( this : access mvc;
                  lvl : in integer32; depth : in integer32 ) is
    begin
      null;
    end info_tuple;

    procedure info_all_dirRed
                ( this : access mvc;
                  depth : in integer32;
                  node : access demics_ftest.class_ftData.ftData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_mFea ( this : access mvc; length : in integer32 ) is
    begin
      null;
    end info_mFea;

    procedure info_firIdx ( this : access mvc; length : in integer32 ) is
    begin
      null;
    end info_firIdx;

    procedure info_fIdx
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end info_fIdx;

    procedure info_candIdx ( this : access mvc ) is
    begin
      null;
    end info_candIdx;

    procedure info_elemNum
                ( this : access mvc;
                  length : in integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  node : demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end info_elemNum;

    procedure info_prop_elemNum
                ( this : access mvc;
                  length : in integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  node : demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end info_prop_elemNum;

    procedure info_table ( this : access mvc ) is
    begin
      null;
    end info_table;

    function new_mvc return mvc is
 
      res : mvc;

    begin
      return res;
    end new_mvc;

    procedure delete_mvc ( this : access mvc ) is
    begin
      null;
    end delete_mvc;

    procedure allocateAndIni
                ( this : access mvc;
                  data : access demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32 ) is
    begin
      null;
    end allocateAndIni;

    procedure initFeasTest ( this : access mvc; depth : in integer32 ) is
    begin
      null;
    end initFeasTest;

    procedure initCheck
                ( this : access mvc;
                  depth : in integer32;
                  data : access demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end initCheck;

    procedure initLP
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  depth : in integer32; idx : in integer32;
                  feaNum : in out integer32 ) is
    begin
      null;
    end initLP;

    function feasTest
                ( this : access mvc;
                  depth : integer32;
                  parent : access demics_ftest.class_theData.theData )
                return integer32 is
    begin
      return 0;
    end feasTest;

    procedure upFeasTest
                ( this : access mvc; depth : in out integer32;
                  flag : out integer32 ) is
    begin
      null;
    end upFeasTest;

    procedure findMixedCell
                ( this : access mvc;
                  depth : in integer32;
                  parent : access demics_ftest.class_theData.theData ) is
    begin
      null;
    end findMixedCell;

    procedure findAllMixedCells ( this : access mvc; depth : in integer32 ) is
    begin
      null;
    end findAllMixedCells;

    function iCheck
                ( this : access mvc;
                  depth : integer32;
                  parent : access demics_ftest.class_theData.theData;
                  data : access demics_ftest.class_ftData.ftData;
                  inifData : access demics_itest.class_inifData.inifData )
                return integer32 is
    begin
      return 0;
    end iCheck;

    procedure iLP ( this : access mvc;
                    parent : access demics_ftest.class_theData.theData;
                    data : access demics_ftest.class_ftData.ftData;
                    depth : in integer32;
                    idx_one : in integer32;
                    fst_pivInIdx : in integer32;
                    sub_fst_pivInIdx : in integer32;
                    preNbN : in integer32;
                    feaNum : in out integer32 ) is
    begin
      null;
    end iLP;

    procedure iLP_Art
                ( this : access mvc;
                  parent : access demics_ftest.class_theData.theData;
                  data : access demics_ftest.class_ftData.ftData;
                  depth : in integer32;
                  idx_one : in integer32;
                  fst_pivInIdx : in integer32;
                  sub_fst_pivInIdx : in integer32;
                  preNbN : in integer32;
                  feaNum : in out integer32 ) is
    begin
      null;
    end iLP_Art;

    procedure findNode
                ( this : access mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  data : access demics_ftest.class_ftData.ftData;
                  flag : out integer32 ) is
    begin
      null;
    end findNode;

    procedure findNextNode
                ( this : access mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  Data : access demics_ftest.class_ftData.ftData;
                  flag : out integer32 ) is
    begin
      null;
    end findNextNode;

    procedure findUpNode
                ( this : access mvc;
                  data : access demics_ftest.class_ftData.ftData;
                  pre : demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  cur : demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  lvl : in out integer32;
                  polyDim : in integer32;
                  depth : in integer32 ) is
    begin
      null;
    end findUpNode;

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
                    length : in integer32; flag : out integer32 ) is
    begin
      null;
    end mLP;

    function checkBasis
                ( this : access mvc;
                  target : access demics_ftest.class_theData.theData;
                  sub_sIdx : integer32 ) return integer32 is
    begin
      return 0;
    end checkBasis;

    function checkAnotherBasis
                ( this : access mvc;
                  repIdx : integer32; dist : integer32;
                  target : demics_ftest.class_theData.Link_to_Array_of_theData
                ) return integer32 is
    begin
      return 0;
    end checkAnotherBasis;

    procedure get_firIdx
                ( this : access mvc;
                  data_a : demics_ftest.class_ftData.ftData;
                  data_b : demics_ftest.class_ftData.ftData;
                  sn : in integer32; lvl : in integer32 ) is
    begin
      null;
    end get_firIdx;

    procedure info_cpuTime
                ( this : access mvc;
                  cpuTime_start : double_float;
                  cpuTime_end : double_float ) is
    begin
      null;
    end info_cpuTime;

    procedure info_final ( this : access mvc ) is
    begin
      null;
    end info_final;

    procedure enum ( this : access mvc ) is
    begin
      null;
    end enum;

  end class_mvc;

end demics_mvc;
