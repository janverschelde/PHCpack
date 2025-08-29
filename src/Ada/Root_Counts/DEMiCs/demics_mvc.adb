package body demics_mvc is

  package body class_mvc is

    procedure getMemory
                ( this : in Link_to_mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32 ) is
    begin
      null;
    end getMemory;

    procedure initMemoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is
    begin
      null;
    end initMemoryCheck;

    procedure memoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is
    begin
      null;
    end memoryCheck;

    procedure get_candIdx
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Link_to_inifData
                ) is
    begin
      null;
    end get_candIdx;

    function chooseSup
                ( this : Link_to_mvc; depth : integer32;
                  curNode : demics_ftest.class_theData.Link_to_theData;
                  curInif : demics_itest.class_inifData.Link_to_inifData;
                  nextInif : demics_itest.class_inifData.Link_to_inifData )
                return integer32 is
    begin
      return 0;
    end chooseSup;

    procedure fUpdateDirRed
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Array_of_inifData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData;
                  curNode : in demics_ftest.Class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) is
    begin
      null;
    end fUpdateDirRed;

    procedure updateDirRed
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Array_of_inifData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData;
                  curNode : in demics_ftest.class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 ) is
    begin
      null;
    end updateDirRed;

    function findUnbDir
                ( this : Link_to_mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : demics_ftest.class_theData.Link_to_theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32 is
    begin
      return 0;
    end findUnbDir;

    function findUnbDir_art
                ( this : Link_to_mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : demics_ftest.class_theData.Link_to_theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32 is
    begin
      return 0;
    end findUnbDir_art;

    function checkDir
                ( this : Link_to_mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : demics_itest.class_uData.Link_to_uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32 is
    begin
      return 0;
    end checkDir;

    function checkDir_art
                ( this : Link_to_mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : demics_itest.class_uData.Link_to_uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32 is
    begin
      return 0;
    end checkDir_art;

    procedure skipPtr
                ( this : in Link_to_mvc;
                  curr : in demics_itest.class_uData.Link_to_Array_of_uData;
                  fHead : in demics_itest.class_uData.Link_to_Array_of_uData
                ) is
    begin
      null;
    end skipPtr;

    procedure get_tuple_index
                ( this : in Link_to_mvc;
                  node : in demics_ftest.class_ftData.Link_to_ftData;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  length : in integer32 ) is
    begin
      null;
    end get_tuple_index;

    procedure dbg_init_transMat
                ( this : in Link_to_mvc;
                  curNode : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end dbg_init_transMat;

    procedure dbg_transMat
                ( this : in Link_to_mvc;
                  preNode : in demics_ftest.class_theData.Link_to_theData;
                  curNode : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end dbg_transMat;

    procedure check_transMat
                ( this : in Link_to_mvc;
                  preNode : in demics_ftest.class_theData.Link_to_theData;
                  curNode : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end check_transMat;

    procedure check_init_transRed
                ( this : in Link_to_mvc;
                  curNode : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end check_init_transRed;

    function checkSign_red
                ( this : Link_to_mvc;
                  curRed : double_float;
                  tarRed : double_float ) return integer32 is
    begin
      return 0;
    end checkSign_red;

    function checkNonNeg_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkNonNeg_dir;

    function checkNonPos_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkNonPos_dir;

    function checkZero_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      return 0;
    end checkZero_dir;

    function table_out
                ( this : in Link_to_mvc;
                  row : integer32; col : integer32 ) return integer32 is
    begin
      return 0;
    end table_out;

    procedure info_neg
                ( this : in Link_to_mvc; termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec ) is
    begin
      null;
    end info_neg;

    procedure info_sp ( this : in Link_to_mvc; depth : in integer32 ) is
    begin
      null;
    end info_sp;

    procedure info_parent_node ( this : in Link_to_mvc;
                                 depth : in integer32 ) is
    begin
      null;
    end info_parent_node;

    procedure info_tuple
                ( this : in Link_to_mvc;
                  lvl : in integer32; depth : in integer32 ) is
    begin
      null;
    end info_tuple;

    procedure info_all_dirRed
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  node : in demics_ftest.class_ftData.Link_to_ftData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData
                ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_mFea ( this : in Link_to_mvc; length : in integer32 ) is
    begin
      null;
    end info_mFea;

    procedure info_firIdx ( this : in Link_to_mvc; length : in integer32 ) is
    begin
      null;
    end info_firIdx;

    procedure info_fIdx
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData ) is
    begin
      null;
    end info_fIdx;

    procedure info_candIdx ( this : in Link_to_mvc ) is
    begin
      null;
    end info_candIdx;

    procedure info_elemNum
                ( this : in Link_to_mvc;
                  length : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  node : in demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end info_elemNum;

    procedure info_prop_elemNum
                ( this : in Link_to_mvc;
                  length : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  node : in demics_ftest.class_ftData.ftData ) is
    begin
      null;
    end info_prop_elemNum;

    procedure info_table ( this : in Link_to_mvc ) is
    begin
      null;
    end info_table;

    function new_mvc return mvc is
 
      res : mvc;

    begin
      return res;
    end new_mvc;

    procedure delete_mvc ( this : in Link_to_mvc ) is
    begin
      null;
    end delete_mvc;

    procedure allocateAndIni
                ( this : in Link_to_mvc;
                  data : in demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32 ) is
    begin
      this.dim := data.dim;
      this.supN := data.supN;
      this.row := data.dim;
      this.termSumNum := data.termSumNum;
      this.termMax := data.termMax;
      this.maxLength := data.typeMax + 1;
      this.col := this.termSumNum - this.supN + this.dim;
      this.termSet := data.termSet;
      this.termStart := data.termStart;
      this.supType := data.supType;
      this.mfNum := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.lvl_1PT := new Standard_Floating_Vectors.Vector(0..this.supN-1);
      this.lvl_2PT := new Standard_Floating_Vectors.Vector(0..this.supN-1);
      this.actnode := new Standard_Floating_Vectors.Vector(0..this.supN-1);
      this.firIdx := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.re_termStart := new Standard_Integer_Vectors.Vector(0..this.supN);
      this.repN := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.sp := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.candIdx := new Standard_Integer_Vectors.Vector(0..this.termMax);
      this.trMat
        := new Standard_Floating_Vectors.Vector(0..this.dim*this.dim-1);
      this.lv := new demics_ftest.class_lvData.Array_of_lvData(0..this.supN-1);
      this.iLv
        := new demics_itest.class_iLvData.Array_of_iLvData(0..this.supN-1);
    end allocateAndIni;

    procedure initFeasTest ( this : in Link_to_mvc; depth : in integer32 ) is
    begin
      null;
    end initFeasTest;

    procedure initCheck
                ( this : in Link_to_mvc; depth : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData ) is
    begin
      null;
    end initCheck;

    procedure initLP
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  depth : in integer32; idx : in integer32;
                  feaNum : in out integer32 ) is
    begin
      null;
    end initLP;

    function feasTest
                ( this : Link_to_mvc; depth : integer32;
                  parent : demics_ftest.class_theData.Link_to_theData )
                return integer32 is
    begin
      return 0;
    end feasTest;

    procedure upFeasTest
                ( this : in Link_to_mvc; depth : in out integer32;
                  flag : out integer32 ) is
    begin
      null;
    end upFeasTest;

    procedure findMixedCell
                ( this : in Link_to_mvc; depth : in integer32;
                  parent : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      null;
    end findMixedCell;

    procedure findAllMixedCells ( this : in Link_to_mvc;
                                  depth : in integer32 ) is
    begin
      null;
    end findAllMixedCells;

    function iCheck
                ( this : Link_to_mvc; depth : integer32;
                  parent : demics_ftest.class_theData.Link_to_theData;
                  data : demics_ftest.class_ftData.Link_to_ftData;
                  inifData : demics_itest.class_inifData.Link_to_inifData )
                return integer32 is
    begin
      return 0;
    end iCheck;

    procedure iLP ( this : in Link_to_mvc;
                    parent : in demics_ftest.class_theData.Link_to_theData;
                    data : in demics_ftest.class_ftData.Link_to_ftData;
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
                ( this : in Link_to_mvc;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
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
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  flag : out integer32 ) is
    begin
      null;
    end findNode;

    procedure findNextNode
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  Data : in demics_ftest.class_ftData.Link_to_ftData;
                  flag : out integer32 ) is
    begin
      null;
    end findNextNode;

    procedure findUpNode
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  pre : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  cur : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  lvl : in out integer32;
                  polyDim : in integer32;
                  depth : in integer32 ) is
    begin
      null;
    end findUpNode;

    procedure mLP ( this : in Link_to_mvc;
                    pre : in demics_ftest.class_ftData.Link_to_ftData;
                    cur : in demics_ftest.class_ftData.Link_to_ftData;
                    data : in demics_ftest.class_ftData.Link_to_ftData;
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
                ( this : Link_to_mvc;
                  target : demics_ftest.class_theData.Link_to_theData;
                  sub_sIdx : integer32 ) return integer32 is
    begin
      return 0;
    end checkBasis;

    function checkAnotherBasis
                ( this : Link_to_mvc;
                  repIdx : integer32; dist : integer32;
                  target : demics_ftest.class_theData.Link_to_Array_of_theData
                ) return integer32 is
    begin
      return 0;
    end checkAnotherBasis;

    procedure get_firIdx
                ( this : in Link_to_mvc;
                  data_a : in demics_ftest.class_ftData.ftData;
                  data_b : in demics_ftest.class_ftData.ftData;
                  sn : in integer32; lvl : in integer32 ) is
    begin
      null;
    end get_firIdx;

    procedure info_cpuTime
                ( this : in Link_to_mvc;
                  cpuTime_start : in double_float;
                  cpuTime_end : in double_float ) is
    begin
      null;
    end info_cpuTime;

    procedure info_final ( this : in Link_to_mvc ) is
    begin
      null;
    end info_final;

    procedure enum ( this : in Link_to_mvc ) is
    begin
      null;
    end enum;

  end class_mvc;

end demics_mvc;
