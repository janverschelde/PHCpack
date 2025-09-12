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
      the_Reltab : demics_reltab.class_reltab.Link_to_reltab;
      the_Simplex : demics_simplex.class_simplex.Link_to_simplex;
      lv : demics_ftest.class_lvData.Link_to_Array_of_lvData;
      iLv : demics_itest.class_iLvData.Link_to_Array_of_iLvData;
    end record;

    type Link_to_mvc is access mvc;

    procedure getMemory
                ( this : in Link_to_mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32; vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates memory for the levels.

    procedure initMemoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32 );

    -- DESCRIPTION :
    --   Adds an element to data if data.cur is null.

    procedure memoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32 );

    -- DESCRIPTION :
    --   Same as initMemoryCheck.

    procedure get_candIdx
                ( this : in Link_to_mvc;
                  curInif : in demics_iTest.class_inifData.Link_to_inifData );

    -- DESCRIPTION :
    --   Defines the elements in this.candIdx.

    function chooseSup
                ( this : Link_to_mvc;
                  depth : integer32;
                  curNode : demics_fTest.class_theData.Link_to_theData;
             curInif : demics_iTest.class_inifData.Link_to_Array_of_inifData;
             nextInif : demics_iTest.class_inifData.Link_to_Array_of_inifData;
                  vrblvl : integer32 := 0 )
                return integer32;

    -- DESCRIPTION :
    --   Invoked in the enum procedure.

    -- NOTE :
    --   Both curInif and nextInif stand out in the specification,
    --   as their definition as inifData* was not obvious from the
    --   original prototype declaration of chooseSup.

    procedure fUpdateDirRed
                ( this : in Link_to_mvc;
                  curInif : in demics_iTest.class_inifData.Array_of_inifData;
                  nextInif : in demics_iTest.class_inifData.Array_of_inifData;
                  curNode : in demics_fTest.Class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 );

    procedure updateDirRed
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Array_of_inifData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData;
                  curNode : in demics_ftest.class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32 );

    function findUnbDir
                ( this : Link_to_mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : demics_ftest.class_theData.Link_to_theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32;

    function findUnbDir_art
                ( this : Link_to_mvc;
                  nextInif : demics_itest.class_inifData.Array_of_inifData;
                  curNode : demics_ftest.class_theData.Link_to_theData;
                  nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                  curRsp : Standard_Integer_Vectors.Link_to_Vector;
                  depth : integer32 ) return integer32;

    function checkDir
                ( this : Link_to_mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : demics_itest.class_uData.Link_to_uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32;

    function checkDir_art
                ( this : Link_to_mvc;
                  corPtr : demics_itest.class_uData.Link_to_Array_of_uData;
                  tarPtr : demics_itest.class_uData.Link_to_uData;
                  tar_dir : Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : double_float;
                  nf_pos : Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : Standard_Integer_Vectors.Link_to_Vector;
                  nfN : integer32 ) return integer32;

    procedure skipPtr
                ( this : in Link_to_mvc;
                  curr : in demics_itest.class_uData.Link_to_Array_of_uData;
                  fHead : in demics_itest.class_uData.Link_to_Array_of_uData );

    procedure get_tuple_index
                ( this : in Link_to_mvc;
                  node : in demics_ftest.class_ftData.Link_to_ftData;
                  data : in demics_ftest.class_ftData.Array_of_ftData;
                  length : in integer32 );

    -- DESCRIPTION :
    --   Sets the value of nodeLabel of node.parent using data.

    procedure dbg_init_transMat
                ( this : in Link_to_mvc;
                  curNode : in demics_ftest.class_theData.Link_to_theData );

    procedure dbg_transMat
                ( this : in Link_to_mvc;
                  preNode : in demics_ftest.class_theData.Link_to_theData;
                  curNode : in demics_ftest.class_theData.Link_to_theData );

    procedure check_transMat
                ( this : in Link_to_mvc;
                  preNode : in demics_ftest.class_theData.Link_to_theData;
                  curNode : in demics_ftest.class_theData.Link_to_theData );

    procedure check_init_transRed
                ( this : in Link_to_mvc;
                  curNode : in demics_ftest.class_theData.Link_to_theData );

    function checkSign_red
                ( this : Link_to_mvc;
                  curRed : double_float;
                  tarRed : double_float ) return integer32;

    -- DESCRIPTION :
    --   Checks the sign of curRed versus tarRed.

    function checkNonNeg_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32;

    function checkNonPos_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float ) return integer32;

    function checkZero_dir
                ( this : Link_to_mvc;
                  curDirElem : double_float;
                  tarDirElem : double_float) return integer32;

    function table_out
                ( this : Link_to_mvc;
                  row : integer32; col : integer32 ) return integer32;

    -- DESCRIPTION :
    --   Returns the element in the table at position defined by row and col.

    procedure info_neg
                ( this : in Link_to_mvc; termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec );

    -- DESCRIPTION :
    --   Writes the numbers in this.trNeg and in negIdx.

    procedure info_sp ( this : in Link_to_mvc; depth : in integer32 );

    -- DESCRIPTION :
    --   Writes the numbers in this.sp.

    procedure info_parent_node ( this : in Link_to_mvc; depth : in integer32 );

    -- DESCRIPTION :
    --   Calls info_parent_node on the elements in this.lv.

    procedure info_tuple
                ( this : in Link_to_mvc;
                  lvl : in integer32 ); -- depth : in integer32 );

    -- DESCRIPTION :
    --   Writes the numbers in this.mFeaIdx.
    --   In the original procedure depth was not referenced.

    procedure info_all_dirRed
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  node : in demics_ftest.class_ftData.Link_to_ftData;
                  nextInif : demics_itest.class_inifData.Array_of_inifData );

    procedure info_mFea ( this : in Link_to_mvc; length : in integer32 );

    -- DESCRIPTION :
    --   Writes the numbers in this.mFea and this.mRepN.

    procedure info_firIdx ( this : in Link_to_mvc; length : in integer32 );

    -- DESCRIPTION :
    --   Writes the numbers in this.firIdx.

    procedure info_fIdx
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData );

   -- DESCRIPTION :
   --   Writes the first index in data.parent.

    procedure info_candIdx ( this : in Link_to_mvc );

    -- DESCRIPTION :
    --   Writes the numbers in this.candIdx.

    procedure info_elemNum
                ( this : in Link_to_mvc;
                  length : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  node : in demics_ftest.class_ftData.ftData );

    procedure info_prop_elemNum
                ( this : in Link_to_mvc;
                  length : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  node : in demics_ftest.class_ftData.ftData );

    procedure info_table ( this : in Link_to_mvc );

    -- DESCRIPTION :
    --   Writes information about the relation table.

    function new_mvc return mvc;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_mvc ( this : in Link_to_mvc );

    -- DESCRIPTION :
    --   Deallocates memory for the vector types.

    procedure allocateAndIni
                ( this : in Link_to_mvc;
                  data : in demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32;
                  vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Allocates and initializes the data.

    procedure initFeasTest
                ( this : in Link_to_mvc; depth : in integer32;
                  vrblvl : in integer32 := 0 );

    procedure initCheck
                ( this : in Link_to_mvc; depth : in integer32;
                  data : in demics_ftest.class_ftData.Link_to_ftData );

    procedure initLP
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  depth : in integer32; idx : in integer32;
                  feaNum : in out integer32 );

    function feasTest
                ( this : Link_to_mvc; depth : integer32;
                  parent : demics_ftest.class_theData.Link_to_theData;
                  vrblvl : integer32 := 0 )
                return integer32;

    procedure upFeasTest
                ( this : in Link_to_mvc; depth : in out integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 );

    -- upFeasTest was declared as a function, assigning to depth
    -- as a side effect, returning a flag

    procedure findMixedCell
                ( this : in Link_to_mvc; depth : in integer32;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  vrblvl : in integer32 := 0 );

    procedure findAllMixedCells
               ( this : in Link_to_mvc; depth : in integer32;
                 vrblvl : in integer32 := 0 );

    function iCheck
                ( this : Link_to_mvc;
                  depth : integer32;
                  parent : demics_ftest.class_theData.Link_to_theData;
                  data : demics_ftest.class_ftData.Link_to_ftData;
                  inifData : demics_itest.class_inifData.Link_to_inifData )
                return integer32;

    procedure iLP ( this : in Link_to_mvc;
                    parent : in demics_ftest.class_theData.Link_to_theData;
                    data : in demics_ftest.class_ftData.Link_to_ftData;
                    depth : in integer32;
                    idx_one : in integer32;
                    fst_pivInIdx : in integer32;
                    sub_fst_pivInIdx : in integer32;
                    preNbN : in integer32;
                    feaNum : in out integer32 );

    procedure iLP_Art
                ( this : in Link_to_mvc;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  depth : in integer32;
                  idx_one : in integer32;
                  fst_pivInIdx : in integer32;
                  sub_fst_pivInIdx : in integer32;
                  preNbN : in integer32;
                  feaNum : in out integer32 );

    procedure findNode
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  data : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  flag : out integer32 );

    -- findNode was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag; and
    -- data was declared as ftData* Data, not obvious that the type had
    -- to be a pointer to an array ...

    procedure findNextNode
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  lvl : in out integer32;
                  feaNum : in out integer32;
                  data : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  flag : out integer32 );

    -- findNode was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag; and
    -- data was declared as ftData* Data, not obvious that the type had
    -- to be a pointer to an array ...

    procedure findUpNode
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  pre : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  cur : in demics_ftest.class_ftData.Link_to_Array_of_ftData;
                  lvl : in out integer32;
                  polyDim : in integer32;
                  depth : in integer32 );

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
                    length : in integer32; flag : out integer32 );

    -- mLP was declared as a function, but with assignments to
    -- its arguments as side effects, returning a flag

    function checkBasis
                ( this : Link_to_mvc;
                  target : demics_ftest.class_theData.Link_to_theData;
                  sub_sIdx : integer32 ) return integer32;

    function checkAnotherBasis
                ( this : Link_to_mvc;
                  repIdx : integer32; dist : integer32;
                  target : demics_ftest.class_theData.Link_to_Array_of_theData
                ) return integer32;

    procedure get_firIdx
                ( this : in Link_to_mvc;
                  data_a : in demics_ftest.class_ftData.ftData;
                  data_b : in demics_ftest.class_ftData.ftData;
                  sn : in integer32; lvl : in integer32 );

    -- DESCRIPTION :
    --   Sets the value of this.firIdx(sn).

    procedure info_cpuTime
                ( this : in Link_to_mvc;
                  cpuTime_start : in double_float;
                  cpuTime_end : in double_float );

    procedure info_final ( this : in Link_to_mvc );

    procedure enum ( this : in Link_to_mvc; vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Enumerates all mixed cells.

  end class_mvc;

end demics_mvc;
