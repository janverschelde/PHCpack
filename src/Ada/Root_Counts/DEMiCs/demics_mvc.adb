with Ada.text_io;                       use Ada.text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with DEMiCs_Global_Constants;

package body demics_mvc is

  package body class_mvc is

    procedure getMemory
                ( this : in Link_to_mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32; vrblvl : in integer32 := 0 ) is

      elemLen,div : integer32;

    begin
      if vrblvl > 0
       then put("-> in demics_mvc.getMemory, ");
      end if;
      if depth /= 0
       then div := 2;
       else div := 1;
      end if;
      if lvl /= length-1 
       then elemLen := this.termSet(depth)/((lvl + 1)*div);
       else elemLen := 1;
      end if;
      if vrblvl > 0 then
        put("depth : "); put(depth,1);
        put(", lvl : "); put(lvl,1);
        put(", elemLen : "); put(elemLen,1); new_line;
      end if;
      for i in 0..this.termSet(depth)-1 loop
        this.lv(depth).ftest(lvl)
          := new demics_fTest.class_ftData.ftData'
                (demics_fTest.class_ftData.new_ftData);
        demics_fTest.class_ftData.create_elem
          (this.lv(depth).ftest(lvl),this.row,this.col,
           this.termSet(depth),this.supType(depth),vrblvl-1);
        demics_fTest.class_ftData.add_elem(this.lv(depth).ftest(lvl),vrblvl-1);
      end loop;
      demics_fTest.class_ftData.mark(this.lv(depth).ftest(lvl),vrblvl-1);
      this.lv(depth).ftest(lvl).cur := this.lv(depth).ftest(lvl).head;
    end getMemory;

    procedure initMemoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is

      sn : constant integer32 := this.sp(depth);

      use demics_fTest.class_theData;

    begin
      if data.cur = null then
        demics_ftest.class_ftData.create_elem
          (data,this.row,this.col,this.termSet(sn),this.supType(sn));
        demics_ftest.class_ftData.add_elem(data);
      end if;
    end initMemoryCheck;

    procedure memoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is
    begin
      initMemoryCheck(this,data,depth); -- same code as initMemoryCheck
    end memoryCheck;

    procedure get_candIdx
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Link_to_inifData
                ) is

      num : integer32 := 0;
      n_curr : demics_iTest.class_uData.Link_to_uData := curInif.fHead;

      use demics_iTest.class_uData;

    begin
      while n_curr /= null loop
        this.candIdx(num+1) := n_curr.supLab;
        n_curr := n_curr.fNext;
        num := num + 1;
      end loop;
      this.candIdx(0) := num;
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
                  data : in demics_ftest.class_ftData.Array_of_ftData;
                  length : in integer32 ) is

      on : constant := DEMiCs_Global_Constants.ON;
      tmpIdx : integer32;

    begin
      for i in 0..length-1 loop
        node.parent.nodeLabel(i) := data(length-1).cur.fIdx;
      end loop;
      node.parent.nodeLabel(length-1) := data(length-1).cur.fIdx;
      if data(1).parent.sw = on then
        tmpIdx := node.parent.nodeLabel(1);
        node.parent.nodeLabel(1) := node.parent.nodeLabel(0);
        node.parent.nodeLabel(0) := tmpIdx;
      end if;
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

      flag : integer32 := 0;
      pluszero : constant := DEMiCs_Global_Constants.PLUSZERO;
      minuszero : constant := DEMiCs_Global_Constants.MINUSZERO;
      positive : constant := DEMiCs_Global_Constants.POSITIVE;
      negative : constant := DEMiCs_Global_Constants.NEGATIVE;

    begin
      if curRed > tarRed + pluszero then
        flag := positive;
      elsif curRed < tarRed + minuszero then
        flag := negative;
      end if;
      return flag;
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
      return this.table(row + col*this.termSumNum);
    end table_out;

    procedure info_neg
                ( this : in Link_to_mvc; termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec ) is
    begin
      put_line("<< trNeg >>");
      for j in 0..termSet-1 loop
        for i in 0..this.row-1 loop
          put(this.trNeg(j)(i),1); put(" ");
        end loop;
        new_line;
      end loop;
      put_line("<< negIdx >>");
      for j in 0..termSet-1 loop
        for i in 0..negIdx(j)(0) loop
          put(negIdx(j)(i+1),1); put(" ");
        end loop;
        new_line;
      end loop;
    end info_neg;

    procedure info_sp ( this : in Link_to_mvc; depth : in integer32 ) is
    begin
      put("sp :");
      for i in 0..depth-1 loop
        put(" "); put(this.sp(i),1);
      end loop;
      new_line;
    end info_sp;

    procedure info_parent_node ( this : in Link_to_mvc;
                                 depth : in integer32 ) is
    begin
      put("Node : ");
      for i in 0..depth-1 loop
        put(this.sp(i),1); put(" : ");
        demics_ftest.class_ftData.info_parent_node(this.lv(this.sp(i)).node);
      end loop;
      new_line;
    end info_parent_node;

    procedure info_tuple
                ( this : in Link_to_mvc;
                  lvl : in integer32 ) is -- ; depth : in integer32 ) is
    begin
      put("( ");
      for i in 0..lvl-1 loop
        put(this.mFeaIdx(i)(this.mRepN(i))+1,1); put(" ");
      end loop;
      put_line(")");
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
      put("mFea :");
      for i in 0..length-1 loop
        put(" "); put(this.mFea(i),1);
      end loop;
      new_line;
      put("mRepN :");
      for i in 0..length-1 loop
        put(" "); put(this.mRepN(i),1);
      end loop;
      new_line;
    end info_mFea;

    procedure info_firIdx ( this : in Link_to_mvc; length : in integer32 ) is
    begin
      put_line("<< firIdx >>");
      for i in 0..length loop
        put(this.firIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_firIdx;

    procedure info_fIdx
                ( this : in Link_to_mvc;
                  data : in demics_ftest.class_ftData.Link_to_ftData ) is
    begin
      put("First Index : "); put(data.parent.fIdx+1,1); new_line;
    end info_fIdx;

    procedure info_candIdx ( this : in Link_to_mvc ) is
    begin
      put("candIdx :");
      for i in 0..this.candIdx(0)-1 loop
        put(this.candIdx(i+1),1); put(" ");
      end loop;
      new_line;
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

    -- NOTE :
    --   This is essentially the same code as the procedure
    --   demics_reltab.class_reltab.info_table.

      u_cnt : integer32 := 0;
      t_cnt : integer32 := 0;

    begin
      put_line("<< Relation table >>");
      for i in 0..this.termSumNum-1 loop
        for k in 0..i-1 loop
          put("  ");
        end loop;
        for j in i+1..this.termSumNum-1 loop
          if table_out(this,j,i) = DEMiCs_Global_Constants.UNBOUNDED
           then u_cnt := u_cnt + 1;
          end if;
          put(table_out(this,j,i),1); put(" ");
          t_cnt := t_cnt + 1;
        end loop;
        new_line;
      end loop;
      new_line;
      put("# Unb. LPs : "); put(u_cnt,1); new_line;
      put("# Elem. : "); put(t_cnt,1); new_line;
      put("Ratio : "); put(double_float(u_cnt)/double_float(t_cnt));
      new_line;
    end info_table;

    function new_mvc return mvc is
 
      res : mvc;

    begin
      res.dim := 0;
      res.supN := 0;
      res.row := 0;
      res.col := 0;
      res.termSumNum := 0;
      res.termMax := 0;
      res.maxLength := 0;
      res.total_iter := 0.0;
      res.total_feasLP := 0.0;
      res.total_LPs := 0.0;
      res.total_1PT := 0.0;
      res.total_2PT := 0.0;
      res.total_triLPs_mLP := 0.0;
      res.total_unbLP_tab := 0.0;
      res.lvl_1PT := null;
      res.lvl_2PT := null;
      res.actNode := null;
      res.mfNum := null;
      res.termSet := null;
      res.termStart := null;
      res.re_termStart := null;
      res.supType := null;
      res.mRepN := null;
      res.mFeaIdx := null;
      res.mFea := null;
      res.trNeg := null;
      res.firIdx := null;
      res.repN := null;
      res.sp := null;
      res.candIdx := null;
      res.trMat := null;
      res.table := null;
      res.lv := null;
      res.iLv := null;
      return res;
    end new_mvc;

    procedure delete_mvc ( this : in Link_to_mvc ) is

      use Standard_Integer_VecVecs;

    begin
      if this.trNeg /= null then
        for i in 0..this.termSet(this.sp(0))-1 loop
          Standard_Integer_Vectors.clear(this.trNeg(i));
        end loop;
        Standard_Integer_VecVecs.shallow_clear(this.trNeg);
        this.trNeg := null;
      end if;
      Standard_Integer_Vectors.clear(this.re_termStart);
      Standard_Integer_Vectors.clear(this.mfNum);
      Standard_Floating_Vectors.clear(this.lvl_1PT);
      Standard_Floating_Vectors.clear(this.lvl_2PT);
      Standard_Floating_Vectors.clear(this.actNode);
      Standard_Integer_Vectors.clear(this.firIdx);
      Standard_Integer_Vectors.clear(this.repN);
      Standard_Integer_Vectors.clear(this.sp);
      Standard_Integer_Vectors.clear(this.candIdx);
      Standard_Floating_Vectors.clear(this.trMat);
     -- delete [] this.lv;
     -- delete [] this.iLv;
    end delete_mvc;

    procedure allocateAndIni
                ( this : in Link_to_mvc;
                  data : in demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32;
                  vrblvl : in integer32 := 0 ) is

      length : integer32;

    begin
      if vrblvl > 0
       then put_line("-> in demics_mvc.allocateAndIni ...");
      end if;
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
      this.lvl_1PT
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.lvl_2PT
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.actnode
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.firIdx
        := new Standard_Integer_Vectors.Vector'(0..this.supN => 0);
      this.re_termStart := new Standard_Integer_Vectors.Vector(0..this.supN);
      this.re_termStart(0) := 0;
      this.repN
        := new Standard_Integer_Vectors.Vector'(0..this.supN-1 => 0);
      this.sp := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.candIdx := new Standard_Integer_Vectors.Vector(0..this.termMax);
      this.trMat
        := new Standard_Floating_Vectors.Vector(0..this.dim*this.dim-1);
      this.lv := new demics_ftest.class_lvData.Array_of_lvData(0..this.supN-1);
      this.iLv
        := new demics_itest.class_iLvData.Array_of_iLvData(0..this.supN-1);
      for i in 0..this.supN-1 loop
        this.re_termStart(i+1) := this.termStart(i+1) - i - 1;
        this.sp(i) := i;
        this.lv(i) := new demics_ftest.class_lvData.lvData'
                         (demics_ftest.class_lvData.new_lvData);
        demics_ftest.class_lvData.create
          (this.lv(i),i,--this.supN,this.dim,
           this.supType(i)+1,this.termMax,vrblvl-1);
        this.iLv(i) := new demics_itest.class_iLvData.ilvData'
                          (demics_itest.class_iLvData.new_iLvData);
        demics_itest.class_iLvData.create
          (this.iLv(i),i,this.supN,this.dim,this.termMax,vrblvl-1);
        length := this.supType(i) + 1;
        for j in 0..length-1 loop
          getMemory(this,i,j,length,vrblvl-1);
        end loop;
      end loop;
      this.the_Simplex := new demics_simplex.class_simplex.simplex'
                             (demics_simplex.class_simplex.new_simplex);
      demics_simplex.class_simplex.allocateAndIni
        (this.the_Simplex,data,this.firIdx,seedNum,output,vrblvl-1);
      this.the_Reltab := new demics_reltab.class_reltab.reltab'
                            (demics_reltab.class_reltab.new_reltab);
      demics_reltab.class_reltab.allocateAndIni
        (this.the_Reltab,this.the_Simplex,this.firIdx,this.dim,this.supN,
         this.termSumNum,this.termSet,this.termStart,this.re_termStart,
         vrblvl-1);
      demics_itest.class_iLvData.getInit
        (this.iLv(0),data,this.the_Simplex.lifting,this.termSet,
         this.termStart,this.dim,this.supN,vrblvl-1);
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

      off : constant := DEMiCs_Global_Constants.OFF;

    begin
      if lvl = 1 then
        this.firIdx(sn) := data_a.parent.fIdx;
      else
        if data_b.parent.sw = off
         then this.firIdx(sn) := data_a.parent.fIdx;
         else this.firIdx(sn) := data_b.parent.fIdx;
        end if;
      end if;
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

    procedure enum ( this : in Link_to_mvc; vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 
       then put_line("-> in demics_mvc.class_mvc.enum ...");
      end if;
      demics_reltab.class_reltab.makeTable
        (this.the_Reltab,this.total_unbLP_tab,vrblvl-1);
    end enum;

  end class_mvc;

end demics_mvc;
