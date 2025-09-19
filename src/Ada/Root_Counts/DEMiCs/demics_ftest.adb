with unchecked_deallocation;
with Ada.text_io;                       use Ada.text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with DEMiCs_Global_Constants;

package body demics_ftest is

  package body class_theData is

    function new_theData return theData is

      res : theData;

    begin
      res.row := 0;
      res.col := 0;
      res.termS := 0;
      res.flag := 0;
      res.polyDim := 0;
      res.nbN := 0;
      res.nfN := 0;
      res.artV := 0;
      res.pivOutNum := 0;
      res.fIdx := 0;
      res.sw := 0;
      res.invB := null;
      res.transMat := null;
      res.transRed := null;
      res.p_sol := null;
      res.d_sol := null;
      res.redVec := null;
      res.basisIdx := null;
      res.nbIdx := null;
      res.nf_pos := null;
      res.rIdx := null;
      res.next := null;
      res.pivOutList := null;
      res.pivOutCheck := null;
      res.invB_ptr := null;  
      res.transMat_ptr := null;
      res.transRed_ptr := null;
      res.p_sol_ptr := null; 
      res.d_sol_ptr := null; 
      res.redVec_ptr := null;
      res.basisIdx_ptr := null; 
      res.nbIdx_ptr := null; 
      res.nf_pos_ptr := null;
      res.nodeLabel := null;
      return res;
    end new_theData;

    procedure delete_theData ( this : in Link_to_theData ) is
    begin
      Standard_Floating_Vectors.clear(this.invB);
      Standard_Floating_Vectors.clear(this.transMat);
      Standard_Floating_Vectors.clear(this.transRed);
      Standard_Floating_Vectors.clear(this.p_sol);
      Standard_Floating_Vectors.clear(this.d_sol);
      Standard_Floating_Vectors.clear(this.redVec);
      Standard_Integer_Vectors.clear(this.basisIdx);
      Standard_Integer_Vectors.clear(this.nbIdx);
      Standard_Integer_Vectors.clear(this.nf_pos);
      Standard_Integer_Vectors.clear(this.rIdx);
      Standard_Integer_Vectors.clear(this.pivOutList);
      Standard_Integer_Vectors.clear(this.pivOutCheck);
      Standard_Integer_Vectors.clear(this.nodeLabel);
    end delete_theData;

    procedure create ( this : in Link_to_theData;
                       ori_row : in integer32; ori_col : in integer32;
                       ori_termS : in integer32;
                       ori_polyDim : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_ftest.class_theData.create, ori_row : ");
        put(ori_row,1); put_line(" ...");
      end if;
      this.row := ori_row;
      this.col := ori_col;
      this.termS := ori_termS;
      this.polyDim := ori_polyDim;
      this.invB
        := new Standard_Floating_Vectors.Vector'(0..this.row*this.row-1 => 0.0);
      this.transMat
        := new Standard_Floating_Vectors.Vector'(0..this.row*this.row-1 => 0.0);
      this.transRed
        := new Standard_Floating_Vectors.Vector'(0..this.row-1 => 0.0);
      this.p_sol
        := new Standard_Floating_Vectors.Vector'(0..this.col-1 => 0.0);
      this.d_sol
        := new Standard_Floating_Vectors.Vector'(0..this.row-1 => 0.0);
      this.basisIdx
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.nf_pos := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.nbIdx := new Standard_Integer_Vectors.Vector'(0..this.col-1 => 0);
      this.redVec
        := new Standard_Floating_Vectors.Vector'(0..this.col-1 => 0.0);
      this.rIdx
        := new Standard_Integer_Vectors.Vector'(0..this.termS-1 => 0);
      this.pivOutList
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.pivOutCheck
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.nodeLabel
        := new Standard_Integer_Vectors.Vector(0..this.polyDim);
    end create;

    procedure joint ( this : in Link_to_theData ) is
    begin
      this.invB_ptr := this.invB;
      this.transMat_ptr := this.transMat;
      this.p_sol_ptr := this.p_sol;
      this.d_sol_ptr := this.d_sol;
      this.redVec_ptr := this.redVec;
      this.basisIdx_ptr := this.basisIdx;
      this.nbIdx_ptr := this.nbIdx;
      this.nf_pos_ptr := this.nf_pos;
    end joint;

    procedure iJoint ( this : in Link_to_theData ) is
    begin
      this.transMat_ptr := this.transMat;
      this.transRed_ptr := this.transRed;
      this.redVec_ptr := this.redVec;
      this.nbIdx_ptr := this.nbIdx;
    end iJoint;

    procedure mJoint ( this : in Link_to_theData ) is
    begin
      this.nf_pos_ptr := this.nf_pos;
    end mJoint;

    procedure clear ( this : in Link_to_theData ) is
    begin
      this.nbN := 0;
      this.nfN := 0;
      this.artV := 0;
      for i in 0..this.row*this.row-1 loop
        this.invB(i) := 0.0;
        this.transMat(i) := 0.0;
      end loop;
      for i in 0..this.col-1 loop
        this.p_sol(i) := 0.0;
        this.nbIdx(i) := 0;
        this.redVec(i) := 0.0;
      end loop;
      for i in 0..this.row-1 loop
        this.d_sol(i) := 0.0;
        this.basisIdx(i) := 0;
        this.nf_pos(i) := 0;
      end loop;
      for i in 0..this.termS-1 loop
        this.rIdx(i) := 0;
      end loop;
    end clear;

    procedure clear_transMat ( this : in Link_to_theData ) is
    begin
      for i in 0..this.row*this.row-1 loop
        this.transMat(i) := 0.0;
      end loop;
    end clear_transMat;

    procedure put_info ( this : in Link_to_theData;
                         repIdx : in integer32; idx2 : out integer32;
                         lNbN : out integer32; lNfN : out integer32 ) is
    begin
      idx2 := this.rIdx(repIdx);
      lNbN := this.nbN;
      lNfN := this.nfN;
    end put_info;

    function invB_out ( this : Link_to_theData;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return this.invB(colIdx + this.row*rowIdx);
    end invB_out;

    function transMat_out ( this : Link_to_theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float is
    begin
      return this.transMat(colIdx + this.row*rowIdx);
    end transMat_out;

    function invB_ptr_out ( this : Link_to_theData;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float is
    begin
      return this.invB_ptr(colIdx + this.row*rowIdx);
    end invB_ptr_out;

    function transMat_ptr_out ( this : Link_to_theData;
                                rowIdx : integer32; colIdx : integer32 )
                              return double_float is
    begin
      return this.transMat_ptr(colIdx + this.row*rowIdx);
    end transMat_ptr_out;

    procedure info_p_sol ( this : in Link_to_theData ) is
    begin
      put_line("<< p_sol >>");
      for i in 0..this.col-1 loop
        put(this.p_sol(i)); put(" ");
      end loop;
      new_line;
    end info_p_sol;

    procedure info_d_sol ( this : in Link_to_theData ) is
    begin
      put_line("<< d_sol >>");
      for i in 0..this.row-1 loop
        put(this.d_sol(i)); put(" ");
      end loop;
      new_line;
    end info_d_sol;

    procedure info_invB ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< invB >>");
      for i in 0..this.row-1 loop
        for j in 0..this.row-1 loop
          val := invB_out(this,i,j);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
    end info_invB;

    procedure info_transMat ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< transMat >>");
      for i in 0..this.row-1 loop
        for j in 0..this.row-1 loop
          val := transMat_out(this,i,j);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
    end info_transMat;

    procedure info_transRed ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< transRed >>");
      for i in 0..this.row-1 loop
        val := this.transRed(i);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
      end loop;
      new_line;
    end info_transRed;

    procedure info_basisIdx ( this : in Link_to_theData ) is
    begin
      put_line("<< basisIdx >>");
      for i in 0..this.row-1 loop
        put(this.basisIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_basisIdx;

    procedure info_nf_pos ( this : in Link_to_theData ) is
    begin
      put_line("<< nf_pos >>");
      for i in 0..this.nfN-1 loop
        put(this.nf_pos(i),1); put(" ");
      end loop;
      new_line;
    end info_nf_pos;

    procedure info_nbIdx ( this : in Link_to_theData ) is
    begin
      put_line("<< nbIdx >>");
      for i in 0..this.col-1 loop
        put(this.nbIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_nbIdx;

    procedure info_redVec ( this : in Link_to_theData ) is
    begin
      put_line("<< redVec >>");
      for i in 0..this.col-1 loop
        put(this.redVec(i)); put(" ");
      end loop;
      new_line;
    end info_redVec;

    procedure info_rIdx ( this : in Link_to_theData ) is
    begin
      put_line("<< rIdx >>");
      for i in 0..this.termS-1 loop
        put(this.rIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_rIdx;

    procedure info_pivOutIdx ( this : in Link_to_theData ) is
    begin
      put("pivOutCheck :");
      for i in 0..this.row-1 loop
        put(" "); put(this.pivOutCheck(i),1);
      end loop;
      new_line;
      put("pivOutList :");
      for i in 0..this.pivOutNum-1 loop
        put(" "); put(this.pivOutList(i),1);
      end loop;
      new_line;
    end info_pivOutIdx;

    procedure info_p_sol_ptr ( this : in Link_to_theData ) is

      use Standard_Floating_Vectors;

    begin
      put_line("<< p_sol_ptr >>");
      put("this.col-1 : "); put(this.col-1,1);
      if this.p_sol_ptr = null
       then put_line("  bug!");
      end if;
      put("  p_sol_ptr'last : "); put(this.p_sol_ptr'last,1); new_line;
      for i in 0..this.col-1 loop
        put(this.p_sol_ptr(i)); put(" ");
      end loop;
      new_line;
    end info_p_sol_ptr;

    procedure info_d_sol_ptr ( this : in Link_to_theData ) is
    begin
      put_line("<< d_sol_ptr >>");
      for i in 0..this.row-1 loop
        put(this.d_sol_ptr(i)); put(" ");
      end loop;
      new_line;
    end info_d_sol_ptr;

    procedure info_invB_ptr ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< invB_ptr >>");
      for i in 0..this.row-1 loop
        for j in 0..this.row-1 loop
          val := invB_ptr_out(this,i,j);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
    end info_invB_ptr;

    procedure info_transMat_ptr ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< transMat_ptr >>");
      for i in 0..this.row-1 loop
        for j in 0..this.row-1 loop
          val := transMat_ptr_out(this,i,j);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
    end info_transMat_ptr;

    procedure info_transRed_ptr ( this : in Link_to_theData ) is

      val : double_float;

    begin
      put_line("<< transRed_ptr >>");
      for i in 0..this.row-1 loop
        val := this.transRed_ptr(i);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO
           then put("0 ");
           else put(val); put(" ");
          end if;
      end loop;
      new_line;
    end info_transRed_ptr;

    procedure info_basisIdx_ptr ( this : in Link_to_theData ) is
    begin
      put_line("<< basisIdx_ptr >>");
      for i in 0..this.row-1 loop
        put(this.basisIdx_ptr(i),1); put(" ");
      end loop;
      new_line;
    end info_basisIdx_ptr;

    procedure info_nf_pos_ptr ( this : in Link_to_theData ) is
    begin
      put_line("<< nf_pos_ptr >>");
      for i in 0..this.nfN-1 loop
        put(this.nf_pos_ptr(i),1); put(" ");
      end loop;
      new_line;
    end info_nf_pos_ptr;

    procedure info_nbIdx_ptr ( this : in Link_to_theData ) is
    begin
      put_line("<< nbIdx_ptr >>");
      for i in 0..this.col-1 loop
        put(this.nbIdx_ptr(i),1); put(" ");
      end loop;
      new_line;
    end info_nbIdx_ptr;

    procedure info_redVec_ptr ( this : in Link_to_theData ) is
    begin
      put_line("<< redVec_ptr >>");
      for i in 0..this.col-1 loop
        put(this.redVec_ptr(i)); put(" ");
      end loop;
      new_line;
    end info_redVec_ptr;

    procedure info_fIdx ( this : in Link_to_theData ) is
    begin
      put(this.fIdx+1,1); new_line;
    end info_fIdx;

    procedure info_node ( this : in Link_to_theData ) is
    begin
      put("( ");
      for i in 0..this.polyDim loop
        put(this.nodeLabel(i)+1,1); put(" ");
      end loop;
      put(") ");
    end info_node;

  end class_theData;

  package body class_ftData is

    function new_ftData return ftData is

      res : ftData;

    begin
      res.elemNum := 0;
      res.cur := null;
      res.parent := null;
      res.limit := null;
      res.head := null;
      res.last := null;
      res.dim := 0;
      return res;
    end new_ftData;

    procedure delete_ftData ( this : in Link_to_ftData ) is
    begin
      null; -- indeed, empty destructor
    end delete_ftData;

    procedure clear ( lftd : in out Link_to_Array_of_ftData ) is

      procedure free is
        new unchecked_deallocation(Array_of_ftData,
                                   Link_to_Array_of_ftData);

    begin
      free(lftd);
    end clear;

    procedure create_elem
                ( this : in Link_to_ftData;
                  row : in integer32; col : in integer32;
                  termS : in integer32; polyDim : in integer32;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_ftest.class_ftData.create_elem, row : ");
        put(row,1); put(", col : "); put(col,1); new_line;
        put("  termS : "); put(termS,1); 
        put(", polyDim : "); put(polyDim,1);put_line(" ...");
      end if;
      declare
        newData : constant Link_to_theData := new theData'(new_theData);
      begin
        class_theData.create(newData,row,col,termS,polyDim,vrblvl-1);
        this.cur := newData;
      end;
      this.dim := row;
    end create_elem;

    procedure add_elem ( this : in Link_to_ftData;
                         vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_ftest.class_ftData.add_elem ...");
      end if;
      if this.last /= null then
        this.last.next := this.cur;
      else
        this.head := this.cur;
        this.parent := this.cur;
      end if;
      this.last := this.cur;
      this.elemNum := this.elemNum + 1;
    end add_elem;

    procedure mark ( this : in Link_to_ftData; vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_ftest.class_ftData.mark ...");
      end if;
      this.limit := this.last;
    end mark;

    procedure clear ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;

    begin
      while curr /= null loop
        class_theData.clear(curr);
        curr := curr.next;
      end loop;
    end clear;

    procedure clear_transMat ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;

    begin
      while curr /= null loop
        class_theData.clear_transMat(curr);
        curr := curr.next;
      end loop;
    end clear_transMat;

    procedure delete_cur ( this : in Link_to_ftData ) is
    begin
      class_theData.delete_theData(this.cur);
    end delete_cur;

    procedure delete_all ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;
      tmp : Link_to_theData;

    begin
      while curr /= null loop
        tmp := curr.next;
        class_theData.delete_theData(curr);
        curr := tmp;
      end loop;
      this.cur := null;
      this.parent := null;
      this.head := null;
      this.last := null;
      this.elemNum := 0;
    end delete_all;

    procedure delete_addedElem ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.limit.next;
      tmp : Link_to_theData;

    begin
      while curr /= null loop
        tmp := curr.next;
        class_theData.delete_theData(curr);
        curr := tmp;
      end loop;
      this.limit.next := null;
      this.last := this.limit;
      this.cur := this.last; 
    end delete_addedElem;

    procedure init_ptr ( this : in Link_to_ftData ) is
    begin
      this.parent := this.head;
      this.cur := this.head;
    end init_ptr;

    procedure make_init_data
                ( this : in Link_to_ftData;
                  termSumNum : in integer32; supN : in integer32;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put_line("-> in demics_fTest.class_ftData.make_init_data ...");
      end if;
      for i in 0..this.dim-1 loop
        this.cur.nf_pos(i) := i;
        this.cur.invB(i*(this.dim+1)) := 1.0;
        this.cur.transMat(i*(this.dim+1)) := 1.0;
        this.cur.basisIdx(i) := termSumNum - supN + i;
        this.cur.d_sol(i) := 1.0;
      end loop;
      for i in 0..termS-2 loop
        this.cur.nbIdx(i) := reTermS + i;
        this.cur.rIdx(i) := -1*(i+1);
      end loop;
    end make_init_data;

    procedure next_data ( this : in Link_to_ftData ) is
    begin
      if this.parent /= null
       then this.parent := this.parent.next;
      end if;
    end next_data;

    procedure copy ( this : in Link_to_ftData;
                     col : in integer32; pre_data : in Link_to_theData ) is
    begin
      for i in 0..col-1 loop
        this.cur.p_sol(i) := pre_data.p_sol_ptr(i);
      end loop;
      for i in 0..this.dim-1 loop
        this.cur.d_sol(i) := pre_data.d_sol_ptr(i);
        this.cur.basisIdx(i) := pre_data.basisIdx_ptr(i);
      end loop;
      for i in 0..col-this.dim-1 loop
        this.cur.nbIdx(i) := pre_data.nbIdx_ptr(i);
      end loop;
    end copy;

    procedure get_ptr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData ) is
    begin
      Standard_Floating_Vectors.clear(this.cur.p_sol);
      Standard_Floating_Vectors.clear(this.cur.d_sol);
      Standard_Integer_Vectors.clear(this.cur.basisIdx);
      Standard_Integer_Vectors.clear(this.cur.nbIdx);
      this.cur.p_sol := pre_data.p_sol;
      this.cur.d_sol := pre_data.d_sol;
      this.cur.basisIdx := pre_data.basisIdx;
      this.cur.nbIdx := pre_data.nbIdx;
    end get_ptr;

    procedure create_rIdx
                ( this : in Link_to_ftData;
                  preNbN : in integer32; repIdx : in integer32;
                  candIdx : in Standard_Integer_Vectors.Link_to_Vector ) is

      cnt : integer32 := 0;
      candNum : constant integer32 := candIdx(0);
      tmp_val : constant integer32 := preNbN - this.dim + 1;
      idx : integer32;

    begin
      for i in 0..candNum-1 loop
        idx := candIdx(i+1);
        if repIdx > idx then
          this.cur.rIdx(idx) := -1*(tmp_val + cnt);
          cnt := cnt + 1;
        elsif repIdx < idx then
          this.cur.rIdx(idx-1) := -1*(tmp_val + cnt);
          cnt := cnt + 1;
        end if;
      end loop;
    end create_rIdx;

    procedure init_info ( this : in Link_to_ftData ) is
    begin
      this.cur.pivOutNum := 0;
      for i in 0..this.dim-1 loop
        this.cur.pivOutCheck(i) := 0;
        this.cur.transRed(i) := 0.0;
      end loop;
    end init_info;

    procedure get_nbIdx_rIdx
                ( this : in Link_to_ftData;
                  preNbN : in integer32; repIdx : in integer32;
                  candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  reTermS : in integer32; pre_data : in Link_to_theData ) is

      pre_length : constant integer32 := preNbN - this.dim;
      candNum : constant integer32 := candIdx(0);
      cnt : integer32 := 0;
      idx,sub_idx : integer32;

    begin
      for i in 0..pre_length-1 loop
        this.cur.nbIdx(i) := pre_data.nbIdx_ptr(i);
      end loop;
      for i in 0..candNum-1 loop
        idx := candIdx(i+1);
        if repIdx > idx then
          sub_idx := pre_length + cnt;
          this.cur.nbIdx(sub_idx) := idx + reTermS;
          this.cur.rIdx(idx) := -1*(sub_idx + 1);
          cnt := cnt + 1;
        elsif repIdx < idx then
          sub_idx := pre_length + cnt;
          this.cur.nbIdx(sub_idx) := idx + reTermS - 1;
          this.cur.rIdx(idx-1) := -1*(sub_idx + 1);
          cnt := cnt + 1;
        end if;
      end loop;
    end get_nbIdx_rIdx;

    procedure iCopy ( this : in Link_to_ftData;
                      preNbN : in integer32; nfN : in integer32;
                      repIdx : in integer32; -- termS : in integer32;
                      reTermS : in integer32;
                      candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                      pre_data : in Link_to_theData ) is

      pre_length : constant integer32 := preNbN - this.dim;
      candNum : constant integer32 := candIdx(0);
      cnt : integer32 := 0;
      idx,sub_idx,ii,nfPos : integer32;

    begin
      for i in 0..nfN-1 loop
        this.cur.nf_pos(i) := pre_data.nf_pos_ptr(i);
      end loop;
      for i in 0..candNum-1 loop
        idx := candIdx(i+1);
        if repIdx > idx then
          sub_idx := pre_length + cnt;
          this.cur.nbIdx(sub_idx) := idx + reTermS;
          this.cur.rIdx(idx) := -1*(sub_idx + 1);
          cnt := cnt + 1;
        elsif repIdx < idx then
          sub_idx := pre_length + cnt;
          this.cur.nbIdx(sub_idx) := idx + reTermS - 1;
          this.cur.rIdx(idx-1) := -1*(sub_idx + 1);
          cnt := cnt + 1;
        end if;
      end loop;
      for j in 0..nfN-1 loop
        nfPos := pre_data.nf_pos_ptr(j);
        ii := nfPos*this.dim;
        for i in 0..this.dim-1 loop
          this.cur.invB(ii+i) := pre_data.invB_ptr(ii+i);
        end loop;
      end loop;
    end iCopy;

    procedure iGetPtr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData ) is
    begin
      this.cur.invB_ptr := pre_data.invB_ptr;
      this.cur.p_sol_ptr := pre_data.p_sol_ptr;
      this.cur.d_sol_ptr := pre_data.d_sol_ptr;
      this.cur.basisIdx_ptr := pre_data.basisIdx_ptr;
      this.cur.nf_pos_ptr := pre_data.nf_pos_ptr;
    end iGetPtr;

    procedure output ( this : in Link_to_ftData;
                       repIdx : in integer32; idx2 : out integer32;
                       nbN : out integer32; nfN : out integer32 ) is
    begin
      idx2 := this.parent.rIdx(repIdx-1);
      nbN := this.parent.nbN;
      nfN := this.parent.nfN;
    end output;

    procedure decrease_nfN ( this : in Link_to_ftData ) is
    begin
      this.cur.nfN := this.cur.nfN - 1;
    end decrease_nfN;

    procedure copy_rIdx
                ( this : in Link_to_ftData; pre_data : in Link_to_theData;
                  termS : in integer32 ) is
    begin
      this.cur.artV := pre_data.artV;
      for i in 0..termS-1 loop
        this.cur.rIdx(i) := pre_data.rIdx(i);
      end loop;
    end copy_rIdx;

    procedure copy_pivOutIdx ( this : in Link_to_ftData;
                               pre_data : in Link_to_theData ) is

      pivOutNum,idx : integer32;

    begin
      for i in 0..this.dim-1 loop
        this.cur.pivOutCheck(i) := 0;
      end loop;
      pivOutNum := pre_data.pivOutNum;
      this.cur.pivOutNum := pivOutNum;
      for i in 0..pivOutNum-1 loop
        idx := pre_data.pivOutList(i);
        this.cur.pivOutCheck(idx) := 1;
        this.cur.pivOutList(i) := idx;
      end loop;
    end copy_pivOutIdx;

    procedure get_nf_pos
                ( this : in Link_to_ftData; pre_data : in Link_to_theData;
                  nfN : in integer32; idx2 : in integer32 ) is

      cnt : integer32 := 0;
     -- nfPos : integer32; -- assigned but never read ...

    begin
      for i in 0..nfN-1 loop
       -- nfPos := pre_data.nf_pos_ptr(i);
        if pre_data.nf_pos_ptr(i) /= idx2 then
          this.cur.nf_pos(cnt) := pre_data.nf_pos_ptr(i);
          cnt := cnt + 1;
        end if;
      end loop;
      this.cur.nfN := this.cur.nfN - 1;
    end get_nf_pos;

    procedure mCopy ( this : in Link_to_ftData;
                      nbN : in integer32; nfN : in integer32;
                      idx2 : in integer32; termS : in integer32;
                      pre_data : in Link_to_theData ) is
    begin
      null;
    end mCopy;

    procedure mGetPtr ( this : in Link_to_ftData;
                        pre_data : in Link_to_theData ) is
    begin
      Standard_Floating_Vectors.clear(this.cur.p_sol);
      Standard_Floating_Vectors.clear(this.cur.d_sol);
      Standard_Integer_Vectors.clear(this.cur.basisIdx);
      Standard_Integer_Vectors.clear(this.cur.nbIdx);
      this.cur.p_sol := pre_data.p_sol;
      this.cur.d_sol := pre_data.d_sol;
      this.cur.basisIdx := pre_data.basisIdx;
      this.cur.nbIdx := pre_data.nbIdx;
    end mGetPtr;

    procedure put_sup ( this : in Link_to_ftData; sup : out integer32 ) is
    begin
      if this.parent /= null
       then sup := this.parent.nodeLabel(0);
      end if;
    end put_sup;

    procedure info_parent_nbN_nfN ( this : in Link_to_ftData ) is
    begin
      put("nbN : "); put(this.parent.nbN,1); new_line;
      put("nfN : "); put(this.parent.nfN,1); new_line;
    end info_parent_nbN_nfN;

    procedure info_parent_p_sol ( this : in Link_to_ftData ) is
    begin
      class_theData.info_p_sol(this.parent);
    end info_parent_p_sol;

    procedure info_parent_d_sol ( this : in Link_to_ftData ) is
    begin
      class_theData.info_d_sol(this.parent);
    end info_parent_d_sol;

    procedure info_parent_invB ( this : in Link_to_ftData ) is
    begin
      class_theData.info_invB(this.parent);
    end info_parent_invB;

    procedure info_parent_transMat ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transMat(this.parent);
    end info_parent_transMat;

    procedure info_parent_transRed ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transRed(this.parent);
    end info_parent_transRed;

    procedure info_parent_basisIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_basisIdx(this.parent);
    end info_parent_basisIdx;

    procedure info_parent_nf_pos ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nf_pos(this.parent);
    end info_parent_nf_pos;

    procedure info_parent_nbIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nbIdx(this.parent);
    end info_parent_nbIdx;

    procedure info_parent_redVec ( this : in Link_to_ftData ) is
    begin
      class_theData.info_redVec(this.parent);
    end info_parent_redVec;

    procedure info_parent_rIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_rIdx(this.parent);
    end info_parent_rIdx;

    procedure info_parent_pivOutIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_pivOutIdx(this.parent);
    end info_parent_pivOutIdx;

    procedure info_parent_p_sol_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_p_sol_ptr(this.parent);
    end info_parent_p_sol_ptr;

    procedure info_parent_d_sol_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_d_sol_ptr(this.parent);
    end info_parent_d_sol_ptr;

    procedure info_parent_invB_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_invB_ptr(this.parent);
    end info_parent_invB_ptr;

    procedure info_parent_transMat_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transMat_ptr(this.parent);
    end info_parent_transMat_ptr;

    procedure info_parent_transRed_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transRed_ptr(this.parent);
    end info_parent_transRed_ptr;

    procedure info_parent_basisIdx_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_basisIdx_ptr(this.parent);
    end info_parent_basisIdx_ptr;

    procedure info_parent_nf_pos_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nf_pos_ptr(this.parent);
    end info_parent_nf_pos_ptr;

    procedure info_parent_nbIdx_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nbIdx_ptr(this.parent);
    end info_parent_nbIdx_ptr;

    procedure info_parent_redVec_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_redVec_ptr(this.parent);
    end info_parent_redVec_ptr;

    procedure info_parent ( this : in Link_to_ftData ) is
    begin
      info_parent_p_sol(this);
      info_parent_d_sol(this);
      info_parent_invB(this);
      info_parent_basisIdx(this);
      info_parent_nf_pos(this);
      info_parent_nbIdx(this);
      info_parent_redVec(this);
      info_parent_rIdx(this);
    end info_parent;

    procedure info_parent_ptr ( this : in Link_to_ftData ) is
    begin
      info_parent_p_sol_ptr(this);
      info_parent_d_sol_ptr(this);
      info_parent_invB_ptr(this);
      info_parent_basisIdx_ptr(this);
      info_parent_nf_pos_ptr(this);
      info_parent_nbIdx_ptr(this);
      info_parent_redVec_ptr(this);
    end info_parent_ptr;

    procedure info_parent_node ( this : in Link_to_ftData ) is
    begin
      class_theData.info_node(this.parent);
    end info_parent_node;

    procedure info_cur_nbN_nfN ( this : in Link_to_ftData ) is
    begin
      put("nbN : "); put(this.cur.nbN,1); new_line;
      put("nfN : "); put(this.cur.nfN,1); new_line;
    end info_cur_nbN_nfN;

    procedure info_cur_p_sol ( this : in Link_to_ftData ) is
    begin
      class_theData.info_p_sol(this.cur);
    end info_cur_p_sol;

    procedure info_cur_d_sol ( this : in Link_to_ftData ) is
    begin
      class_theData.info_d_sol(this.cur);
    end info_cur_d_sol;

    procedure info_cur_invB ( this : in Link_to_ftData ) is
    begin
      class_theData.info_invB(this.cur);
    end info_cur_invB;

    procedure info_cur_transMat ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transMat(this.cur);
    end info_cur_transMat;

    procedure info_cur_transRed ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transRed(this.cur);
    end info_cur_transRed;

    procedure info_cur_basisIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_basisIdx(this.cur);
    end info_cur_basisIdx;

    procedure info_cur_nf_pos ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nf_pos(this.cur);
    end info_cur_nf_pos;

    procedure info_cur_nbIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nbIdx(this.cur);
    end info_cur_nbIdx;

    procedure info_cur_redVec ( this : in Link_to_ftData ) is
    begin
      class_theData.info_redVec(this.cur);
    end info_cur_redVec;

    procedure info_cur_rIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_rIdx(this.cur);
    end info_cur_rIdx;

    procedure info_cur_pivOutIdx ( this : in Link_to_ftData ) is
    begin
      class_theData.info_pivOutIdx(this.cur);
    end info_cur_pivOutIdx;

    procedure info_cur_p_sol_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_p_sol_ptr(this.cur);
    end info_cur_p_sol_ptr;

    procedure info_cur_d_sol_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_d_sol_ptr(this.cur);
    end info_cur_d_sol_ptr;

    procedure info_cur_invB_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_invB_ptr(this.cur);
    end info_cur_invB_ptr;

    procedure info_cur_transMat_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transMat_ptr(this.cur);
    end info_cur_transMat_ptr;

    procedure info_cur_transRed_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_transRed_ptr(this.cur);
    end info_cur_transRed_ptr;

    procedure info_cur_basisIdx_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_basisIdx_ptr(this.cur);
    end info_cur_basisIdx_ptr;

    procedure info_cur_nf_pos_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nf_pos_ptr(this.cur);
    end info_cur_nf_pos_ptr;

    procedure info_cur_nbIdx_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_nbIdx_ptr(this.cur);
    end info_cur_nbIdx_ptr;

    procedure info_cur_redVec_ptr ( this : in Link_to_ftData ) is
    begin
      class_theData.info_redVec_ptr(this.cur);
    end info_cur_redVec_ptr;

    procedure info_cur ( this : in Link_to_ftData ) is
    begin
      info_cur_p_sol(this);
      info_cur_d_sol(this);
      info_cur_invB(this);
      info_cur_transMat(this);
      info_cur_basisIdx(this);
      info_cur_nf_pos(this);
      info_cur_nbIdx(this);
      info_cur_redVec(this);
      info_cur_rIdx(this);
    end info_cur;

    procedure info_cur_ptr ( this : in Link_to_ftData ) is
    begin
      info_cur_p_sol_ptr(this);
      info_cur_d_sol_ptr(this);
      info_cur_invB_ptr(this);
      info_cur_transMat_ptr(this);
      info_cur_basisIdx_ptr(this);
      info_cur_nf_pos_ptr(this);
      info_cur_nbIdx_ptr(this);
      info_cur_redVec_ptr(this);
    end info_cur_ptr;

    procedure info_cur_node ( this : in Link_to_ftData ) is
    begin
      class_theData.info_node(this.cur);
    end info_cur_node;

    procedure info_all_node ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;
      i : integer32 := 0;

    begin
      put_line("<< info_all_node >>");
      while curr /= null loop
        put("# "); put(i+1,1); new_line;
        info_node(curr);
        curr := curr.next;
        i := i+1;
      end loop;
      new_line;
    end info_all_node;

    procedure info_all_cur ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;
      num : integer32 := 0;

    begin
      put_line("<< info_all_cur >>");
      while curr /= null loop
        put("# "); put(num+1,1); new_line;
        info_p_sol(curr);
        info_d_sol(curr);
        info_invB(curr);
        info_basisIdx(curr);
        info_nf_pos(curr);
        info_nbIdx(curr);
        info_redVec(curr);
        info_rIdx(curr);
        curr := curr.next;
        num := num+1;
      end loop;
    end info_all_cur;

    procedure info_all_nodeNum ( this : in Link_to_ftData ) is

      curr : Link_to_theData := this.head;
      i : integer32 := 0;

    begin
      put_line("<< info_all_nodeNum >>");
      while curr /= null loop
        put("# "); put(i+1,1); new_line;
        curr := curr.next;
        i := i+1;
      end loop;
      new_line;
    end info_all_nodeNum;

    procedure info_numElem ( this : in Link_to_ftData ) is

      num : integer32 := 0;
      curr : Link_to_theData := this.head;

    begin
      while curr /= null loop
        curr := curr.next;
        num := num + 1;
      end loop;
      put(num,1); put(" ");
    end info_numElem;

  end class_ftData;

  package body class_lvData is

    function new_lvData return lvData is

      res : lvData;

    begin
      res.fTest := null;
      res.node := null;
      res.mRepN := null;
      res.mFeaIdx := null;
      res.mFea := null;
      res.dim := 0;
      res.length := 0;
      res.termMax := 0;
      return res;
    end new_lvData;

    procedure delete_lvData ( this : in Link_to_lvData ) is

      use Standard_Integer_VecVecs;

    begin
      for i in 0..this.length-1 loop
        class_ftData.delete_all(this.fTest(i));
      end loop;
      class_ftData.clear(this.fTest);
      Standard_Integer_Vectors.clear(this.mRepN);
      Standard_Integer_Vectors.clear(this.mFea);
      if this.mFeaIdx /= null then
        Standard_Integer_VecVecs.Deep_Clear(this.mFeaIdx);
        this.mFeaIdx := null;
      end if;
    end delete_lvData;

    procedure clear ( lvd : in out Link_to_Array_of_lvData ) is

      procedure free is
        new unchecked_deallocation(Array_of_lvData,
                                   Link_to_Array_of_lvData);

    begin
      free(lvd);
    end clear;

    procedure create ( this : in Link_to_lvData; depth : in integer32;
                      -- supN : in integer32; dim : in integer32;
                       ori_length : in integer32; ori_termMax : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_ftest.class_lvData.create, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      this.length := ori_length;
      this.termMax := ori_termMax;
      this.ftest := new Array_of_ftData(0..this.length-1);
      this.node := this.ftest(this.length-1);
     -- this.node := this.ftest(this.length-1)'access;
      this.mRepN
        := new Standard_Integer_Vectors.Vector'(0..this.length-1 => 0);
      this.mFeaIdx := new Standard_Integer_VecVecs.VecVec(0..this.length-1);
      for i in 0..this.length-1 loop
        this.mFeaIdx(i)
          := new Standard_Integer_Vectors.Vector(0..this.termMax-1);
      end loop;
      this.mFea := new Standard_Integer_Vectors.Vector'(0..this.length-1 => 0);
    end create;

    procedure get_info
                ( this : in Link_to_lvData;
                  g_mRepN : out Standard_Integer_Vectors.Link_to_Vector;
                  g_mFeaIdx : out Standard_Integer_VecVecs.Link_to_VecVec;
                  g_mFea : out Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      g_mRepN := this.mRepN;
      g_mFeaIdx := this.mFeaIdx;
      g_mFea := this.mFea;
    end get_info;

    procedure init_ptr ( this : in Link_to_lvData ) is
    begin
      for i in 0..this.length-1 loop
        class_ftData.init_ptr(this.fTest(i));
      end loop;
    end init_ptr;

    procedure info_mFea ( this : in Link_to_lvData ) is
    begin
      put_line("mFea :");
      for i in 0..this.length-1 loop
        put(" "); put(this.mFea(i),1);
      end loop;
      new_line;
      put_line("mRepN :");
      for i in 0..this.length-1 loop
        put(" "); put(this.mRepN(i),1);
      end loop;
      new_line;
    end info_mFea;

  end class_lvData;

end demics_ftest;
