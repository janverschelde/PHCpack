with Ada.text_io;                       use Ada.text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with DEMiCs_Global_Constants;

package body demics_simplex is

  package body class_supportSet is

    function new_supportSet return supportSet is

      res : supportSet;

    begin
      res.row := 0;
      res.col := 0;
      res.supMat := null;
      res.costVec := null;
      return res;
    end new_supportSet;

    procedure delete_supportSet ( this : in Link_to_supportSet ) is
    begin
      Standard_Floating_Vectors.clear(this.supMat);
      Standard_Floating_Vectors.clear(this.costVec);
    end delete_supportSet;

    procedure allocSupp
                ( this : in Link_to_supportSet;
                  data : demics_input_data.class_dataSet.dataSet;
                  level : in integer32; num : in integer32;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;

      use demics_input_data.class_dataSet; -- for support_out

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_supportSet.allocSupp, level : ");
        put(level,1); put(", num : "); put(num,1); put_line(" ...");      
      end if;
      this.row := data.dim;
      this.col := data.termSet(level);
      this.supMat := new Standard_Floating_Vectors.Vector'
                           (0..this.row*this.col-1 => 0.0);
      this.costVec
        := new Standard_Floating_Vectors.Vector'(0..this.col-2 => 0.0);
      for j in 0..this.col-1 loop
        if j /= num then
          for i in 0..this.row-1 loop
            supMat_in(this,i,cnt,
                      support_out(data,data.termStart(level) + j,i)
                    - support_out(data,data.termStart(level) + num,i));
          end loop;
          this.costVec(cnt) := lifting(data.termStart(level) + j)
                             - lifting(data.termStart(level) + num);
          cnt := cnt + 1;
        end if;
      end loop;
      this.col := this.col - 1;
    end allocSupp;

    procedure allocAux
                ( this : in Link_to_supportSet;
                  data : in demics_input_data.class_dataSet.dataSet;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_supportSet.allocAux ...");
      end if;
      this.row := data.dim;
      this.col := data.dim;
      this.supMat := new Standard_Floating_Vectors.Vector'
                           (0..this.row*this.col-1 => 0.0);
      this.costVec := new Standard_Floating_Vectors.Vector(0..this.col-1);
      for i in 0..this.col-1 loop
        supMat_in(this,i,i,1.0);
        this.costVec(i) := 1.0;
      end loop;
    end allocAux;

    procedure supMat_in ( this : in Link_to_supportSet;
                          rowIdx : in integer32; colIdx : in integer32;
                          elem : in double_float ) is
    begin
      this.supMat(rowIdx + colIdx*this.row) := elem;
    end supMat_in;

    procedure supMat_neg ( this : in Link_to_supportSet;
                           rowIdx : in integer32; colIdx : in integer32 ) is
    begin
      this.supMat(rowIdx + colIdx*this.row) 
        := -this.supMat(rowIdx + colIdx*this.row);
    end supMat_neg;

    function supMat_out ( this : Link_to_supportSet;
                          rowIdx : integer32; colIdx : integer32 )
                        return double_float is

     -- idx : constant integer32 := rowIdx + colIdx*this.row;

    begin
     -- put("idx : "); put(idx,1);
     -- put(", this.supMat'last : "); put(this.supMat'last,1); new_line;
      return this.supMat(rowIdx + colIdx*this.row);
    end supMat_out;

    function redVal ( this : Link_to_supportSet;
                      d_sol : Standard_Floating_Vectors.Link_to_Vector;
                      idx : integer32; ii : integer32 )
                    return double_float is

      val : double_float := 0.0;

    begin
      for i in 0..this.row-1 loop
        val := val + d_sol(i)*this.supMat(i+ii);
      end loop;
      return (this.costVec(idx) - val);
    end redVal;

    procedure info_sup ( this : in Link_to_supportSet ) is
    begin
      for j in 0..this.row-1 loop
        for i in 0..this.col-1 loop
          put(supMat_out(this,j,i)); put(" ");
        end loop;
        new_line;
      end loop;
    end info_sup;

    procedure info_costVec ( this : in Link_to_supportSet ) is
    begin
      for i in 0..this.col-1 loop
        put(this.costVec(i)); put(" ");
      end loop;
      new_line;
    end info_costVec;

  end class_supportSet;

  package body class_simplex is

-- relation table

    function checkFrIdx ( this : Link_to_simplex;
                          vrblvl : integer32 := 0 ) return integer32 is

      flag : integer32 := DEMiCs_Global_Constants.CONTINUE;
      redFlag,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx : integer32;
      redCost,theta : double_float;

    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_simplex.checkFrIdx ...");
      end if;
      for i in 0..this.nbN - this.dim - 1 loop
        if this.nbIdx(i) = this.frIdx then
          pivInIdx := this.frIdx;
          sub_pivInIdx := i;
          flag := DEMiCs_Global_Constants.PIVOT_IN;
          exit;
        end if;
      end loop;
      if flag = DEMiCs_Global_Constants.PIVOT_IN then
        calRedCost(this,pivInIdx,redCost,vrblvl-1);
        if redCost > DEMiCs_Global_Constants.PLUSZERO
         then redFlag := DEMiCs_Global_Constants.NEGTHETA;
         else redFlag := DEMiCs_Global_Constants.POSTHETA;
        end if;
        ratioTest_art(this,redFlag,pivInIdx,--sub_pivInIdx,
                      pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        if flag = DEMiCs_Global_Constants.CONTINUE then
          elimFrIdx(this,sub_pivOutIdx,vrblvl-1);
          createNewBandN_tab
            (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx, 
             theta, redCost,vrblvl-1);
        else
          if vrblvl > 0 
           then put_line("----- UNB.checkFrIdx -----");
          end if;
        end if;
      else
        for i in 0..this.dim-1 loop
          if this.basisIdx(i) = this.frIdx then
            elimFrIdx(this,i,vrblvl-1);
            exit;
          end if;
        end loop;
      end if;
      return flag;
    end checkFrIdx;

    procedure elimFrIdx ( this : in Link_to_simplex;
                          sub_pivOutIdx : in integer32;
                          vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.elimFrIdx, sub_pivOutIdx : ");
        put(sub_pivOutIdx,1); put_line(" ...");
      end if;
      for i in 0..this.nfN-1 loop
        if i /= sub_pivOutIdx then
          this.nf_pos(cnt) := this.nf_pos(i);
          cnt := cnt + 1;
        end if;
      end loop;
      this.nfN := this.nfN - 1;
    end elimFrIdx;

-- phase 1

    procedure reMakeNonBasisIdx
                ( this : in Link_to_simplex; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;
      artNbIdxN : integer32 := 0;
      tmp_nbIdx : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.reMakeNonBasisIdx, ");
        put("reTermS : "); put(reTermS,1); put_line(" ...");
      end if;
      for i in 0..this.nbN - this.dim - 1 loop
        if this.nbIdx(i) < this.termSumNum - this.supN then
          tmp_nbIdx := this.nbIdx(i);
          this.nbIdx(cnt) := tmp_nbIdx;
          this.rIdx(tmp_nbIdx - reTermS) := -(cnt+1);
          cnt := cnt + 1;
        else
          this.nbIdx(i) := 0;
          artNbIdxN := artNbIdxN + 1;
        end if;
      end loop;
      if artNbIdxN = this.dim
       then this.artV := 0;
       else this.artV := 1;
      end if;
      IP_vec_mat(this,vrblvl-1);
      this.nbN := this.dim + cnt;
    end reMakeNonBasisIdx;

    procedure reMakeNonBasisIdx_tab
                ( this : in Link_to_simplex; vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;
      tmp_nbIdx : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.");
        put_line("reMakeNonBasisIdx_tab ...");
      end if;
      for i in 0..this.nbN - this.dim - 1 loop
        if this.nbIdx(i) < this.termSumNum - this.supN then
          tmp_nbIdx := this.nbIdx(i);
          this.nbIdx(cnt) := tmp_nbIdx;
          cnt := cnt + 1;
        else
          this.nbIdx(i) := 0;
        end if;
      end loop;
      IP_vec_mat(this,vrblvl-1);
      this.nbN := this.dim + cnt;
    end reMakeNonBasisIdx_tab;

    procedure elimArt ( this : in Link_to_simplex;
                        depth : in integer32; preNbN : in integer32;
                        termS : in integer32; reTermS : in integer32;
                        iter : in out integer32;
                        vrblvl : in integer32 := 0 ) is

      sub_pivInIdx,flag : integer32;
      cnt : integer32 := 0;
      cnt_t : integer32 := 0;
      redCost : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.elimArt, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      for i in 0..this.dim-1 loop
        if this.basisIdx(i) >= this.termSumNum - this.supN then
          isZeroDirEle(this,termS,i,preNbN,sub_pivInIdx,flag,vrblvl-1);
          cnt := cnt + 1;
         -- True  -- val > 0
         -- False -- val =0
          if flag = DEMiCs_Global_Constants.TRUE then
	    calRedCost(this,this.nbIdx(sub_pivInIdx),redCost,vrblvl-1);
            IP_mat_vec(this,this.nbIdx(sub_pivInIdx),vrblvl-1);
            createNewBandN_art
              (this,this.nbIdx(sub_pivInIdx),sub_pivInIdx,--this.basisIdx(i),
               i,redCost,termS,reTermS,vrblvl-1);
            iter := iter + 1;
            cnt_t := cnt_t + 1;
          end if;
        end if;
      end loop;
      if cnt = cnt_t
       then this.artV := 0;
       else this.artV := 1;
      end if;
    end elimArt;

    procedure calRedCost ( this : in Link_to_simplex;
                           pivInIdx : in integer32;
                           redCost : out double_float;
                           vrblvl : in integer32 := 0 ) is

      ii,level,idx,idx2 : integer32;
      val : double_float := 0.0;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.calRedCost, pivInIdx : ");
        put(pivInIdx,1); put_line(" ...");
      end if;
      ii := 2*pivInIdx;
      level := this.nIdx(ii);
      idx := this.nIdx(ii + 1);
      idx2 := this.firIdx(level);
      ii := this.dim * idx;
      for i in 0..this.dim-1 loop  
        val := val + this.d_sol(i) * this.supp(level)(idx2).supMat(i + ii);
      end loop;
      redCost := this.supp(level)(idx2).costVec(idx) - val;
    end calRedCost;

    procedure isZeroDirEle
                ( this : in Link_to_simplex;
                  termS : in integer32; idx : in integer32;
                  preNbN : in integer32; sub_pivInIdx : out integer32;
                  result : out integer32; vrblvl : in integer32 := 0 ) is

      ii,iii,idx2,level : integer32;
      val : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.isZeroDirEle, termS : ");
        put(termS,1); put(", idx : "); put(idx,1); put_line(" ...");
      end if;
      for j in 0..termS-2 loop
        ii := 2*this.nbIdx(preNbN - this.dim + j);
        level := this.nIdx(ii);
        idx2 := this.firIdx(level);
        iii := this.dim*this.nIdx(ii + 1);
        val := 0.0;
        ii := this.dim*idx;
        for i in 0..this.dim-1 loop
          val := val
               + this.invB(ii + i)*this.supp(level)(idx2).supMat(i+iii);
        end loop;
        if (val < DEMiCs_Global_Constants.MINUSZERO) or
           (val > DEMiCs_Global_Constants.PLUSZERO) then
          sub_pivInIdx := preNbN - this.dim + j;
          result := DEMiCs_Global_Constants.TRUE; return;
        end if;
      end loop;
      result := DEMiCs_Global_Constants.FALSE;
    end isZeroDirEle;

    procedure IP_vec_mat ( this : in Link_to_simplex;
                           vrblvl : in integer32 := 0 ) is

      ii,level,idx,idx2 : integer32;

    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_simplex.IP_vec_mat ...");
      end if;
      for j in 0..this.dim-1 loop
        if this.basisIdx(j) < this.termSumNum - this.supN then
          ii := 2*this.basisIdx(j);
          level := this.nIdx(ii);
          idx := this.nIdx(ii+1);
          idx2 := this.firIdx(level);
          ii := this.dim*j;
          for i in 0..this.dim-1 loop
            this.d_sol(i) := this.d_sol(i)
              + this.invB(ii+i)*this.supp(level)(idx2).costVec(idx);
          end loop;
        end if;
      end loop;
    end IP_vec_mat;

-- reduced cost

    procedure reducedCost_tab_p1
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      val,tmp_redCost : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_tab_p1 ...");
      end if;
      flag := opt;
      redCost := 1.0E-8;
      for j in 0..this.nbN - this.dim - 1 loop
        val := 0.0;
        tmp_non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        for i in 0..this.dim-1 loop
          val := val + this.d_sol(i)*this.supp(level)(idx2).supMat(i+ii);
        end loop;
        tmp_redCost := this.aux_cvec(tmp_non_basisIdx) - val;
        if (tmp_redcost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_tab_p1;

    procedure reducedCost_tab
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      tmp_redCost : double_float;
     -- val : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_tab ...");
      end if;
      flag := opt;
      redCost := 1.0E-8;
      for j in 0..this.nbN - this.dim - 1 loop
       -- val := 0.0;
        tmp_non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        tmp_redCost := class_supportSet.redVal
          (this.supp(level)(idx2),this.d_sol,idx,ii);
        if (tmp_redcost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_tab;

    procedure reducedCost_p1
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      val,tmp_redCost : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_p1 ...");
      end if;
      flag := opt;
      redCost := 1.0E-8;
      for j in 0..this.nbN - this.dim - 1 loop
        val := 0.0;
        tmp_non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        for i in 0..this.dim-1 loop
           val := val + this.d_sol(i)*this.Supp(level)(idx2).supMat(i+ii);
        end loop;
        tmp_redCost := this.aux_cvec(tmp_non_basisIdx) - val;
        this.redVec(tmp_non_basisIdx) := tmp_redCost;
        if (tmp_redCost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_p1;

    procedure reducedCost
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      tmp_redCost : double_float;
     -- val : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost ...");
      end if;
      flag := opt;
      redCost := 1.0E-8;
      for j in 0..this.nbN - this.dim - 1 loop
       -- val := 0.0;
        tmp_non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        tmp_redCost := class_supportSet.redVal
                         (this.supp(level)(idx2),this.d_sol,idx,ii);
        this.redVec(tmp_non_basisIdx) := tmp_redCost;

        if (tmp_redCost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost;

    procedure reducedCost_Bland
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      tmp_redCost : double_float;
     -- val : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_Bland ...");
      end if;
      flag := opt;
      pivInIdx := DEMiCs_Global_Constants.BIGINT;
      for j in 0..this.nbN - this.dim - 1 loop
       -- val := 0.0;
        tmp_non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        tmp_redCost := class_supportSet.redVal
                         (this.supp(level)(idx2),this.d_sol,idx,ii);
        this.redVec(tmp_non_basisIdx) := tmp_redCost;
        if (tmp_redCost < DEMiCs_Global_Constants.MINUSZERO) and
           (tmp_non_basisIdx < pivInIdx) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_Bland;

    procedure reducedCost_mFst
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32; sub_pivInIdx : in out integer32;
                  pivOutIdx : in integer32;
                 -- sub_pivOutIdx : in integer32;
                  redCost : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,tmp_non_basisIdx,level,idx,idx2 : integer32;
      pre_pivOutIdx,pre_sub_pivInIdx : integer32;
      tmp_redCost : double_float;
     -- val : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_mFst ...");
      end if;
      flag := opt;
      pre_pivOutIdx := pivOutIdx;
      pre_sub_pivInIdx := sub_pivInIdx;
      redCost := 1.0E-8;
      for j in 0..this.nbN - this.Dim - 1 loop
        if j = pre_sub_pivInIdx then
          this.nbIdx(j) := pre_pivOutIdx;
          tmp_non_basisIdx := pre_pivOutIdx;
        else
          this.nbIdx(j) := this.pre_nbIdx(j);
          tmp_non_basisIdx := this.pre_nbIdx(j);
        end if;
       -- val := 0.0;
        getIdx(this,level,idx,idx2,ii,2*tmp_non_basisIdx,vrblvl-1);
        tmp_redCost := class_supportSet.redVal
          (this.supp(level)(idx2),this.d_sol,idx,ii);
        this.redVec(tmp_non_basisIdx) := tmp_redCost;
        if (tmp_redCost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := tmp_non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_mFst;

    procedure reducedCost_iFst
                ( this : in Link_to_simplex;
                  pivInIdx : in out integer32;
                  sub_pivInIdx : in out integer32;
                  pivOutIdx : in integer32; -- sub_pivOutIdx : in integer32;
                  redCost : out double_float; -- termS : in integer32;
                  reTermS : in integer32; preNbN : in integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      ii,non_basisIdx,level,idx,idx2 : integer32;
      pre_pivOutIdx,pre_pivInIdx : integer32;
     -- pre_sub_pivInIdx : integer32;
      length,pre_length,candNum,cnt : integer32;
      tmp_redCost : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_iFst ...");
      end if;
      flag := opt;
      length := this.nbN - this.dim;
      pre_length := preNbN - this.dim;
      pre_pivInIdx := pivInIdx;
     -- pre_sub_pivInIdx := sub_pivInIdx;
      pre_pivOutIdx := pivOutIdx;
      for i in 0..pre_length loop
        this.nbIdx(i) := this.pre_nbIdx(i);
      end loop;
      candNum := this.nbN - preNbN + 1;
      cnt := 0;
      for i in 0..candNum-1 loop
        extend_nbIdx
          (this,this.candIdx(i+1),pre_pivInIdx,pre_pivOutIdx,
           pre_length,reTermS,cnt,vrblvl-1);
      end loop;
      redCost := 1.0E-8;
      for j in 0..length-1 loop
        non_basisIdx := this.nbIdx(j);
        getIdx(this,level,idx,idx2,ii,2*non_basisIdx,vrblvl-1);
        tmp_redCost := class_supportSet.redVal
          (this.supp(level)(idx2),this.d_sol,idx,ii);
        this.redVec(non_basisIdx) := tmp_redCost;
        if (tmp_redCost < DEMiCs_Global_Constants.MINUSZERO) and
           (abs(tmp_redCost) > abs(redCost)) then
          redCost := tmp_redCost;
          pivInIdx := non_basisIdx;
          sub_pivInIdx := j;
          flag := DEMiCs_Global_Constants.POSTHETA;
        end if;
      end loop;
    end reducedCost_iFst;

    procedure extend_nbIdx
                ( this : in Link_to_simplex; cIdx : in integer32;
                  pre_pivInIdx : in integer32; pre_pivOutIdx : in integer32;
                  pre_length : in integer32; reTermS : in integer32;
                  cnt : in out integer32; vrblvl : in integer32 := 0 ) is

      idx : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.extend_nbIdx, cIdx : ");
        put(cIdx,1); put_line(" ...");
      end if;
      if this.repIdx > cIdx then
        idx := cIdx + reTermS;
        if idx = pre_pivInIdx
         then this.nbIdx(pre_length + cnt) := pre_pivOutIdx;
         else this.nbIdx(pre_length + cnt) := cIdx + reTermS;
        end if;
        cnt := cnt + 1;
      elsif this.repIdx < cIdx then
        idx := cIdx + reTermS - 1;
        if idx = pre_pivInIdx
         then this.nbIdx(pre_length + cnt) := pre_pivOutIdx;
         else this.nbIdx(pre_length + cnt) := cIdx + reTermS - 1;
        end if;
        cnt := cnt + 1;
      end if;
    end extend_nbIdx;

    procedure extend_nbIdx_comp
                ( this : in Link_to_simplex;
                  non_basisIdx : out integer32; cIdx : in integer32;
                  pre_pivInIdx : in integer32; pre_pivOutIdx : in integer32;
                  pre_length : in integer32; reTermS : in integer32;
                  cnt : in out integer32; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      idx : integer32;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.extend_nbIdx_comp ... ");
      end if;
      flag := DEMiCs_Global_Constants.TRUE;
      if this.repIdx > cIdx then
         idx := cIdx + reTermS;
         if idx = pre_pivInIdx then
           this.nbIdx(pre_length + cnt) := pre_pivOutIdx;
           non_basisIdx := pre_pivOutIdx;
         else
           this.nbIdx(pre_length + cnt) := cIdx + reTermS;
           non_basisIdx := this.nbIdx(pre_length + cnt);
         end if;
         cnt := cnt + 1;
      elsif this.repIdx < cIdx then
         idx := cIdx + reTermS - 1;
         if idx = pre_pivInIdx then
           this.nbIdx(pre_length + cnt) := pre_pivOutIdx;
           non_basisIdx := pre_pivOutIdx;
         else
           this.nbIdx(pre_length + cnt) := cIdx + reTermS - 1;
           non_basisIdx := this.nbIdx(pre_length + cnt);
         end if;
         cnt := cnt + 1;
      else
        flag := DEMiCs_Global_Constants.FALSE;
      end if;
    end extend_nbIdx_comp;

    procedure getIdx ( this : in Link_to_simplex; level : out integer32;
                       idx : out integer32; idx2 : out integer32;
                       ii : out integer32; d_nbIdx : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_simplex.getIdx ... ");
      end if;
      level := this.nIdx(d_nbIdx);
      idx := this.nIdx(d_nbIdx+1);
      idx2 := this.firIdx(level);
      ii := this.dim*idx;
    end getIdx;

-- ratio test

    procedure ratioTest
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; -- sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      nfPos,checker,nonNegVarNum : integer32;
      tmp_theta : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest, redFlag : ");
        put(redFlag,1); put_line(" ...");
      end if;
      flag := 0;
      IP_mat_vec(this,pivInIdx,vrblvl-1);
      nonNegVarNum := 0;
      checker := 0;
      case redFlag is
        when DEMiCs_Global_Constants.POSTHETA =>
          tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop
            nfPos := this.nf_pos(i);
            nonNegVarNum := nonNegVarNum + 1;
            if this.dir(nfPos) < DEMiCs_Global_Constants.MINUSZERO then
              tmp_theta := this.p_sol(this.basisIdx(nfPos))/this.dir(nfPos);
            else
              checker := checker + 1;
              tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
            end if;
            if theta < tmp_theta then
              theta := tmp_theta;
              pivOutIdx := this.basisIdx(nfPos);
              sub_pivOutIdx := nfPos;
              flag := DEMiCs_Global_Constants.CONTINUE;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
          end if;
        when DEMiCs_Global_Constants.NEGTHETA =>
          tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop
            nfPos := this.nf_pos(i);
            nonNegVarNum := nonNegVarNum + 1;
            if this.dir(nfPos) > DEMiCs_Global_Constants.PLUSZERO then
              tmp_theta := this.p_sol(this.basisIdx(nfPos))/this.dir(nfPos);
            else
              checker := checker + 1;
              tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
            end if;
            if theta > tmp_theta then
              theta := tmp_theta;
              pivOutIdx := this.basisIdx(nfPos);
              sub_pivOutIdx := nfPos;
              flag := DEMiCs_Global_Constants.CONTINUE;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
          end if;
        when others => null;
      end case;
      theta := -theta;
    end ratioTest;

    procedure ratioTest_artFst
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; -- sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      nfPos,checker,nonNegVarNum : integer32;
      tmp_theta : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest_artFst, ");
        put("redFlag : "); put(redFlag,1); put_line(" ...");
      end if;
      flag := 0;
      IP_mat_vec_fst(this,pivInIdx,vrblvl-1);
      nonNegVarNum := 0;
      checker := 0;
      case redFlag is
        when DEMiCs_Global_Constants.POSTHETA =>
          tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop
            nfPos := this.pre_nf_pos(i);
            if this.pre_basisIdx(nfPos) < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) < DEMiCs_Global_Constants.MINUSZERO then
                tmp_theta
                  := this.pre_p_sol(this.pre_basisIdx(nfPos))/this.dir(nfPos);
              else
                checker := checker + 1;
                tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
              end if;
              if theta < tmp_theta then
                theta := tmp_theta;
                pivOutIdx := this.pre_basisIdx(nfPos);
                sub_pivOutIdx := nfPos;
                flag := DEMiCs_Global_Constants.CONTINUE;
              end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when DEMiCs_Global_Constants.NEGTHETA =>
          tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop
            nfPos := this.pre_nf_pos(i);
            if this.pre_basisIdx(nfPos) < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) > DEMiCs_Global_Constants.PLUSZERO then
                tmp_theta
                  := this.pre_p_sol(this.pre_basisIdx(nfPos))/this.dir(nfPos);
              else
                checker := checker + 1;
                tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
              end if;
              if theta > tmp_theta then
                theta := tmp_theta;
                pivOutIdx := this.pre_basisIdx(nfPos);
                sub_pivOutIdx := nfPos;
                flag := DEMiCs_Global_Constants.CONTINUE;
              end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when others => null;
      end case;
      theta := -theta;
    end ratioTest_artFst;

    procedure ratioTest_art
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; -- sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      nfPos,checker,nonNegVarNum : integer32;
      tmp_theta : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest_art, redFlag : ");
        put(redFlag,1); put_line(" ...");
      end if;
      flag := 0;
      IP_mat_vec(this,pivInIdx,vrblvl-1);
      nonNegVarNum := 0;
      checker := 0;
      case redFlag is
        when DEMiCs_Global_Constants.POSTHETA =>
          tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop   
            nfPos := this.nf_pos(i);
            if this.basisIdx(nfPos) < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) < DEMiCs_Global_Constants.MINUSZERO then
                tmp_theta := this.p_sol(this.basisIdx(nfPos))/this.dir(nfPos);
              else
                checker := checker + 1;
                tmp_theta := DEMiCs_Global_Constants.SMALLDOUBLE;
              end if;
              if theta < tmp_theta then
                theta := tmp_theta;
                pivOutIdx := this.basisIdx(nfPos);
                sub_pivOutIdx := nfPos;
                flag := DEMiCs_Global_Constants.CONTINUE;
	      end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when DEMiCs_Global_Constants.NEGTHETA =>
          tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
          theta := tmp_theta;
          for i in 0..this.nfN-1 loop
            nfPos := this.nf_pos(i);
            if this.basisIdx(nfPos) < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) > DEMiCs_Global_Constants.PLUSZERO then
                tmp_theta := this.p_sol(this.basisIdx(nfPos))/this.dir(nfPos);
              else
                checker := checker + 1;
                tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
	      end if;
              if theta > tmp_theta then
                theta := tmp_theta;
                pivOutIdx := this.basisIdx(nfPos);
                sub_pivOutIdx := nfPos;
                flag := DEMiCs_Global_Constants.CONTINUE;
              end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when others => null;
      end case;
      theta := -theta;
    end ratioTest_art;

    procedure ratioTest_art_Bland
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; -- sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

     nfPos,checker,nonNegVarNum,tmp_basisIdx : integer32;
     tmp_theta : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest_art_Bland, ");
        put("redFlag : "); put(redFlag,1); put_line(" ...");
      end if;
      flag := 0;
      IP_mat_vec(this,pivInIdx,vrblvl-1);
      nonNegVarNum := 0;
      checker := 0;
      case redFlag is
        when DEMiCs_Global_Constants.POSTHETA =>
          theta := 0.0;
          pivOutIdx := DEMiCs_Global_Constants.BIGINT;
          for i in 0..this.nfN-1 loop
            nfPos := this.nf_pos(i);
            tmp_basisIdx := this.basisIdx(nfPos);
            if tmp_basisIdx < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) < DEMiCs_Global_Constants.MINUSZERO then
                tmp_theta := this.p_sol(tmp_basisIdx)/this.dir(nfPos);
                if tmp_basisIdx < pivOutIdx then
                  theta := tmp_theta;
                  pivOutIdx := tmp_basisIdx;
                  sub_pivOutIdx := nfPos;
                  flag := DEMiCs_Global_Constants.CONTINUE;
                end if;
              else
                checker := checker + 1;
              end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when DEMiCs_Global_Constants.NEGTHETA =>
          theta := 0.0;
          pivOutIdx := DEMiCs_Global_Constants.BIGINT;
          for i in 0..this.nfN-1 loop
            nfPos := this.nf_pos(i);
            tmp_basisIdx := this.basisIdx(nfPos);
            if tmp_basisIdx < this.termSumNum - this.supN then
              nonNegVarNum := nonNegVarNum + 1;
              if this.dir(nfPos) > DEMiCs_Global_Constants.PLUSZERO then
                tmp_theta := this.p_sol(tmp_basisIdx)/this.dir(nfPos);
                if tmp_basisIdx < pivOutIdx then
                  theta := tmp_theta;
                  pivOutIdx := tmp_basisIdx;
                  sub_pivOutIdx := nfPos;
                  flag := DEMiCs_Global_Constants.CONTINUE;
                end if;
              else
                checker := checker + 1;
              end if;
            end if;
          end loop;
          if checker = nonNegVarNum
           then flag := DEMiCs_Global_Constants.UNBOUNDED;
           else flag := DEMiCs_Global_Constants.CONTINUE;
          end if;
        when others => null;
      end case;
      theta := -theta;
    end ratioTest_art_Bland;

    function ratioTest_frIdx
               ( this : Link_to_simplex; pivInIdx : integer32;
                 vrblvl : integer32 := 0 ) return integer32 is

      flag,nfPos,checker,nonNegVarNum : integer32;
      theta,tmp_theta : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest_frIdx, pivInIdx : ");
        put(pivInIdx,1); put_line(" ...");
      end if;
      flag := 0;
      IP_mat_vec(this,pivInIdx,vrblvl-1);
      nonNegVarNum := 0;
      checker := 0;
      tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
      theta := tmp_theta;
      for i in 0..this.nfN-1 loop
        nfPos := this.nf_pos(i);
        if this.basisIdx(nfPos) < this.termSumNum - this.supN then
          nonNegVarNum := nonNegVarNum + 1;
          if this.dir(nfPos) > DEMiCs_Global_Constants.PLUSZERO then
            tmp_theta := this.p_sol(this.basisIdx(nfPos))/this.dir(nfPos);
          else
            checker := checker + 1;
            tmp_theta := DEMiCs_Global_Constants.BIGDOUBLE;
          end if; 
          if theta > tmp_theta
           then theta := tmp_theta;
          end if;
        end if; 
      end loop;   
      if checker = nonNegVarNum
       then flag := DEMiCs_Global_Constants.UNBOUNDED;
       else flag := DEMiCs_Global_Constants.OPT;
      end if;
      return flag;
    end ratioTest_frIdx;

    procedure IP_mat_vec ( this : in Link_to_simplex;
                           pivInIdx : in integer32;
                           vrblvl : in integer32 := 0 ) is

      ii,iii,level,idx2,nfPos : integer32;
      val : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.IP_mat_vec, pivInIdx : ");
        put(pivInIdx,1); put_line(" ...");
      end if;
      ii := 2 * pivInIdx;
      level := this.nIdx(ii);
      idx2 := this.firIdx(level);
      iii := this.dim * this.nIdx(ii + 1);
      for j in 0..this.nfN-1 loop 
        nfPos := this.nf_pos(j);
        val := 0.0;
        ii := this.dim*nfPos;
        for i in 0..this.dim-1 loop
          val := val
               + this.invB(ii + i)*this.supp(level)(idx2).supMat(i + iii);
        end loop;
        this.dir(nfPos) := -val;
      end loop;
    end IP_mat_vec;

    procedure IP_mat_vec_fst ( this : in Link_to_simplex;
                               pivInIdx : in integer32;
                               vrblvl : in integer32 := 0 ) is

      ii,iii,level,idx2,nfPos : integer32;
      val : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.IP_mat_vec_fst, pivInIdx : ");
        put(pivInIdx,1); put_line(" ...");
      end if;
      ii := 2 * pivInIdx;
      level := this.nIdx(ii);
      idx2 := this.firIdx(level);
      iii := this.dim*this.nIdx(ii + 1);
      for j in 0..this.nfN-1 loop
        nfPos := this.pre_nf_pos(j);
        val := 0.0;
        ii := this.dim * nfPos;
        for i in 0..this.dim-1 loop
          val := val + this.pre_invB(ii + i)
                     * this.supp(level)(idx2).supMat(i + iii);
        end loop;
        this.dir(nfPos) := -val;
      end loop;
    end IP_mat_vec_fst;

    procedure update_p1_d_sol
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivOutIdx : in integer32;
                  vrblvl : in integer32 := 0 ) is

      ii, idx, idx2, level : integer32;
      redCost, val : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.update_p1_d_sol, ");
        put("pivInIdx : "); put(pivInIdx,1);
        put(", sub_pivOutIdx : ");
        put(sub_pivOutIdx,1); put_line(" ...");
      end if;
      getIdx(this,level,idx,idx2,ii,2*pivInIdx,vrblvl-1);
      redCost := class_supportSet.redVal
                   (this.supp(level)(idx2),this.p1_d_sol,idx, ii);
      ii := sub_pivOutIdx*this.dim;
      val := redCost/this.dir(sub_pivOutIdx);
     -- make new d_sol
      for i in 0..this.dim-1 loop
        this.p1_d_sol(i) := this.p1_d_sol(i) - val*this.invB(ii+i);
        this.transRed(i) := this.transRed(i) - val*this.transMat(ii+i);
      end loop;
    end update_p1_d_sol;

    procedure modify_p_sol ( this : in Link_to_simplex;
                             pivInIdx : in integer32;
                             vrblvl : in integer32 := 0 ) is

      level : constant integer32 := this.nIdx(2*pivInIdx);
      idx : constant integer32 := this.nIdx(2*pivInIdx+1);
      idx2 : constant integer32 := this.firIdx(level);

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.modify_p_sol, pivInIdx : ");
        put(pivInIdx,1); put_line(" ...");
      end if;
      if vrblvl > 0 then -- #if DBG_MODIFY
        put_line("------- << modify_p_sol >> -------");
        info_p_sol(this);
        info_basisIdx(this);
        put("pivInIdx : "); put(pivInIdx,1); new_line;
        put("supVec : ");
        for i in 0..this.dim-1 loop
          put(class_supportSet.supMat_out(this.supp(level)(idx2),i,idx));
          put(" ");
        end loop;
        new_line;
      end if;
      for i in 0..this.dim-1 loop
        if (isZero(this.p_sol(this.basisIdx(i)))
             = DEMiCs_Global_Constants.TRUE) and
           (isZero(class_supportSet.supMat_out(this.supp(level)(idx2),i,idx))
             = DEMiCs_Global_Constants.FALSE)
         then calElem(this,i,vrblvl-1);
        end if;
      end loop;
    end modify_p_sol;

    procedure calElem ( this : in Link_to_simplex; idx : in integer32;
                        vrblvl : in integer32 := 0 ) is

      val,ddelta,randNum : double_float;
      -- The original delta variable is reserved in Ada,
      -- replaced original delta variable by ddelta.
      bigdouble : constant := DEMiCs_Global_Constants.BIGDOUBLE;
      smalldouble : constant := DEMiCs_Global_Constants.SMALLDOUBLE;
      upper : double_float := bigdouble;
      lower : double_float := smalldouble;
      tmp_upper,tmp_lower : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.calElem, idx : ");
        put(idx,1); put_line(" ...");
      end if;
      if vrblvl > 0 then -- #if DBG_MODIFY
       -- put_line("<< calElem >>");
       -- put("idx : "); put(idx,1); new_line;
        info_invB(this);
        put("p_sol : ");
        for i in 0..this.dim-1 loop
          put(this.p_sol(this.basisIdx(i))); put(" ");
        end loop;
        new_line;
        put("colum-invB : ");
        for i in 0..this.dim-1 loop
          put(invB_out(this,i,idx)); put(" ");
        end loop;
        new_line;
      end if;
      for i in 0..this.dim-1 loop
        val := invB_out(this,i,idx);
        if val > DEMiCs_Global_Constants.PLUSZERO then
          tmp_lower := -this.p_sol(this.basisIdx(i))/val;
          if tmp_lower > lower
           then lower := tmp_lower;
          end if;
        end if;
        if val < DEMiCs_Global_Constants.MINUSZERO then
          tmp_upper := -this.p_sol(this.basisIdx(i))/val;
          if tmp_upper < upper
           then upper := tmp_upper;
          end if;
        end if;
      end loop;
      if vrblvl > 0 then -- #if DBG_MODIFY
        put("lower : "); put(lower); new_line;
        put("upper : "); put(upper); new_line;
      end if;
      if (lower = 0.0) and (upper = 0.0) then
        ddelta := 0.0;
      else
        if (lower = smalldouble) and (upper = bigdouble) then
          Standard_Random_Numbers.Set_Seed(2);
          ddelta := abs(Standard_Random_Numbers.Random);
        elsif (lower = smalldouble) and (upper < bigdouble) then
          Standard_Random_Numbers.Set_Seed(3);        
          randNum := abs(Standard_Random_Numbers.Random);
          ddelta := upper - randNum/2.0;
        elsif (lower > smalldouble) and (upper = bigdouble) then
          Standard_Random_Numbers.Set_Seed(4);
          randNum := abs(Standard_Random_Numbers.Random);
          ddelta := lower + randNum/2.0;
        else
          ddelta := (lower + upper)/2.0;
        end if;
      end if;
      if vrblvl > 0 then -- #if DBG_MODIFY
        put("delta : "); put(ddelta); new_line;
      end if;
      for i in 0..this.dim-1 loop
        this.p_sol(this.basisIdx(i))
          := this.p_sol(this.basisIdx(i)) + ddelta*invB_out(this,i,idx);
      end loop;
      if vrblvl > 0 then -- #if DBG_MODIFY
        put("-- modified --");
        info_p_sol(this);
      end if;
    end calElem;

-- create new basis and nonbasis

    procedure createNewBandN_tab
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  vrblvl : in integer32 := 0 ) is

      ii,nfPos : integer32;
      elem,vval : double_float;
     -- val : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_tab ...");
      end if;
      ii := sub_pivOutIdx * this.dim;
      elem := this.dir(sub_pivOutIdx);
     -- val := (-1.0 - elem) / elem;
      vval := redCost / elem;
      for i in 0..this.dim-1 loop -- make new d_sol
        this.d_sol(i) := this.d_sol(i) - vval * this.invB(ii + i);
      end loop;
      for i in 0..this.nfN-1 loop -- make new p-sol
        nfPos := this.nf_pos(i);
        if nfPos /= sub_pivOutIdx then
          this.p_sol(this.basisIdx(nfPos))
            := this.p_sol(this.basisIdx(nfPos)) + theta * this.dir(nfPos);
        else
          this.p_sol(this.basisIdx(nfPos)) := 0.0;
        end if;
      end loop;
      this.p_sol(pivInIdx) := theta;
     -- make new basisIdx and non_basisIdx
      this.basisIdx(sub_pivOutIdx) := pivInIdx;
      this.nbIdx(sub_pivInIdx) := pivOutIdx;
      for i in 0..this.dim-1 loop -- make the new basis matrix
        this.invB(ii + i) := this.invB(ii + i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii + i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim * nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii + i)
              := this.invB(ii + i) + this.dir(nfPos) * this.tmp_newInvB(i);
          end loop;
        end if;
      end loop;
    end createNewBandN_tab;

    procedure createNewBandN_p1
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is

      ii,nfPos : integer32;
      elem,val,vval : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_p1 ...");
      end if;
      ii := sub_pivOutIdx*this.dim;
      elem := this.dir(sub_pivOutIdx);
      val := (-1.0 - elem)/elem;
      vval := redCost / elem;
     -- make new d_sol
      for i in 0..this.dim-1 loop
        this.d_sol(i):= this.d_sol(i) - vval*this.invB(ii+i);
      end loop;
     -- make new p_sol
      for i in 0..this.nfN-1 loop
        nfPos := this.nf_pos(i);
        if nfPos /= sub_pivOutIdx then
          this.p_sol(this.basisIdx(nfPos))
            := this.p_sol(this.basisIdx(nfPos)) + theta*this.dir(nfPos);
        else
          this.p_sol(this.basisIdx(nfPos)) := 0.0;
        end if;
      end loop;
      this.p_sol(pivInIdx) := theta;
     -- make new basisIdx and non_basisIdx
      this.basisIdx(sub_pivOutIdx) := pivInIdx;
      this.nbIdx(sub_pivInIdx) := pivOutIdx;
      if this.pivOutCheck(sub_pivOutIdx) = 0 then
        this.pivOutCheck(sub_pivOutIdx) := 1;
        this.pivOutList(this.pivOutNum) := sub_pivOutIdx;
        this.pivOutNum := this.pivOutNum + 1;
      end if;
      if (reTermS <= pivInIdx) and (pivInIdx < reTermS + termS - 1) then
        this.rIdx(pivInIdx - reTermS) := sub_pivOutIdx;
      end if;
      if(reTermS <= pivOutIdx) and (pivOutIdx < reTermS + termS - 1) then
        this.rIdx(pivOutIdx - reTermS) := -(1 + sub_pivInIdx);
      end if;
     -- make the new basis matrix   
      for i in 0..this.dim-1 loop
       -- tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);
        this.invB(ii+i) := this.invB(ii+i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii+i);
        this.tmp_transMat(i) := this.transMat(ii+i);
        this.transMat(ii+i) := this.transMat(ii+i) + val*this.tmp_transMat(i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim*nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii+i)
              := this.invB(ii+i) + this.dir(nfPos)*this.tmp_newInvB(i);
            this.transMat(ii+i) := this.transMat(ii+i)
               - this.tmp_transMat(i)*(this.dir(nfPos)/elem);
          end loop;
        end if;
      end loop;
    end createNewBandN_p1;

    procedure createNewBandN
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is

      ii,nfPos : integer32;
      elem,val,vval : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN ...");
      end if;
      ii := sub_pivOutIdx*this.dim;
      elem := this.dir(sub_pivOutIdx);
      val := (-1.0 - elem)/elem;
      vval := redCost/elem;
     -- make new d_sol
      for i in 0..this.dim-1 loop
        this.d_sol(i) := this.d_sol(i) - vval*this.invB(ii+i);
        this.transRed(i) := this.transRed(i) - vval*this.transMat(ii+i);
      end loop;
     -- make new p_sol
      for i in 0..this.nfN-1 loop
        nfPos := this.nf_pos(i);
        if nfPos /= sub_pivOutIdx then
          this.p_sol(this.basisIdx(nfPos))
            := this.p_sol(this.basisIdx(nfPos)) + theta*this.dir(nfPos);
        else
          this.p_sol(this.basisIdx(nfPos)) := 0.0;
        end if;
      end loop;
      this.p_sol(pivInIdx) := theta;
     -- make new basisIdx and non_basisIdx
      this.basisIdx(sub_pivOutIdx) := pivInIdx;
      this.nbIdx(sub_pivInIdx) := pivOutIdx;
      if this.pivOutCheck(sub_pivOutIdx) = 0 then
        this.pivOutCheck(sub_pivOutIdx) := 1;
        this.pivOutList(this.pivOutNum) := sub_pivOutIdx;
        this.pivOutNum := this.pivOutNum + 1;
      end if;
      if(reTermS <= pivInIdx) and (pivInIdx < reTermS + termS - 1) then
        this.rIdx(pivInIdx - reTermS) := sub_pivOutIdx;
      end if;
      if(reTermS <= pivOutIdx) and (pivOutIdx < reTermS + termS - 1) then
        this.rIdx(pivOutIdx - reTermS) := -(1 + sub_pivInIdx);
      end if;
     -- make the new basis matrix   
      for i in 0..this.dim-1 loop
        this.invB(ii + i) := this.invB(ii + i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii+i);
        this.tmp_transMat(i) := this.transMat(ii+i);
        this.transMat(ii+i)
          := this.transMat(ii+i) +  val*this.tmp_transMat(i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim*nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii+i)
              := this.invB(ii+i) + this.dir(nfPos)*this.tmp_newInvB(i);
            this.transMat(ii+i) := this.transMat(ii+i)
              -  this.tmp_transMat(i)*(this.dir(nfPos)/elem);
          end loop;
        end if;
      end loop;
    end createNewBandN;

    procedure createNewBandN_iFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is

      ii,nfPos : integer32;
     -- cnt : integer32 := 0;
      elem,val,vval : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_iFst ...");
      end if;
      ii := sub_pivOutIdx*this.dim;
      elem := this.dir(sub_pivOutIdx);
      val := (-1.0 - elem)/elem;
      vval := redCost/elem;
    -- make new basisIdx and non_basisIdx
      for i in 0..this.dim-1 loop
        this.d_sol(i) := this.pre_d_sol(i) - vval*this.pre_invB(ii+i);
        this.transRed(i) := this.transRed(i) - vval*this.eye(ii+i);
        this.nf_pos(i) := this.pre_nf_pos(i);
        if i /= sub_pivOutIdx
         then this.basisIdx(i) := this.pre_basisIdx(i);
         else this.basisIdx(i) := pivInIdx;
        end if;
      end loop;
      if this.pivOutCheck(sub_pivOutIdx) = 0 then
        this.pivOutCheck(sub_pivOutIdx) := 1;
        this.pivOutList(this.pivOutNum) := sub_pivOutIdx;
        this.pivOutNum := this.pivOutNum + 1;
      end if;
      if(reTermS <= pivInIdx) and (pivInIdx < reTermS + termS - 1) then
        this.rIdx(pivInIdx - reTermS) := sub_pivOutIdx;
      end if;
      if(reTermS <= pivOutIdx) and (pivOutIdx < reTermS + termS - 1) then
        this.rIdx(pivOutIdx - reTermS) := -(1 + sub_pivInIdx);
      end if;
     -- make new p_sol
      for i in 0..this.nfN-1 loop
        nfPos := this.nf_pos(i);
        if nfPos /= sub_pivOutIdx then
          this.p_sol(this.basisIdx(nfPos))
            := this.pre_p_sol(this.basisIdx(nfPos)) + theta*this.dir(nfPos);
        else
          this.p_sol(this.basisIdx(nfPos)) := 0.0;
        end if;
      end loop;
      this.p_sol(pivInIdx) := theta;
     -- make the new basis matrix   
      for i in 0..this.dim-1 loop
        this.invB(ii+i) := this.pre_invB(ii+i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii+i);
        this.tmp_transMat(i) := this.eye(ii+i);
        this.transMat(ii+i) := this.eye(ii+i) + val*this.tmp_transMat(i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim*nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii+i)
              := this.pre_invB(ii+i) + this.dir(nfPos)*this.tmp_newInvB(i);
            this.transMat(ii+i)
              := this.eye(ii+i) - this.tmp_transMat(i)*(this.dir(nfPos)/elem);
          end loop;
        end if;
      end loop;
    end createNewBandN_iFst;

    procedure createNewBandN_mFst
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  termS : in integer32; reTermS : in integer32;
                  vrblvl : in integer32 := 0 ) is

      ii,cnt,nfPos : integer32;
      elem,vval : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_mFst ...");
      end if;
      ii := sub_pivOutIdx*this.dim;
      elem := this.dir(sub_pivOutIdx);
      vval := redCost/elem;
     -- make new basisIdx and non_basisIdx
      cnt := 0;
      for i in 0..this.dim-1 loop
        this.d_sol(i) := this.pre_d_sol(i) - vval*this.pre_invB(ii+i);
        this.transRed(i)
          := this.pre_transRed(i) - vval*this.pre_transMat(ii+i);
        if this.pre_nf_pos(i) /= sub_pivOutIdx then
          this.nf_pos(cnt) := this.pre_nf_pos(i);
          cnt := cnt + 1;
        end if;
        if i /= sub_pivOutIdx
         then this.basisIdx(i) := this.pre_basisIdx(i);
         else this.basisIdx(i) := pivInIdx;
        end if;
      end loop;
      this.nfN := this.nfN - 1;
      if this.pivOutCheck(sub_pivOutIdx) = 0 then
        this.pivOutCheck(sub_pivOutIdx) := 1;
        this.pivOutList(this.pivOutNum) := sub_pivOutIdx;
        this.pivOutNum := this.pivOutNum + 1;
      end if;
      if (reTermS <= pivInIdx) and (pivInIdx < reTermS + termS - 1) then
        this.rIdx(pivInIdx - reTermS) := sub_pivOutIdx;
      end if;
      if (reTermS <= pivOutIdx) and (pivOutIdx < reTermS + termS - 1) then
        this.rIdx(pivOutIdx - reTermS) := -(1 + sub_pivInIdx);
      end if;
     -- make new p_sol
      for i in 0..this.nfN-1 loop
        nfPos := this.nf_pos(i);
        if nfPos /= sub_pivOutIdx then
          this.p_sol(this.basisIdx(nfPos))
            := this.pre_p_sol(this.basisIdx(nfPos)) + theta*this.dir(nfPos);
        else
          this.p_sol(this.basisIdx(nfPos)) := 0.0;
        end if;
      end loop;
      this.p_sol(pivInIdx) := theta;
     -- make the new basis matrix   
      for i in 0..this.dim-1 loop
        this.invB(ii+i) := this.pre_invB(ii+i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii+i);
        this.tmp_transMat(i) := this.pre_transMat(ii+i);
        this.transMat(ii+i)
          := this.pre_transMat(ii+i) + (-1.0-elem)/elem*this.tmp_transMat(i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim*nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii+i) :=
              this.pre_invB(ii+i) + this.dir(nfPos)*this.tmp_newInvB(i);
            this.transMat(ii+i) := this.pre_transMat(ii+i)
              - this.tmp_transMat(i)*(this.dir(nfPos)/elem);
          end loop;
        end if;
      end loop;
    end createNewBandN_mFst;

    procedure createNewBandN_art
                ( this : in Link_to_simplex;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                 -- pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : in double_float; termS : in integer32;
                  reTermS : in integer32; vrblvl : in integer32 := 0 ) is

      ii,nfPos : integer32;
      elem,val,vval : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_art ...");
      end if;
      ii := sub_pivOutIdx*this.dim;
      elem := this.dir(sub_pivOutIdx);
      val := (-1.0 - elem)/elem;
      vval := redCost/elem;
     -- make new d_sol
      for i in 0..this.dim-1 loop
        this.d_sol(i) := this.d_sol(i) - vval*this.invB(ii+i);
        this.transRed(i) := this.transRed(i) - vval*this.transMat(ii+i);
      end loop;
     -- make new p_sol
      modify_p_sol(this,pivInIdx,vrblvl-1);
     -- make new basisIdx and non_basisIdx
      this.basisIdx(sub_pivOutIdx) := pivInIdx;
      if this.pivOutCheck(sub_pivOutIdx) = 0 then
        this.pivOutCheck(sub_pivOutIdx) := 1;
        this.pivOutList(this.pivOutNum) := sub_pivOutIdx;
        this.pivOutNum := this.pivOutNum + 1;
      end if;
      if (reTermS <= pivInIdx) and (pivInIdx < reTermS + termS - 1) then
        this.rIdx(pivInIdx - reTermS) := sub_pivOutIdx;
      end if;
      for i in sub_pivInIdx..this.nbN - this.dim - 1 loop
        this.nbIdx(i) := this.nbIdx(i+1);
        if (reTermS <= this.nbIdx(i)) and
           (this.nbIdx(i) < reTermS + termS - 1) then
          this.rIdx(this.nbIdx(i) - reTermS) := -(i+1);
        end if;
      end loop;
      this.nbN := this.nbN - 1;
     -- make the new basis matrix   
      for i in 0..this.dim-1 loop
       -- tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);
        this.invB(ii+i) := this.invB(ii+i)/(-elem);
        this.tmp_newInvB(i) := this.invB(ii+i);
        this.tmp_transMat(i) := this.transMat(ii+i);
        this.transMat(ii+i) := this.transMat(ii+i) + val*this.tmp_transMat(i);
      end loop;
      for j in 0..this.nfN-1 loop
        nfPos := this.nf_pos(j);
        if nfPos /= sub_pivOutIdx then
          ii := this.dim*nfPos;
          for i in 0..this.dim-1 loop
            this.invB(ii+i)
              := this.invB(ii+i) + this.dir(nfPos)*this.tmp_newInvB(i);
            this.transMat(ii+i) := this.transMat(ii+i)
              - this.tmp_transMat(i)*(this.dir(nfPos)/elem);
          end loop; 
        end if;
      end loop;
    end createNewBandN_art;

    procedure invB_in ( this : in Link_to_simplex;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ) is
    begin
      this.invB(colIdx + this.dim*rowIdx) := elem;
    end invB_in;

    function invB_out ( this : Link_to_simplex;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return this.invB(colIdx + this.dim*rowIdx);
    end invB_out;

    function transMat_out ( this : Link_to_simplex;
                            rowIdx : integer32; colIdx : integer32 )
                          return double_float is
    begin
      return this.transMat(colIdx + this.dim*rowIdx);
    end transMat_out;

    procedure supp_in ( this : in Link_to_simplex; lvl : in integer32;
                        rowIdx : in integer32; colIdx : in integer32;
                        elem : in double_float ) is
    begin
      this.oriSupp(lvl)(rowIdx + colIdx*this.dim) := elem;
    end supp_in;

    function supp_out ( this : Link_to_simplex; lvl : integer32;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return this.oriSupp(lvl)(rowIdx + colIdx*this.dim);
    end supp_out;

    function isZero ( val : double_float ) return integer32 is
    begin
      if val > DEMiCs_Global_Constants.MINUSZERO and
         val < DEMiCs_Global_Constants.PLUSZERO
       then return DEMiCs_Global_Constants.TRUE;
       else return DEMiCs_Global_Constants.FALSE;
      end if;
    end isZero;

    procedure info_p_sol ( this : in Link_to_simplex ) is
    begin
      put_line("<< p_sol >>");
      for i in 0..this.termSumNum - this.supN + this.dim - 1 loop
        put(this.p_sol(i)); put(" ");
      end loop;
      new_line;
    end info_p_sol;

    procedure info_d_sol ( this : in Link_to_simplex ) is
    begin
      put_line("<< d_sol >>");
      for i in 0..this.dim-1 loop
        put(this.d_sol(i)); put(" ");
      end loop;
      new_line;
    end info_d_sol;

    procedure info_p1_d_sol ( this : in Link_to_simplex ) is
    begin
      put_line("<< p1_d_sol >>");
      for i in 0..this.dim-1 loop
        put(this.p1_d_sol(i)); put(" ");
      end loop;
      new_line;
    end info_p1_d_sol;

    procedure info_invB ( this : in Link_to_simplex ) is
    begin
      put_line("<< invB >>");
      for i in 0..this.dim-1 loop
        for j in 0..this.dim-1 loop
          put(invB_out(this,i,j)); put(" ");
        end loop;
        new_line;
      end loop;
      new_line;
    end info_invB;

    procedure info_transMat ( this : in Link_to_simplex ) is

      val : double_float;

    begin
      put_line("<< transMat >>");
      for i in 0..this.dim-1 loop
        for j in 0..this.dim-1 loop
          val := transMat_out(this,i,j);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO then
            put("0 ");
          else
            put(transMat_out(this,i,j)); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
      new_line;
    end info_transMat;

    procedure info_transRed ( this : in Link_to_simplex ) is

      val : double_float;

    begin
      put_line("<< transRed >>");
      for i in 0..this.dim-1 loop
        for j in 0..this.dim-1 loop
          val := this.transRed(i);
          if val < DEMiCs_Global_Constants.PLUSZERO and
             val > DEMiCs_Global_Constants.MINUSZERO then
            put("0 ");
          else
            put(val); put(" ");
          end if;
        end loop;
        new_line;
      end loop;
      new_line;
    end info_transRed;

    procedure info_basisIdx ( this : in Link_to_simplex ) is
    begin
      put_line("<< basisIdx >>");
      for i in 0..this.dim-1 loop
        put(this.basisIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_basisIdx;

    procedure info_nf_pos ( this : in Link_to_simplex ) is
    begin
      put_line("<< nf_pos >>");
      for i in 0..this.nfN-1 loop
        put(this.nf_pos(i),1); put(" ");
      end loop;
      new_line;
    end info_nf_pos;

    procedure info_nbIdx ( this : in Link_to_simplex ) is
    begin
      put_line("<< nbIdx >>");
      for i in 0..this.nbN - this.dim - 1 loop
        put(this.nbIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_nbIdx;

    procedure info_rIdx ( this : in Link_to_simplex ) is
    begin
      put_line("<< rIdx >>");
      put("this.rIdx'last : "); put(this.rIdx'last,1); new_line;
      put("this.nbN-1 : "); put(this.nbN-1,1);
      if this.rIdx'last = this.nbN
       then put_line(" OK.");
       else put_line(" Bug??!");
      end if;
     -- for i in 0..this.nbN-1 loop
      for i in 0..this.rIdx'last loop
        put(this.rIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_rIdx;

    procedure info_redVec ( this : in Link_to_simplex ) is
    begin
      put_line("<< redVec >>");
      for i in 0..this.nbN - this.dim - 1 loop
        put(this.redVec(i)); put(" ");
      end loop;
      new_line;
    end info_redVec;

    procedure info_dir ( this : in Link_to_simplex ) is
    begin
      put_line("<< dir >>");
      for i in 0..this.dim-1 loop
        put(this.dir(i)); put(" ");
      end loop;
      new_line;
    end info_dir;

    procedure info_frIdx ( this : in Link_to_simplex ) is
    begin
      put("frIdx : "); put(this.frIdx,1); new_line;
    end info_frIdx;

    procedure info_candIdx ( this : in Link_to_simplex ) is
    begin
      put_line("<< candIdx >>");
      for i in 0..this.candIdx(0)-1 loop
        put(this.candIdx(i + 1),1); put(" ");
      end loop;
      new_line;
    end info_candIdx;

    procedure info_repIdx ( this : in Link_to_simplex ) is
    begin
      put("repIdx : "); put(this.repIdx,1); new_line;
    end info_repIdx;

    procedure info_oriSup ( this : in Link_to_simplex ) is
    begin
      put_line("<< oriSup >>");
      for k in 0..this.supN-1 loop
        for j in 0..this.dim-1 loop
          for i in 0..this.termSet(k)-1 loop
            put(supp_out(this,k,j,i)); put(" ");
          end loop;
          new_line;
        end loop;
        new_line;
      end loop;
      new_line;
    end info_oriSup;

    function new_simplex return simplex is

      res : simplex;

    begin
      res.dim := 0;
      res.supN := 0;
      res.termSumNum := 0;
      res.repIdx := 0;
      res.output := 0;
      res.mixedVol := 0.0;
      res.mixedCell := 0;
      res.ip := null;
      res.weight := null;
      res.vol := null;
      res.eye := null;
      res.nbN := 0;
      res.nfN := 0;
      res.artV := 0;
      res.pivOutNum := 0;
      res.frIdx := 0;
      res.Supp := null;
      res.lifting := null;
      res.oriSupp := null;
      res.firIdx := null;
      res.termSet := null;
      res.termStart := null;
      res.re_termStart := null;
      res.invB := null;
      res.transMat := null;
      res.transRed := null;
      res.p_sol := null;
      res.d_sol := null;
      res.p1_d_sol := null;
      res.fst_d_sol := null;
      res.aux_cvec := null;
      res.dir := null;
      res.fst_redVec := null;
      res.redVec := null;
      res.basisIdx := null;
      res.nf_pos := null;
      res.nbIdx := null;
      res.rIdx := null;
      res.pivOutList := null;
      res.pivOutCheck := null;
      res.tmp_newInvB := null;
      res.tmp_transMat := null;
      res.nIdx := null;
      res.pre_p_sol := null;
      res.pre_d_sol := null;
      res.pre_redVec := null;
      res.pre_basisIdx := null;
      res.pre_nbIdx := null;
      res.pre_nf_pos := null;
      res.pre_invB := null;
      res.pre_transMat := null;
      res.pre_transRed := null;
      return res;
    end new_simplex;

    procedure delete_simplex ( this : in Link_to_simplex ) is

      use Standard_Floating_VecVecs;

    begin
      Standard_Integer_Vectors.Clear(this.re_termStart);
      Standard_Floating_Vectors.Clear(this.p1_d_sol);
      Standard_Integer_Vectors.Clear(this.ip);
      Standard_Floating_Vectors.Clear(this.weight);
      Standard_Floating_Vectors.Clear(this.vol);
      Standard_Floating_Vectors.Clear(this.eye);
      Standard_Floating_Vectors.Clear(this.lifting);
      Standard_Floating_Vectors.Clear(this.fst_d_sol);
      Standard_Floating_Vectors.Clear(this.aux_cvec);
      Standard_Floating_Vectors.Clear(this.dir);
      Standard_Floating_Vectors.Clear(this.fst_redVec);
      Standard_Floating_Vectors.Clear(this.tmp_newInvB);
      Standard_Floating_Vectors.Clear(this.tmp_transMat);
      Standard_Integer_Vectors.Clear(this.nIdx);
      if this.Supp /= null then
        for i in 0..this.supN loop
          for j in this.Supp(i)'range loop
            class_supportSet.delete_supportSet(this.Supp(i)(j));
          end loop;
        end loop;
        this.supp := null;
      end if;
      if this.oriSupp /= null then
        for i in 0..this.supN-1 loop
          Standard_Floating_Vectors.Clear(this.oriSupp(i));
        end loop;
        Standard_Floating_VecVecs.Shallow_Clear(this.oriSupp);
        this.oriSupp := null;
      end if;
    end delete_simplex;

    procedure get_iNbN_nfN
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData;
                  lNbN : in integer32; lNfN : in integer32 ) is
    begin
      cur.nbN := lNbN; this.nbN := lNbN;
      cur.nfN := lNfN; this.nfN := lNfN;
    end get_iNbN_nfN;

    procedure get_mNbN_nfN
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  cur : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      this.nbN := parent.nbN; cur.nbN := parent.nbN;
      this.nfN := parent.nfN; cur.nfN := parent.nfN;
    end get_mNbN_nfN;

    procedure get_repIdx_candIdx
                ( this : in Link_to_simplex;
                  ori_candIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_repIdx : in integer32 ) is
    begin
      this.repIdx := ori_repIdx;
      this.candIdx := ori_candIdx;
    end get_repIdx_candIdx;

    procedure get_parent
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      this.pre_p_sol := parent.p_sol_ptr;
      this.pre_d_sol := parent.d_sol_ptr;
      this.pre_redVec := parent.redVec_ptr;
      this.pre_basisIdx := parent.basisIdx_ptr;
      this.pre_nbIdx := parent.nbIdx_ptr;
      this.pre_nf_pos := parent.nf_pos_ptr;
      this.pre_invB := parent.invB_ptr;
      this.pre_transMat := parent.transMat_ptr;
      this.pre_transRed := parent.transRed_ptr;
    end get_parent;

    procedure get_cur
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      this.p_sol := cur.p_sol;
      this.d_sol := cur.d_sol;
      this.redVec := cur.redVec;
      this.basisIdx := cur.basisIdx;
      this.nf_pos := cur.nf_pos;
      this.nbIdx := cur.nbIdx;
      this.rIdx := cur.rIdx;
      this.invB := cur.invB;
      this.transMat := cur.transMat;
      this.transRed := cur.transRed;
      this.pivOutList := cur.pivOutList;
      this.pivOutCheck := cur.pivOutCheck;
      this.pivOutNum := cur.pivOutNum;
    end get_cur;

    procedure get_res ( this : in Link_to_simplex; 
                        iData : in demics_ftest.class_ftData.Link_to_ftData
                      ) is
    begin
      iData.cur.artV := this.artV;
      this.artV := 0;
      iData.cur.nbN := this.nbN;
    end get_res;

    procedure get_pivOutNum
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData) is
    begin
      cur.pivOutNum := this.pivOutNum;
    end get_pivOutNum;

    procedure get_nbN_nfN ( this : in Link_to_simplex;
                            ori_nbN : in integer32; ori_nfN : in integer32 ) is
    begin
      this.nbN := ori_nbN;
      this.nfN := ori_nfN;
    end get_nbN_nfN;

    procedure get_p_sol
                ( this : in Link_to_simplex;
                  ori_p_sol : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      this.p_sol := ori_p_sol;
    end get_p_sol;

    procedure get_d_sol
                ( this : in Link_to_simplex;
                  ori_d_sol : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      this.d_sol := ori_d_sol;
    end get_d_sol;

    procedure get_basisIdx
                ( this : in Link_to_simplex;
                  ori_basisIdx : in Standard_Integer_Vectors.Link_to_Vector
                ) is
    begin
      this.basisIdx := ori_basisIdx;
    end get_basisIdx; 

    procedure get_nf_pos
                ( this : in Link_to_simplex;
                  ori_nf_pos : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      this.nf_pos := ori_nf_pos;
    end get_nf_pos;

    procedure get_nbIdx
                ( this : in Link_to_simplex;
                  ori_nbIdx : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      this.nbIdx := ori_nbIdx;
    end get_nbIdx;

    procedure get_invB
                ( this : in Link_to_simplex;
                  ori_invB : in Standard_Floating_Vectors.Link_to_Vector ) is
    begin
      this.invB := ori_invB;
    end get_invB;

    procedure get_frIdx ( this : in Link_to_simplex;
                          ori_frIdx : in integer32 ) is
    begin
      this.frIdx := ori_frIdx;
    end get_frIdx;

    procedure copy_p1_d_sol
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      for i in 0..this.dim-1 loop
        this.p1_d_sol(i) := cur.d_sol(i);
      end loop;
    end copy_p1_d_sol;

    procedure copy_eye
                ( this : in Link_to_simplex;
                  cur : in demics_ftest.class_theData.Link_to_theData ) is
    begin
      for i in 0..this.dim*this.dim-1 loop
        cur.transMat(i) := this.eye(i);
      end loop;
    end copy_eye;

    procedure allocateAndIni
                ( this : in Link_to_simplex;
                  data : demics_input_data.class_dataSet.dataSet;
                  ori_firIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  seedNum : in integer32; ori_output : in integer32;
                  vrblvl : in integer32 := 0 ) is

      cnt,tmp_nDim,tmp_mDim : integer32;
     -- rand_max : double_float;

    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_simplex.allocateAndIni ...");
      end if;
      tmp_nDim := data.dim;
      tmp_mDim := data.dim;
      this.dim := data.dim;
      this.supN := data.supN;
      this.termSumNum := data.termSumNum;
      tmp_nDim := tmp_nDim + this.termSumNum - this.supN;
      this.firIdx := ori_firIdx;
      this.output := ori_output;
      this.termSet := data.termSet;    
      this.termStart := data.termStart;
      this.re_termStart := new Standard_Integer_Vectors.Vector(0..this.supN);
      this.re_termStart(0) := 0;
      this.p1_d_sol
        := new Standard_Floating_Vectors.Vector'(0..this.dim-1 => 0.0);
      this.ip  := new Standard_Integer_Vectors.Vector'(0..this.dim-1 => 0);
      this.weight := new Standard_Floating_Vectors.Vector(0..this.Dim-1);
      this.vol
        := new Standard_Floating_Vectors.Vector'(0..this.dim*this.dim-1 => 0.0);
      this.eye
        := new Standard_Floating_Vectors.Vector'(0..this.dim*this.dim-1 => 0.0);
      this.lifting
        := new Standard_Floating_Vectors.Vector(0..this.termSumNum-1);
      this.oriSupp := new Standard_Floating_VecVecs.VecVec(0..this.supN-1);
      for i in 0..this.supN-1 loop
        this.oriSupp(i)
          := new Standard_Floating_Vectors.Vector'
                   (0..this.dim*this.termSet(i)-1 => 0.0);
        for j in 0..this.dim-1 loop
          for k in this.termStart(i)..this.termStart(i)+this.termSet(i)-1 loop
            supp_in(this,i,j,k - this.termStart(i),
                    demics_input_data.class_dataSet.support_out(data,k, j));
          end loop;
        end loop;
      end loop;
      this.fst_d_sol
        := new Standard_Floating_Vectors.Vector'(0..tmp_mDim-1 => 0.0);
      this.aux_cvec
        := new Standard_Floating_Vectors.Vector'(0..tmp_nDim-1 => 0.0);
      this.dir := new Standard_Floating_Vectors.Vector'(0..tmp_mDim-1 => 0.0);
      this.fst_redVec
        := new Standard_Floating_Vectors.Vector'(0..data.termMax-1 => 0.0);
      this.tmp_newInvB := new Standard_Floating_Vectors.Vector(0..this.dim-1);
      this.tmp_transMat := new Standard_Floating_Vectors.Vector(0..this.dim-1);
      this.nIdx := new Standard_Integer_Vectors.Vector(0..2*tmp_nDim-1);
      Standard_Random_Numbers.Set_Seed(natural32(seedNum));
     -- rand_max was needed because of rand()/rand_max*10 was used ...
     -- rand_max := 2.0;  
     -- for i in 1..29 loop
     --   rand_max := 2.0*rand_max;
     -- end loop;
     -- rand_max := (rand_max - 1.0)*2.0 + 1.0;
      for i in 0..this.termSumNum-1 loop
        this.lifting(i) := 10.0*abs(Standard_Random_Numbers.Random);
      end loop;
      this.supp := new VecVec_of_supportSets(0..this.supN);
      for i in 0..this.supN-1 loop
        this.re_termStart(i+1) := data.termStart(i+1) - i - 1;
        this.supp(i) := new Array_of_supportSets(0..this.termSet(i)-1);
      end loop;
      for i in 0..this.supN-1 loop
        for j in 0..this.termSet(i)-1 loop
          this.supp(i)(j) := new supportSet'(new_supportSet);
          class_supportSet.allocSupp
            (this.supp(i)(j),data,i,j,this.lifting,vrblvl-1);
        end loop;
      end loop;
      this.supp(this.supN) := new Array_of_supportSets(0..0);
      this.supp(this.supN)(0) := new supportSet'(new_supportSet);
      class_supportSet.allocAux(this.supp(this.supN)(0),data,vrblvl-1);
      cnt := 0;
      for i in 0..this.supN-1 loop
        for j in 0..this.termSet(i)-2 loop
          this.nIdx(2*cnt) := i;
          this.nIdx(2*cnt+1) := j;
          cnt := cnt + 1;
        end loop;
      end loop;
      for i in 0..this.dim-1 loop
        this.aux_cvec(this.termSumNum - this.supN + i) := 1.0;
        this.nIdx(2*cnt) := this.supN;
        this.nIdx(2*cnt+1) := 0;
        this.eye(i*(this.dim+1)) := 1.0;
        cnt := cnt + 1;
      end loop;
      if this.output > 0 then
        put_line("----------------------------------");
        put("* Seed number = "); put(seedNum,1); new_line;
        put_line("* Lifting values for elements in each support set");
        cnt := 0;
        for i in 0..this.supN-1 loop
          put("S"); put(i+1,1); put(" : ");
          for j in 
              this.termStart(i)..this.termStart(i) + this.termSet(i)-1 loop
            put(this.lifting(cnt)); put(" ");
            cnt := cnt + 1;
          end loop;
          new_line;
        end loop;
        put_line("----------------------------------");
      end if;
    end allocateAndIni;

-- for relation table

    procedure tSolLP ( this : in Link_to_simplex;
                       iter : in out integer32;
                       mode : in integer32; flag : out integer32;
                       vrblvl : in integer32 := 0 ) is

      iter_1 : integer32 := 0;
      iter_2 : integer32 := 0;
      pivInIdx,sub_pivInIdx,redFlag,pivOutIdx,sub_pivOutIdx : integer32;
      redCost,theta : double_float;
      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.tSolLP, iter : ");
        put(iter,1); put_line(" ...");      
      end if;
     -- phase 1
      if vrblvl > 0 then
        info_basisIdx(this);
        info_nbIdx(this);
        info_invB(this);
        info_p_sol(this);
        info_d_sol(this);
        info_frIdx(this);
      end if;
      loop
        if vrblvl > 0 then
          put("----- Phase-1. Iter : "); put(iter_1,1);
          put_line(" -----");
        end if;
        reducedCost_tab_p1
          (this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
        if redflag = opt then
          reMakeNonBasisIdx_tab(this,vrblvl-1);
          if vrblvl > 0 then
            put_line("----- OPT.Phase-1 -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_dir(this);
            info_nf_pos(this);
          end if;
          exit;
        end if;
        ratioTest(this,redFlag,pivInIdx,--sub_pivInIdx,
                  pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        createNewBandN_tab
          (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
           theta,redCost,vrblvl-1);
        if vrblvl > 0 then
          info_basisIdx(this);
          info_nbIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_dir(this);
          info_nf_pos(this);
        end if;
        iter_1 := iter_1 + 1;
        if iter_1 > DEMiCs_Global_Constants.ITER then
          flag := DEMiCs_Global_Constants.ERROR_ITER;
          exit;
        end if;
      end loop;
     -- phase 2
      if (mode = DEMiCs_Global_Constants.TRIANGLE) and
         (flag = DEMiCs_Global_Constants.CONTINUE) then
        flag := checkFrIdx(this,vrblvl-1);
        if vrblvl > 0 then
          put_line("----- checkFrIdx -----");
          info_basisIdx(this);
          info_nbIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_dir(this);
          info_nf_pos(this);
        end if;
      end if;
      if flag = DEMiCs_Global_Constants.CONTINUE then
        loop
          if vrblvl > 0 then   
            put("----- Phase-2. Iter : "); put(iter_2,1);
            put_line(" -----");
          end if;
          reducedCost_tab(this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
          if redFlag = opt then
            flag := opt;
            if vrblvl > 0 then
              put_line("----- OPT.Phase-2 -----");
              info_basisIdx(this);
              info_nbIdx(this);
              info_invB(this);
              info_p_sol(this);
              info_d_sol(this);
              info_dir(this);
              info_nf_pos(this);
            end if;
            exit;
          end if;
          ratioTest_art(this,redFlag,pivInIdx,--sub_pivInIdx,
                        pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNBOUNDED then
            if vrblvl > 0 then
              put_line("----- UNB.Phase-2 -----");
	    end if;
            exit;
          end if;
          createNewBandN_tab
            (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,theta,redCost,
             vrblvl-1);
          iter_2 := iter_2 + 1;
          if vrblvl > 0 then
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_dir(this);
            info_nf_pos(this);
          end if;
          if iter_2 > DEMiCs_Global_Constants.ITER then
            flag := DEMiCs_Global_Constants.ERROR_ITER;
            exit;
          end if;
        end loop;
      end if;
      iter := iter_1 + iter_2;
    end tSolLP;

-- for phase 1 and 2

    procedure fSolLP ( this : in Link_to_simplex;
                       termS : in integer32; reTermS : in integer32;
                       iter : in out integer32; flag : out integer32;
                       vrblvl : in integer32 := 0 ) is

      iter_1 : integer32 := 0;
      iter_2 : integer32 := 0;
      pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx : integer32;
      redFlag : integer32;
      theta,redCost : double_float;
      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.fSolLP, iter : ");
        put(iter,1); put_line(" ...");      
      end if;
     -- phase 1
      if vrblvl > 0 then
        info_basisIdx(this);
        info_nbIdx(this);
        info_rIdx(this);
        info_invB(this);
        info_p_sol(this);
        info_d_sol(this);
        info_dir(this);
      end if;
      loop
        if vrblvl > 0 then
          put("----- Phase-1. Iter : "); put(iter_1,1);
          put_line(" -----");
        end if;
        reducedCost_p1(this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
        if redFlag = opt then
          reMakeNonBasisIdx(this,reTermS,vrblvl-1);
          if vrblvl > 0 then
            put_line("----- OPT.Phase-1 -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_rIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_dir(this);
            info_nf_pos(this);
          end if;
          exit;
        end if;
        ratioTest(this,redFlag,pivInIdx,--sub_pivInIdx,
                  pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        update_p1_d_sol(this,pivInIdx,sub_pivOutIdx,vrblvl-1);
        createNewBandN_p1(this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx, 
                          theta,redCost,termS,reTermS,vrblvl-1);
        if vrblvl > 0 then
          info_basisIdx(this);
          info_nbIdx(this);
          info_rIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_dir(this);
          info_nf_pos(this);
        end if;
        iter_1 := iter_1 + 1;
        if iter_1 > DEMiCs_Global_Constants.ITER then
          flag := DEMiCs_Global_Constants.ERROR_ITER;
          exit;
        end if;
      end loop;
     -- phase 2
      if flag = DEMiCs_Global_Constants.CONTINUE then
        loop
          if vrblvl > 0 then
            put("----- Phase-2. Iter : "); put(iter_2,1);
            put_line(" -----");
          end if;
          reducedCost(this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
          if redFlag = opt then
            flag := opt;
            if vrblvl > 0 then
              put_line("----- OPT.Phase-2 -----");
              info_basisIdx(this);
              info_nbIdx(this);
              info_rIdx(this);
              info_invB(this);
              info_p_sol(this);
              info_d_sol(this);
              info_dir(this);
              info_nf_pos(this);
            end if;
            exit;
          end if;
          ratioTest_art(this,redFlag,pivInIdx,--sub_pivInIdx,
                        pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNBOUNDED
           then exit;
          end if;
          createNewBandN(this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx, 
                         theta,redCost,termS,reTermS,vrblvl-1);
          iter_2 := iter_2 + 1;
          if vrblvl > 0 then
            info_basisIdx(this);
            info_nbIdx(this);
            info_rIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_dir(this);
            info_nf_pos(this);
          end if;
          if iter > DEMiCs_Global_Constants.ITER_BLAND then
            solLP_art_Bland
              (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,redFlag,
               theta,redCost,termS,reTermS,iter_2,flag,vrblvl-1);
	    exit;
          end if;
        end loop;
      end if;
      iter := iter_1 + iter_2;
    end fSolLP;

-- iCheck

    procedure fstRed_candIdx
                ( this : in Link_to_simplex;
                  curInif : in demics_iTest.class_inifData.Link_to_inifData;
                  mCandIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  vrblvl : in integer32 := 0 ) is

      use demics_itest.class_uData;

      num : integer32 := 0;
      idx : integer32;
      tmp_redCost : double_float;
      redCost : double_float := 1.0E+8;
      n_curr : Link_to_uData := curInif.fHead;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.fstRed_candIdx ...");
      end if;
      while n_curr /= null loop
        tmp_redCost := n_curr.red;
        idx := n_curr.supLab;
        this.fst_redVec(idx) := tmp_redCost;
        mCandIdx(num+1) := idx;
        if tmp_redCost < redCost then
          redCost := tmp_redCost;
          pivInIdx := idx;
          sub_pivInIdx := num;
        end if;
        n_curr := n_curr.fNext;
        num := num + 1;
      end loop;
      mCandIdx(0) := num;
    end fstRed_candIdx;

    procedure cal_redVec
                ( this : in Link_to_simplex;
                  termS : in integer32; reTermS : in integer32;
                  fst_pivInIdx : in integer32;
                  cur : in demics_ftest.class_theData.Link_to_theData;
                  vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.cal_redVec, termS : ");
        put(termS,1); put_line(" ...");
      end if;
      for i in 0..termS-1 loop
        if i /= fst_pivInIdx then
          cur.redVec(cnt+reTermS)
            := this.fst_redVec(i) - this.fst_redVec(fst_pivInIdx);
          cnt := cnt + 1;
        end if;
      end loop;
    end cal_redVec;

    function put_redCost
               ( this : Link_to_simplex; fst_pivInIdx : integer32 )
               return double_float is
    begin
      return (this.fst_redVec(fst_pivInIdx) - this.fst_redVec(this.repIdx));
    end put_redCost;

-- iCheck_art

    procedure solLP_art
                ( this : in Link_to_simplex; depth : in integer32;
                 -- idx_one : in integer32; fst_pivIn : in integer32;
                  preNbN : in integer32; termS : in integer32;
                  reTermS : in integer32; iter : in out integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,redFlag : integer32;
      redCost,theta : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.solLP_art ...");
      end if;
     -- phase 2
      if vrblvl > 0 then
        put_line("<< iSolLP_Art >>");
        info_basisIdx(this);
        info_nbIdx(this);
        info_invB(this);
        info_p_sol(this);
        info_d_sol(this);
        info_nf_pos(this);
        info_rIdx(this);
        info_dir(this);
      end if;
      elimArt(this,depth,preNbN,termS,reTermS,iter,vrblvl-1);
      loop
        if vrblvl > 0 then
          put("----- Iter : "); put(iter,1); put_line(" -----");
          info_basisIdx(this);
          info_nbIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_rIdx(this);
          info_dir(this);
        end if;
        reducedCost(this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
        if redFlag = opt then
          if vrblvl > 0 then
            put_line("----- OPT -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_rIdx(this);
            info_dir(this);
          end if;
          flag := opt;
          exit;
        end if;
        ratioTest_art(this,redFlag,pivInIdx,--sub_pivInIdx,
                      pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        if flag = DEMiCs_Global_Constants.UNBOUNDED
         then exit;
        end if;
        createNewBandN(this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx, 
                       theta,redCost,termS,reTermS,vrblvl-1);
        iter := iter + 1;
        if iter > DEMiCs_Global_Constants.ITER_BLAND then
          solLP_art_Bland
            (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
             redFlag,theta,redCost,termS,reTermS,iter,flag,vrblvl-1);
          exit;
        end if;
      end loop;
    end solLP_art;

    procedure solLP_art_Bland
                ( this : in Link_to_simplex;
                  pivInIdx : in out integer32; sub_pivInIdx : in out integer32;
                  pivOutIdx : in out integer32;
                  sub_pivOutIdx : in out integer32;
                  redFlag : out integer32; theta : in out double_float;
                  redCost : in out double_float; termS : in integer32;
                  reTermS : in integer32; iter : in out integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.solLP_art_Bland ...");
      end if;
      loop
        if vrblvl > 0 then -- #if DBG_MSOLLP
          put("----- Iter : "); put(iter,1); put_line(" -----");
          info_basisIdx(this);
          info_nbIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_rIdx(this);
          info_dir(this);
        end if;
        reducedCost_Bland
          (this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
        if redFlag = opt then
          if vrblvl > 0 then -- #if DBG_MSOLLP
            put_line("----- OPT -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_rIdx(this);
            info_dir(this);
          end if;
          flag := opt;
          exit;
        end if;
        ratioTest_art_Bland(this,redFlag,pivInIdx,--sub_pivInIdx,
                            pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        exit when (flag = DEMiCs_Global_Constants.UNBOUNDED);
        createNewBandN
          (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,theta,
           redCost,termS,reTermS,vrblvl-1);
        iter := iter + 1;
        if iter > DEMiCs_Global_Constants.ITER then
          flag := DEMiCs_Global_Constants.ERROR_ITER;
          info_redVec(this);
          info_dir(this);
          info_basisIdx(this);
          info_nbIdx(this);
          info_nf_pos(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          exit;
        end if;
      end loop;
    end solLP_art_Bland;

-- for mLP

    procedure solLP ( this : in Link_to_simplex; depth : in integer32;
                      fst_pivInIdx : in out integer32;
                      fst_sub_pivInIdx : in out integer32;
                      fst_redCost : in double_float;
                      mode : in integer32; termS : in integer32;
                      reTermS : in integer32; preNbN : in integer32;
                      iter : in out integer32; flag : out integer32;
                      vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;
      pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,redFlag : integer32;
      theta,redCost : double_float;
  
    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.solLP, depth : ");
        put(depth,1); put_line(" ...");      
      end if;
     -- phase 2
      initIter(this,mode,fst_pivInIdx,fst_sub_pivInIdx,fst_redCost,
               redFlag,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
               theta,redCost,termS,reTermS,preNbN,flag,vrblvl-1);
      iter := iter + 1;
     -- if(flag != CONTINUE)
     --    return (flag);
     -- else
      if flag = DEMiCs_Global_Constants.CONTINUE then
        loop
          if vrblvl > 0 then -- #if DBG_MSOLLP
            put("----- Iter : "); put(iter,1); put_line(" -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_rIdx(this);
            info_dir(this);
          end if;
          ratioTest_art(this,redFlag,pivInIdx,--sub_pivInIdx,
                        pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
          exit when (flag = DEMiCs_Global_Constants.UNBOUNDED);
          createNewBandN
            (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
             theta,redCost,termS,reTermS,vrblvl-1);
          reducedCost(this,pivInIdx, sub_pivInIdx, redCost,redFlag,vrblvl-1);
          if redFlag = opt then
            if vrblvl > 0 then -- #if DBG_MSOLLP
              put_line("----- OPT -----");
              info_basisIdx(this);
              info_nbIdx(this);
              info_invB(this);
              info_p_sol(this);
              info_d_sol(this);
              info_rIdx(this);
              info_dir(this);
            end if;
            flag := opt;
            exit;
          end if;
          iter := iter + 1;
          if iter > DEMiCs_Global_Constants.ITER_BLAND then
            solLP_Bland
              (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
               redFlag,theta,redCost,termS,reTermS,iter,flag,vrblvl-1);
            exit;
          end if;
        end loop;
      end if;
    end solLP;

    procedure solLP_Bland
                ( this : in Link_to_simplex;
                  pivInIdx : in out integer32;
                  sub_pivInIdx : in out integer32;
                  pivOutIdx : in out integer32;
                  sub_pivOutIdx : in out integer32;
                  redFlag : in out integer32; theta : in out double_float;
                  redCost : in out double_float; termS : in integer32;
                  reTermS : in integer32; iter : in out integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.solLP_Bland, iter : ");
        put(iter,1); put_line(" ...");      
      end if;
      loop
        if vrblvl > 0 then -- #if DBG_MSOLLP
          put("----- Iter : "); put(iter,1); put_line(" -----");
          info_basisIdx(this);
          info_nbIdx(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          info_rIdx(this);
          info_dir(this);
        end if;
        ratioTest_art_Bland(this,redFlag,pivInIdx,--sub_pivInIdx,
                            pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
        exit when (flag = DEMiCs_Global_Constants.UNBOUNDED);
        createNewBandN
          (this,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,theta,
           redCost,termS,reTermS,vrblvl-1);
        reducedCost_Bland(this,pivInIdx,sub_pivInIdx,redCost,redFlag,vrblvl-1);
        if redFlag = opt then
          if vrblvl > 0 then -- #if DBG_MSOLLP
            put_line("----- OPT -----");
            info_basisIdx(this);
            info_nbIdx(this);
            info_invB(this);
            info_p_sol(this);
            info_d_sol(this);
            info_rIdx(this);
            info_dir(this);
          end if;
          flag := opt;
          exit;
        end if;
        iter := iter + 1;
        if iter > DEMiCs_Global_Constants.ITER then
          flag := DEMiCs_Global_Constants.ERROR_ITER;
          info_redVec(this);
          info_dir(this);
          info_basisIdx(this);
          info_nbIdx(this);
          info_nf_pos(this);
          info_invB(this);
          info_p_sol(this);
          info_d_sol(this);
          exit;
        end if;
      end loop;
    end solLP_Bland;

    procedure initIter
                ( this : in Link_to_simplex; mode : in integer32;
                  fst_pivInIdx : in out integer32;
                  fst_sub_pivInIdx : in out integer32;
                  fst_redCost : in double_float; redFlag : out integer32;
                  pivInIdx : out integer32; sub_pivInIdx : out integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; redCost : out double_float;
                  termS : in integer32; reTermS : in integer32;
                  preNbN : in integer32; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.initIter, mode : ");
        put(mode,1); put_line(" ...");      
      end if;
      flag := -1;
      case mode is
        when DEMiCs_Global_Constants.ICHECK =>
          ratioTest_artFst
            (this,DEMiCs_Global_Constants.POSTHETA,fst_pivInIdx,
             --fst_sub_pivInIdx,
             pivOutIdx,sub_pivOutIdx,theta,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNBOUNDED then
            return;
          else
            redCost := fst_redCost;
            createNewBandN_iFst
              (this,fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
               theta,redCost,termS,reTermS,vrblvl-1);
            reducedCost_iFst
              (this,fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,--sub_pivOutIdx,
               redCost,--termS,
               reTermS,preNbN,redFlag,vrblvl-1);
            if redFlag = opt then
              if vrblvl < 0 then -- #if DBG_MSOLLP
                put_line("----- OPT -----");
                info_basisIdx(this);
                info_nbIdx(this);
                info_invB(this);
                info_p_sol(this);
                info_d_sol(this);
                info_dir(this);
              end if;
              flag := redFlag;
              return;
            else
              pivInIdx := fst_pivInIdx;
              sub_pivInIdx := fst_sub_pivInIdx;
              flag := DEMiCs_Global_Constants.CONTINUE;
              return;
            end if;
          end if;
        when DEMiCs_Global_Constants.MCHECK =>
          ratioTest_artFst
            (this,DEMiCs_Global_Constants.NEGTHETA,fst_pivInIdx,
             -- fst_sub_pivInIdx,
             pivOutIdx, sub_pivOutIdx, theta,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNBOUNDED then
            return;
          else
            redCost := fst_redCost;
            createNewBandN_mFst
              (this,fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
               theta,redCost,termS,reTermS,vrblvl-1);
            reducedCost_mFst
              (this,fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,--sub_pivOutIdx,
               redCost,redFlag,vrblvl-1);
            if redFlag = opt then
              if vrblvl > 0 then -- #if DBG_MSOLLP
                put_line("----- OPT -----");
                info_basisIdx(this);
                info_nbIdx(this);
                info_invB(this);
                info_p_sol(this);
                info_d_sol(this);
                info_dir(this);
              end if;
              flag := redFlag;
              return;
            else
              pivInIdx := fst_pivInIdx;
              sub_pivInIdx := fst_sub_pivInIdx;
              flag := DEMiCs_Global_Constants.CONTINUE;
              return;
            end if;
          end if;
        when others => null;
      end case;
    end initIter;

    procedure calMixedVol
                ( this : in Link_to_simplex;
                  lv : in demics_fTest.class_lvData.Array_of_lvData;
                  sp : in Standard_Integer_Vectors.Link_to_Vector;
                  supN : in integer32; vrblvl : in integer32 := 0 ) is

      ii,jj,idx,cnt,polyDim,fIdx : integer32;
      det : double_float;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.calMixedVol ...");
      end if;
      this.mixedCell := this.mixedCell + 1;
      if this.output /= 0
       then put("# "); put(this.mixedCell,1); put(" : ");
      end if;
      cnt := 0;
      for i in 0..supN-1 loop
        polyDim := lv(sp(i)).Node.parent.polyDim; 
        fIdx := lv(sp(i)).Node.parent.nodeLabel(0);
        jj := fIdx*this.dim;
        if this.output /= 0 then
          put(sp(i)+1,1); put(" : ( ");
          put(fIdx+1,1); put(" ");
        end if;
        for j in 0..polyDim-1 loop
          idx := lv(sp(i)).Node.parent.nodeLabel(j+1);
          ii := idx*this.dim;
          if this.output /= 0 then
            put(idx+1,1); put(" ");
          end if;
          for k in 0..this.dim-1 loop
            this.vol(this.dim*cnt+k)
              := this.oriSupp(sp(i))(k+ii) - this.oriSupp(sp(i))(k+jj);
          end loop;
          cnt := cnt + 1;
        end loop; 
        if this.output /= 0
         then put(") ");
        end if;
      end loop;
      lu(this,this.dim,this.vol,det);
      det := abs(det);
      this.mixedVol := this.mixedVol + det;
      if this.output /= 0 then
        new_line;
        put("Volume : "); put(det);
        new_line;
      end if;
    end calMixedVol;

    procedure lu ( this : in Link_to_simplex; n : in integer32;
                   a : in Standard_Floating_Vectors.Link_to_Vector;
                   det : out double_float ) is

      ii,ik,piv : integer32;
      t,u : double_float;

    begin
      det := 0.0;
      for k in 0..n-1 loop
        this.ip(k) := k;
        u := 0.0;                 
        for j in 0..n-1 loop
          t := abs(a(n*k+j));
          if t > u
           then u := t;
          end if;
        end loop;
        if u = 0.0
         then return;
        end if;
        this.weight(k) := 1.0/u;     
      end loop;
      det := 1.0;
      for k in 0..n-1 loop
        u := -1.0;
        piv := k;
        for i in k..n-1 loop
          ii := this.ip(i);            
          t := abs(a(n*ii+k))*this.weight(ii);
          if t > u 
           then u := t; piv := i;
          end if;
        end loop;
        ik := this.ip(piv);
        if piv /= k then
          this.ip(piv) := this.ip(k);
          this.ip(k) := ik;
          det := -det; 
        end if;
        u := a(n*ik+k);
        det := det*u; 
        if u = 0.0
         then return;
        end if;
        for i in k+1..n-1 loop
          ii := this.ip(i);
          a(n*ii+k) := a(n*ii+k)/u;
          t := a(n*ii+k);
          for j in k+1..n-1 loop
            a(n*ii+j) := a(n*ii+j) - t*a(n*ik+j);
          end loop;
        end loop;
      end loop;
    end lu;

    procedure matinv ( this : in Link_to_simplex; n : in integer32;
                       a : in Standard_Floating_Vectors.Link_to_Vector;
                       a_inv : in Standard_Floating_Vectors.Link_to_Vector;
                       det : out double_float ) is

      ii : integer32;
      t : double_float;

    begin
      lu(this,n,a,det);
      if det /= 0.0 then
        for k in 0..n-1 loop
          for i in 0..n-1 loop
            ii := this.ip(i);
           -- t := double_float(ii = k);
            if ii = k
             then t := 1.0;
             else t := 0.0;
            end if;
            for j in 0..i-1 loop
              t := t - a(n*ii+j)*a_inv(n*j+k);
            end loop;
            a_inv(n*i+k) := t;
          end loop;
          for i in reverse 0..n-1 loop
            t := a_inv(n*i+k);
            ii := this.ip(i);
            for j in i+1..n-1 loop
              t := t - a(n*ii+j)*a_inv(n*j+k);
            end loop;
            a_inv(n*i+k) := t/a(n*ii+i);
          end loop;
        end loop;
      end if;
    end matinv;

    function put_elem_supp ( this : Link_to_simplex;
                             lvl : integer32; idx : integer32;
                             row : integer32; col : integer32;
                             vrblvl : integer32 := 0 )
                           return double_float is
    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.put_elem_supp, lvl : ");
        put(lvl,1); put(" idx : "); put(idx,1); new_line;
        put("row : "); put(row,1);
        put(", col : "); put(col,1); put_line(" ...");
      end if;
      return class_supportSet.supMat_Out(this.supp(lvl)(idx),row,col);
    end put_elem_supp;

    procedure mult_elem_supp ( this : in Link_to_simplex;
                               lvl : in integer32; idx : in integer32;
                               row : in integer32; col : in integer32 ) is
    begin
      class_supportSet.supMat_neg(this.supp(lvl)(idx),row,col);
    end mult_elem_supp;

    procedure check_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  depth : in integer32 ) is

      nfPos,nfN,cnt : integer32;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      val,red : double_float;
      invB,d_sol : Standard_Floating_Vectors.Link_to_Vector;

    begin  
      put_line("----- << check_dirRed >> -----");
      invB := parent.invB;
      d_sol := parent.d_sol;
      nf_pos := parent.nf_pos;
      nfN := parent.nfN;
      put_line("[ Direction and Reduced Cost ]");
      for ii in depth+1..this.supN-1 loop
        cnt := 0;
        put("--- Support : "); put(ii+1,1); put_line(" ---");
        for k in 0..this.termSet(ii)-1 loop
          put(cnt+1,3); put(" : ");
          for j in 0..nfN-1 loop
            nfPos := nf_pos(j);
            val := 0.0;
            for i in 0..this.dim-1 loop
              val := val + invB(i+nfPos*this.dim)*supp_out(this,ii,i,k);
            end loop;
            if (val < DEMiCs_Global_Constants.PLUSZERO) and
               (val > DEMiCs_Global_Constants.MINUSZERO) then
              put("0 ");
            else
              put(val,4); put(" ");
            end if;
          end loop;
          val := 0.0;
          for i in 0..this.dim-1 loop
            val := val + d_sol(i)*supp_out(this,ii,i,k);
          end loop;
          red := this.lifting(this.termStart(ii)+k) - val;
          put(" : "); put(red,4); new_line;
          cnt := cnt + 1;
        end loop;
        new_line;
      end loop;
    end check_dirRed;

    procedure dbg_dirRed
                ( this : in Link_to_simplex;
                  parent : in demics_ftest.class_theData.Link_to_theData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData;
                  depth : in integer32 ) is

      nfPos,nfN,cnt : integer32;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      val,red,diff : double_float;
      invB,d_sol : Standard_Floating_Vectors.Link_to_Vector;
      n_curr : demics_iTest.class_uData.Link_to_uData;

    begin
      invB := parent.invB;
      d_sol := parent.d_sol;
      nf_pos := parent.nf_pos;
      nfN := parent.nfN;
      for ii in depth+1..this.supN-1 loop
        cnt := 0;
        n_curr := nextInif(ii).fHead;
        for k in 0..this.termSet(ii)-1 loop
          for j in 0..nfN-1 loop
            nfPos := nf_pos(j);
            val := 0.0;
            for i in 0..this.dim-1 loop
              val := val + invB(i+nfPos*this.dim)*supp_out(this,ii,i,k);
            end loop;
            diff := val - n_curr.dir(nfPos);
            if (diff > DEMiCs_Global_Constants.PLUSZERO) or
               (diff < DEMiCs_Global_Constants.MINUSZERO) then
              put_line("dbg_dirRed:  ERROR -- Direction!!");
            end if;
          end loop;
          val := 0.0;
          for i in 0..this.dim-1 loop
            val := val + d_sol(i)*supp_out(this,ii,i,k);
          end loop;
          red := this.lifting(this.termStart(ii)+k) - val;
          diff := red - n_curr.red;
          if (diff > DEMiCs_Global_Constants.PLUSZERO) or
             (diff < DEMiCs_Global_Constants.MINUSZERO) then
            put_line("dbg_dirRed:  ERROR -- Reduced Cost!!");
          end if;
          n_curr := n_curr.fNext;
          cnt := cnt + 1;
        end loop; 
      end loop;
    end dbg_dirRed;

    procedure info_mv ( this : in Link_to_simplex ) is
    begin
      put("# Mixed Cells : "); put(this.mixedCell); new_line;
      put("Mixed Volume : "); put(this.mixedVol); new_line;
    end info_mv;

    procedure info_allSup ( this : in Link_to_simplex ) is
    begin
      put_line("<< Support Set >>");
      for i in 0..this.supN-1 loop
        put("---- Level : "); put(i,1); put_line(" ----");
        for j in 0..this.termSet(i)-1 loop
          put("* FrIdx : "); put(j,1); new_line;
          class_supportSet.info_sup(this.supp(i)(j));
          new_line;
        end loop;
      end loop;
      put_line("-- AuxMat --");
      class_supportSet.info_sup(this.supp(this.supN)(0));
      new_line;
    end info_allSup;

    procedure info_allCostVec ( this : in Link_to_simplex ) is
    begin
      put_line("<< Cost Vector >>");
      for i in 0..this.supN-1 loop
        put("---- Level : "); put(i,1); put_line(" ----");
        for j in 0..this.termSet(i)-1 loop
          put("* FrIdx : "); put(j,1); new_line;
          class_supportSet.info_costVec(this.supp(i)(j));
          new_line;
        end loop;
      end loop;
      new_line;
    end info_allCostVec;

    procedure info_lifting ( this : in Link_to_simplex ) is

      counter : integer32 := 0;

    begin
      new_line;
      put_line("Lifting :");
      for i in 0..this.supN-1 loop
        for j in this.termStart(i)..this.termStart(i)+this.termSet(i)-1 loop
          put(this.lifting(counter)); put(" ");
          counter := counter + 1;
        end loop;
        new_line;
      end loop;
      new_line;
      for k in 0..this.supN-1 loop
        put("level : "); put(k,1); new_line;
        for j in 0..this.termSet(k)-1 loop
          put("free index : "); put(j,1); new_line;
          for i in this.termStart(k)..this.termStart(k)+this.termSet(k)-1 loop
            if i /= this.termStart(k) + j then
              put(this.lifting(i) - this.lifting(j + this.termStart(k)));
              put(" ");
            end if;
          end loop;
          new_line;
        end loop;
        new_line;
      end loop;
    end info_lifting;

    procedure info_simplexData ( this : in Link_to_simplex ) is
    begin
      info_invB(this);
      info_p_sol(this);
      info_d_sol(this);
      info_basisIdx(this);
      info_nbIdx(this); 
      info_nf_pos(this);
    end info_simplexData;

  end class_simplex;

end demics_simplex;
