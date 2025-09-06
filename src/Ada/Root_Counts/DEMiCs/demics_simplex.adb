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
                  level : in integer32;
                  num : in integer32;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  vrblvl : in integer32 := 0 ) is

      cnt : integer32 := 0;

      use demics_input_data.class_dataSet; -- for support_out

    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_supportSet.allocSupp, level : ");
        put(level,1); put_line(" ...");      
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
                          rowIdx : in integer32;
                          colIdx : in integer32;
                          elem : in double_float ) is
    begin
      this.supMat(rowIdx + colIdx*this.row) := elem;
    end supMat_in;

    procedure supMat_neg ( this : in Link_to_supportSet;
                           rowIdx : in integer32;
                           colIdx : in integer32 ) is
    begin
      this.supMat(rowIdx + colIdx*this.row) 
        := -this.supMat(rowIdx + colIdx*this.row);
    end supMat_neg;

    function supMat_out ( this : Link_to_supportSet;
                          rowIdx : integer32;
                          colIdx : integer32 ) return double_float is
    begin
      return this.supMat(rowIdx + colIdx*this.row);
    end supMat_out;

    function redVal ( this : Link_to_supportSet;
                      d_sol : Standard_Floating_Vectors.Link_to_Vector;
                      idx : integer32;
                      ii : integer32 ) return double_float is

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
    begin
      if vrblvl > 0
       then put_line("-> in demics_simplex.class_simplex.checkFrIdx ...");
      end if;
      return 0;
    end checkFrIdx;

    procedure elimFrIdx ( this : in Link_to_simplex;
                          sub_pivOutIdx : in integer32 ) is
    begin
      null;
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
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
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
    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.reducedCost_tab ...");
      end if;
    end reducedCost_tab;

    procedure reducedCost_p1
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_p1;

    procedure reducedCost
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost;

    procedure reducedCost_Bland
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_Bland;

    procedure reducedCost_mFst
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
                  pivOutIdx : in integer32;
                  sub_pivOutIdx : in integer32;
                  redCost : out double_float; flag : out integer32 ) is
    begin
      null;
    end reducedCost_mFst;

    procedure reducedCost_iFst
                ( this : in Link_to_simplex;
                  pivInIdx : out integer32;
                  sub_pivInIdx : out integer32;
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

    procedure getIdx ( this : in Link_to_simplex; level : out integer32;
                       idx : out integer32; idx2 : out integer32;
                       ii : out integer32; d_nbIdx : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.getIdx, level : ");
        put(level,1); put_line(" ...");
      end if;
      level := this.nIdx(d_nbIdx);
      idx := this.nIdx(d_nbIdx+1);
      idx2 := this.firIdx(level);
      ii := this.dim*idx;
    end getIdx;

-- ratio test

    procedure ratioTest
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest, redFlag : ");
        put(redFlag,1); put_line(" ...");
      end if;
    end ratioTest;

    procedure ratioTest_artFst
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32 ) is
    begin
      null;
    end ratioTest_artFst;

    procedure ratioTest_art
                ( this : in Link_to_simplex; redFlag : in integer32;
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : out integer32; sub_pivOutIdx : out integer32;
                  theta : out double_float; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_simplex.class_simplex.ratioTest_art, redFlag : ");
        put(redFlag,1); put_line(" ...");
      end if;
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
                  pivInIdx : in integer32; sub_pivInIdx : in integer32;
                  pivOutIdx : in integer32; sub_pivOutIdx : in integer32;
                  theta : in double_float; redCost : in double_float;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put_line("-> in demics_simplex.class_simplex.createNewBandN_tab ...");
      end if;
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
      for i in 0..this.nbN-1 loop
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
    begin
      null;
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
                   (0..this.dim*this.termSet(i) => 0.0);
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
        put("-> in demics_simplex.class_supportSet.tSolLP, iter : ");
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
        ratioTest(this,redFlag,pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
                  theta,flag,vrblvl-1);
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
          ratioTest_art(this,redFlag,pivInIdx,sub_pivInIdx,
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
