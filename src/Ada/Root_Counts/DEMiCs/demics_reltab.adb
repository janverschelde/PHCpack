with Ada.text_io;                       use Ada.text_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with DEMiCs_Global_Constants;

package body demics_reltab is

  package body class_reltab is

    procedure get_init_triData
                ( this : in Link_to_reltab;
                  lab : in integer32; idx : in integer32 ) is
    begin
      null;
    end get_init_triData;

    procedure get_init_squData
                ( this : in Link_to_reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 ) is
    begin
      null;
    end get_init_squData;

    procedure init_data ( this : in Link_to_reltab ) is
    begin
      for i in 0..this.row*this.row-1 loop
        this.invB(i) := 0.0;
      end loop;
      for i in 0..this.col-1 loop
        this.p_sol(i) := 0.0;
      end loop;
    end init_data;

    procedure init_tri ( this : in Link_to_reltab;
                         lab : in integer32; idx : in integer32 ) is

      constNum : constant integer32 := this.termSet(lab) - 1;

    begin
      for j in 0..this.negIdx(0)-1 loop
        for i in 0..constNum-1 loop
          demics_simplex.class_simplex.mult_elem_supp
            (this.the_Simplex,lab,idx,this.negIdx(j+1),i);
        end loop;
      end loop;
    end init_tri;

    procedure init_squ ( this : in Link_to_reltab;
                         lab_a : in integer32; lab_b : in integer32;
                         idx_a : in integer32; idx_b : in integer32 ) is

      constNum_a : constant integer32 := this.termSet(lab_a) - 1;
      constNum_b : constant integer32 := this.termSet(lab_b) - 1;

    begin
      for j in 0..this.negIdx(0)-1 loop
        for i in 0..constNum_a-1 loop
          demics_simplex.class_simplex.mult_elem_supp
            (this.the_Simplex,lab_a,idx_a,this.negIdx(j+1),i);
        end loop;
        for i in 0..constNum_b-1 loop
          demics_simplex.class_simplex.mult_elem_supp
            (this.the_Simplex,lab_b,idx_b,this.negIdx(j+1),i);
        end loop;
      end loop;
    end init_squ;

    procedure put_data ( this : in Link_to_reltab ) is
    begin
      demics_simplex.class_simplex.get_nbN_nfN
        (this.the_Simplex,this.nbN,this.nfN);
      demics_simplex.class_simplex.get_p_sol(this.the_Simplex,this.p_sol);
      demics_simplex.class_simplex.get_d_sol(this.the_Simplex,this.d_sol);
      demics_simplex.class_simplex.get_basisIdx
        (this.the_Simplex,this.basisIdx);
      demics_simplex.class_simplex.get_nf_pos(this.the_Simplex,this.nf_pos);
      demics_simplex.class_simplex.get_nbIdx(this.the_Simplex,this.nbIdx);
      demics_simplex.class_simplex.get_invB(this.the_Simplex,this.invB);
    end put_data;

    procedure put_frIdx ( this : in Link_to_reltab; frIdx : in integer32 ) is
    begin
      demics_simplex.class_simplex.get_frIdx(this.the_Simplex,frIdx);
    end put_frIdx;

    procedure makeTri ( this : in Link_to_reltab;
                        vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_reltab.class_reltab.makeTri ...");
      end if;
    end makeTri;

    procedure makeSqu ( this : in Link_to_reltab;
                        vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_reltab.class_reltab.makeSqu ...");
      end if;
    end makeSqu;

    procedure findAllFeasLPs_tri
                ( this : in Link_to_reltab;
                  lab : in integer32; idx : in integer32 ) is
                 -- frIdx : in integer32 ) is

      reTermS : constant integer32 := this.re_termStart(lab);
      termSta : constant integer32 := this.termStart(lab);
      bIdx,tmp_idx : integer32;
      fIdx,fIdx_a,fIdx_b : integer32;
      num : integer32 := 0;
      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      for i in 0..this.dim-1 loop
        bIdx := this.basisIdx(i);
        if bIdx < this.termSumNum - this.supN then
          tmp_idx := bIdx - reTermS;
          if tmp_idx >= idx
           then fIdx := tmp_idx + 1;
           else fIdx := tmp_idx;
          end if;
          table_in(this,termSta +  idx,termSta + fIdx,opt);
          table_in(this,termSta + fidx,termSta +  idx,opt);
          this.feasIdx_a(num) := fIdx;
          num := num + 1;
        end if;
      end loop;
      for j in 0..num-1 loop
        fIdx_b := this.feasIdx_a(j);
        for i in j+1..num-1 loop
          fIdx_a := this.feasIdx_a(i);
          table_in(this,termSta + fIdx_b,termSta + fIdx_a,opt);
          table_in(this,termSta + fIdx_a,termSta + fIdx_b,opt);
        end loop;
      end loop;
    end findAllFeasLPs_tri;

    procedure findAllFeasLPs_squ
                ( this : in Link_to_reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 ) is

      constNum_a : constant integer32 := this.termSet(lab_a) - 1;
     -- constNum_b : constant integer32 := this.termSet(lab_b) - 1;
      reTermS_a : constant integer32 := this.re_termStart(lab_a);
      reTermS_b : constant integer32 := this.re_termStart(lab_b);
      bIdx,tmp_idx : integer32;
      fIdx,fIdx_a,fIdx_b : integer32;
      num_a : integer32 := 0;
      num_b : integer32 := 0;
      opt : constant := DEMiCs_Global_Constants.OPT;

    begin
      table_in(this,colPos + idx_a,rowPos + idx_b,opt);
      table_in(this,rowPos + idx_b,colPos + idx_a,opt);
      for i in 0..this.dim-1 loop
        bIdx := this.basisIdx(i);
        if bIdx < this.termSumNum - this.supN then
          if (reTermS_a <= bIdx) and (bIdx < reTermS_a + constNum_a) then
            tmp_idx := bIdx - reTermS_a;
            if tmp_idx >= idx_a
             then fIdx := tmp_idx + 1;
             else fIdx := tmp_idx;
            end if;
            table_in(this,colPos + fIdx,  rowPos + idx_b,opt);
            table_in(this,rowPos + idx_b, colPos + fIdx, opt);
            this.feasIdx_a(num_a) := fIdx;
            num_a := num_a + 1;
          else
            tmp_idx := bIdx - reTermS_b;
            if tmp_idx >= idx_b
             then fIdx := tmp_idx + 1;
             else fIdx := tmp_idx;
            end if;
            table_in(this,colPos + idx_a, rowPos + fIdx, opt);
            table_in(this,rowPos + fIdx,  colPos + idx_a,opt);
            this.feasIdx_b(num_b) := fIdx;
            num_b := num_b + 1;
          end if;
        end if;
      end loop;
      for j in 0..num_b-1 loop
        fIdx_b := this.feasIdx_b(j);
        for i in 0..num_a-1 loop
          fIdx_a := this.feasIdx_a(i);
          table_in(this,colPos + fIdx_a, rowPos + fIdx_b, opt);
          table_in(this,rowPos + fIdx_b, colPos + fIdx_a, opt);
        end loop;
      end loop;
    end findAllFeasLPs_squ;

    procedure table_in ( this : in Link_to_reltab;
                         row : in integer32; col : in integer32;
                         elem : in integer32 ) is
    begin
      this.table(row + col*this.termSumNum) := elem;
    end table_in;

    function table_out ( this : in Link_to_reltab;
                         row : integer32; col : integer32 )
                       return integer32 is
    begin
      return this.table(row + col*this.termSumNum);
    end table_out;

    procedure info_invB ( this : in Link_to_reltab ) is
    begin
      put_line("<< invB >>");
      for i in 0..this.dim-1 loop
        for j in 0..this.dim-1 loop
          put(invB_out(this,i,j)); put(" ");
        end loop;
        new_line;
      end loop;
    end info_invB;

    procedure info_p_sol ( this : in Link_to_reltab ) is
    begin
      put_line("<< p_sol >>");
      for i in 0..this.maxConst + this.dim-1 loop
        put(this.p_sol(i)); put(" ");
      end loop;
      new_line;
    end info_p_sol;

    procedure info_d_sol ( this : in Link_to_reltab ) is
    begin
      put_line("<< d_sol >>");
      for i in 0..this.dim-1 loop
        put(this.d_sol(i)); put(" ");
      end loop;
      new_line;
    end info_d_sol;

    procedure info_basisIdx ( this : in Link_to_reltab ) is
    begin
      put_line("<< basisIdx >>");
      for i in 0..this.dim-1 loop
        put(this.basisIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_basisIdx;

    procedure info_nbIdx ( this : in Link_to_reltab ) is
    begin
      put_line("<< nbIdx >>");
      for i in 0..this.dim-1 loop
        put(this.nbIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_nbIdx;

    procedure info_nf_pos ( this : in Link_to_reltab ) is
    begin
      put_line("<< nf_pos >>");
      for i in 0..this.dim-1 loop
        put(this.nf_pos(i),1); put(" ");
      end loop;
      new_line;
    end info_nf_pos;

    procedure info_feasIdx_tri ( this : in Link_to_reltab;
                                 num : in integer32 ) is
    begin
      put("feasIdx :");
      for i in 0..num-1 loop
        put(" "); put(this.feasIdx_a(i),1);
      end loop;
      new_line;
    end info_feasIdx_tri;

    procedure info_feasIdx_squ
                ( this : in Link_to_reltab;
                  num_a : in integer32; num_b : in integer32 ) is
    begin
      put("feasIdx_a :");
      for i in 0..num_a-1 loop
        put(" "); put(this.feasIdx_a(i),1);
      end loop;
      new_line;
      put("feasIdx_b :");
      for i in 0..num_b-1 loop
        put(" "); put(this.feasIdx_b(i),1);
      end loop;
      new_line;
    end info_feasIdx_squ;

    procedure info_allTable ( this : in Link_to_reltab ) is

      unbNum : integer32 := 0;

    begin
      put_line("<< all elements in Relation Table >>");
      for j in 0..this.termSumNum-1 loop
        for i in 0..this.termSumNum-1 loop
          put(table_out(this,j,i),1); put(" ");
          if table_out(this,j,i) = DEMiCs_Global_Constants.UNBOUNDED
           then unbNum := unbNum + 1;
          end if;
        end loop;
        new_line;
      end loop;
      new_line;
      put("# Unb. LPs : "); put(unbNum/2,1); new_line;
    end info_allTable;

    procedure info_table ( this : in Link_to_reltab ) is

      u_cnt : integer32 := 0;
      t_cnt : integer32 := 0;

    begin
      put_line("<< Relation Table >>");
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

    function invB_out ( this : Link_to_reltab;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return this.invB(colIdx + this.dim*rowIdx);
    end invB_out;

    function new_reltab return reltab is

      res : reltab;

    begin
      res.dim := 0;
      res.supN := 0;
      res.maxConst := 0;
      res.termSumNum := 0;
      res.row := 0;
      res.col := 0;
      res.unbLP := 0.0;
      res.totalLP := 0.0;
      res.nbN := 0;
      res.nfN := 0;
      res.termSet := null;
      res.termStart := null;
      res.firIdx := null;
      res.invB := null;
      res.p_sol := null;
      res.d_sol := null;
      res.basisIdx := null;
      res.nbIdx := null;
      res.nf_pos := null;
      res.negIdx := null;
      res.val := null;
      res.feasIdx_a := null;
      res.feasIdx_b := null;
      res.table := null;
      return res;
    end new_reltab;

    procedure delete_reltab ( this : in Link_to_reltab ) is
    begin
      Standard_Floating_Vectors.clear(this.invB);
      Standard_Floating_Vectors.clear(this.p_sol);
      Standard_Floating_Vectors.clear(this.d_sol);
      Standard_Integer_Vectors.clear(this.basisIdx);
      Standard_Integer_Vectors.clear(this.nbIdx);
      Standard_Integer_Vectors.clear(this.nf_pos);
      Standard_Integer_Vectors.clear(this.negIdx);
      Standard_Floating_Vectors.clear(this.val);
      Standard_Integer_Vectors.clear(this.feasIdx_a);
      Standard_Integer_Vectors.clear(this.feasIdx_b);
      Standard_Integer_Vectors.clear(this.table);
    end delete_reltab;

    procedure allocateAndIni
                ( this : in Link_to_reltab;
                  ori_Simplex
                    : in demics_simplex.class_simplex.Link_to_simplex;
                  ori_firIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_dim : in integer32;
                  ori_supN : in integer32;
                  ori_termSumNum : in integer32;
                  ori_termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_re_termStart
                    : in Standard_Integer_Vectors.Link_to_Vector;
                  vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_reltab.class_reltab.allocateAndIni ...");
      end if;
      this.dim := ori_Dim;
      this.supN := ori_supN;
      this.termSumNum := ori_termSumNum;
      this.maxConst := this.termSumNum - this.supN;
      this.termSet := ori_termSet;
      this.termStart := ori_termStart;
      this.re_termStart := ori_re_termStart;
      this.firIdx := ori_firIdx;
      this.the_Simplex := ori_Simplex;
      this.row := ori_Dim;
      this.col := this.maxConst + this.dim;
      this.invB := new Standard_Floating_Vectors.Vector'
                         (0..this.row*this.row-1 => 0.0);
      this.p_sol := new Standard_Floating_Vectors.Vector'
                          (0..this.col-1 => 0.0);
      this.d_sol := new Standard_Floating_Vectors.Vector'
                          (0..this.row-1 => 0.0);
      this.basisIdx
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.nbIdx
        := new Standard_Integer_Vectors.Vector'(0..this.col-1 => 0);
      this.nf_pos
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.val := new Standard_Floating_Vectors.Vector(0..this.col-1);
      this.negIdx := new Standard_Integer_Vectors.Vector(0..this.row);
      this.feasIdx_a
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.feasIdx_b
        := new Standard_Integer_Vectors.Vector'(0..this.row-1 => 0);
      this.table := new Standard_Integer_Vectors.Vector'
                          (0..this.termSumNum*this.termSumNum-1 => 0);
    end allocateAndIni;

    procedure makeTable ( this : in Link_to_reltab;
                          total_unbLP_tab : out double_float;
                          vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0
       then put_line("-> in demics_reltab.class_reltab.makeTable ...");
      end if;
      makeTri(this,vrblvl-1);
      makeSqu(this,vrblvl-1);
      total_unbLP_tab := this.unbLP;
    end makeTable;

  end class_reltab;

end demics_reltab;
