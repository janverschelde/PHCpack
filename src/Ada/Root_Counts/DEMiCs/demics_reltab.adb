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
      null;
    end init_data;

    procedure init_tri ( this : in Link_to_reltab;
                         lab : in integer32; idx : in integer32 ) is
    begin
      null;
    end init_tri;

    procedure init_squ ( this : in Link_to_reltab;
                         lab_a : in integer32; lab_b : in integer32;
                         idx_a : in integer32; idx_b : in integer32 ) is
    begin
      null;
    end init_squ;

    procedure put_data ( this : in Link_to_reltab ) is
    begin
      null;
    end put_data;

    procedure put_frIdx ( this : in Link_to_reltab; frIdx : in integer32 ) is
    begin
      null;
    end put_frIdx;

    procedure makeTri ( this : in Link_to_reltab ) is
    begin
      null;
    end makeTri;

    procedure makeSqu ( this : in Link_to_reltab ) is
    begin
      null;
    end makeSqu;

    procedure findAllFeasLPs_tri
                ( this : in Link_to_reltab;
                  lab : in integer32; idx : in integer32;
                  frIdx : in integer32 ) is
    begin
      null;
    end findAllFeasLPs_tri;

    procedure findAllFeasLPs_squ
                ( this : in Link_to_reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 ) is
    begin
      null;
    end findAllFeasLPs_squ;

    procedure table_in ( this : in Link_to_reltab;
                         row : in integer32; col : in integer32;
                         elem : in integer32 ) is
    begin
      null;
    end table_in;

    function table_out ( this : in Link_to_reltab;
                         row : integer32; col : integer32 )
                       return integer32 is
    begin
      return 0;
    end table_out;

    procedure info_invB ( this : in Link_to_reltab ) is
    begin
      null;
    end info_invB;

    procedure info_p_sol ( this : in Link_to_reltab ) is
    begin
      null;
    end info_p_sol;

    procedure info_d_sol ( this : in Link_to_reltab ) is
    begin
      null;
    end info_d_sol;

    procedure info_basisIdx ( this : in Link_to_reltab ) is
    begin
      null;
    end info_basisIdx;

    procedure info_nbIdx ( this : in Link_to_reltab ) is
    begin
      null;
    end info_nbIdx;

    procedure info_nf_pos ( this : in Link_to_reltab ) is
    begin
      null;
    end info_nf_pos;

    procedure info_feasIdx_tri ( this : in Link_to_reltab;
                                 num : in integer32 ) is
    begin
      null;
    end info_feasIdx_tri;

    procedure info_feasIdx_squ
                ( this : in Link_to_reltab;
                  num_a : in integer32; num_b : in integer32 ) is
    begin
      null;
    end info_feasIdx_squ;

    procedure info_allTable ( this : in Link_to_reltab ) is
    begin
      null;
    end info_allTable;

    procedure info_table ( this : in Link_to_reltab ) is
    begin
      null;
    end info_table;

    function invB_out ( this : Link_to_reltab;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float is
    begin
      return 0.0;
    end invB_out;

    function new_reltab return reltab is

      res : reltab;

    begin
      return res;
    end new_reltab;

    procedure delete_reltab ( this : in Link_to_reltab ) is
    begin
      null;
    end delete_reltab;

    procedure allocateAndIni
                ( this : in Link_to_reltab;
                  ori_Simplex
                    : in demics_simplex.class_simplex.Link_to_simplex;
                  ori_firIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  ori_dim : in integer32;
                  ori_supN : in integer32;
                  ori_termSumNum : in integer32;
                  ori_termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_re_termStart
                    : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      null;
    end allocateAndIni;

    procedure makeTable ( this : in Link_to_reltab;
                          total_unbLP_tab : out double_float ) is
    begin
      null;
    end makeTable;

  end class_reltab;

end demics_reltab;
