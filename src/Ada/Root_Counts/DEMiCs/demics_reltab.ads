with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with demics_simplex;

package demics_reltab is

-- DESCRIPTION :
--   Defines the relation table.
--   Translation initiated by g++ -c -fdump-ada-spec reltab.h.

  package class_reltab is

    type reltab is record
      dim : integer32;
      supN : integer32;
      maxConst : integer32;
      termSumNum : integer32;
      termSet : Standard_Integer_Vectors.Link_to_Vector;
      termStart : Standard_Integer_Vectors.Link_to_Vector;
      re_termStart : Standard_Integer_Vectors.Link_to_Vector;
      firIdx : Standard_Integer_Vectors.Link_to_Vector;
      unbLP : double_float;
      totalLP : double_float;
      row : integer32;
      col : integer32;
      nbN : integer32;
      nfN : integer32;
      invB : Standard_Floating_Vectors.Link_to_Vector;
      p_sol : Standard_Floating_Vectors.Link_to_Vector;
      d_sol : Standard_Floating_Vectors.Link_to_Vector;
      basisIdx : Standard_Integer_Vectors.Link_to_Vector;
      nbIdx : Standard_Integer_Vectors.Link_to_Vector;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      negIdx : Standard_Integer_Vectors.Link_to_Vector;
      val : Standard_Floating_Vectors.Link_to_Vector;
      feasIdx_a : Standard_Integer_Vectors.Link_to_Vector;
      feasIdx_b : Standard_Integer_Vectors.Link_to_Vector;
      the_Simplex : access demics_simplex.class_simplex.simplex;
      table : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    procedure get_init_triData
                ( this : access reltab;
                  lab : in integer32; idx : in integer32 );

    procedure get_init_squData
                ( this : access reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 );

    procedure init_data ( this : access reltab );

    procedure init_tri
                ( this : access reltab;
                  lab : in integer32; idx : in integer32 );

    procedure init_squ
                ( this : access reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32 );

    procedure put_data ( this : access reltab );

    procedure put_frIdx ( this : access reltab; frIdx : in integer32 );

    procedure makeTri ( this : access reltab );

    procedure makeSqu ( this : access reltab );

    procedure findAllFeasLPs_tri
                ( this : access reltab;
                  lab : in integer32; idx : in integer32;
                  frIdx : in integer32 );

    procedure findAllFeasLPs_squ
                ( this : access reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 );

    procedure table_in
                ( this : access reltab;
                  row : in integer32; col : in integer32;
                  elem : in integer32 );

    function table_out
                ( this : access reltab;
                  row : integer32; col : integer32 ) return integer32;

    procedure info_invB ( this : access reltab );

    procedure info_p_sol ( this : access reltab );

    procedure info_d_sol ( this : access reltab );

    procedure info_basisIdx ( this : access reltab );

    procedure info_nbIdx ( this : access reltab );

    procedure info_nf_pos ( this : access reltab );

    procedure info_feasIdx_tri ( this : access reltab; num : in integer32 );

    procedure info_feasIdx_squ
                ( this : access reltab;
                  num_a : in integer32; num_b : in integer32 );

    procedure info_allTable ( this : access reltab );

    procedure info_table ( this : access reltab );

    function invB_out
               ( this : access reltab;
                 rowIdx : integer32; colIdx : integer32 )
               return double_float; 

    function new_reltab return reltab;

    procedure delete_reltab ( this : access reltab );

    procedure allocateAndIni
                ( this : access reltab;
                  ori_Simplex : access demics_simplex.class_simplex.simplex;
                  ori_firIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  ori_dim : in integer32;
                  ori_supN : in integer32;
                  ori_termSumNum : in integer32;
                  ori_termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  ori_re_termStart
                    : in Standard_Integer_Vectors.Link_to_Vector );

    procedure makeTable ( this : access reltab;
                          total_unbLP_tab : out double_float );

  end class_reltab;

end demics_reltab;
