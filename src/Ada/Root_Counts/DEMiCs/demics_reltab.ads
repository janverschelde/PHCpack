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
      the_Simplex : demics_simplex.class_simplex.Link_to_simplex;
      table : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    type Link_to_reltab is access reltab;

    procedure get_init_triData
                ( this : in Link_to_reltab;
                  lab : in integer32; idx : in integer32 );

    procedure get_init_squData
                ( this : in Link_to_reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 );

    procedure init_data ( this : in Link_to_reltab );

    procedure init_tri ( this : in Link_to_reltab;
                         lab : in integer32; idx : in integer32 );

    procedure init_squ ( this : in Link_to_reltab;
                         lab_a : in integer32; lab_b : in integer32;
                         idx_a : in integer32; idx_b : in integer32 );

    procedure put_data ( this : in Link_to_reltab );

    procedure put_frIdx ( this : in Link_to_reltab; frIdx : in integer32 );

    procedure makeTri ( this : in Link_to_reltab );

    procedure makeSqu ( this : in Link_to_reltab );

    procedure findAllFeasLPs_tri
                ( this : in Link_to_reltab;
                  lab : in integer32; idx : in integer32;
                  frIdx : in integer32 );

    procedure findAllFeasLPs_squ
                ( this : in Link_to_reltab;
                  lab_a : in integer32; lab_b : in integer32;
                  idx_a : in integer32; idx_b : in integer32;
                  colPos : in integer32; rowPos : in integer32 );

    procedure table_in ( this : in Link_to_reltab;
                         row : in integer32; col : in integer32;
                         elem : in integer32 );

    function table_out ( this : in Link_to_reltab;
                         row : integer32; col : integer32 ) return integer32;

    procedure info_invB ( this : in Link_to_reltab );

    procedure info_p_sol ( this : in Link_to_reltab );

    procedure info_d_sol ( this : in Link_to_reltab );

    procedure info_basisIdx ( this : in Link_to_reltab );

    procedure info_nbIdx ( this : in Link_to_reltab );

    procedure info_nf_pos ( this : in Link_to_reltab );

    procedure info_feasIdx_tri ( this : in Link_to_reltab;
                                 num : in integer32 );

    procedure info_feasIdx_squ
                ( this : in Link_to_reltab;
                  num_a : in integer32; num_b : in integer32 );

    procedure info_allTable ( this : in Link_to_reltab );

    procedure info_table ( this : in Link_to_reltab );

    function invB_out ( this : Link_to_reltab;
                        rowIdx : integer32; colIdx : integer32 )
                      return double_float; 

    function new_reltab return reltab;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_reltab ( this : in Link_to_reltab );

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
                  vrblvl : in integer32 := 0 );

    -- NOTE :
    --   The ori_firIdx was declare as int** which could indicate a
    --   vector of vectors data structure, but then in the allocateAndIni
    --   of mvc.cpp the argument was called as &firIdx.
    --   Thus, as the int** was immediately dereferenced by the &,
    --   the argument refers to a one dimensional data type.

    procedure makeTable ( this : in Link_to_reltab;
                          total_unbLP_tab : out double_float );

  end class_reltab;

end demics_reltab;
