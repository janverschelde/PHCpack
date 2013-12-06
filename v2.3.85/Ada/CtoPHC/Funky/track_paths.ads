with C_Integer_Arrays,C_Double_Arrays;  use C_Integer_Arrays,C_Double_Arrays;

procedure Track_Paths ( n,m : in integer;
                        p_mc : in C_intarrs.Pointer;
                        p_ns : in integer; p_s : in C_intarrs.Pointer;
                        p_nc : in integer; p_c : in C_dblarrs.Pointer;
                        q_mc : in C_intarrs.Pointer;
                        q_ns : in integer; q_s : in C_intarrs.Pointer;
                        q_nc : in integer; q_c : in C_dblarrs.Pointer;
                        nbmul : in integer; mul : in C_intarrs.Pointer;
                        nbcfs : in integer; cfs : in C_dblarrs.Pointer );

-- DESCRIPTION :
--   Calls the Ada path tracking routines.