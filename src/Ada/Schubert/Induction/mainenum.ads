with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure mainenum ( nt : in natural32; verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the routine for the homotopies in enumerative geometry,
--   as called by the option handlers of the phc executable.
--
-- ON ENTRY :
--   nt       the number of tasks, 0 for no multitasking;
--   vrblvl   is the verbose level.
