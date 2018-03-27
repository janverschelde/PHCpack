with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

procedure Write_Number_of_Tasks ( file : in file_type; nt : in natural32 );

-- DESCRIPTION :
--   This procedure defines the layout of the message to write the
--   number of tasks used in a run of phc.
--   If nt equals zero, then no multitasking was used.
