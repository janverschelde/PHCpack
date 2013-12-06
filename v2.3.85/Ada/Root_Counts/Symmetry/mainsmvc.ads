with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

procedure mainsmvc ( nt : in natural32; infilename,outfilename : in string );

-- DESCRIPTION :
--   This procedure implements the dispatcher for the mvc tool
--   by choosing the menu and calling the selected driver.

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     the name of the input file;
--   outfilename    the name of the output file.
