with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;

procedure ts_mtpermtrol is

-- DESCRIPTION :
--   Development of controlling the permanent computation with multitasking.

  procedure Main is

  -- DESCRPTION :
  --   Prompts the user for the dimension of a matrix.

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
  end Main;

begin
  Main;
end ts_mtpermtrol;
