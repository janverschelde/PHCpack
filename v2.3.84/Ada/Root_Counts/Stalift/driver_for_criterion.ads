with text_io;                            use text_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

procedure Driver_for_Criterion
             ( file : in file_type; points : in out Array_of_Lists );

-- DESCRIPTION :
--   Allows to apply a criterion to sweep out the points that do not
--   contribute to the mixed volume of the tuple of point sets.

-- ON ENTRY :
--   file      must be opened for output;
--   points    tuple of point sets.

-- ON RETURN :
--   points    possibly reduced tuple of point sets.
