with String_Splitters;
with Generic_Lists;

package Lists_of_Strings is
  new Generic_Lists(String_Splitters.Link_to_String);

-- DESCRIPTION :
--   A list of strings is represented as a list of pointers,
--   each pointer in the list points to a string.
