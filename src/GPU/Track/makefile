# The code has been tested on Linux and Windows.
# Swap the # symbol in the definition of the MAKEFILE below
# to use the makefile_windows or the makefile_unix,
# respectively for a Windows or GNU Linux/Unix computer.

# MAKEFILE = makefile_windows
MAKEFILE = makefile_unix

all: test_dbl_track test_dbl2_track test_dbl4_track test_dbl8_track

test_dbl_track:
	make -f $(MAKEFILE) test_dbl_track

test_dbl2_track:
	make -f $(MAKEFILE) test_dbl2_track

test_dbl4_track:
	make -f $(MAKEFILE) test_dbl4_track

test_dbl8_track:
	make -f $(MAKEFILE) test_dbl8_track

clean:
	make -f $(MAKEFILE) clean

cleanall:
	make -f $(MAKEFILE) cleanall
