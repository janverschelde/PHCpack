# The code has been tested on Linux and Windows.
# Swap the # symbol in the definition of the MAKEFILE below
# to use the makefile_windows or the makefile_unix,
# respectively for a Windows or GNU Linux/Unix computer.

# MAKEFILE = makefile_windows2
# MAKEFILE = makefile_windows
MAKEFILE = makefile_unix

all: test_dbl_newton test_dbl2_newton test_dbl4_newton test_dbl8_newton

test_unimodular:
	make -f $(MAKEFILE) test_unimodular

test_columns:
	make -f $(MAKEFILE) test_columns

test_monomial_systems:
	make -f $(MAKEFILE) test_monomial_systems

test2_monomial_systems:
	make -f $(MAKEFILE) test2_monomial_systems

test4_monomial_systems:
	make -f $(MAKEFILE) test4_monomial_systems

test8_monomial_systems:
	make -f $(MAKEFILE) test8_monomial_systems

test_monomial_jobs:
	make -f $(MAKEFILE) test_monomial_jobs

test_dbl_cyclic:
	make -f $(MAKEFILE) test_dbl_cyclic

test_dbl2_cyclic:
	make -f $(MAKEFILE) test_dbl2_cyclic

test_dbl4_cyclic:
	make -f $(MAKEFILE) test_dbl4_cyclic

test_dbl8_cyclic:
	make -f $(MAKEFILE) test_dbl8_cyclic

test_dbl_flopbals:
	make -f $(MAKEFILE) test_dbl_flopbals

test_dbl2_flopbals:
	make -f $(MAKEFILE) test_dbl2_flopbals

test_dbl4_flopbals:
	make -f $(MAKEFILE) test_dbl4_flopbals

test_dbl8_flopbals:
	make -f $(MAKEFILE) test_dbl8_flopbals

test_dbl_newton:
	make -f $(MAKEFILE) test_dbl_newton

test_dbl2_newton:
	make -f $(MAKEFILE) test_dbl2_newton

test_dbl4_newton:
	make -f $(MAKEFILE) test_dbl4_newton

test_dbl8_newton:
	make -f $(MAKEFILE) test_dbl8_newton

clean:
	make -f $(MAKEFILE) clean

cleanall:
	make -f $(MAKEFILE) cleanall
