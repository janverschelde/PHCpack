#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string.h>

#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"

#include "poly.h"

extern "C" void adainit(void);
extern "C" void adafinal(void);

void ada_read_homotopy(char* start_file, char* target_file, \
		PolySys& start_sys, PolySys& target_sys, PolySolSet& sols);
