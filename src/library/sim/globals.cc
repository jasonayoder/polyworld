/********************************************************************/
/* PolyWorld:  An Artificial Life Ecological Simulator              */
/* by Larry Yaeger                                                  */
/* Copyright Apple Computer 1990,1991,1992                          */
/********************************************************************/
// globals.cp - global variables for pw (PolyWorld)

// Self
#include "globals.h"

float   globals::worldsize;		// agent, food
float   globals::worldwidth;
float   globals::worlddepth;
float   globals::camerax;
float   globals::cameray;
float   globals::cameraz;

bool	globals::wraparound;	// agent, sceneview, simulation
bool	globals::blockedEdges;	// agent, simulation
bool	globals::stickyEdges;
int     globals::numEnergyTypes;
AbstractFile::ConcreteFileType globals::recordFileType;

