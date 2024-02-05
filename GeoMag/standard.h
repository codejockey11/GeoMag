#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI    ((2)*(acos(0.0)))
#endif

#define RAD2DEG(rad)    ((rad)*(180.0L/M_PI))
#define DEG2RAD(deg)    ((deg)*(M_PI/180.0L))