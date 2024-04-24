#ifndef AETHER_FIELD_UTILITIES_H
#define AETHER_FIELD_UTILITIES_H

#include "symbol_export.h"

SHO_PUBLIC void scale_flow_to_inlet(double INLET_FLOW,const char *VORQ);
SHO_PUBLIC void calculate_blood_volume(const char *OUTPUT);

#endif /* AETHER_FIELD_UTILITIES_H */
