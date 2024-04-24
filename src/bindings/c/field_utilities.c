
#include "field_utilities.h"
#include "utils.h"

#include <string.h>

void scale_flow_to_inlet_c(double *INLET_FLOW,const char *VORQ,int *VQLEN);
void calculate_blood_volume_c(const char *OUTPUT, int *OUTLEN);

void scale_flow_to_inlet(double INLET_FLOW,const char *VORQ)
{
	int vq_len = strlen(VORQ);
	scale_flow_to_inlet_c(&INLET_FLOW, VORQ, &vq_len);
}

void calculate_blood_volume(const char *OUTPUT)
{
  int out_len = (int)strlen(OUTPUT);
  calculate_blood_volume_c(OUTPUT, &out_len);
}


