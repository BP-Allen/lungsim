#include "lymphatics.h"
#include "string.h"

void alveolar_capillary_flux_c(int *num_nodes);
void lymphatic_transport_c(const char *filename, int *filename_len);


void alveolar_capillary_flux(int num_nodes)
{
  alveolar_capillary_flux_c(&num_nodes);
}

void lymphatic_transport(const char *filename)
{
  int filename_len = (int)strlen(filename);
  lymphatic_transport_c(filename, &filename_len);
}
