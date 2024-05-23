#include "calpgm.h"
#include <math.h>
#include <stdlib.h>

int comparator(const void *line_1, const void *line_2) {
  double freq_1 = ((SXLINE *)line_1)->xfrq;
  double freq_2 = ((SXLINE *)line_2)->xfrq;
  if (freq_1 < freq_2) {
    return -1;
  } else if (freq_1 > freq_2) {
    return 1;
  } else {
    return 0;
  }
}

void sort_lines(SXLINE *lines, int num_lines) {
  qsort(lines, num_lines, sizeof(SXLINE), comparator);
}

int nearest_calc_freq(SXLINE *lines, int num_lines, double freq) {
  int i;
  int nearest;
  double min_dif = 1e10;
  for (i = 0; i < num_lines; i++) {
    double dif = fabs(lines[i].cfrq - freq);
    if (dif < min_dif) {
      min_dif = dif;
      nearest = i;
    }
  }
  return nearest;
}
