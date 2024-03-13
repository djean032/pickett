#include <stdio.h>
struct cat_entry {
  float freq;
  float err;
  float lgint;
  int dr;
  float elo;
  int gup;
  int tag;
  int qnfmt;
  char *qn;
};

int main(int argc, char *argv[]) {
  FILE *cat;
  char contents[100];
  char *cat_name = argv[1];
  struct cat_entry cat_entry;
  sprintf(cat_name, "%s.cat", argv[1]);
  cat = fopen(cat_name, "r");
  if (!cat) {
    printf("Error: cannot open %s\n", cat_name);
    return 1;
  }
  while (fgets(contents, 100, cat)) {
    printf("%12s", &contents[55]);
  }
}
