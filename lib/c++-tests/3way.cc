#include <stdio.h>
#include <stdlib.h>
#include "paircomp.hh"

paircomp::ImmutableComparison * load_paircomp(char * filename,
					      int top, int bot, int win)
{
  FILE * fp = fopen(filename, "r");
  fseek(fp, 0, SEEK_END);
  unsigned int filesize = ftell(fp);
  rewind(fp);

  char buf[filesize + 1];
  fread(buf, 1, filesize, fp);
  buf[filesize] = 0;

  paircomp::MutableComparison c(top, bot, win);
  c.parse_paircomp_format(buf);

  return new paircomp::ImmutableComparison(&c);
}


int main(int argc, char * argv[])
{
  paircomp::ImmutableComparison * ab, * bc, * ac;

  ab = load_paircomp("tmp.cmp", 40, 40, 21);
  bc = load_paircomp("tmp.cmp", 40, 40, 21);
  ac = load_paircomp("tmp.cmp", 40, 40, 21);

  printf("loaded.\n");

  paircomp::ImmutableComparison * new_ab, * new_bc, * new_ac;
  ab->filter_transitively(*bc, *ac, &new_ab, &new_bc, &new_ac);
}
