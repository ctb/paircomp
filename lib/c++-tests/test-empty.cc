#include "paircomp.hh"

using namespace paircomp;

int main(int argc, char * argv[])
{
  system("../../bin/paircomp ../../tests/el.txt ../../tests/br.txt 20 .7 a");

  MutableComparison * c1 = load_paircomp_file(352, 357, 20, "a");
  ImmutableComparison * c1b = new ImmutableComparison(c1);
  ImmutableComparison * c2b = c1b->filter_by_threshold(.9);

  assert(c1b->contains(*c2b));
  assert(!c2b->contains(*c1b));
  assert(!c1b->equals(*c2b));
  
  ImmutableComparison * c3b = c2b->subtract(*c1b);
  assert(c3b->is_empty());

  ImmutableComparison * c4b = c1b->subtract(*c2b);
  assert(!c4b->is_empty());

  printf("All tests OK.\n");
}
