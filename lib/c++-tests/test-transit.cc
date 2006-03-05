#include "paircomp.hh"

using namespace paircomp;

void test_3way(std::string A, std::string B, std::string C)
{
  ImmutableComparison * ab = rolling_nxn_comparison(A, B, 20, .9);
  ImmutableComparison * bc = rolling_nxn_comparison(B, C, 20, .9);
  ImmutableComparison * ac = rolling_nxn_comparison(A, C, 20, .9);

  ImmutableComparison *new_ab, *new_bc, *new_ac;

  ab->filter_transitively(*bc, *ac, &new_ab, &new_bc, &new_ac);

  assert(!new_ab->is_empty());
  assert(!new_bc->is_empty());
  assert(!new_ac->is_empty());

  assert(ab->contains(*new_ab));
  assert(bc->contains(*new_bc));
  assert(ac->contains(*new_ac));

}

int main(int argc, char * argv[])
{
  std::string all_A = "AAAAAAAAAAAAAAAAAAAA";
  std::string all_T = "TTTTTTTTTTTTTTTTTTTT";

  printf("%s x %s x %s\n", all_A.c_str(), all_T.c_str(), all_T.c_str());
  test_3way(all_A, all_T, all_T);
  printf("%s x %s x %s\n", all_A.c_str(), all_T.c_str(), all_A.c_str());
  test_3way(all_A, all_T, all_A);
  printf("%s x %s x %s\n", all_A.c_str(), all_A.c_str(), all_T.c_str());
  test_3way(all_A, all_A, all_T);

  printf("All tests passed ok.\n");
}
