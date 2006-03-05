//
// window-extract extracts windows.
//
// basically, it's spaghetti code for selecting out "interesting" windows
// from paircomp/find_patch-style 10bp comparisons.
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <MutableComparison.hh>
#include <ImmutableComparison.hh>
using namespace paircomp;

#define WINDOWSIZE 100
#define CUTOFF 15
#define DIAGONAL_TOLERANCE 10
#define OFFDIAG_CUTOFF 3

int cmpfn(const void * ap, const void * bp)
{
  int a = *((int *)ap);
  int b = *((int *)bp);

  return (a - b);
}

//
// main
//

int main(int argc, char ** argv)
{
  if (argc != 5) {
    fprintf(stderr, "Usage:\n\t%s top_len bot_len cmpfile outfile\n\n", argv[0]);
    exit(-1);
  }
  unsigned int top_len = atoi(argv[1]);
  unsigned int bot_len = atoi(argv[2]);
  char * cmpfile = argv[3];
  char * outfile = argv[4];

  MutableComparison * c = load_paircomp_file(top_len, bot_len, 10, cmpfile);
  ImmutableComparison * cmp = new ImmutableComparison(c);
  MutableComparison * out_cmp = new MutableComparison(top_len, bot_len, 10);

  //
  // Build a list indexed by top position that, for each position,
  // contains a list of matching windows sorted by bottom position.
  //

  int * bot_lists[top_len + 1];

  unsigned int top;
  matches_iterator iter;
  for (top = 0; top < top_len; top++) {
    const _MatchContainer * cont = cmp->get_matches(top);

    if (cont == NULL) {
      bot_lists[top] = NULL;
      continue;
    }

    int * bot_pos_list = (int *) malloc(2*sizeof(int) * (cont->num + 1));

    unsigned int i;
    for (i = 0; i < cont->num; i++) {
      const Match * m = &cont->block[i];
      bot_pos_list[2*i] = m->get_bot_pos();
      bot_pos_list[2*i + 1] = m->get_n_matching();
    }
    bot_pos_list[2*i] = -1;

    qsort(bot_pos_list, cont->num, 2*sizeof(int), cmpfn);

    bot_lists[top] = bot_pos_list;
  }

  //
  // now, iterate over the whole kit & caboodle and just count.
  //

  printf("starting...\n");

  int top_start, bot_start, bot_i;
  for (top_start = 0; top_start < top_len - WINDOWSIZE; top_start++) {
    if (bot_lists[top_start] != NULL) {
      int * bot_list = bot_lists[top_start];
      for (bot_i = 0; bot_list[2*bot_i] != -1; bot_i++) {
	bot_start = bot_list[2*bot_i];

	//
	// Now, for this top_start, bot_start pair count how many
	// matches are within a square of size WINDOWSIZE starting
	// at those coordinates.
	//
	int window_count = 0, offdiag_count = 0;
	int i, j;

	for (i = top_start; i < top_start + WINDOWSIZE; i++) {
	  if (bot_lists[i] != NULL) {
	    int * b = bot_lists[i];
	    for (j = 0; b[2*j] != -1; j++) {
	      if (b[2*j] >= bot_start && b[2*j] < bot_start + WINDOWSIZE) {
		window_count++;
	      }
	    }
	  }
	}

	// If it's over the cutoff for significance, then count how
	// many of the matches are within/without the given diagonal
	// tolerance.  This allows for gaps etc.
	
	if (window_count >= CUTOFF) {
	  window_count = offdiag_count = 0;
	  for (i = 0; i < WINDOWSIZE; i++) {
	    if (bot_lists[i + top_start] != NULL) {
	      int * b = bot_lists[i + top_start];
	      for (j = 0; b[2*j] != -1; j++) {
		if (b[2*j] >= bot_start && b[2*j] < bot_start + WINDOWSIZE) {
		  if (abs(b[2*j] - bot_start - i) <= DIAGONAL_TOLERANCE) {
		    window_count++;
		  } else {
		    offdiag_count++;
		  }
		}
	      }
	    }
	  }

	  // Now, finally, if we have a good window without too many
	  // off-diagonal matches, save.

	  if (window_count >= CUTOFF && offdiag_count < OFFDIAG_CUTOFF) {
	    for (i = 0; i < WINDOWSIZE; i++) {
	      if (bot_lists[i + top_start] != NULL) {
		int * b = bot_lists[i + top_start];
		for (j = 0; b[2*j] != -1; j++) {
		  if (abs(b[2*j] - bot_start - i) <= DIAGONAL_TOLERANCE) {
		    Match * m = new Match(i + top_start, b[2*j],
					  10, b[2*j + 1], 1);
		    out_cmp->add_match(m, false);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  ImmutableComparison * out_c = new ImmutableComparison(out_cmp);
  out_c->save_as_paircomp(outfile);

  exit(0);
}
