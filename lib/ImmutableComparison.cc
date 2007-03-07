/*!
*   \file ImmutableComparison.cc
*/
#include "ImmutableComparison.hh"
#include "algorithms.hh"

using namespace paircomp;

// utility fn to check parameters.

/*!
*   \fn static void _check_comparison_compatibility(const ImmutableComparison * c1,
					    const ImmutableComparison * c2)
*   \brief A utility function to check that the comparison parameters match, throws an exception when incompatible.
*   \param c1 A pointer to the first ImmutableComparison object.
*   \param c2 A pointer to the second ImmutableComparison object
*   \exception paircomp_exception top lengths do not match; were these comparisons done on the same top sequence?
*   \exception paircomp_exception bottom lengths do not match; were these comparisons done on the same bottom sequence?
*   \exception paircomp_exception window sizes do not match
*/
static void _check_comparison_compatibility(const ImmutableComparison * c1,
					    const ImmutableComparison * c2)
{
  if (c1->get_top_length() != c2->get_top_length()) {
    throw paircomp_exception("top lengths do not match; were these comparisons done on the same top sequence?");
  }
  if (c1->get_bottom_length() != c2->get_bottom_length()) {
    throw paircomp_exception("bottom lengths do not match; were these comparisons done on the same bottom sequence?");
  }
  if (c1->get_windowsize() != c2->get_windowsize()) {
    throw paircomp_exception("windowsizes do not match");
  }
}





/*!
*   \fn static int _compare_matches_by_strength(const void * p1, const void * p2)
*   \brief Comparse the number of matches between the two Match objects.
*   \param p1 A void pointer to the first Match object for comparison.
*   \param p2 A void pointer to the second Match object for comparison.
*   \return An integer representing the difference in the number of matches between p1 and p2.
*/
static int _compare_matches_by_strength(const void * p1, const void * p2)
{
  const Match * m1 = (const Match *) p1;
  const Match * m2 = (const Match *) p2;

  int i1 = m1->get_n_matching();
  int i2 = m2->get_n_matching();

  return i2 - i1;		// sort in descending.
}



/*!
*   \fn static void _sort_matches(_MatchContainer * cont)
*   \brief Performs a quicksort on the given _MatchContainer object.
*   \param cont A pointer to a _MatchContainer object to be sorted.
*   
*/
static void _sort_matches(_MatchContainer * cont)
{
  qsort(cont->block, cont->num, sizeof(Match), _compare_matches_by_strength);
}

/*!
*   \var const std::string forward_match("ACGTN")
*   \brief A string object with an ordering of characters for a forward match.
*/
const std::string forward_match("ACGTN");

/*!
*   \var const std::string reverse_complement_match("TGCAN")
*   \brief A string object with an ordering of characters for a reverse match.
*/
const std::string reverse_complement_match("TGCAN");


/*!
*   \fn bool test_forward_match(char ch1, char ch2)
*   \brief Function for testing a forward match between two characters.
*   \param ch1 The first character representing a nucelotide.
*   \param ch2 The second character representing a nucelotide.
*   \return True if ch1 and ch2 match and ch1 is not 'N', false otherwise.
*   \exception paircomp_exception test_forward_match failed
*/
bool test_forward_match(char ch1, char ch2)
{
  std::string::size_type i = forward_match.find(ch1);
  if (i == std::string::npos) {
    throw paircomp_exception("test_forward_match failed");
  }

  std::string::size_type j = forward_match.find(ch2);
  if (j == std::string::npos)  {
    throw paircomp_exception("test_forward_match failed");
  }

  if (i == j && ch1 != 'N') return true;
  return false;
}


/*!
*   \fn bool test_reverse_complement_match(char ch1, char ch2)
*   \brief Function for testing a reverse match between two characters.
*   \param ch1 The first character representing a nucelotide.
*   \param ch2 The second character representing a nucelotide.
*   \return True if ch1 and ch2 match and ch1 is not 'N', false otherwise.
*   \exception pair_comp test_reverse_complement_match failed
*    
*
*/
bool test_reverse_complement_match(char ch1, char ch2)
{
  std::string::size_type i = forward_match.find(ch1);
  if (i == std::string::npos)  {
    throw paircomp_exception("test_reverse_complement_match failed");
  }

  std::string::size_type j = reverse_complement_match.find(ch2);
  if (i == std::string::npos)  {
    throw paircomp_exception("test_reverse_complement_match failed");
  }

  if (i == j && ch1 != 'N') return true;
  return false;
}



//
// ImmutableComparison: Class implementation.
//

// public constructor: take a MutableComparison & transform it.

ImmutableComparison::ImmutableComparison(MutableComparison * cmp)
  : Comparison(cmp->get_top_length(),
	       cmp->get_bottom_length(),
	       cmp->get_windowsize())
{
  std::map<unsigned int, std::vector<Match*> >::const_iterator iter;

  // Loop through & transfer all of the matches into a fixed-size
  // container.

  for (iter = cmp->_matches.begin(); iter != cmp->_matches.end(); iter++) {
    unsigned int pos = iter->first;
    std::vector<Match *> v = iter->second;

    _MatchContainer * cont = new _MatchContainer;

    cont->block = (Match *) calloc(sizeof(Match), v.size());
    cont->num = v.size();

    for (unsigned int j = 0; j < v.size(); j++) {
      Match * m = v[j];
      memcpy(&cont->block[j], m, sizeof(Match));
    }

    // order by match strength, inverted.
    _sort_matches(cont);

    _matches[pos] = cont;
  }
}

ImmutableComparison::~ImmutableComparison()
{
  matches_iterator iter = _matches.begin();

  while (iter != _matches.end()) {
    _MatchContainer * cont = iter->second;

    free(cont->block);
    delete cont;

    // be sure not to touch cont after this...

    iter++;
  }
}

bool ImmutableComparison::contains_match(const Match &m) const
{
  matches_iterator iter;
  iter = _matches.find(m.get_top_pos());

  if (iter == _matches.end()) {	// no matches at that pos.
    return false;
  }

  const _MatchContainer * cont = iter->second;

  // Iterate over all matches at this position, checking to see if they're
  // equal.

  for (unsigned int i = 0; i < cont->num; i++) {
    const Match * this_m = (const Match *) &cont->block[i];
    if (this_m->equals(&m)) {
      return true;
    }
  }

  return false;
}


/******************************************************************************/
ImmutableComparison * ImmutableComparison::filter_by_threshold(float threshold)
  const
{
  matches_iterator iter;
  unsigned int ithreshold = (unsigned int)
    (threshold * (float) _windowsize + 0.5); // round up...

  ImmutableComparison * new_cmp = new ImmutableComparison(_top_length,
							  _bottom_length,
							  _windowsize);

  // Loop through & transfer all matches above threshold into a new
  // container.

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    unsigned int pos = iter->first;
    _MatchContainer * old = iter->second;

    // find the watermark where n_matches drops below threshold:
    unsigned int j;
    for (j = 0; j < old->num; j++) {
      Match * m = &old->block[j];
      if (m->get_n_matching() < ithreshold) break;
    }

    if (j == 0) {		// none matching!
      continue;
    }

    // build a _MatchContainer and populate it with matches.
    _MatchContainer * young = new _MatchContainer;

    young->num = j;
    young->block = (Match *) calloc(sizeof(Match), j);
    memcpy(&young->block[0], &old->block[0], j * sizeof(Match));

    new_cmp->_matches[pos] = young;
  }

  return new_cmp;
}


/******************************************************************************/
ImmutableComparison * ImmutableComparison::reverse_top_matches() const
{
  MutableComparison new_cmp(_top_length, _bottom_length, _windowsize);

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    _MatchContainer * cont = iter->second;

    for (unsigned int i = 0; i < cont->num; i++) {
      Match * m = &cont->block[i];
      Match * new_m = m->copy();
      new_m->reverse_top(_top_length, _bottom_length);

      new_cmp.add_match(new_m);
    }
  }
  
  return new ImmutableComparison(&new_cmp);
}



/******************************************************************************/
ImmutableComparison * ImmutableComparison::reverse_bot_matches() const
{
  MutableComparison new_cmp(_top_length, _bottom_length, _windowsize);

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    _MatchContainer * cont = iter->second;

    for (unsigned int i = 0; i < cont->num; i++) {
      Match * m = &cont->block[i];
      Match * new_m = m->copy();
      new_m->reverse_bottom(_bottom_length);

      new_cmp.add_match(new_m);
    }
  }
  
  return new ImmutableComparison(&new_cmp);
}


/******************************************************************************/
//
// 
//
/******************************************************************************/
ImmutableComparison * ImmutableComparison::isolate_matching_bases(
					     std::string top,
					     std::string bot)
  const
{
  MutableComparison new_cmp(_top_length, _bottom_length, 1);

  // do some quick checking... remember to account for windowsize!
  if (top.length() != _top_length) {
    throw paircomp_exception("lengths do not match");
  }

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    _MatchContainer * cont = iter->second;

    for (unsigned int i = 0; i < cont->num; i++) {
      Match * m = &cont->block[i];

      unsigned int top_pos, bot_pos, len, n_orig;
      int o;

      top_pos = m->get_top_pos();
      bot_pos = m->get_bot_pos();
      len = m->get_length();
      o = m->get_orientation();
      n_orig = m->get_n_matching();

      for (unsigned int j = 0; j < len; j++) {

	//
	// test for single-bp matches, and create a new match obj for each one.
	//
	
	if ((o == 1 && test_forward_match(top[top_pos + j], bot[bot_pos + j]))
	    ||
    (o == -1 &&
          test_reverse_complement_match(top[top_pos + j], bot[bot_pos - j]))) {

	  Match * new_m = new Match(top_pos + j, bot_pos + o*j, 1, 1, o);
	  new_cmp.add_match(new_m, false);
	}
      }
    }
  }

  return new ImmutableComparison(&new_cmp);
}



/******************************************************************************/
ImmutableComparison * ImmutableComparison::invert() const
{
  MutableComparison inverted(_bottom_length, _top_length,
			     _windowsize);

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    const _MatchContainer * cont = iter->second;

    for (unsigned int i = 0; i < cont->num; i++) {
      const Match * old_m = &cont->block[i];
      Match * new_m = old_m->copy();
      new_m->invert(_top_length, _bottom_length);
      inverted.add_match(new_m);
    }
  }

  return new ImmutableComparison(&inverted);
}



/******************************************************************************/
ImmutableComparison * ImmutableComparison::subtract(const 
						    ImmutableComparison &other)
  const
{
  _check_comparison_compatibility(this, &other);

  MutableComparison subtracted(_top_length, _bottom_length,
			       _windowsize);

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    const _MatchContainer * this_cont = iter->second;

    // iterate over all matches in this_cont
    for (unsigned int i = 0; i < this_cont->num; i++) {
      Match * this_m = &this_cont->block[i];

      // If other doesn't contain this match, copy it.
      if (!other.contains_match(*this_m)) {
	     Match * copy = this_m->copy();
	     subtracted.add_match(copy);
      }
    }
  }
  
  return new ImmutableComparison(&subtracted);
}



/******************************************************************************/
bool ImmutableComparison::is_empty() const
{
  if (_matches.begin() == _matches.end()) {
    return true;
  }
  return false;
}



/******************************************************************************/
bool ImmutableComparison::contains(const ImmutableComparison &other)
  const
{
  _check_comparison_compatibility(this, &other);

  ImmutableComparison * sub = other.subtract(*this);
  bool answer;
  
  if (sub->is_empty()) {
    answer = true;
  } else {
    answer = false;
  }
 
  delete sub;
  return answer;
}



/******************************************************************************/
bool ImmutableComparison::equals(const ImmutableComparison &other)
  const
{
  _check_comparison_compatibility(this, &other);

  return (this->contains(other) && other.contains(*this));
}



/******************************************************************************/
//
// save_as_seqcomp -- save to a given file in Tristan-style b3.5 output.
//

void ImmutableComparison::save_as_seqcomp(char * filename)
  const
{
  FILE * fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    throw ComparisonFileException("file cannot be created");
  }

  fprintf(fp, "b3.5\n%d\n", _windowsize);

  matches_iterator iter;

  unsigned int pos = 0;
  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    unsigned int cur = iter->first;
    _MatchContainer * cont = iter->second;

    while (pos < cur) {
      fprintf(fp, "0.0\n");
      pos ++;
    }
    pos ++;

    Match * m;
    for (unsigned int j = 0; j < cont->num; j++) {
      m = &cont->block[j];

      fprintf(fp, "%d%c %d ",
	      m->get_bot_pos(),
	      m->get_orientation() == 1 ? '+' : '-',
	      m->get_n_matching());
    }
    fprintf(fp, "0.0\n");
  }

  while (pos < _top_length - _windowsize + 1) {
    fprintf(fp, "0.0\n");
    pos ++;
  }

  fclose(fp);
}




/******************************************************************************/
//
// save_as_paircomp -- saves a comparison in the (simple) paircomp format.
//

void ImmutableComparison::save_as_paircomp(char * filename)
  const
{
  FILE * fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    throw ComparisonFileException("file cannot be created");
  }

  matches_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    _MatchContainer * cont = iter->second;

    Match * m;
    for (unsigned int j = 0; j < cont->num; j++) {
      m = &cont->block[j];

      fprintf(fp, "%d\t%d\t%d\t%d\n",
	      m->get_top_pos(), m->get_bot_pos(), m->get_n_matching(),
	      m->get_orientation());
    }
  }

  fclose(fp);
}




/******************************************************************************/
ImmutableComparison *
ImmutableComparison::intersect(const ImmutableComparison &other)
  const
{
  _check_comparison_compatibility(this, &other);

  MutableComparison cmp(_top_length, _bottom_length, _windowsize);

  // easy to do: just go through & grab.
  matches_iterator iter1, iter2;

  for (iter1 = _matches.begin(); iter1 != _matches.end(); iter1++) {
    const unsigned int top_pos = iter1->first;
    const _MatchContainer * cont1 = iter1->second;

    iter2 = other._matches.find(top_pos);

    if (iter2 == other._matches.end()) { // no matches at this pos.
      continue;
    }

    // ok, now take & save all intersecting matches.
    const _MatchContainer * cont2 = iter2->second;
    for (unsigned int i = 0; i < cont1->num; i++) {
      const Match * this_m = (const Match *) &cont1->block[i];

      for (unsigned int j = 0; j < cont2->num; j++) {
	const Match * other_m = (const Match *) &cont2->block[j];

	if (other_m->equals(this_m)) {
	  Match * new_m = other_m->copy();
	  cmp.add_match(new_m);
	  break;
	}
      }
    }
  }

  return new ImmutableComparison(&cmp);
}



/*******************************************************************************/
ImmutableComparison * ImmutableComparison::build_transitive(const
	    ImmutableComparison &bc,
            const std::string a_seq,
	    const std::string c_seq, float threshold) const
{
  if (_windowsize != bc._windowsize ||
      _bottom_length != bc._top_length ||
      a_seq.length() != _top_length ||
      c_seq.length() != bc._bottom_length) {
      throw paircomp_exception("AB/BC are incompatible comparisons, OR the sequences passed in have different lengths.");
  }

  std::string top_seq = a_seq;
  std::string bot_seq = c_seq;

  _preprocess(top_seq);
  _preprocess(bot_seq);
  _check_parameters(top_seq, bot_seq, _windowsize, threshold);
  std::string seq2_rev = _reverse_complement(bot_seq);

  const char * top = top_seq.c_str();
  const char * bot = bot_seq.c_str();
  const char * bot_rev = seq2_rev.c_str();

  MutableComparison new_ac(_top_length, bc._bottom_length, _windowsize);

  unsigned int min_matches = (unsigned int) (threshold * (float)_windowsize + .5);
  this->_build_forward(bc, top, bot, bot_rev, min_matches, new_ac);

  ImmutableComparison * ab_r = this->reverse_bot_matches();
  ImmutableComparison * bc_r = bc.reverse_top_matches();

  ab_r->_build_forward(*bc_r, top, bot, bot_rev, min_matches, new_ac);

  delete ab_r;
  delete bc_r;

  return new ImmutableComparison(&new_ac);
}
/*******************************************************************************/

/*******************************************************************************/
void ImmutableComparison::_build_forward(const ImmutableComparison &bc,
		    const char * top, const char * bot, const char * bot_rev,
		    unsigned int ithreshold, MutableComparison &new_ac)
  const
{
  matches_iterator ab_iter, bc_iter;
  const unsigned int bot_len = strlen(bot);

  for (ab_iter = _matches.begin(); ab_iter != _matches.end(); ab_iter++) {
    const _MatchContainer * ab_cont = ab_iter->second;

    for (unsigned int i = 0; i < ab_cont->num; i++) {
      const Match * ab_match = &ab_cont->block[i];

      // only consider forward matches from a-->b...
      if (ab_match->get_orientation() == -1) { continue; }

      // (...because this line only makes sense if you consider only
      // forward matches from a-->b)
      bc_iter = bc._matches.find(ab_match->get_bot_pos());
      if (bc_iter == bc._matches.end()) { // no b-->c matches for this b.
	continue;
      }

      const _MatchContainer * bc_cont = bc_iter->second;

      for (unsigned int j = 0; j < bc_cont->num; j++) {
	const Match * bc_match = &bc_cont->block[j];
	unsigned int m = 0;

	const char * top_w = ab_match->get_top_window(top);
	const char * bot_w = bc_match->get_bot_window(bot, bot_rev, bot_len);
	m = _n_matching(top_w, 0, bot_w, 0, _windowsize);

	if (m >= ithreshold) {
	  Match * new_match = new Match(ab_match->get_top_pos(),
					bc_match->get_bot_pos(),
					_windowsize, m,
					bc_match->get_orientation());

	  const char * new_bot_w = new_match->get_bot_window(bot, bot_rev, bot_len);
	  assert (bot_w == new_bot_w);

	  // add this new match into AC.
	  new_ac.add_match(new_match, false);
	}
      }
    }
  }

  return;
}
/*******************************************************************************/

/*******************************************************************************/
void ImmutableComparison::filter_transitively(const ImmutableComparison &bc,
					      const ImmutableComparison &ac,
					      ImmutableComparison **new_ab_i,
					      ImmutableComparison **new_bc_i,
					      ImmutableComparison **new_ac_i)
  const
{
  // check to make sure that at least some of the parameters are appropriate...
  // in the absence of sequence, we can't know if these are *actually*
  // ab/bc/ac comparisons.

  if (_windowsize != bc._windowsize ||
      _windowsize != ac._windowsize ||
      _top_length != ac._top_length ||
      _bottom_length != bc._top_length ||
      ac._bottom_length != bc._bottom_length) {
    throw paircomp_exception("must pass AB/BC/AC comparisons with the same windowsize into filter_transitively");
  }
      
  MutableComparison new_ab(_top_length, _bottom_length, _windowsize);

  MutableComparison new_bc(bc.get_top_length(),
			   bc.get_bottom_length(),
			   bc.get_windowsize());
  MutableComparison new_ac(ac.get_top_length(),
			   ac.get_bottom_length(),
			   ac.get_windowsize());

  this->_filter_forward(bc, false, ac, false, &new_ab, &new_bc, &new_ac);

  ImmutableComparison * ab_r = this->reverse_bot_matches();
  ImmutableComparison * bc_r = bc.reverse_top_matches();
  ab_r->_filter_forward(*bc_r, true, ac, false, &new_ab, &new_bc, &new_ac);

  ImmutableComparison * bc_rr = bc_r->reverse_bot_matches();
  ImmutableComparison * ac_r = ac.reverse_bot_matches();
  ab_r->_filter_forward(*bc_rr, true, *ac_r, true, &new_ab, &new_bc, &new_ac);

  ImmutableComparison * bc_r2 = bc.reverse_bot_matches();
  this->_filter_forward(*bc_r2, false, *ac_r, true, &new_ab, &new_bc, &new_ac);

  delete ab_r;
  delete bc_r;
  delete bc_rr;
  delete ac_r;
  delete bc_r2;

  *new_ab_i = new ImmutableComparison(&new_ab);
  *new_bc_i = new ImmutableComparison(&new_bc);
  *new_ac_i = new ImmutableComparison(&new_ac);
}

// #define _TEST_TRANSITIVE







/******************************************************************************/
void ImmutableComparison::_filter_forward(
		  const ImmutableComparison &bc, bool rev_b,
		  const ImmutableComparison &ac, bool rev_c,
		  MutableComparison * new_ab,
		  MutableComparison * new_bc,
		  MutableComparison * new_ac)
  const
{
  matches_iterator iter_ab, iter_bc, iter_ac;

  for (iter_ab = _matches.begin(); iter_ab != _matches.end(); iter_ab++) {
    const unsigned int a_pos = iter_ab->first;
    const _MatchContainer * ab_cont = iter_ab->second;

    iter_ac = ac._matches.find(a_pos);

    if (iter_ac == ac._matches.end()) { // no a-->c matches for this a, exeunt.
#ifdef _TEST_TRANSITIVE
      printf("FAIL at ac contains %d\n", a_pos);
#endif // TEST_TRANSITIVE
      continue;
    }

    const _MatchContainer * ac_cont = iter_ac->second;

    for (unsigned int i = 0; i < ab_cont->num; i++) { // check a-->b
      const Match * ab_match = &ab_cont->block[i];

      if (ab_match->get_orientation() == -1) {
	     continue;
      }

      iter_bc = bc._matches.find(ab_match->get_bot_pos());

      if (iter_bc == bc._matches.end()) { // no b-->c matches for this b, exit.
#ifdef _TEST_TRANSITIVE
	printf("FAIL at bc contains %d\n", ab_match->get_bot_pos());
#endif // TEST_TRANSITIVE
	continue;
      }

      const _MatchContainer * bc_cont = iter_bc->second;

      for (unsigned int j = 0; j < ac_cont->num; j++) {
	     const Match * ac_match = &ac_cont->block[j];

	for (unsigned int k = 0; k < bc_cont->num; k++) {
	   const Match * bc_match = &bc_cont->block[k];

#ifdef _TEST_TRANSITIVE
	  printf("(%d %d) --> (%d %d) --> (%d %d)??\n",
		 ab_match->get_top_pos(), ab_match->get_bot_pos(),
		 bc_match->get_top_pos(), bc_match->get_bot_pos(),
		 ac_match->get_top_pos(), ac_match->get_bot_pos());
#endif // TEST_TRANSITIVE
      
	  /*******************************************************************
	  * Check if there is  match at the same positions between bc and ac,
	  * if the position and orientation match then it is transitive. 
	  ********************************************************************/
	  if (bc_match->get_bot_pos() == ac_match->get_bot_pos() &&
	      (ab_match->get_orientation() * bc_match->get_orientation() ==
	       ac_match->get_orientation())) {

#ifdef _TEST_TRANSITIVE
	    printf("(%d %d) --> (%d %d) --> (%d %d)\n",
		   ab_match->get_top_pos(), ab_match->get_bot_pos(),
		   bc_match->get_top_pos(), bc_match->get_bot_pos(),
		   ac_match->get_top_pos(), ac_match->get_bot_pos());
	    
#endif // TEST_TRANSITIVE

	    // found a transitive match!
	    Match * new_ab_match = ab_match->copy();
	    if (rev_b) {
	      new_ab_match->reverse_bottom(this->_bottom_length);
	    }
	    new_ab->add_match(new_ab_match, false);

	    Match * new_bc_match = bc_match->copy();
	    if (rev_b) {
	      new_bc_match->reverse_top(bc._top_length,
					bc._bottom_length);
	    }
	    if (rev_c) {
	      new_bc_match->reverse_bottom(ac._bottom_length);
	    }
	    new_bc->add_match(new_bc_match, false);

	    Match * new_ac_match = ac_match->copy();
	    if (rev_c) {
	      new_ac_match->reverse_bottom(ac._bottom_length);
	    }
	    new_ac->add_match(new_ac_match, false);
	  }
	}
      }
    }
  }
}
/*******************************************************************************/