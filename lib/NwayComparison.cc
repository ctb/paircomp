#include "NwayComparison.hh"
#include "algorithms.hh"

using namespace paircomp;

void NwayComparison::do_comparisons()
{
  std::map<MapIndexPair, ImmutableComparison *>::iterator iter;

  for (iter = _comparisons.begin(); iter != _comparisons.end(); iter++) {
    MapIndexPair coord = iter->first;
    ImmutableComparison * cmp = iter->second;

    if (cmp == NULL) {
      unsigned int a = coord.x;
      unsigned int b = coord.y;

      std::string seq1 = _sequences[a];
      std::string seq2 = _sequences[b];

      cmp = rolling_nxn_comparison(seq1, seq2, windowsize, threshold);
      _comparisons[coord] = cmp;
    }
  }
}

std::vector<NwayPath> NwayComparison::make_paths(unsigned int start_seq, unsigned int pos, int max_paths)
{
  std::vector<NwayPath> paths;

  unsigned int next_seq = start_seq + 1;
  if (next_seq == _sequences.size()) {
    NwayPath p;
    p.push_back(PosAndO(pos, 1));
    paths.push_back(p);
    return paths;
  }
  
  ImmutableComparison * cmp = this->get_comparison(start_seq, next_seq);
  const _MatchContainer * cont = cmp->get_matches(pos);

  if (cont == NULL) { return paths; }

  max_paths -= cont->num;

  if (max_paths < 0) { return paths; }
  
  for (unsigned int j = 0; j < cont->num; j++) {
    const Match * m = &cont->block[j];
    unsigned int next_pos;
    if (m->get_orientation() == 1) {
      next_pos = m->get_bot_pos();
    } else {
      next_pos = m->get_bot_pos() - windowsize + 1;
    }

    // recurse
    std::vector<NwayPath> sub_paths = make_paths(next_seq, next_pos, max_paths);
    
    // insert this position at the beginning.
    for (unsigned int k = 0; k < sub_paths.size(); k++) {
      NwayPath sub_path = sub_paths[k];
      sub_path.insert(sub_path.begin(), PosAndO(pos, m->get_orientation()));
      paths.push_back(sub_path);
    }
  }

  return paths;
}

bool NwayComparison::check_match(unsigned int from, unsigned int to,
				 const NwayPath& path)
{
  unsigned int top_pos = path[from].pos;
  unsigned int bot_pos = path[to].pos;

  int o = 1;
  for (unsigned int i = from; i < to; i++) {
    o = o * path[i].orient;
  }

  if (o == -1) {		// correct for orientation
    bot_pos = bot_pos + windowsize - 1;
  }

  ImmutableComparison * cmp = get_comparison(from, to);
  const _MatchContainer * cont = cmp->get_matches(top_pos);

  if (cont == NULL) { return false; }

  for (unsigned int i = 0; i < cont->num; i++) {
    const Match * m = &cont->block[i];
    if (m->get_bot_pos() == bot_pos) {
      return true;
    }
  }
  
  return false;
}

std::vector<NwayPath> NwayComparison::filter()
{
  std::vector<NwayPath> ret_paths;
  if (_sequences.size() < 2) {
    return ret_paths;
  }

  do_comparisons();

  std::string seq0 = _sequences[0];

  for (unsigned int i = 0; i < seq0.length(); i++) {
    std::vector<NwayPath> paths = make_paths(0, i, _max_paths);

    for (unsigned int z = 0; z < paths.size(); z++) {
      bool keep_path = true;
      NwayPath path = paths[z];
      for (unsigned int j = 2; j < _sequences.size() && keep_path; j++) {
	for (unsigned int k = 0; k < _sequences.size() - j; k++) {
	  if (!check_match(k, k + j, path)) {
	    keep_path = false;
	    break;
	  }
	}
      }
      if (keep_path) {
	ret_paths.push_back(path);
      }
    }
  }
  return ret_paths; 
}
