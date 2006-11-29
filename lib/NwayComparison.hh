#include "ImmutableComparison.hh"

#define MAX_SEQUENCES 1000

#include <vector>
#include <map>

namespace paircomp {
  class PosAndO {
  public:
    const unsigned int pos;
    const int orient;
    PosAndO(unsigned int p, int o) : pos(p), orient(o) {};
    bool operator==(PosAndO other) const {
      return (other.pos == pos) && (other.orient == orient);
    }
    PosAndO operator=(const PosAndO other) {
      return PosAndO(other.pos, other.orient);
    }
  };

  typedef std::vector<PosAndO> NwayPath;

  class MapIndexPair {
  public:
    unsigned int x, y;
    MapIndexPair(unsigned int a, unsigned int b) : x(a), y(b) {};
    bool operator==(MapIndexPair other) const { 
      return (other.x == x && other.y == y);
    }
    bool operator<(MapIndexPair other) const { 
      return (((other.x*MAX_SEQUENCES)+ other.y) < (x*MAX_SEQUENCES)+y);
    }
    void operator=(MapIndexPair other) {
      x = other.x; y = other.y;
    }
  };

  class NwayComparison {
  protected:
    int _max_paths;
    std::vector<std::string> _sequences;
    std::map<MapIndexPair, ImmutableComparison *> _comparisons;
  public:
    const float threshold;
    const unsigned int windowsize;

    NwayComparison(unsigned int w, float t) : threshold(t), windowsize(w) {
      _max_paths = 50;
    };

    void add_sequence(const std::string seq) {
      for (unsigned int i = 0; i < _sequences.size(); i++) {
	_comparisons[MapIndexPair(i, _sequences.size())] = NULL;
	;
      }

      _sequences.push_back(seq);
    }
    std::string get_sequence(unsigned int i) { return _sequences[i]; }
    unsigned int n_sequences() { return _sequences.size(); }

    ImmutableComparison * get_comparison(unsigned int a, unsigned int b) {
      return _comparisons[MapIndexPair(a, b)];
    }

    void set_comparison(unsigned int a, unsigned int b, \
			ImmutableComparison *cmp) {
      _comparisons[MapIndexPair(a, b)] = cmp;
    }

    void do_comparisons();

    std::vector<NwayPath> make_paths(unsigned int start_seq, unsigned int pos, int max_paths);
    bool check_match(unsigned int from, unsigned int to, const NwayPath& path);
    std::vector<NwayPath> filter();

    int max_paths_at_position() { return _max_paths; }
    void max_paths_at_position(int m) { _max_paths = m; }
  };
};
