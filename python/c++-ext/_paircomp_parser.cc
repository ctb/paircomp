// Python interface functions for Comparison objects.
//
// See README.txt for license and copyright information.

#include "Python.h"

#include "paircomp.hh"

using namespace paircomp;

// a class to automatically handle saving of thread state.
class py_thread_saver {
protected:
  PyThreadState * _tstate;
public:
  py_thread_saver() { _tstate = PyEval_SaveThread(); }
  ~py_thread_saver() { PyEval_RestoreThread(_tstate); }
};

//
// Function necessary for Python loading:
//

extern "C" {
  void init_paircomp_parser();
}

static PyObject * PaircompParserError;

//
// cleanup function for ImmutableComparison objects.
//

static void cleanup_ImmComparison(void * p)
{
  ImmutableComparison * c = (ImmutableComparison *) p;
  delete c;
}

//
// cleanup function for NwayComparison objects.
//

static void cleanup_NwayComparison(void * p)
{
  NwayComparison * c = (NwayComparison *) p;
  delete c;
}

//
// Various exported functions.
//

// read in a seqcomp comparison.

static PyObject * parse_seqcomp_comparison(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  char * buf;
  int windowsize, top_len, bot_len;
  ImmutableComparison * c = NULL;

  if (!PyArg_ParseTuple(args, "siii", &buf,
			&top_len, &bot_len,
			&windowsize)) {
    return NULL;
  }

  try {
    { py_thread_saver save;

    MutableComparison cmp(top_len, bot_len, windowsize);
    cmp.parse_seqcomp_format(buf);

    c = new ImmutableComparison(&cmp);
    }

    ret = PyCObject_FromVoidPtr(c, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

// save a seqcomp comparison in seqcomp format.

static PyObject * save_seqcomp_comparison(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  char * filename;
  PyObject * p;
  ImmutableComparison * c;

  if (!PyArg_ParseTuple(args, "sO", &filename, &p)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  try {
    c->save_as_seqcomp(filename);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

// read in a paircomp comparison.

static PyObject * parse_paircomp_comparison(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  char * buf;
  int windowsize, top_len, bot_len;
  ImmutableComparison * c;

  if (!PyArg_ParseTuple(args, "siii", &buf,
			&top_len, &bot_len, &windowsize)) {
    return NULL;
  }

  try {
    { py_thread_saver save;
    MutableComparison cmp(top_len, bot_len, windowsize);
    cmp.parse_paircomp_format(buf);

    c = new ImmutableComparison(&cmp);
    }
    ret = PyCObject_FromVoidPtr(c, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

// save a seqcomp comparison in shoudan's format.

static PyObject * save_paircomp_comparison(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  char * filename;
  PyObject * p;
  ImmutableComparison * c;

  if (!PyArg_ParseTuple(args, "sO", &filename, &p)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  try {
    c->save_as_paircomp(filename);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

// get the top-sequence length from a comparison.

static PyObject * get_top_length(PyObject * self, PyObject * args)
{
  Comparison * c;
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  return Py_BuildValue("i", c->get_top_length());
}

// get the bottom-sequence length from a comparison.

static PyObject * get_bottom_length(PyObject * self, PyObject * args)
{
  Comparison * c;
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  return Py_BuildValue("i", c->get_bottom_length());
}

// get the windowsize from a comparison.

static PyObject * get_windowsize(PyObject * self, PyObject * args)
{
  Comparison * c;
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  return Py_BuildValue("i", c->get_windowsize());
}

// get the matches at a particular position

static PyObject * get_matches(PyObject * self, PyObject * args)
{
  int pos;
  ImmutableComparison * c;
  PyObject * p;

  if (!PyArg_ParseTuple(args, "Oi", &p, &pos)) {
    return NULL;
  }

  c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  const _MatchContainer * cont = c->get_matches((unsigned int) pos);

  if (cont == NULL) {
    return Py_BuildValue("()");
  }

  PyObject * t = PyTuple_New(cont->num);
  for (unsigned int i = 0; i < cont->num; i++) {
    const Match * m = &cont->block[i];

    PyObject * item = Py_BuildValue("(iiiii)",
				    m->get_top_pos(),
				    m->get_bot_pos(),
				    m->get_length(),
				    m->get_n_matching(),
				    m->get_orientation());
    PyTuple_SET_ITEM(t, i, item);
  }

  return t;
}

// reverse the top sequence coordinates.

static PyObject * reverse_top(PyObject * self, PyObject * args)
{
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * r = c->reverse_top_matches();

  return PyCObject_FromVoidPtr(r, cleanup_ImmComparison);
}

// reverse the bottom sequence coordinates.

static PyObject * reverse_bottom(PyObject * self, PyObject * args)
{
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }
  
  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * r = c->reverse_bot_matches();

  return PyCObject_FromVoidPtr(r, cleanup_ImmComparison);
}

// filter the matches on a floating-point threshold.

static PyObject * filter_matches(PyObject * self, PyObject * args)
{
  PyObject * p;
  float threshold;

  if (!PyArg_ParseTuple(args, "Of", &p, &threshold)) {
    return NULL;
  }

  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * f = c->filter_by_threshold(threshold);

  return PyCObject_FromVoidPtr(f, cleanup_ImmComparison);
}

// filter the matches by orientation

static PyObject * filter_orientation(PyObject * self, PyObject * args)
{
  PyObject * p;
  int forward_filter, reverse_filter;

  if (!PyArg_ParseTuple(args, "Oii", &p, &forward_filter, &reverse_filter)) {
    return NULL;
  }

  bool fwd = (forward_filter != 0);
  bool rev = (reverse_filter != 0);

  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * f = c->filter_by_orientation(fwd, rev);

  return PyCObject_FromVoidPtr(f, cleanup_ImmComparison);
}

// invert the matches.

static PyObject * invert(PyObject * self, PyObject * args)
{
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * r = c->invert();

  return PyCObject_FromVoidPtr(r, cleanup_ImmComparison);
}

// create a set of 1-bp matches

static PyObject * isolate_matching_bases(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * p;
  char * top_sp, * bot_sp;

  if (!PyArg_ParseTuple(args, "Oss", &p, &top_sp, &bot_sp)) {
    return NULL;
  }

  std::string top_seq(top_sp), bot_seq(bot_sp);
  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  try {
    ImmutableComparison * f;
    { py_thread_saver save;

    f = c->isolate_matching_bases(top_seq, bot_seq);

    }
    ret = PyCObject_FromVoidPtr(f, cleanup_ImmComparison);
  } catch (paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

// build transitive bridge ab+bc ==> ac.

static PyObject * build_transitive(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * ab, * bc;
  char * seq1;
  char * seq2;
  float threshold;

  if (!PyArg_ParseTuple(args, "OOssf", &ab, &bc, &seq1, &seq2, &threshold)) {
    return NULL;
  }

  const ImmutableComparison * map_ab, * map_bc;
  map_ab = (const ImmutableComparison *) PyCObject_AsVoidPtr(ab);
  map_bc = (const ImmutableComparison *) PyCObject_AsVoidPtr(bc);
  
  try {
    ImmutableComparison * map_ac;
    { py_thread_saver save;

    map_ac = map_ab->build_transitive(*map_bc, seq1, seq2, threshold);

    }
    ret = PyCObject_FromVoidPtr(map_ac, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

// filter three (ab, bc, ac) to two (ab, bc) requiring transitivity.

static PyObject * filter_transitively(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * ab, * bc, *ac;
  
  if (!PyArg_ParseTuple(args, "OOO", &ab, &bc, &ac)) {
    return NULL;
  }

  const ImmutableComparison * map_ab, * map_bc, * map_ac;
  ImmutableComparison * new_ab, * new_bc, * new_ac;

  map_ab = (ImmutableComparison *) PyCObject_AsVoidPtr(ab);
  map_bc = (ImmutableComparison *) PyCObject_AsVoidPtr(bc);
  map_ac = (ImmutableComparison *) PyCObject_AsVoidPtr(ac);

  try {
    { py_thread_saver save;
    map_ab->filter_transitively(*map_bc, *map_ac, &new_ab, &new_bc, &new_ac);
    }

    PyObject * p1, *p2, *p3;

    p1 = PyCObject_FromVoidPtr(new_ab, cleanup_ImmComparison);
    p2 = PyCObject_FromVoidPtr(new_bc, cleanup_ImmComparison);
    p3 = PyCObject_FromVoidPtr(new_ac, cleanup_ImmComparison);

    ret = Py_BuildValue("OOO", p1, p2, p3);
    Py_DECREF(p1);
    Py_DECREF(p2);
    Py_DECREF(p3);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * contains(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * big_p, * small_p;
  
  if (!PyArg_ParseTuple(args, "OO", &big_p, &small_p)) {
    return NULL;
  }

  const ImmutableComparison * big, * small;
  
  big = (ImmutableComparison *) PyCObject_AsVoidPtr(big_p);
  small = (ImmutableComparison *) PyCObject_AsVoidPtr(small_p);

  int answer = 0;
  try {
    if (big->contains(*small)) {
      answer = 1;
    }
    ret = Py_BuildValue("i", answer);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * is_empty(PyObject * self, PyObject * args)
{
  PyObject * p;
  
  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  const ImmutableComparison * cmp;
  
  cmp = (ImmutableComparison *) PyCObject_AsVoidPtr(p);

  int answer = 0;
  if (cmp->is_empty()) {
    answer = 1;
  }

  return Py_BuildValue("i", answer);
}

// intersect the matches.

static PyObject * intersect(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * p, * o;

  if (!PyArg_ParseTuple(args, "OO", &p, &o)) {
    return NULL;
  }

  ImmutableComparison * c = (ImmutableComparison *) PyCObject_AsVoidPtr(p);
  ImmutableComparison * other = (ImmutableComparison *) PyCObject_AsVoidPtr(o);

  try {
    ImmutableComparison * r;

    { py_thread_saver save;
    r = c->intersect(*other);
    }
    ret = PyCObject_FromVoidPtr(r, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * subtract(PyObject * self, PyObject * args)
{
  PyObject * ret = NULL;
  PyObject * big_p, * small_p;
  
  if (!PyArg_ParseTuple(args, "OO", &big_p, &small_p)) {
    return NULL;
  }

  const ImmutableComparison * big, * small;
  
  big = (ImmutableComparison *) PyCObject_AsVoidPtr(big_p);
  small = (ImmutableComparison *) PyCObject_AsVoidPtr(small_p);

  ImmutableComparison * diff = NULL;

  try {
    { 
      py_thread_saver save;

      diff = big->subtract(*small);
    }

    ret = PyCObject_FromVoidPtr(diff, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {

    PyErr_SetString(PaircompParserError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * create_nway(PyObject * self, PyObject * args)
{
  unsigned int windowsize;
  float threshold;

  if (!PyArg_ParseTuple(args, "If", &windowsize, &threshold)) {
    return NULL;
  }

  NwayComparison * nway = new NwayComparison(windowsize, threshold);

  return PyCObject_FromVoidPtr(nway, cleanup_NwayComparison); // @CTB
}

static PyObject * add_sequence_to_nway(PyObject * self, PyObject * args)
{
  char * seq;
  PyObject * p;

  if (!PyArg_ParseTuple(args, "Os", &p, &seq)) {
    return NULL;
  }

  NwayComparison * nway = (NwayComparison *) PyCObject_AsVoidPtr(p);
  nway->add_sequence(seq);

  Py_INCREF(Py_None);
  return Py_None;
}

std::string print_poso_v(NwayPath v)
{
  char buf[50];
  std::string ret;
  for (unsigned int i = 0; i < v.size(); i++) {
    PosAndO p = v[i];
    sprintf(buf, "%d(%c) ", p.pos, p.orient > 0 ? '+' : '-');
    ret += buf;
  }
  ret += "\n";

  return ret;
}

static PyObject * get_nway_filtered_paths(PyObject * self, PyObject * args)
{
  PyObject * p;

  if (!PyArg_ParseTuple(args, "O", &p)) {
    return NULL;
  }

  NwayComparison * nway = (NwayComparison *) PyCObject_AsVoidPtr(p);

  std::vector<NwayPath> paths = nway->filter();

  PyObject * t = PyTuple_New(paths.size());
  for (unsigned int i = 0; i < paths.size(); i++) {
    NwayPath p = paths[i];

    PyObject * t2 = PyTuple_New(p.size());

    for (unsigned int j = 0; j < p.size(); j++) {
      PosAndO po = p[j];
      PyObject * item = Py_BuildValue("(ii)", po.pos, po.orient);
      PyTuple_SET_ITEM(t2, j, item);
    }
    PyTuple_SET_ITEM(t, i, t2);
  }
  
  return t;
}

//
// Module machinery.
//

static PyMethodDef SeqcompParserMethods[] = {
  { "parse_seqcomp_comparison", parse_seqcomp_comparison, METH_VARARGS },
  { "save_seqcomp_comparison", save_seqcomp_comparison, METH_VARARGS },
  { "parse_paircomp_comparison", parse_paircomp_comparison, METH_VARARGS },
  { "save_paircomp_comparison", save_paircomp_comparison, METH_VARARGS },
  { "get_top_length", get_top_length, METH_VARARGS },
  { "get_bottom_length", get_bottom_length, METH_VARARGS },
  { "get_windowsize", get_windowsize, METH_VARARGS },
  { "get_matches", get_matches, METH_VARARGS },
  { "reverse_top", reverse_top, METH_VARARGS },
  { "reverse_bot", reverse_bottom, METH_VARARGS },
  { "invert", invert, METH_VARARGS },
  { "filter_matches", filter_matches, METH_VARARGS },
  { "filter_orientation", filter_orientation, METH_VARARGS },
  { "build_transitive", build_transitive, METH_VARARGS },
  { "filter_transitively", filter_transitively, METH_VARARGS },
  { "isolate_matching_bases", isolate_matching_bases, METH_VARARGS },
  { "intersect", intersect, METH_VARARGS },
  { "contains", contains, METH_VARARGS },
  { "is_empty", is_empty, METH_VARARGS },
  { "subtract", subtract, METH_VARARGS },
  { "create_nway", create_nway, METH_VARARGS },
  { "add_sequence_to_nway", add_sequence_to_nway, METH_VARARGS },
  { "get_nway_filtered_paths", get_nway_filtered_paths, METH_VARARGS },
  { NULL, NULL }
};

void init_paircomp_parser()
{
  PyObject * m;

  m = Py_InitModule("_paircomp_parser", SeqcompParserMethods);

  PaircompParserError = PyErr_NewException("_paircomp_parser.error", NULL, NULL);
  Py_INCREF(PaircompParserError);

  PyModule_AddObject(m, "error", PaircompParserError);
}
