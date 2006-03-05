// Python interface functions for the various paircomp algorithms.
//
// See README.txt for license and copyright information.

#include "Python.h"

#include "Comparison.hh"
#include "algorithms.hh"
#include "algorithms2.hh"

using namespace paircomp;

//
// Function necessary for Python loading:
//

extern "C" {
  void init_paircomp_algorithms();
}

// a class to automatically handle saving of thread state.
class py_thread_saver {
protected:
  PyThreadState * _tstate;
public:
  py_thread_saver() { _tstate = PyEval_SaveThread(); }
  ~py_thread_saver() { PyEval_RestoreThread(_tstate); }
};

static PyObject * PaircompAlgorithmsError;

//
// cleanup functions for ImmutableComparison objects.
//

static void cleanup_ImmComparison(void * p)
{
  ImmutableComparison * c = (ImmutableComparison *) p;
  delete c;
}

static PyObject * do_simple_nxn_compare(PyObject * self, PyObject * args)
{
  PyObject * ret;
  char * top_seq, * bot_seq;
  unsigned int windowsize;
  float threshold;

  if (!PyArg_ParseTuple(args, "ssIf", &top_seq, &bot_seq, &windowsize,
			&threshold)) {
    return NULL;
  }

  ImmutableComparison * cmp;
  try {
    { py_thread_saver save;
    cmp = simple_nxn_comparison(top_seq, bot_seq, windowsize, threshold);
    }

    ret = PyCObject_FromVoidPtr(cmp, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompAlgorithmsError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * do_rolling_nxn_compare(PyObject * self, PyObject * args)
{
  PyObject * ret;
  char * top_seq, * bot_seq;
  unsigned int windowsize;
  float threshold;

  if (!PyArg_ParseTuple(args, "ssIf", &top_seq, &bot_seq, &windowsize,
			&threshold)) {
    return NULL;
  }

  ImmutableComparison * cmp;

  try {
    { py_thread_saver save;
    cmp = rolling_nxn_comparison(top_seq, bot_seq, windowsize, threshold);
    }
    ret = PyCObject_FromVoidPtr(cmp, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompAlgorithmsError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

static PyObject * do_hashed_n_compare(PyObject * self, PyObject * args)
{
  PyObject * ret;
  char * top_seq, * bot_seq;
  unsigned int windowsize;
  float threshold;

  if (!PyArg_ParseTuple(args, "ssIf", &top_seq, &bot_seq, &windowsize,
			&threshold)) {
    return NULL;
  }

  ImmutableComparison * cmp;
  try {
    { py_thread_saver save;
    cmp = hashed_n_comparison(top_seq, bot_seq, windowsize, threshold);
    }

    ret = PyCObject_FromVoidPtr(cmp, cleanup_ImmComparison);
  } catch(paircomp_exception & e) {
    PyErr_SetString(PaircompAlgorithmsError, e.error_msg.c_str());

    ret = NULL;
  }

  return ret;
}

//
// Module machinery.
//

static PyMethodDef PaircompAlgorithmMethods[] = {
  { "do_simple_nxn_compare", do_simple_nxn_compare, METH_VARARGS },
  { "do_rolling_nxn_compare", do_rolling_nxn_compare, METH_VARARGS },
  { "do_hashed_n_compare", do_hashed_n_compare, METH_VARARGS },
  { NULL, NULL }
};

void init_paircomp_algorithms()
{
  PyObject * m;

  m = Py_InitModule("_paircomp_algorithms", PaircompAlgorithmMethods);
  
  PaircompAlgorithmsError = PyErr_NewException("_paircomp_parser.error", NULL, NULL);
  Py_INCREF(PaircompAlgorithmsError);

  PyModule_AddObject(m, "error", PaircompAlgorithmsError);
}
