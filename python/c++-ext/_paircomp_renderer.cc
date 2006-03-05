// Python interface functions for a seqcomp renderer.
//
// See README.txt for license and copyright information.

#include <stdlib.h>
#include <math.h>

#include "Python.h"

#include "MatchList.hh"
#include "Comparison.hh"

using namespace paircomp;

//
// Function necessary for Python loading:
//

extern "C" {
  void init_paircomp_renderer();
}

// render a comparison object.

static PyObject * render_comparison(PyObject * self, PyObject * args)
{
  int top_start, top_end, bot_start, bot_end;
  PyObject * p, * fn, * passthru;

  if (!PyArg_ParseTuple(args, "OiiiiOO",
			&p,
			&top_start, &top_end,
			&bot_start, &bot_end,
			&fn, &passthru)) {
    return NULL;
  }

  if (top_start < 0 || top_start > top_end ||
      bot_start < 0 || bot_start > bot_end) {
    PyErr_SetString(PyExc_Exception, "invalid arguments -- check your numbers");
    return NULL;
  }

  // Extract the comparison object.
  Comparison * c = (Comparison *) PyCObject_AsVoidPtr(p);

  // Make sure that the function passed in is indeed callable.
  if (!PyCallable_Check(fn)) {
    return NULL;
  }

  if (top_end > c->getTopLength()) {
    PyErr_SetString(PyExc_Exception, "top_len as given is too long");
    return NULL;
  }

  int start = top_start - c->getWindowsize();
  if (start < 0) { start = 0; }

  int end = top_end + c->getWindowsize();
  if (end >= c->getTopLength()) { end = c->getTopLength(); }

  // Run through the comparisons and matches...
  for (int i = start; i < end; i++) {
    MatchList * list = c->getMatchList(i);
    if (list != NULL) {
      MatchListIterator * iter = list->getIterator();

      const Match * m = iter->first();
      while(m != NULL) {

	//
	// For each match, check to see if it overlaps with the visible
	// portion of things.  If it does, call the passed-in Python function
	// to plot it.
	// 

	int o = m->get_orientation(), top_pos = m->get_top_pos(),
	  bot_pos = m->get_bot_pos(), length = m->get_length();

	int bot_a, bot_b;
	if (o == -1) {
	  bot_a = bot_pos - length + 1;
	  bot_b = bot_pos;
	} else {
	  bot_a = bot_pos;
	  bot_b = bot_pos + length;
	}

	if (bot_a >= bot_start && bot_b < bot_end) {
	  // plot it.
	  PyObject * r = PyObject_CallFunction(fn, "Oiiiii",
					       passthru,
					       top_pos, bot_pos, o,
					       length, m->get_n_matching());
	  if (r == NULL) {
	    delete iter;
	    return NULL;
	  } 
	  Py_DECREF(r);
	}

	m = iter->next();
      }
      delete iter;
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}

// pixelize and then render a comparison object.

// a class to keep track of pixel positions & matches
class _pixel_pos_link {
public:
  _pixel_pos_link(int p, int m) { pixel = p; max_match = m; next = NULL; }

  int pixel;
  int max_match;
  _pixel_pos_link * next;

  void set_match(int m) {
    if (m > max_match) { max_match = m; }
  }
};

// a class to keep track of the lists of pixel positions & matches.
class _pixel_pos_head {
public:
  _pixel_pos_link * head;

  _pixel_pos_head() { head = NULL; }

  // Destructor -- de-allocate linked list.
  ~_pixel_pos_head() {
    _pixel_pos_link * tmp = head;
    while(head) {
      tmp = head; head = head->next;
      delete tmp;
    }
  }

  // Add in a new link, keeping it sorted by pixel pos.
  void add(int p, int m) {
    _pixel_pos_link * l = NULL;

    if (head == NULL) {
      head = new _pixel_pos_link(p, m);
      return;
    }

    // Insert at head?
    if (p < head->pixel) {
      l = new _pixel_pos_link(p, m);
      l->next = head; head = l;
    } else if (p == head->pixel) { // update head?
      head->set_match(m);
    } else {			// add in somewhere later: traverse.
      _pixel_pos_link * c = head;

      while(c->next && c->next->pixel >= l->pixel) { c = c->next; }

      if (c->next && c->next->pixel == p) { // update next
	c->next->set_match(m);
      } else {			// insert new link
	l = new _pixel_pos_link(p, m);
	l->next = c->next; c->next = l;
      }
    }
  }
};

// actually render things.

static PyObject * render_pixelized_comparison(PyObject * self, PyObject * args)
{
  int top_start, top_end, bot_start, bot_end, top_pixels, bot_pixels;
  PyObject * p, * fn, * passthru;

  if (!PyArg_ParseTuple(args, "OiiiiiiOO",
			&p,
			&top_start, &top_end, &top_pixels,
			&bot_start, &bot_end, &bot_pixels,
			&fn, &passthru)) {
    return NULL;
  }

  if (top_start < 0 || top_start > top_end || top_pixels <= 0 ||
      bot_start < 0 || bot_start > bot_end || bot_pixels <= 0) {
    PyErr_SetString(PyExc_Exception, "invalid arguments -- check your numbers");
    return NULL;
  }

  // Make sure that the function passed in is indeed callable.
  if (!PyCallable_Check(fn)) {
    PyErr_SetString(PyExc_Exception, "fn must be callable!");
    return NULL;
  }

  // Extract the comparison object.
  Comparison * c = (Comparison *) PyCObject_AsVoidPtr(p);

  if (top_end > c->getTopLength()) {
    PyErr_SetString(PyExc_Exception, "top_len as given is too long");
    return NULL;
  }

  // Figure out by how much to scale things:
  int top_width = top_end - top_start;
  int bot_width = bot_end - bot_start;
  float top_scale = (float)top_pixels / (float)top_width;
  float bot_scale = (float)bot_pixels / (float)bot_width;

  // Allocate the array of pixelization holders.
  _pixel_pos_head ** top_pixel_list;
  top_pixel_list = (_pixel_pos_head **) calloc(sizeof(_pixel_pos_head *),
					       top_pixels);

  Py_BEGIN_ALLOW_THREADS

  // Run through the comparisons & matches:

  for (int i = top_start; i < top_end; i++) {
    int top_pixel = (int)((i - top_start) * top_scale + .5);
    MatchList * list = c->getMatchList(i);

    if (list != NULL) {
      MatchListIterator * iter = list->getIterator();

      const Match * m = iter->first();
      while(m != NULL) {

	//
	// For each match, insert it into the pixelization lists.
	// 

	int bot_pos = m->get_bot_pos(), n_matching = m->get_n_matching();

	if (bot_pos >= bot_start && bot_pos < bot_end) {
	  int bot_pixel = (int)((bot_pos - bot_start) * bot_scale + .5);

	  if (top_pixel_list[top_pixel] == NULL) {
	    top_pixel_list[top_pixel] = new _pixel_pos_head();
	  }
	  top_pixel_list[top_pixel]->add(bot_pixel, n_matching);
	}

	m = iter->next();
      }
      delete iter;
    }
  }

  Py_END_ALLOW_THREADS

  // OK: now, run through all of the various pixel pairs, and plot 'em.
  int break_through = 0;
  for (int i = 0; i < top_pixels; i++) {
    if (top_pixel_list[i] != NULL) {
      _pixel_pos_link * c;

      c = top_pixel_list[i]->head;
      while(c != NULL) {
	PyObject * r = PyObject_CallFunction(fn, "Oiii", passthru,
					     i, c->pixel, c->max_match);

	if (r == NULL) {
	  break_through = 1;	// still gotta clean up...
	  break;
	}
	Py_DECREF(r);

	c = c->next;
      }
    }
    if (break_through) {
      break;
    }
  }

  // Clean things up.

  for (int i = 0; i < top_pixels; i++) {
    if (top_pixel_list[i] != NULL) {
      delete top_pixel_list[i];
    }
  }

  if (break_through) { return NULL; }

  Py_INCREF(Py_None);
  return Py_None;
}

//
// Module machinery.
//

static PyMethodDef SeqcompRendererMethods[] = {
  { "render_comparison", render_comparison, METH_VARARGS },
  { "render_pixelized_comparison", render_pixelized_comparison, METH_VARARGS },
  { NULL, NULL }
};

void init_paircomp_renderer()
{
  (void) Py_InitModule("_paircomp_renderer", SeqcompRendererMethods);
  ;
}
