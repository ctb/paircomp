#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
/*!
*   \file MutableComparison.cc
*/
#include "MutableComparison.hh"

using namespace paircomp;

MutableComparison::MutableComparison(unsigned int top, unsigned int bot,
				     unsigned int w)
  : Comparison(top, bot, w)
{
  ;
}

MutableComparison::~MutableComparison()
{
  std::map<unsigned int, std::vector<Match*> >::const_iterator iter;

  for (iter = _matches.begin(); iter != _matches.end(); iter++) {
    std::vector<Match *> v = iter->second;
    
    for (unsigned int i = 0; i < v.size(); i++) {
      delete v[i];
    }
  }
}

//
// add_match -- adds Match to pos, consuming it in the process.
//

void MutableComparison::add_match(Match * &match, bool require_uniq)
{
  unsigned int pos = match->get_top_pos();

  std::vector<Match*> v = _matches[pos];

  bool do_add = true;

  for (unsigned int i = 0; i < v.size(); i++) {
    Match * m = v[i];
    if (m->equals(match)) {
      if (require_uniq) {
	throw paircomp_exception("match is not unique, as promised!");
      } else {			// eliminate redundancy
	delete match; match = NULL;
	do_add = false;
	break;
      }
    }
  }

  if (do_add) {
    v.push_back(match);

    _matches[pos] = v;
  }

  match = NULL;

}

// Load from bufs in various formats.

void MutableComparison::parse_paircomp_format(char * buf)
{
  char * line, * eol;

  //
  // Now run through all the lines & figure it out.
  //

  line = buf;

  while (1) {
    char * top_p, * bot_p, * n_p, * o_p;
    unsigned int top, bot, n;
    int o;

    eol = index(line, '\n');
    if (eol == NULL) {
      break;
    } else if (eol - line == 0) {
      line++;
      continue;
    }
    
    top_p = line;
    top = atoi(top_p);

    bot_p = index(top_p, '\t') + 1;
    bot = atoi(bot_p);

    n_p = index(bot_p, '\t') + 1;
    n = atoi(n_p);

    o_p = index(n_p, '\t');
    if (o_p != NULL && o_p < eol) {
      o = atoi(o_p);
    } else {
      o = 1;			// by default, orientation = +1
    }

    if (!(top < this->_top_length)) {
      throw ComparisonParserException("top length as given is incorrect");
    }

    if (bot > this->_bottom_length) {
      throw ComparisonParserException("bottom match out of bounds");
    }

    if (n > this->_windowsize) {
      throw ComparisonParserException("n matches larger than windowsize!");
    }

    // Construct the match, and add it:
    Match * m = new Match(top, bot, _windowsize, n, o);
    add_match(m);

    line = eol + 1;
  }
}

void MutableComparison::parse_seqcomp_format(const char * inbuf)
{
  char * buf = (char *) malloc(strlen(inbuf) + 5); // @CTB keep 5.
  char * freebuf = buf;	// @CTB hack hack hack.
  memcpy(buf, inbuf, strlen(inbuf) + 1);

  while (isspace(*buf)) buf++; // iterate past first whitespace

  char * line = &buf[strlen(buf) - 1];
  while (isspace(*line)) line--;

  *(line + 1) = '\n';		// terminate last line.
  *(line + 2) = 0;
  line = buf;

  char * eol = index(line, '\n'); // find first real line
  int line_length = (int) (eol - line);

  // Check to make sure the right magic number is there.
  if (!line_length || strncmp(line, "b3.5", line_length)) {
    free(freebuf);

    throw ComparisonParserException("magic number 'b3.5' not present in file");
  }

  // Get the windowsize:
  line = eol + 1;
  while (isspace(*line)) line++;
  eol = index(line, '\n');

  char windowsize_str[50];
  strncpy(windowsize_str, line, (eol - line));
  windowsize_str[(eol-line)] = 0;
  
  unsigned int windowsize = (unsigned int) atoi(windowsize_str);
  if (windowsize != this->_windowsize) {
    std::string exc = "windowsize (in buf) does not match expected windowsize";

    free(freebuf);
    throw ComparisonParserException(exc);
  }

  //
  // Find out how many lines there are; that's the top length.
  //

  int i = 0;
  while(eol) {
    line = eol;
    eol = index(eol + 1, '\n');
    i++;
  }

  i -= 1;			// for the last newline.

  // top_len is set to be the number of lines.
  unsigned int top_len = (unsigned int) i;

  if (top_len != _top_length - _windowsize + 1) {
    std::string exc = "seqcomp lines in buf do not match expected top length.";

    // exc += "seqcomp lines " + top_len;
    //exc += " (in buf) does not match expected top length " +
    //  (this->_top_length - _windowsize + 1);
    // printf("error here is: %s\n", exc.c_str());

    free(freebuf);
    throw ComparisonParserException(exc);
  }

  //
  // OK: now iterate through all the lines, collecting stuff.
  //
  eol = index(buf, '\n');	// end of first line.
  eol = index(eol + 1, '\n');	// end of second line.
  line = eol + 1;		// beginning of third line.
  eol = index(line, '\n');	// end of third line.


  char * match_pos, * n_matches_pos, * orientation_pos;
  unsigned int match, n_matches;
  int orientation;
  bool done;

  unsigned int line_n;

  line_n = 0;
  while(line_n < top_len) {
    match_pos = line;
    done = false;

    //
    // Run through the file, and for each line, read 0 or more pairs of
    // matches.
    //

    while(!done) {
      n_matches_pos = index(match_pos, ' ');

      //
      // If there are no more spaces on the line, then we've no more pairs
      // of matches to read.  Finish the line.
      //
      
      if (n_matches_pos > eol || n_matches_pos == NULL) {
	done = true;

	// The average number of matches is the last number on each line.
	// averages[line_n] = atof(match_pos);
      }
      else {			// Read another pair

	orientation_pos = n_matches_pos - 1;
	n_matches_pos++;

	if (orientation_pos[0] == '+') {
	  orientation = 1;
	} else if (orientation_pos[0] == '-') {
	  orientation = -1;
	} else {
	  free(freebuf);
	  throw ComparisonParserException("unable to understand orientation");
	}

	match = atoi(match_pos);
	if (match > this->_bottom_length) {
	  free(freebuf);
	  throw ComparisonParserException("bottom match out of bounds");
	}

	n_matches = atoi(n_matches_pos);
	if (n_matches > this->_windowsize) {
	  free(freebuf);
	  throw ComparisonParserException("n matches larger than windowsize!");
	}


	//
	// Create a match object & save it.
	//

	Match * matchObj = new Match(line_n, match, windowsize,
				     n_matches, orientation);

	add_match(matchObj);

	match_pos = index(n_matches_pos, ' ');
	match_pos++;
      }
    }

    // Increment stuff.

    line = eol + 1;
    eol = index(line, '\n');
    line_n ++;
  }

  free(freebuf);
}

///// File loading stuff
/*!
*   \fn char * _load_file(char * filename)
*   \brief A helper function for loading a file and returning the character text.
*   \param filename A valid filename.
*/
char * _load_file(char * filename)
{
  FILE * fp;
  long filesize;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    throw ComparisonParserException("file does not exist");
  }

  // Determine the file size.
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  rewind(fp);

  // Allocate a buffer & slurp it all into the buffer.
  char * buf = (char *) malloc(filesize + 1);
  long items_read = fread(buf, filesize, 1, fp);
  if (items_read != 1) {
    throw ComparisonParserException("error reading from file");
  }

  buf[filesize] = 0;		// null terminate string.

  return buf;
}

MutableComparison * paircomp::load_paircomp_file(unsigned int top_len,
						 unsigned int bot_len,
						 unsigned int windowsize,
						 char * filename)
{
  MutableComparison * cmp = new MutableComparison(top_len, bot_len,
						  windowsize);

  char * buf = _load_file(filename);
  cmp->parse_paircomp_format(buf);

  free(buf);

  return cmp;
}

MutableComparison * paircomp::load_seqcomp_file(unsigned int top_len,
						unsigned int bot_len,
						unsigned int windowsize,
						char * filename)
{
  MutableComparison * cmp = new MutableComparison(top_len, bot_len,
						  windowsize);

  char * buf = _load_file(filename);
  cmp->parse_seqcomp_format(buf);

  free(buf);

  return cmp;
}

