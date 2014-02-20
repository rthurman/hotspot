/**
 * File: ByLine.hpp
 * Author: Shane Neph
 * Date: Autumn 2006
 * Version: $Id$
 *
 * Comments:
 *   see below
 */

// Macro guard
#ifndef BYLINE_H
#define BYLINE_H

#include <string>
#include <iostream>

//==============================================================================
// ByLine structure:
//  Simple extension of the std::string class to allow reading input by lines
//   rather than by whitespace by default.
//==============================================================================
struct ByLine : public std::string {
  friend std::istream& operator>>(std::istream& is, ByLine& b) {
    std::getline(is, b);
    return(is);
  }
};

#endif // BYLINE_H

