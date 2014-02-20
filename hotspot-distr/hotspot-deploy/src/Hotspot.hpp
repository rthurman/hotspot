/**
 * File: Hotspot.hpp
 * Author: Scott Kuehn
 * Date: July 25, 2009
 * Version: $Id$
 *
 * Comments:
 *   This class represents information pertaining to an
 *   individual Hotspot, such as windowing variables and statistics.
 */

#ifndef HOTSPOT_HPP_
#define HOTSPOT_HPP_

#include <cstdio>

struct Hotspot
{
  // Windowing
  int filterDist;
  int filterIndexLeft;
  int filterIndexRight;
  int filterDensIndexLeft;
  int filterDensIndexRight;
  int filterSize;
  double filterWidth;
  int minSite;
  int maxSite;
  int averagePos;
  int maxWindow;

  // Statistics
  double filteredZScore;
  double filteredZScoreAdjusted;
  double weightedAvgSD;

  // Background data
  double densCount;

  /**
   * Initialize a Hotspot
   */
  Hotspot( )
    : filterDist( -1 ),
      filterIndexLeft( -1 ),
      filterIndexRight( -1 ),
      filterDensIndexLeft( -1 ),
      filterDensIndexRight( -1 ),
      filterSize( 0 ),
      filterWidth( 0.0 ),
      minSite( -1 ),
      maxSite( -1 ),
      averagePos( -1 ),
      maxWindow( -1 ),
      filteredZScore( 0.0 ),
      filteredZScoreAdjusted( 0.0 ),
      weightedAvgSD( 0.0 ),
      densCount( 0.0 )
  { /* */ }

  /**
   * Write a detailed header to <fp>
   */
  static void printHeaderVerbose( std::FILE *fp );

  /**
   * Write detailed information about this object to <fp>,
   *  with a column associating the hotspot with <chrom>
   */
  void printOutVerbose( const char* chrom, std::FILE* fp );

  /**
   * Write a header to <fp>
   */
  static void printHeader( std::FILE *fp );

  /**
   * Write information about this object to <fp>,
   *  with a column associating the hotspot with <chrom>
   */
  void printOut( const char* chrom, std::FILE* fp );
};

#endif /* HOTSPOT_HPP_ */
