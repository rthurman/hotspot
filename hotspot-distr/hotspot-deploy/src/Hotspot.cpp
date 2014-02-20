/**
 * File: Hotspot.cpp
 * Author: Scott Kuehn
 * Date: July 25, 2009
 * Version: $Id$
 *
 * Comments:
 *  this file contains an implementation of Hotspot.hpp
 */

#include <cstdio>
#include <iostream>
#include "Hotspot.hpp"

void Hotspot::printHeaderVerbose( std::FILE *outputFile )
{
  if( outputFile == NULL )
    {
      return;
    }
  std::fprintf( outputFile,
		"Chrome\tPosition\tIntensity\tAvgSd\tClusterSize\tInterDist\tZScore\tWindowWidth\tMinSite\tMaxSite\tZScore2\tMinSiteInd\tMaxSiteInd\tMinDensWinInd\tMaxDensWinInd\n" );
}


void Hotspot::printOutVerbose( const char *chrom, std::FILE *fp )
{
  if( fp == NULL || chrom == NULL )
    {
      return;
    }
  std::fprintf( fp, "%s\t%d\t%5.2f\t%5.2f\t%d\t%d\t%f\t%5f\t%d\t%d\t%f\t%d\t%d\t%d\t%d\n",
		chrom, averagePos, densCount,
		weightedAvgSD, filterSize, filterDist,
		filteredZScore, filterWidth, minSite, maxSite, filteredZScoreAdjusted,
		filterIndexLeft, filterIndexRight, filterDensIndexLeft,
		filterDensIndexRight );
}

void Hotspot::printHeader( std::FILE *outputFile )
{
  if( outputFile == NULL )
    {
      return;
    }
  std::fprintf( outputFile, "Chrome\tPosition\tClusterSize\tInterDist\tWindowWidth\tMinSite\tMaxSite\tZScore2\n" );
}

void Hotspot::printOut( const char *chrom, FILE *fp )
{
  if( fp == NULL || chrom == NULL )
    {
      return;
    }
  std::fprintf( fp, "%s\t%d\t%d\t%d\t%5f\t%d\t%d\t%f\n",
		chrom, averagePos, filterSize, filterDist, filterWidth, minSite, maxSite, filteredZScoreAdjusted );
}
