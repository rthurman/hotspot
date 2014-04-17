/**
 * File: Cluster.hpp
 * Original Authors: Mike Hawrylycz, Bob Thurman
 * Contributors: Scott Kuehn, Eric Haugen
 * Date: 2004-2009
 * Version: $Id$
 *
 * Comments:
 *   This file contains the declaration of functions and data residing in the
 *   hotspot namespace.  The data is hotspot global run-time paramaters, and is necessary
 *   for running an instance of the hotspot program.  The functions provide ability
 *   to read input arguments, collect background data, and perform the hotspot clustering
 *   calculations.
 */

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include <map>
#include <vector>
#include "HotspotDefaults.hpp"
#include "Hotspot.hpp"

namespace hotspot
{
    // Print version statement
    void Version();

    // Print usage statement
    void Usage(std::ostream& outputStream);

	// Process program Input
	void GetArgs (int argc, char **argv);

	// Background data Input
	int countMappableSites( int base, int densityWin, int densityWinSmall,
							const std::vector< int >& mappableCounts );
	// Clustering calculations
	double calculateZScore( int basesSpannedByCluster, int numSitesInCluster,
							int densityWindowSize, int numSitesInDensityWindow, int mappableSites );
	int countDensity2( int base, const std::vector< int >& inputData );
	double ComputeHotSpots( const std::vector<int>& inputData, int winLow, int winHigh,
							int winInc, std::map< int, Hotspot* >& hotspots );
	void FilterHotspots( const std::map<int, Hotspot* >& hotspots,
						 std::map< int, Hotspot* >& filteredHotspots );
	void ClusterSize( const std::vector<int>& inputData, int densityWin, std::map< int, Hotspot* >& filteredHotspots,
					  const std::vector< int >& mappableCounts );

	// Input arguments
	int totaltagcount, densityWin, fuzzySeed;
	int lowInt, highInt, incInt;
	double numSD, genomeSize;
	std::string outputFileName;
	FILE *fpout; // file associated with outputFileName
}

#endif // __CLUSTER_H__
