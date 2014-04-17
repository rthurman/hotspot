/**
 * File: Cluster.cpp
 * Original Authors: Mike Hawrylycz, Bob Thurman
 * Contributors: Scott Kuehn, Eric Haugen
 * Date: 2003-2009
 * Version: $Id$
 *
 * Comments:  Implementation of Cluster.hpp, as well as the main( )
 *   program entry point.
 *   Most of the statements in this file belong to the declarations
 *   of the clustering calculation functions.  Consider Hotspot.hpp
 *   for details on the global variables.  All global variables are
 *   storage for program input parameters.
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <map>
#include <exception>

extern "C"
{
	#include <gsl/gsl_sf_gamma.h>
	#include <gsl/gsl_errno.h>
}

#include "Cluster.hpp"
#include "HotspotDefaults.hpp"
#include "Hotspot.hpp"
#include "InputDataReader.hpp"
#include "MappableCountsDataReader.hpp"

namespace hotspot
{
	std::string densitypath = HotspotDefaults::DENSITY_PATH;
	std::string libpath = hotspot::HotspotDefaults::LIB_PATH;
	int densityWinSmall = HotspotDefaults::DENSITY_WIN_SMALL;
	bool useGenomeDensWin = HotspotDefaults::USE_GENOME_DENS_WIN; // flag to use alternate density window if it gives a lower z-score
	int numGenomeDens, numLocalDens;
	double genomeDensZ, localDensZ;
	double mpblGenomeSize = HotspotDefaults::MAPPABLE_GENOME_SIZE;
	bool useDefaultBackgroundTags = HotspotDefaults::USE_DEFAULT_BACKGROUND_TAGS; // flag to determine what tag count to use for z-score genome-wide background calculations
	bool useFuzzyThreshold = HotspotDefaults::USE_FUZZY_THRESHOLD;
	int backgroundTotalTagCount;

    struct NotEnoughArgsException { /* */ };
    struct HelpException { /* */ };
    struct VersionException { /* */ };

    const std::string revision = "4.1.0";
    const std::string authors = "Mike Hawrylycz, Bob Thurman";
    const std::string contributors = "Scott Kuehn, Eric Haugen";
    const std::string citation = "Sam John et al., Chromatin accessibility pre-determines glucocorticoid receptor binding patterns, Nature Genetics 43, 264-268";

} // namespace

int main( int argc, char **argv )
{
    // We use a GSL special function to compute binomial cdfs.  
    // Turn off the GSL error handler so we can trap those errors ourself.
    gsl_set_error_handler_off();

    // Set defaults program input values
    hotspot::numSD = hotspot::HotspotDefaults::MINSD;  // minimum sd and intensity for detection
    hotspot::lowInt = hotspot::HotspotDefaults::LOW_INTERVAL_WIDTH;   // range and interval for scanning
    hotspot::highInt = hotspot::HotspotDefaults::HIGH_INTEVAL_WIDTH;
    hotspot::incInt = hotspot::HotspotDefaults::INTERVAL_INCREMENT;
    hotspot::fuzzySeed = hotspot::HotspotDefaults::FUZZY_SEED;
    hotspot::totaltagcount = 0;
    hotspot::densityWin = hotspot::HotspotDefaults::DENSITY_WIN;

    // Read in arguments
    try {
	hotspot::GetArgs( argc, argv );
    } catch (hotspot::HelpException& h) {
        hotspot::Usage( std::cout );
    } catch (hotspot::NotEnoughArgsException& n) {
        hotspot::Usage( std::cerr );
    } catch (hotspot::VersionException& v) {
        hotspot::Version();
    }

	// Fetch input data
	hotspot::InputDataReader inputDataReader( hotspot::libpath );
	hotspot::totaltagcount = inputDataReader.numLines( );
	std::cout << "TotalTagCount: " << hotspot::totaltagcount << std::endl;
	hotspot::MappableCountsDataReader mappableCountsDataReader( hotspot::densitypath );
	if( hotspot::useDefaultBackgroundTags )
	{
		hotspot::backgroundTotalTagCount = hotspot::totaltagcount;
	}
    if( hotspot::useFuzzyThreshold )
    {
    	std::srand( hotspot::fuzzySeed );
    }

	bool headerPrinted = false;
	std::vector< int > inputData;
	std::vector< int > mappableCounts;

	// Main processing loop: each pass considers each chromosome in the input data set
	while( inputDataReader.readNextChrom( inputData ) > 0 )
    {
		std::cerr << "Processing chrom: " << inputDataReader.currentChromName( ) << std::endl;

		// get counts of 'background' mappable K-mers on each 50kb interval on that chromosome
    	int numRead = mappableCountsDataReader.readChrom(inputDataReader.currentChromName( ), mappableCounts );
    	if( numRead < 0 )
    	{
    		std::cerr << "Error reading background file. Aborting" << std::endl;
			exit( EXIT_FAILURE );
    	}

    	// Compute the hot spots and filter them
   		std::cerr << "Compute Hot Spots " << std::endl;
       	std::map< int, Hotspot* > hotspots;
  		hotspot::ComputeHotSpots( inputData, hotspot::lowInt, hotspot::highInt, hotspot::incInt, hotspots );
   		std::cerr << "Filter Hot Spots " << std::endl;
       	std::map< int, Hotspot* > filteredHotspots;
   		hotspot::FilterHotspots( hotspots, filteredHotspots );

   		// Calculate cluster size, and other hotspot statistics
   		if( hotspot::useGenomeDensWin)
   		{
   			// These vars are updated (side-effect) of the ClusterSize routine
   			hotspot::numGenomeDens = 0; hotspot::numLocalDens = 0;
   			hotspot::genomeDensZ = 0.0; hotspot::localDensZ = 0.0;
   		}
   		std::cerr << "Cluster Size" << std::endl;
   		hotspot::ClusterSize( inputData, hotspot::densityWin, filteredHotspots, mappableCounts );

   		if( hotspot::useGenomeDensWin )
   		{
   			std::cerr << hotspot::numGenomeDens
						<< " clusters scored using genome-wide density, avg. z = "
						<< hotspot::genomeDensZ / hotspot::numGenomeDens
						<< "; " << hotspot::numLocalDens
						<< " scored using local density, avg. z = "
						<< hotspot::localDensZ / hotspot::numLocalDens
						<< std::endl;
   		}

   		// Summarize the results of this chromosome. Reset temp data structures
   		std::map< int, Hotspot* >::iterator iter;
   		for( iter = hotspots.begin( ); iter != hotspots.end( ); ++iter )
   		{
   			delete iter->second;
   		}
   		if(! headerPrinted )
   		{
   			Hotspot::printHeader( hotspot::fpout );
   			headerPrinted = true;
   		}
   		for( iter = filteredHotspots.begin(); iter != filteredHotspots.end(); ++iter )
   		{
   			//TODO remove costly char lookup to one-time lookup outside loop
   			iter->second->printOut( inputDataReader.currentChromName().c_str( ), hotspot::fpout );
   			delete iter->second;
   		}

   		std::cerr << "Chrom summary: " << filteredHotspots.size( ) << std::endl;
   		std::fflush( hotspot::fpout );
   		mappableCounts.clear( );
   		inputData.clear( );

    }  // end loop over all chromosomes

	// Release open resources
    if( hotspot::fpout )
    {
    	std::fclose( hotspot::fpout );
	}

    std::exit( EXIT_SUCCESS );
}

namespace hotspot
{

	double ComputeHotSpots( const std::vector< int>& inputData, int winLow,	int winHigh,
							int winInc, std::map< int, Hotspot* >& hotspots )
	{
		/* computes an estimate of the discrepancy using the class
		   of 1-dimensional intervals of width winLow to winHigh  */

		genomeSize = mpblGenomeSize; // RET:  changed from EDH's value of 3.0E9

		double disc = 0.0;
		int wincount = 0;
		for (int winsize = winLow; winsize <= winHigh; winsize += winInc, wincount++ )
		{
			unsigned int startmarker=0; // market to the first point that can be in window
			double prob = winsize / genomeSize;
			double mean = prob * totaltagcount;  // RET: adjust for sampling fraction

			double sd = std::sqrt(prob*(1-prob)*totaltagcount);  // RET: adjust for sampling fraction
			double detectThresh = 1 + mean + numSD * sd;  // includes offset of 1 as we are centering on clones

			//    cout << "Comparative Probs " << endl;
			//    cout << "Window size " << winsize << endl;
			//    cout << "Mean " << mean <<  endl;
			//    cout << "sd " << sd <<  endl;
			//    cout << "detect " << detectThresh  << endl << endl;

			for( unsigned int i = 0; i < inputData.size( ); i++ )
			{
				int contained = 0;
				double clonePosAvg = 0.0;

				// center interval on lib point, assume sorted library points
				double leftEnd  = inputData[i] - winsize/2.0;
				double rightEnd = inputData[i] + winsize/2.0;
				// std::printf("leftEnd %f rightEnd %f\n",leftEnd,rightEnd);
				if ((leftEnd > inputData[startmarker]) && (startmarker<inputData.size( ))) { startmarker++; }
				unsigned int walkmarker = startmarker;

				while ((inputData[walkmarker] <= rightEnd)&&(walkmarker<inputData.size( ))) {
					if ((inputData[walkmarker] >= leftEnd) && (inputData[walkmarker] <= rightEnd)) {
						contained++;
						clonePosAvg += inputData[walkmarker];
					}
					walkmarker++;
				}
				if (contained > 0)
				{
					clonePosAvg = clonePosAvg / contained;
					if ( useFuzzyThreshold )
					  {
					    double detectThresh = 1 + mean + numSD * sd;  // includes offset of 1 as we are centering on clones
					    if (std::fabs(contained - detectThresh) <= 0.5)
					      {
						double randVal = std::rand( ) /  static_cast< double >( RAND_MAX );
						detectThresh += ( randVal - 0.5 );
					      }
					  }
				}

				double contFrac = contained /(double)totaltagcount; // RET: adjust for sampling fraction
				double diff = std::fabs(contFrac-prob);

				if (contained > detectThresh)
				{
					// determine number of sd's corresponding to this intensity.
					// algorithm finds largest window over range containing anomaly.
					double currSD = (contained - 1 - mean) / sd;

					// Make a new hotspot for this index
					if( hotspots.count( i ) == 0 )
					{
						Hotspot* h = new Hotspot;
						hotspots[ i ] = h;
					}
					hotspots[i]->densCount += 1;
					hotspots[i]->weightedAvgSD += currSD;
					hotspots[i]->averagePos = static_cast< int >(clonePosAvg + 0.5);
					hotspots[i]->maxWindow = winsize;
				}
				if (diff > disc) {disc=diff;}
			}  // over all clones
		}  // over all window sizes

		std::map< int, Hotspot* >::iterator iter;
		for( iter = hotspots.begin( ); iter != hotspots.end( ); ++iter )
		{
			iter->second->weightedAvgSD /= iter->second->densCount;
		}

		std::cerr << "Completing HotSpot Identification" << std::endl;
		return disc;
	}

	void FilterHotspots( const std::map<int, Hotspot* >& hotspots,
				std::map< int, Hotspot* >& filteredHotspots )
	{
		// Iteration/windowing helper vars
	        // TODO: filterCluster and lastCenter might not be properly initialized
		int filterCluster = 0, lastCenter = 0;
		double averageCount = 0.0;
		double averageSd = 0.0;
		double adjustedCenter = 0.0;
		int numFilteredHotspots = 0;
		bool firstPassInit = true;

		// Considering each known hotspot, create filteredHotspots
		std::map<int, Hotspot*>::const_iterator iter;
		for( iter = hotspots.begin( ); iter != hotspots.end( ); ++iter )
		{
			Hotspot* currHotspot = iter->second;
			int tagNum = iter->first;

			// Init first cluster
			if( firstPassInit )
			{
				Hotspot* h = new Hotspot;
				h->densCount = currHotspot->densCount;
				h->averagePos = currHotspot->averagePos;
				h->weightedAvgSD = currHotspot->weightedAvgSD;
				h->filterIndexLeft = tagNum;
				h->filterIndexRight = tagNum;
				h->filterWidth = currHotspot->maxWindow;
				filteredHotspots[ numFilteredHotspots ] = h;

				lastCenter = currHotspot->averagePos;
				adjustedCenter = static_cast< double >( lastCenter );
				averageCount = currHotspot->densCount;
				averageSd = currHotspot->weightedAvgSD;
				filterCluster = 1;
				firstPassInit = false;
				continue;
			}

			if( std::abs( lastCenter - currHotspot->averagePos ) < currHotspot->maxWindow )
			{
				// same cluster
				filteredHotspots[ numFilteredHotspots ]->filterIndexRight = tagNum; // shift over right boundary
				filteredHotspots[ numFilteredHotspots ]->filterWidth += currHotspot->maxWindow;

				// Adjust windowing and filtering helper vars
				adjustedCenter += static_cast< double >( currHotspot->averagePos );
				averageCount += currHotspot->densCount;
				averageSd += currHotspot->weightedAvgSD;
				filterCluster++;
			}
			else // new cluster
			{
				// Note that this clause takes into effect if a continuing cluster gets too big, OR if a new
				// cluster has started with intervening denscount = 0 and the new cluster centroid is greater than maxwindow[i]
				// from the old centroid.
				// Before starting a new cluster move the centroid to the average position of
				// assigned clones.
				filteredHotspots[ numFilteredHotspots ]->averagePos = static_cast< int >( adjustedCenter / filterCluster );
				filteredHotspots[ numFilteredHotspots ]->filterWidth /= filterCluster;
				filteredHotspots[ numFilteredHotspots ]->densCount = averageCount / filterCluster;
				filteredHotspots[ numFilteredHotspots ]->weightedAvgSD = averageSd / filterCluster;
				if(filteredHotspots[ numFilteredHotspots ]->filterWidth > highInt*2)
				  {
				    std::cerr << "else: filterwidth=" << filteredHotspots[ numFilteredHotspots ]->filterWidth << ", filterCluster=" << filterCluster << std::endl;;
				  }

				// then increment and initialize a new cluster
				numFilteredHotspots++;

				Hotspot* h = new Hotspot;
				filteredHotspots[ numFilteredHotspots ] = h;
				h->filterIndexLeft = tagNum;
				h->filterIndexRight = tagNum;
				h->filterWidth = currHotspot->maxWindow;
				lastCenter =  currHotspot->averagePos;
				adjustedCenter = static_cast< double >( lastCenter );
				averageCount = currHotspot->densCount;
				averageSd = currHotspot->weightedAvgSD;
				filterCluster = 1;  // reset to 1 item in cluster
			}
		}
		// finish last cluster
		if( filterCluster > 0 ) // need this test in case there were no hotspots
		  {
		    filteredHotspots[ numFilteredHotspots ]->averagePos = static_cast< int >( adjustedCenter / filterCluster );
		    filteredHotspots[ numFilteredHotspots ]->densCount = averageCount / filterCluster;
		    filteredHotspots[ numFilteredHotspots ]->weightedAvgSD = averageSd / filterCluster;
		    filteredHotspots[ numFilteredHotspots ]->filterWidth /= filterCluster;
		  }
	}

  void ClusterSize( const std::vector<int>& inputData, int densityWin, std::map< int, Hotspot* >& filteredHotspots,
		    const std::vector< int >& mappableCounts )
  {
    // finally go through and determine the number of library clones contained
    // in filterwidth, also get the maximum inter-cluster width

    // TODO: leftindex might not be properly initialized
    int contcount,leftindex = -1, rightindex = -1,leftdens,rightdens,dencount;
    double leftcent, rightcent;
    int halfDensityWin = densityWin / 2;
    /* double probz,meanz,sdz; */

    std::map<int, Hotspot* >::iterator iter;
    for( iter = filteredHotspots.begin( ); iter != filteredHotspots.end( ); ++iter)
      {
	Hotspot *currHotspot = iter->second;

	leftcent  = currHotspot->averagePos - currHotspot->filterWidth / 2;
	rightcent = currHotspot->averagePos + currHotspot->filterWidth / 2;
	leftdens = currHotspot->averagePos - halfDensityWin;
	rightdens = currHotspot->averagePos + halfDensityWin;

	contcount=0;
	dencount=0;
	unsigned int libIdx;
	for (libIdx = 0; libIdx < inputData.size(); libIdx++)
	  {   // for every tag on that chrome
	    if ((inputData[libIdx] >= leftcent) && (inputData[libIdx] <= rightcent))
	      {
		if (contcount==0) 
		  {
		    leftindex=libIdx;
		  }
		rightindex=libIdx;
		contcount++;
	      }
	    if ((inputData[libIdx] >= leftdens) && (inputData[libIdx] <= rightdens))
	      {
		if(dencount == 0)
		  {
		    currHotspot->filterDensIndexLeft = libIdx;
		  }
		dencount++;
	      }
	    if (inputData[libIdx] > rightdens && inputData[libIdx] > rightcent) 
	      {
		break;
	      }
	  }
	currHotspot->filterDensIndexRight = libIdx - 1;
	currHotspot->filterSize = contcount;
	currHotspot->filterDist = inputData[rightindex] - inputData[leftindex] + 1; // changed to add 1 -- RET
	currHotspot->minSite = inputData[currHotspot->filterIndexLeft];
	currHotspot->maxSite = inputData[currHotspot->filterIndexRight];

	int uniquelyMappableSitesInWindow = countMappableSites( currHotspot->averagePos, densityWin, densityWinSmall,
								mappableCounts );


	// Following counts tags in the nearest 50kb window starting on a 10kb boundary, to
	// match the windows used in countMappableSites.
	int dencount2 = countDensity2( currHotspot->averagePos, inputData );
	currHotspot->filteredZScoreAdjusted = calculateZScore( lround( currHotspot->filterWidth), currHotspot->filterSize,
							       densityWin, dencount2, uniquelyMappableSitesInWindow );
      }
    return;
  }

	int countMappableSites( int base, int densityWin, int densityWinSmall,
							const std::vector< int>& mappableCounts )
	{
		int sum = 0;
		int subWindows = densityWin / densityWinSmall;
		int start = ( base / densityWinSmall ) - ( subWindows / 2 );
		if ( start < 0 ) start = 0;
		for( int i = 0; i < subWindows; ++i)
		{
			int winCount = 0;

			//TODO fix storage size discrepancy
			if( (start + i ) >= static_cast< int >( mappableCounts.size( ) ) )
			{
				winCount = densityWinSmall;
			}
			else
			{
				winCount = mappableCounts[start + i];
			}

			if ( winCount > densityWinSmall )
			{
				std::cerr <<  "Warning: " << winCount << " > " << densityWinSmall << std::endl;
				winCount = densityWinSmall;
			}
			sum += winCount;
		}

		if( sum > densityWin )
		{
			std::cerr << "Error: " << sum << " > " << densityWin << std::endl;
			sum = densityWin;
		}
		return sum;
	}

	// Called on each Hotspot,during ClusterSize() operation.
	// Counting the number of tags in the window marked by [leftdens <====> rightdens]
	// Optimize by starting at current obs, look left until out-of-bounds, look right until out-of-bounds

	int countDensity2( int base, const std::vector< int >& inputData ) {
		int subWindows = densityWin / densityWinSmall;
		int start = (base / densityWinSmall) - (subWindows / 2);
		int leftdens = start * densityWinSmall;
		int rightdens = leftdens + densityWin - 1;
		int dencount = 0;
		for (unsigned int j = 0; j < inputData.size( ); j++) {   // for every tag on that chrome
			if ((inputData[ j ] >= leftdens) && (inputData[ j ] <= rightdens)) {
				dencount++;
			}
			if (inputData[ j ] > rightdens) {break;}
		}
		return dencount;
	}

	/* To adjust for local mappable K-mer density, just reduce the densityWindowSize
	   (from 50,000 etc) to however many sites in that window are uniquely mappable
	   by K-mers.   But that throws off the probZ calc in a bad way.
	   I use adjustedNumSitesInCluster as an estimate of how many sites I would observe
	   in the density window if all the sites were actually mappable.
	*/
	double calculateZScore( int basesSpannedByCluster, int numSitesInCluster,
			int densityWindowSize, int numSitesInDensityWindow, int numMappableSites )
	{

		int fewEnoughSites = 2; //densityWinSmall; //densityWindowSize / 2;
		// Now using exact counts from the interval, so we shouldn't need this anymore:
		if (numMappableSites < fewEnoughSites)
		{
			// let's put a limit on the adjustment
			std::cerr << "Warning: only " << numMappableSites
						<< " of " << densityWindowSize
						<< " are mappable, increasing to "
						<< fewEnoughSites << "..." << std::endl;
			numMappableSites = fewEnoughSites;
		}

		// Adjust probZ using mappable sites.
		double probZ = ((double)basesSpannedByCluster) / (double)numMappableSites;
		double meanZ = numSitesInDensityWindow * probZ;
		double sdZ = std::sqrt(numSitesInDensityWindow * probZ * (1-probZ));
		double zScore = (numSitesInCluster - meanZ) / sdZ;
		// p-value = P(X >= k), where k = number tags in the cluster
		// = 1 - P(X < k) = 1 - P(X <= k - 1).  For the binomial distribution, the cdf is given by
		// the normalized (or regularized) incomplete beta function, I_x(a,b) (see wikipedia), and
		// P(X <= k-1) = I_(1-p)(n - (k-1), k), where p = probability, n = number of tags in the
		// density window.

		// Number mappable bases genome-wide, from ~rthurman/proj/dhs-peaks/results/fdr/fdr.R
		if ( !useGenomeDensWin )
		  return zScore;
		else
		{
			double probZgw = ((double)basesSpannedByCluster) / mpblGenomeSize;
			double meanZgw = backgroundTotalTagCount * probZgw;
			double sdZgw = std::sqrt(backgroundTotalTagCount * probZgw * (1-probZgw));
			double zScoregw = (numSitesInCluster - meanZgw) / sdZgw;
			//std::printf("zScoregw = %f, zScore = %f\n", zScoregw, zScore);
			gsl_sf_result beta_result;
			double pVal, pValgw;
			pVal = 0;
			pValgw = 0;
			int status = gsl_sf_beta_inc_e(numSitesInCluster, numSitesInDensityWindow - numSitesInCluster + 1, probZ, &beta_result);
			if(status == 15){
				//std::printf("Warning: underflow in pval computation, setting pVal = 0.\n");
				pVal = 0.0;
			}
			else if(status){
				//std::printf("gsl_sf_beta_inc_e error: status = %d, error = %s, setting pVal = 1.\n", status, gsl_strerror(status));
				pVal = 1.0;
			}
			else
			{
				pVal = beta_result.val;
			}
			status = gsl_sf_beta_inc_e(numSitesInCluster, backgroundTotalTagCount - numSitesInCluster + 1, probZgw, &beta_result);
			if(status == 15){
				//std::printf("Warning: underflow in pval computation, setting pValgw = 0.\n");
				pValgw = 0.0;
			}
			else if(status){
				//std::printf("gsl_sf_beta_inc_e error: status = %d, error = %s, setting pValgw = 1.\n", status, gsl_strerror(status));
				pValgw = 1.0;
			}
			else
			{
				pValgw = beta_result.val;
			}
			//std::printf("zScoregw = %f, zScore = %f, pVal = %g, pValgw = %g\n", zScoregw, zScore, pVal, pValgw);
			if (zScoregw < zScore){
				//std::printf("Genome-wide density used.\n");
				genomeDensZ += zScoregw;
				numGenomeDens++;
				return zScoregw;
			}else{
				numLocalDens++;
				localDensZ += zScore;
				return zScore;
			}
		}
	}

    void Usage(std::ostream& outputStream) 
    {
        std::string msg  = "HotSpot5 Usage:";
        msg += "\n    -h | --help (show this usage message)";
        msg += "\n    -v | --version (show application version and other details)";
        msg += "\n    -range <int> <int> <int> (lower upper increment windows)";
        msg += "\n    -densWin <int> (background window)";
        msg += "\n    -minsd <float> (minimum for anomaly)";
        msg += "\n    -fuzzy (flag to randomly adjust the hotspot selection threshold by 0-0.5)";
        msg += "\n    -fuzzy-seed <int> (for use with fuzzy, seed the random number gen. Default = 1 )";
        msg += "\n    -i <file-name> (input library file, must be in lexicographical sorted order)";
        msg += "\n    -k <file-name> (input K-mer density file, must be in lexicographical sorted order)";
        msg += "\n    -o <file-name> (output file for results)";
        msg += "\n    -gendw (flag to use genome-wide density window if it gives lower z-score)";
        msg += "\n    -bckgnmsize <float> (for computing background - default=2.55E9)";
        msg += "\n    -bckntags <float> (for computing background - default=number of tags in library)";
        msg += "\n";

        outputStream << msg << std::endl;
        
        if (&outputStream == &std::cout)
            std::exit(EXIT_SUCCESS);
        
        std::exit(EXIT_FAILURE);
    }

    void Version()
    {
        std::cout << "Hotspot5 About:" << std::endl;
        std::cout << "  citation:      " << hotspot::citation << std::endl;
        std::cout << "  version:       " << hotspot::revision << std::endl;
        std::cout << "  authors:       " << hotspot::authors << std::endl;
        std::cout << "  contributors:  " << hotspot::contributors << std::endl;
        std::exit(EXIT_SUCCESS);
    }

	/**
	 * Process program input and initialize this instance of Hotspot
	 */
	void GetArgs (int argc, char **argv)
	{
            if (argc < 2) {
                throw(NotEnoughArgsException());
            }

	  for( int i = 1 ; i < argc; i++ )
	  {
              if (( std::strcmp( argv[ i ], "-h" ) == 0 ) || ( std::strcmp( argv[ i ], "--help" ) == 0 )) {
                  throw(HelpException());
              }
              else if (( std::strcmp( argv[ i ], "-v" ) == 0 ) || ( std::strcmp( argv[ i ], "--version" ) == 0 )) {
                  throw(VersionException());
              }
              else if( std::strcmp( argv[ i ], "-range" ) == 0 )
		{
		  lowInt = std::atoi( argv[ i + 1 ] );
		  highInt = std::atoi( argv[ i + 2 ] );
		  incInt = std::atoi( argv[ i + 3 ] );
		  i += 3;
		}
		else if( std::strcmp( argv[ i ],"-minsd" ) == 0 )
		{
		  numSD = std::atof( argv[ i + 1 ] );
		  i++;
		}
		else if( std::strcmp( argv[ i ], "-o" ) == 0 )
		{
		  std::string outfile = argv[ i + 1 ];
		  fpout = std::fopen( outfile.c_str( ), "w" );
		  if( fpout == NULL )
		  {
			  std::cerr << "Error: unable to access " << outfile << std::endl;
			  std::exit( EXIT_FAILURE );
		  }
		  i++;
		}
		else if( std::strcmp(argv[ i ], "-i" ) == 0 )
		{
		  libpath = argv[ i + 1 ];
		  if( access( libpath.c_str( ), R_OK ) )
		  {
			  std::cerr << "Error: unable to access " << libpath << std::endl;
			  std::exit( EXIT_FAILURE);
		  }
		  i++;
		}
		else if( std::strcmp( argv[ i ], "-k" ) == 0 )
		{
			densitypath = argv[ i + 1 ];
			if( access( densitypath.c_str( ), R_OK ) )
			{
				std::cerr << "Error: unable to access " << densitypath << std::endl;
				std::exit( EXIT_FAILURE );
			}
			i++;
		}
		else if( std::strcmp( argv[ i ], "-fuzzy") == 0 )
		{
			useFuzzyThreshold = true;
		}
		else if( std::strcmp( argv[ i ], "-fuzzy-seed") == 0 )
		{
			fuzzySeed = std::atoi( argv[ i + 1 ] );
			i++;
		}
		else if( std::strcmp( argv[ i ], "-gendw" ) == 0 )
		{
		  useGenomeDensWin = true;
		}
		else if( std::strcmp( argv[ i ], "-densWin" ) == 0 )
		{
		  densityWin = std::atoi( argv[ i + 1 ] );
		  i++;
		}
		else if( std::strcmp( argv[ i ], "-bckgnmsize" ) == 0 )
		{
		  mpblGenomeSize = std::atof( argv[ i + 1 ] );
		  i++;
		}
		else if( std::strcmp( argv[ i ], "-bckntags" ) == 0 )
		{
		  useDefaultBackgroundTags = false;
		  backgroundTotalTagCount = std::atoi( argv[ i + 1 ] );
		  i++;
		}
		else
		{
			std::cerr << "Unrecognized option: " << argv[ i ] << ". Aborting." << std::endl;
			std::exit( EXIT_FAILURE );
		}
	  }

	  hotspot::outputFileName = libpath.c_str( );

	  if (!fpout)
	  {
		  std::cerr << "Output file required" << std::endl;
		  std::exit( EXIT_FAILURE);
	  }
	  return;
	}
} // namespace
