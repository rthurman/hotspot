/**
 * File: HotspotDefaults.hpp
 * Author: Scott Kuehn
 * Date: July 25, 2009
 * Version: $Id$
 *
 * Comments:
 *  This file is intended to serve as a centralized location for configurable
 *   constants available to the hotspot namespace.  These constants
 *   typically provide program parameter default values.
 */

#ifndef HOTSPOTDEFAULTS_HPP_
#define HOTSPOTDEFAULTS_HPP_

namespace hotspot
{
	class HotspotDefaults
	{
	public:
		// Behavior
		static const bool USE_DEFAULT_BACKGROUND_TAGS = true;
		static const bool USE_GENOME_DENS_WIN = false;
		static const bool USE_FUZZY_THRESHOLD = false;
		static const int FUZZY_SEED = 1;

		// Windowing
		static const float MINSD;
		static const int LOW_INTERVAL_WIDTH = 250;
		static const int HIGH_INTEVAL_WIDTH = 1000;
		static const int INTERVAL_INCREMENT = 100;
		static const int DENSITY_WIN_SMALL = 10000;
		static const int DENSITY_WIN = 50000;

		// Input
		static const char *LIB_PATH;
		static const char *DENSITY_PATH;

		// Limits
		static const int MAXLINE = 100000;
		static const float MAPPABLE_GENOME_SIZE;
		static const int MAX_CHROM_NAME_LEN = 127;
	};

} // namespace

#endif /* HOTSPOTDEFAULTS_HPP_ */

