/**
 * File: HotspotDefaults.cpp
 * Author: Scott Kuehn
 * Date: July 25, 2009
 * Version: $Id$
 *
 * Comments:
 *  This file contains the initialization of non-int variables from
 *   HotspotDefaults.hpp
 */

#include "HotspotDefaults.hpp"

namespace hotspot
{
	const float HotspotDefaults::MAPPABLE_GENOME_SIZE = 2.55E9;
	const float HotspotDefaults::MINSD = 3.0;
	const char *HotspotDefaults::LIB_PATH = "input.lib";
	const char *HotspotDefaults::DENSITY_PATH = "mappable_site.counts";
}
