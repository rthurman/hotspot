/**
 * File: MappableCountsDataReader.hpp
 * Author: Scott Kuehn
 * Date: August 3, 2009
 * Version: $Id$
 *
 * Comments:
 *  Read background data.  Input is fetched a chromosome at a time, and
 *  stored in a vector.  An input data record is formatted as: <string> <int> <int>
 *  (a single space delimits fields).
 */

#ifndef MAPPABLE_COUNTS_DATA_READER_HPP_
#define MAPPABLE_COUNTS_DATA_READER_HPP_

#include <string>
#include <vector>
#include <stack>
#include <fstream>

#include "ByLine.hpp"

namespace hotspot
{

	class MappableCountsDataReader
	{
	public:
		/**
		 * Init the mappable counts data reader with a valid file name
		 */
		MappableCountsDataReader( const std::string& inputFileName );
		/**
		 * Release system resources
		 */
		~MappableCountsDataReader();

		/**
		 *  Obtain the next chrom of information from file.
		 *    results are appended to <tags>.  This method is
		 *    a forward-only read operation, so each call returns
		 *    results from the next chromosome in the input file.
		 */
		int readChrom( const std::string& chromName,  std::vector< int>& tags );

		/**
		 * Count and return the number of lines in the input file. This method
		 * does not affect the underlying position used by readChrom( )
		 */
		int numLines( ) const;

	private:
		std::string _inputFileName;
		int _recordNum;
		std::ifstream _inputDataStream;
		std::stack< const std::string* > _nextLine;
};

} // namespace hotspot

#endif /* MAPPABLE_COUNTS_DATA_READER_HPP_ */
