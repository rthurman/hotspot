/**
 * File: InputDataReader.hpp
 * Author: Scott Kuehn
 * Date: July 28, 2009
 * Version: $Id$
 *
 * Comments:
 *  Read input data, formatted as one or more lines
 *  of <string> <int> (a single space delimits the fields).
 *  Input is fetched a chromosome at a time, and stored in a vector
 */

#ifndef INPUTDATAREADER_HPP_
#define INPUTDATAREADER_HPP_

#include <string>
#include <vector>
#include <fstream>

#include "ByLine.hpp"

namespace hotspot
{

	class InputDataReader
	{
	public:

		/**
		 * Init the input data reader with a valid file name
		 */
		InputDataReader( const std::string& inputFileName );

		/**
		 *  Obtain the next chrom of information from file.
		 *    results are appended to <tags>.  This method is
		 *    a forward-only read operation, so each call returns
		 *    results from the next chromosome in the input file.
		 */
		int readNextChrom( std::vector< int >& tags );

		/**
		 * Returns the name of the most recently processed chromosome
		 */
		std::string currentChromName( ) const;

		/**
		 * Release system resources
		 */
		virtual ~InputDataReader( );

		/**
		 * Count and return the number of lines in the input file. This method
		 * does not affect the underlying position used by readChrom( )
		 */
		int numLines( ) const;

	private:
		std::string _inputFileName;
		std::ifstream _inputDataStream;
		std::string _currentChromName;

		// Iteration helpers
		int _nextTag;
		std::string _nextChromName;
		bool _hasMoreData;
		int _numChromsProcessed;
	};

} // namespace hotspot

#endif /* INPUTDATAREADER_HPP_ */
