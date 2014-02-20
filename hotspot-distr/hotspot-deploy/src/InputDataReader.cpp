/**
 * File: InputDataReader.cpp
 * Author: Scott Kuehn
 * Date: July 28, 2009
 * Version: $Id$
 *
 * Comments:
 *  this file contains an implementation of InputDataReader.hpp
 */

#include "InputDataReader.hpp"
#include "HotspotDefaults.hpp"
#include "ByLine.hpp"

#include <string>
#include <cstdio>
#include <istream>
#include <exception>
#include <vector>

namespace hotspot
{

	InputDataReader::InputDataReader( const std::string& inputFileName )
					: _inputFileName( inputFileName ), _currentChromName( "" )
					, _nextTag( -1 ), _nextChromName( "" ), _hasMoreData( false ),
					_numChromsProcessed( 0 )

	{
		_inputDataStream.open( _inputFileName.c_str( ) );
	}

	InputDataReader::~InputDataReader()
	{
		if( _inputDataStream )
		{
			_inputDataStream.close( );
		}
	}

	int InputDataReader::numLines( ) const
	{
		std::ifstream inf;
		inf.open( _inputFileName.c_str( ) );

		int numLines = 0;
		ByLine inputDataRecord;
		while( inf >> inputDataRecord ) numLines++;
		inf.close( );
		return numLines;
	}

	std::string InputDataReader::currentChromName( ) const
	{
		return _currentChromName;
	}

	int InputDataReader::readNextChrom( std::vector< int >& tags )
	{
		std::string temp;
		int tagLoc = -1;
		int numLines = 0;
		char scannedChromName[ HotspotDefaults::MAX_CHROM_NAME_LEN + 1 ];
		std::string chromName;

		ByLine inputDataRecord;

		// Record residual data from last record scan
		if( _hasMoreData )
		{
			tags.push_back( _nextTag );
			_currentChromName = _nextChromName;
			_numChromsProcessed++;
			_hasMoreData = false;
		}

		while( _inputDataStream >> inputDataRecord  )
		{
			// Scan and validate an input record
			int numScanned = std::sscanf( inputDataRecord.c_str( ), "%s %d",
					scannedChromName, &tagLoc );
			if( numScanned != 2 )
			{
				std::fprintf( stderr, "Error: input file %s contains a malformed entry on line %d\n",
						_inputFileName.c_str( ), numLines + 1 );
				return -1;
			}
			chromName = scannedChromName;

			if( _numChromsProcessed == 0 )
			{
				_currentChromName = chromName;
				_numChromsProcessed++;
			}
			else
			{
				// This record belongs to the next chrom.  Save the data for the next caller
				if( chromName.compare( _currentChromName ) != 0 )
				{
					_nextTag = tagLoc;
					_nextChromName = chromName;
					_hasMoreData = true;
					return numLines;
				}
			}

			tags.push_back( tagLoc );
			numLines++;
		}
		return numLines;
	}
}
