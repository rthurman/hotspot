/**
 * File: MappableCountsDataReader.cpp
 * Author: Scott Kuehn
 * Date: August 3, 2009
 * Version: $Id$
 *
 * Comments:
 *  An implementation of MappableCountsDataReader.hpp
 */

#include "MappableCountsDataReader.hpp"
#include "HotspotDefaults.hpp"
#include "ByLine.hpp"

#include <string>
#include <cstdio>
#include <istream>
#include <exception>
#include <vector>

namespace hotspot
{


	MappableCountsDataReader::MappableCountsDataReader( const std::string& inputFileName )
				: _inputFileName( inputFileName ), _recordNum( 0 )
	{
		_inputDataStream.open( _inputFileName.c_str( ) );
	}

	MappableCountsDataReader::~MappableCountsDataReader()
	{
		if( _inputDataStream )
		{
			_inputDataStream.close( );
		}
	}

	int MappableCountsDataReader::numLines( ) const
	{
		std::ifstream inf;
		inf.open( _inputFileName.c_str( ) );

		int numLines = 0;
		ByLine inputDataRecord;
		while( inf >> inputDataRecord ) numLines++;
		inf.close( );
		return numLines;
	}

	int MappableCountsDataReader::readChrom( const std::string& matchChromName, std::vector< int >& mappableCounts )
	{

		// Storage for input reads
		int start;
		char scannedChromName[ HotspotDefaults::MAX_CHROM_NAME_LEN + 1 ];
		int densityCount;

		// Iteration helpers
		int numChromLines = 0;
		ByLine inputDataRecord;

		while( _inputDataStream >> inputDataRecord  )
		{
			_recordNum++;
			// Process residual record, if present
			if( ! _nextLine.empty( ) )
			{
				const std::string *scanResidual = _nextLine.top( );
				_nextLine.pop( );
				// Scan and validate the current record
				int numScanned = std::sscanf( scanResidual->c_str( ), "%s %d %d",
						scannedChromName, &start, &densityCount );
				if( numScanned != 3 )
				{
					std::fprintf( stderr, "Error processing residual: input file %s contains a malformed entry on line %d\n",
							_inputFileName.c_str( ), _recordNum );
					std::fprintf( stderr, "Details: input line: %s, num records scanned: %d, num residual items: %ld\n", scanResidual->c_str( ), numScanned, _nextLine.size( )  );
					return -1;
				}
				if(! matchChromName.compare( scannedChromName ) )
				{
						// Record matches chrom
					mappableCounts.push_back( densityCount );
					numChromLines++;
				}
				delete scanResidual;
			}

			// Scan and validate the current record
			int numScanned = std::sscanf( inputDataRecord.c_str( ), "%s %d %d",
					scannedChromName, &start, &densityCount );
			if( numScanned != 3 )
			{
				std::fprintf( stderr, "Error: input file %s contains a malformed entry on line %d\n",
						_inputFileName.c_str( ), _recordNum );
				std::fprintf( stderr, "Details: input line: %s, num records scanned: %d\n", inputDataRecord.c_str( ), numScanned );
				return -1;
			}

			if( matchChromName.compare( scannedChromName ) )
			{
				// Record does not match chrom
				if (numChromLines == 0 )
				{
					continue;
				}
				else
				{
					//std::cerr << "Pushing: " << inputDataRecord.c_str( ) << std::endl;
					_nextLine.push( new std::string( inputDataRecord.c_str( ) ) );
					numChromLines++;
					break;
				}
			}
			else
			{
				// Record matches chrom
				mappableCounts.push_back( densityCount );
				numChromLines++;
			}
		}
		return numChromLines;
	}
}
