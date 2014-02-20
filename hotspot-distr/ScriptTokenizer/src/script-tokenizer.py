#!/usr/bin/env python

# FILE: script-tokenizer.py
# AUTHOR: Scott Kuehn
# CREATE DATE: Wed, Sep 30 2009

usage = """USAGE: script-tokenizer.py [options] <token-file> <script-name...>
      Options:
        --clobber :  overwrite any existing files
        --output-dir=<dirName> :  place generated files in this location
        --execute-scripts : execute the newly generated files
        --break-on-error : for use with --execute-scripts, if a script returns
                           unsuccessfully, stop execution
        --help :  print this message and exit
        --version :  print the program version and exit
        """
version = "0.1"

import sys, os, tokenize
import ConfigParser 
from getopt import getopt
import shutil
    
def readCommandLine( args, pvars ):
    """Parse and evaluate program arguments in <args>.  Stores
        results in the dict backing <pvars> 
    """
                
    longOptions = [ 'clobber', 'execute-scripts', 'output-dir=', 'help', 'version', 'break-on-error' ]
    # Parse command line
    opts, pargs = getopt( args, [], longOptions )

    for optName, optValue in opts:

        if optName[ 0:2 ] != "--":
            optName = "--" + optName                                

        if optName == '--clobber':
            pvars[ 'clobber' ] = True
        elif optName == '--execute-scripts':
            pvars[ 'executeScripts' ] = True
        elif optName == '--output-dir':
            pvars[ 'outputDir' ] = optValue
        elif optName == '--break-on-error':
            pvars[ 'breakOnError' ] = True
        elif optName == '--help':
            print >> sys.stdout, usage
            sys.exit( 0 )
        elif optName == '--version':
            print >> sys.stdout, version
            sys.exit( 0 )

    if len( pargs ) < 2:
        print >> sys.stderr, usage
        sys.exit( 1 )
            
    pvars[ 'tokenFile' ] = pargs[ 0 ]
    if not os.access( pvars[ 'tokenFile' ], os.R_OK ):
        raise Exception, "Unable to access: " + str( pvars[ 'tokenFile' ] )
    pvars[ 'inputScripts' ].extend( pargs[ 1: ] )
    
def readTokenFile( tokenFileName ):
    """
        Parse and load a token file, return a dictionary
        of the parsed name/val pairs
    """

    config = ConfigParser.SafeConfigParser( )
    config.optionxform = str # Make config parser case sensitive
    config.readfp( open( tokenFileName ) )
    tokenSet = config.items( "script-tokenizer" )
    return dict( tokenSet )

def makeTokenizedFile( scriptName,  doClobber, outputDir, tokSuffix ):
    """ 
        Create the file that will store the tokenized output. Returns
        a tuple (file-name, file-handle).  The file created will be named:
            <outputDir>/<scriptName>.<tokSuffix>.  If outputDir is not specified
            the cwd will be used.  If doClobber is False, any existing output
            file matching the new name will be rewritten.
    """

    outputProcDir = outputDir
    if outputProcDir == None :
        outputProcDir = "."
    tokenizedScriptName = outputProcDir + "/" + os.path.basename(scriptName) + tokSuffix
    if not doClobber and os.path.exists(tokenizedScriptName ):
        raise Exception( "Error: " + tokenizedScriptName + " already exists" )
    if not os.path.exists(outputProcDir):
        os.makedirs(outputProcDir)
    tokenizedScript = open(tokenizedScriptName,  'w')
    shutil.copymode( scriptName, tokenizedScriptName )
    
    return ( tokenizedScriptName, tokenizedScript )

def tokenize( tokenizedFile, tokens, scriptName ):
    """
        Replace any matching tokens from <tokens> in <scriptName>, 
        writing the result to <tokenizedFile>.  Returns void.
    """
    
    result = []
    scriptFile = open( scriptName )
    for line in scriptFile:
        tokenizedLine = line
        for tokenName in tokens.keys( ):
            tokenizedLine = tokenizedLine.replace( tokenName, tokens[ tokenName ] )
        tokenizedFile.write( tokenizedLine )
    scriptFile.close( )

def main( ):

    # Program data
    pvars = { 
             'outputDir' : None,
             'clobber' :  False,
             'executeScripts' : False,
             'tokenFile' : None,
             'inputScripts' : [],
             'tokSuffix' : '.tok',
             'breakOnError' : False
             }
    tokenizedFileNames = []
    
    # Read input
    #try:
    readCommandLine( sys.argv[1:], pvars )
    tokens = readTokenFile( pvars['tokenFile'] )
    
    # Apply tokens to scripts
    for scriptName in pvars[ 'inputScripts' ]:    
        tokenizedFileName, tokenizedFile = makeTokenizedFile( scriptName, pvars['clobber'], \
                                                                  pvars['outputDir'], pvars['tokSuffix'] )
        tokenizedFileNames.append( tokenizedFileName )
        tokenize( tokenizedFile, tokens, scriptName )      
        tokenizedFile.close( )

    # Run tokenized files
    if( pvars['executeScripts'] ):
        for tokenizedFileName in tokenizedFileNames:
            ret = os.system( tokenizedFileName )
            if ret and pvars[ 'breakOnError' ]:
                print >> sys.stderr, "Script " + tokenizedFileName + \
                    " failed with exit code " + str( ret ) +". Aborting."
                sys.exit( 1 )
    
    #except Exception as e:
    #    print >> sys.stderr, str(e) + ". Aborting"
   #     sys.exit( 1 )

if __name__ == '__main__':
    main( )
    
