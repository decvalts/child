/****************************************************************************/
/**
**  @file tInputFile.cpp
**  @brief Member functions for class tInputFile.
**
**  (see tInputFile.h for a description of this class)
**
**  Greg Tucker, November 1997
**  Re-written, AD, July 2003
**
**  $Id: tInputFile.cpp,v 1.28 2003-07-31 13:13:00 childcvs Exp $
*/
/****************************************************************************/

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "tInputFile.h"
#include "../errors/errors.h"
#include "../tList/tList.h"

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
# include <fstream>
using namespace std;
#else
# include <iostream.h>
# include <fstream.h>
#endif

#define kCommentMark '#'


/****************************************************************************\
**
**  tKeyPair
**
**  a simple key-value class.
**  - 11/07/2003 AD
\****************************************************************************/

tKeyPair::tKeyPair(const char *key, const char *value) :
  key_(0), value_(0)
{
  setKey(key);
  setValue(value);
}

tKeyPair::tKeyPair(tKeyPair const &k) :
  key_(0), value_(0)
{
  setKey(k.key_);
  setValue(k.value_);
}

tKeyPair& tKeyPair::operator=(tKeyPair const&k){
  if (&k != this) {
    setKey(k.key_);
    setValue(k.value_);
  }
  return *this;
}

tKeyPair::~tKeyPair(){
  clear();
}

void tKeyPair::setKey(const char *s){
  clearKey();
  if (s == NULL) return;
  key_ = new char[strlen(s)+1];
  strcpy(key_,s);
}

void tKeyPair::setValue(const char *s){
  clearValue();
  if (s == NULL) return;
  value_ = new char[strlen(s)+1];
  strcpy(value_,s);
}

void tKeyPair::clear(){
  clearKey();
  clearValue();
}

void tKeyPair::clearKey(){
  delete [] key_; key_ = NULL;
}

void tKeyPair::clearValue(){
  delete [] value_; value_ = NULL;
}

/****************************************************************************\
**
**  setKeyWordTable
**
**  Set KeyWordTable
**  - 11/07/2003 AD
\****************************************************************************/
static
void setKeyWordTable(tArray< tKeyPair > &KeyWordTable,
		     tList< tKeyPair > &KeyWordList) {
  KeyWordTable.setSize(KeyWordList.getSize());
  int i = 0;
  tListIter< tKeyPair > ki(KeyWordList);
  for( tKeyPair *k=ki.FirstP(); !(ki.AtEnd()); k=ki.NextP(), ++i ) {
    KeyWordTable[i] = *k;
  }
}

/****************************************************************************\
**
**  various functions used for reading the input file
**
**  - 11/07/2003 AD
\****************************************************************************/
// read a line in headerLine and discards the any remaining characters
static
void readLine(char *headerLine, ifstream& infile){
  infile.getline( headerLine, kMaxNameLength );
  // still some characters left on the current line ?
  if ( !infile.eof() && !infile.bad() && infile.rdstate() & ios::failbit){
    // clear failbit
    infile.clear(infile.rdstate() & ~ios::failbit);
    // discard characters.
    char c;
    while( infile.get(c) && c != '\n');
  }
}

inline bool isComment(const char *headerLine){
  assert(headerLine != NULL);
  return BOOL(headerLine[0]==kCommentMark);
}

static
void skipCommentsAndReadValue(char *headerLine, ifstream& infile){
  do {
    readLine(headerLine, infile);
  } while( !infile.eof() && !infile.bad() &&
	   isComment(headerLine) );
}

// strip a string at the first ' ' or ':' found. 
static
void stripKey(char *key){
  char c;
  int i = 0;
  while((c = key[i]) != '\n'){
    if (c == ':' || c == ' '){
      key[i] = '\0';
      break;
    }
    ++i;
  }
}

static
void stripTrailingBlanks(char *s){
  const size_t len = strlen(s);
  int i = len - 1;
  while(i>=0){
    if (s[i] != ' '){
      s[i+1] = '\0';
      break;
    }
    --i;
  }
}


/* Read input file and build a list of key,value pair */
static
void readInputFile(ifstream& infile, tList< tKeyPair > &KeyWordList){
  char headerLine[kMaxNameLength];
  enum{ KEY, VALUE} state = KEY;
  tKeyPair aPair;

  for(;;){
    skipCommentsAndReadValue(headerLine, infile);
    if (infile.eof() || infile.bad())
      break;
    if (isComment(headerLine))
      goto fail;
    stripTrailingBlanks(headerLine);
    // if a blank line is reached, we stop.
    if (headerLine[0] == '\0')
      break;
    switch (state){
    case KEY:
      stripKey(headerLine);
      aPair.setKey(headerLine);
      state = VALUE;
      break;
    case VALUE:
      aPair.setValue(headerLine);
      assert(aPair.key() != NULL);
      KeyWordList.insertAtBack(aPair);
      aPair.clear();
      state = KEY;
      break;
    }
  }

  return;
 fail:
  cerr
    << "I expected to read a parameter or a value"
    "', but reached EOF first" << endl;
  ReportFatalError( "Error in input file" );
  /*NOTREACHED*/
}

/****************************************************************************\
**
**  tInputFile Constructor
**
**  Looks for a file called filename, opens it if found or generates an
**  error if not. Then reads the base name for output files and creates
**  a file called <filename>.inputs which will act as a log file for
**  parameters read. (This is often useful when the original input file
**  gets lost or modified).
**
**  Modifications:
**    - 2/02: error check for inoutfile added (GT)
**    - rewritten 11/07/2003 AD
**
\****************************************************************************/
tInputFile::tInputFile( const char *filename )
{
   ifstream infile;     // the input file

   // Open file
   infile.open( filename );
   if( !infile.good() )
   {
      cerr << "tInputFile::tInputFile: Unable to open '" << filename
	   << "'." << endl;
      ReportFatalError( "The file may not exist or may be mis-named." );
   }

   // Set KeyWordTable
   {
     tList< tKeyPair > KeyWordList;
     readInputFile(infile, KeyWordList);
     setKeyWordTable(KeyWordTable, KeyWordList);
   }
   infile.close();

   // write log File
   writeLogFile();
}


/****************************************************************************\
**
**  tInputFile::writeLogFile
**
**  Write log file.
**  - 11/07/2003 AD
\****************************************************************************/
void tInputFile::writeLogFile() const
{
  ofstream inoutfile;  // output file in which items are recorded
  char inoutname[kMaxNameLength];
  // Create log file for inputs
  ReadItem( inoutname, sizeof(kMaxNameLength), "OUTFILENAME" );
  strcat( inoutname, ".inputs" );
  inoutfile.open( inoutname );
  if( !inoutfile.good() )
    {
      cerr << "Unable to open '" << inoutname << "'.\n"
	"(Error generated in module tInputFile,"
	" function tInputFile::tInputFile( const char * ) )" << endl;
      ReportFatalError( "The specified path name may not exist.\n" );
    }
  // write header
  {
    const time_t now = time(NULL);
    char *p = ctime(&now);
    *strchr(p, '\n') = '\0';
    inoutfile << "# Created by CHILD on " << p << "." << endl;
  }
  // dump content of KeyWordTable
  {
    const int len = KeyWordTable.getSize();
     for( int i = 0; i < len; ++i ) {
       inoutfile << KeyWordTable[i].key() << '\n'
		 << KeyWordTable[i].value() << '\n';
     }
  }
  inoutfile.close();
}

/****************************************************************************\
**
**  tInputFile::findKeyWord
**
**  Find a tKeyPair in KeyWordTable. A match is found when the beginning of a
**  key is equal to the argument 'key'.
**  Returns notFound in case of failure.
**  - 11/07/2003 AD
\****************************************************************************/
int tInputFile::findKeyWord( const char *key ) const
{
  assert(key != NULL);
  const int len = KeyWordTable.getSize();
  const size_t sizeKey = strlen(key);

  for(int i=0; i < len; ++i){
    if (0 == strncmp(key, KeyWordTable[i].key(), sizeKey))
      return i;
  }
  return notFound;
}

/****************************************************************************\
**
**  tInputFile::ReadItem
**
**  Reads one parameter from the file. The format is assumed to be a line
**  of text that begins with the code "itemCode", followed by a line containing
**  the parameter to be read. The function is overloaded according to the
**  type of data desired (datType simply governs which overloaded function
**  will be called; it is not used by the routines).
**
**  Inputs:  datType -- dummy variable indicating the data type to be read
**                      (in the case of the string version, the string read
**                      is placed here)
**           itemCode -- string that describes the parameter to be read
**  Returns:  the item read (except in the case of the string version)
**  Modifications:
**    - revised to allow arbitrary ordering of items in infile and/or
**      ReadItem calls in code; routine searches through
**      list until it either finds the right itemCode or reaches EOF.
**      12/23/97 SL
**    - rewritten 11/07/2003 AD
**
\****************************************************************************/
static
void ReportNonExistingKeyWord(const char *itemCode){
  cerr << "Cannot find  '" << itemCode
       << "' in the input file." << endl;
  ReportFatalError( "Missing parameter in input file" );
}

int tInputFile::ReadItem( const int & /*datType*/, const char *itemCode ) const
{
  const int i = findKeyWord( itemCode );
  if (i == notFound)
    ReportNonExistingKeyWord( itemCode );
  return atoi(KeyWordTable[i].value());
}

long tInputFile::ReadItem( const long & /*datType*/, const char *itemCode ) const
{
  const int i = findKeyWord( itemCode );
  if (i == notFound)
    ReportNonExistingKeyWord( itemCode );
  return atol(KeyWordTable[i].value());
}

double tInputFile::ReadItem( const double & /*datType*/, const char *itemCode ) const
{
  const int i = findKeyWord( itemCode );
  if (i == notFound)
    ReportNonExistingKeyWord( itemCode );
  return atof(KeyWordTable[i].value());
}

// The size of 'theString' is 'len' including the trailing '\0'
void tInputFile::ReadItem( char * theString, size_t len,
			   const char *itemCode ) const
{
  assert(len>0);
  const int i = findKeyWord( itemCode );
  if (i == notFound)
    ReportNonExistingKeyWord( itemCode );
  strncpy(theString, KeyWordTable[i].value(), len);
  theString[len-1] = '\0';

  const size_t llen = strlen(theString);
  if (llen == 0) return;
  // strip trailing '\r' if we are dealing with a windows CR/LF text file
  if (theString[llen-1] == '\r')
    theString[llen-1] = '\0';
}
