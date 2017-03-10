/*
 * Unix wc command, but for biological sequences.
 *
 */

#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include <seqan/seq_io.h>

#include <tclap/CmdLine.h>

#include <SeqFileInWrapper.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::queue;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

int main(int argc, char * argv[]) {
  
  TCLAP::CmdLine cmd("Equivalent of `wc' for sequence files", ' ', "0.0");
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");

  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  SeqFileInWrapper seq_handle;

  for(string& infile: infiles) {

    try {
        seq_handle.open(infile);
    } catch(Exception const &e) {
      cerr << "Could not open " << infile << endl;
      seq_handle.close();
      return 1;
    }
    
    int nrecs_read = 0;

    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq, seq_handle.sqh);
        nrecs_read++; 

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.
    } // End single file reading loop

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        return 1;
    }

    cout << infile << "\t" << nrecs_read << endl;

  } // End loop over files

  return 0;
}
