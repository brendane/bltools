/*
 * Unix head command, but for biological sequence files
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
  
  TCLAP::CmdLine cmd("Equivalent of `head' for sequence files", ' ', "0.0");
  TCLAP::ValueArg<string> format_arg("f", "format",
                                     "Output format: fasta or fastq; fasta is default",
                                     false, "fasta", "fast[aq]", cmd);
  TCLAP::ValueArg<int> nlines_arg("n", "lines",
                                  "print the first n lines of each file",
                                  false, 10, "int", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  string format = format_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");
  int nlines = nlines_arg.getValue();
  unsigned look_ahead = 0;
  if(nlines < 0) {
    look_ahead = -1 * nlines;
  }
  
  SeqFileOut out_handle(cout, Fasta());
  //open(out_handle, cout);
  if(format == "fasta") {
    setFormat(out_handle, Fasta());
  } else if(format == "fastq") {
    setFormat(out_handle, Fastq());
  } else {
    cerr << "Unrecognized output format";
    return 1;
  }

  CharString id;
  queue<CharString> ids;
  CharString seq;              // CharString more flexible than Dna5String
  queue<CharString> seqs;
  CharString qual;
  queue<CharString> quals;
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
    // Fill up seqs, quals, ids until look_ahead is reached, then for
    // every additional record, pop one off of seqs, quals, and ids, and
    // push the new one on until the end of the file is reached.
    while(!seq_handle.atEnd() && ((look_ahead == 0 && nrecs_read < nlines) || look_ahead > 0)) {

      try {

        readRecord(id, seq, qual, seq_handle.sqh);

      } catch (Exception const &e) {

        cout << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.


      if(look_ahead > 0) {
        seqs.push(seq); ids.push(id); quals.push(qual);
        if(seqs.size() > look_ahead) {
          id = ids.front(); seq = seqs.front(); qual = quals.front();
          ids.pop(); seqs.pop(); quals.pop();
        } else {
          continue;
        }
      }

      // Write output
      try {
        writeRecord(out_handle, id, seq, qual);
        nrecs_read++;
      } catch (Exception const &e) {
        cerr << "Error writing output";
        seq_handle.close();
        return 1;
      }
      
    } // End single file reading loop

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
    }

  } // End loop over files
  close(out_handle);

  return 0;
}
