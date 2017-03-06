/*
 * Equivalent of tail for biological sequences.
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
using std::stoi;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

int main(int argc, char * argv[]) {
  
  TCLAP::CmdLine cmd("Equivalent of `tail' for sequence files", ' ', "0.0");
  TCLAP::ValueArg<string> format_arg("o", "output-format",
                                     "Output format: fasta or fastq; fasta is default",
                                     false, "fasta", "fast[aq]", cmd);
  TCLAP::ValueArg<string> nlines_arg("n", "lines",
                                     "print the last n lines of each file or all lines but the first +n",
                                     false, "10", "[+]int", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  string format = format_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");
  string nlines_string = nlines_arg.getValue();
  int nskip = 0;
  int nlines = 0;
  if(nlines_string[0] == '+') {
    nlines_string.erase(0, 1);
    nskip = stoi(nlines_string);
    nlines = 0;
  } else {
    nlines = stoi(nlines_string);
    nskip = 0;
  }
  if(nlines < 0 || nskip < 0) {
    cerr << "Can't have a negative number of lines" << endl;
    return 1;
  }
  
  SeqFileOut out_handle(cout, Fasta());
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
    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq, qual, seq_handle.sqh);
        nrecs_read++;

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        close(out_handle);
        return 1;

      } // End try-catch for record reading.

      // If nskip > 0, just continue until nrecs_read > nskip, then write
      // output as file is read.
      //
      // Otherwise, keep pushing to the queue (after queue.size() == nlines,
      // also pop a record each time). Then, after the while loop, write all
      // the records in the queue.
      
      if(nskip > 0) {
        if(nrecs_read >= nskip) {
          try {
            writeRecord(out_handle, id, seq, qual);
          } catch (Exception const &e) {
            cerr << "Error writing output";
            seq_handle.close();
            return 1;
          }
        } else {
          continue;
        }
      } // End if(nskip > 0)
      else if(nlines > 0) {
        seqs.push(seq); ids.push(id); quals.push(qual);
        if(seqs.size() > (unsigned) nlines) {
          ids.pop(); seqs.pop(); quals.pop();
        }
      } // End if(nlines > 0)

    } // End single file reading loop
    
    // Write output if nlines > 0
    // Can we do for(StringChar id: ids; StringChar seq:seqs...)?
    if(nlines > 0) {
      for(int i=0; i < nlines; i++) {
        try {
          writeRecord(out_handle, ids.front(), seqs.front(),
                      quals.front());
          ids.pop(); seqs.pop(); quals.pop();
        } catch (Exception const &e) {
          cerr << "Error writing output";
          seq_handle.close();
          return 1;
        }
      }
    }

    if(!seq_handle.close()) {
      cerr << "Problem closing " << infile << endl;
      close(out_handle);
      return 1;
    }

  } // End loop over files
  close(out_handle);

  return 0;
}
