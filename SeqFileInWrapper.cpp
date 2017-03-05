/*
 * Class to wrap SeqFileIn
 *
 * Might be better to just inherit from SeqFileIn, add input
 * and input_stream as private variables, and make a new
 * definition of open and close, but that's a more involved
 * project, and the method I use below will work for my own
 * programs.
 *
 * FormattedFile is defined in seqan/stream/formatted_file.h
 *
 * Note that Seqan is supposed to be able to handle gzipped files,
 * but only is SEQAN_HAS_ZLIB is defined as 1 and it is compiled
 * with zlib support. This seems to be a somewhat sketchy feature.
 *
 */

#include <string>
#include <iostream>
#include <seqan/seq_io.h>
#include <SeqFileInWrapper.h>

using std::string;
using std::cin;
using std::ifstream;
using std::istream;
using namespace seqan;

namespace bltools {

    void SeqFileInWrapper::open(char * infile) {
        string inf(infile);
        open(inf);
    }

    void SeqFileInWrapper::open(string &infile) {

        bool file_ok = false;

        if(infile == "-") {
            input_stream = &cin;
            file_ok = true;
        } else {
            input.open(infile.c_str(), ifstream::in);
            file_ok = input.is_open() && input.good();
            input_stream = &input;
        }
        file_ok &= seqan::open(sqh, *input_stream);
        if(!file_ok) {
            throw "problem opening file";
        }
    }

    bool SeqFileInWrapper::close() {
        bool close_ok = seqan::close(sqh);
        input.close();
        return close_ok;
    }

    bool SeqFileInWrapper::atEnd() {
        return seqan::atEnd(sqh);
    }
}
