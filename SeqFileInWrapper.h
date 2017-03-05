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

using std::string;
using std::cin;
using std::ifstream;
using std::istream;
using namespace seqan;

namespace bltools {

    struct SeqFileInWrapper {

        private:
            ifstream input;
            istream * input_stream;

        public:
            SeqFileIn sqh;

            void open(char * infile);
            void open(string &infile); 
            bool close();
            bool atEnd();
    };
}
