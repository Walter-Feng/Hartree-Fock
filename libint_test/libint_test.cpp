#include <libint2.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace libint2;
using namespace std;

int main(int argc, char* argv[]) {
    libint2::initialize();  // safe to use libint now .. do `libint2::initialize(true)` to produce diagnostic messages

    // all other code snippets go here

    string xyzfilename = "/path/to/input_dot_xyz_file"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
    ifstream input_file(xyzfilename);
    vector<Atom> atoms = read_dotxyz(input_file);

    BasisSet obs("cc-pVDZ", atoms);

    // #include <algorithm> and <iterator>

    // print out the original basis
    std::copy(begin(obs), end(obs),
            std::ostream_iterator<Shell>(std::cout, "\n"));

    // stable sort preserves the original order where possible
    stable_sort(begin(obs), end(obs),
                [](const Shell& a, const Shell& b){
                return a.contr[0].l < b.contr[0].l;
                });

    // print out the resorted basis
    std::copy(begin(obs), end(obs),
            std::ostream_iterator<Shell>(std::cout, "\n"));


    libint2::finalize();  // do not use libint after this

    // can repeat the libint2::initialize() ... finalize() cycle as many times as
    // necessary

  return 0;
}