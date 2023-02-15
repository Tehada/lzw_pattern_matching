#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

using CodeType = std::uint16_t;

namespace globals {

const CodeType dms{std::numeric_limits<CodeType>::max()};

}  // namespace globals

std::vector<char> operator+(std::vector<char> vc, char c) {
  vc.push_back(c);
  return vc;
}

void compress(std::istream &is, std::ostream &os) {
  std::map<std::vector<char>, CodeType> dictionary;

  // "named" lambda function, used to reset the dictionary to its initial
  // contents
  const auto reset_dictionary = [&dictionary] {
    dictionary.clear();

    const long int minc = std::numeric_limits<char>::min();
    const long int maxc = std::numeric_limits<char>::max();

    for (long int c = minc; c <= maxc; ++c) {
      // to prevent Undefined Behavior, resulting from reading and modifying
      // the dictionary object at the same time
      const CodeType dictionary_size = dictionary.size();

      dictionary[{static_cast<char>(c)}] = dictionary_size;
    }
  };

  reset_dictionary();

  std::vector<char> s;  // String
  char c;

  while (is.get(c)) {
    // dictionary's maximum size was reached
    if (dictionary.size() == globals::dms) reset_dictionary();

    s.push_back(c);

    if (dictionary.count(s) == 0) {
      // to prevent Undefined Behavior, resulting from reading and modifying
      // the dictionary object at the same time
      const CodeType dictionary_size = dictionary.size();

      dictionary[s] = dictionary_size;
      s.pop_back();
      os.write(reinterpret_cast<const char *>(&dictionary.at(s)),
               sizeof(CodeType));
      s = {c};
    }
  }

  if (!s.empty())
    os.write(reinterpret_cast<const char *>(&dictionary.at(s)),
             sizeof(CodeType));
}

void decompress(std::istream &is, std::ostream &os) {
  std::vector<std::vector<char>> dictionary;

  // "named" lambda function, used to reset the dictionary to its initial
  // contents
  const auto reset_dictionary = [&dictionary] {
    dictionary.clear();
    dictionary.reserve(globals::dms);

    const long int minc = std::numeric_limits<char>::min();
    const long int maxc = std::numeric_limits<char>::max();

    for (long int c = minc; c <= maxc; ++c)
      dictionary.push_back({static_cast<char>(c)});
  };

  reset_dictionary();

  std::vector<char> s;  // String
  CodeType k;           // Key

  while (is.read(reinterpret_cast<char *>(&k), sizeof(CodeType))) {
    // dictionary's maximum size was reached
    if (dictionary.size() == globals::dms) reset_dictionary();

    if (k > dictionary.size())
      throw std::runtime_error("invalid compressed code");

    if (k == dictionary.size())
      dictionary.push_back(s + s.front());
    else if (!s.empty())
      dictionary.push_back(s + dictionary.at(k).front());

    os.write(&dictionary.at(k).front(), dictionary.at(k).size());
    s = dictionary.at(k);
  }

  if (!is.eof() || is.gcount() != 0)
    throw std::runtime_error("corrupted compressed file");
}

void print_usage(const std::string &s = "", bool su = true) {
  if (!s.empty()) std::cerr << "\nERROR: " << s << '\n';

  if (su) {
    std::cerr << "\nUsage:\n";
    std::cerr << "\tprogram -flag input_file output_file\n\n";
    std::cerr << "Where `flag' is either `c' for compressing, or `d' for "
                 "decompressing, and\n";
    std::cerr << "`input_file' and `output_file' are distinct files.\n\n";
    std::cerr << "Examples:\n";
    std::cerr << "\tlzw_v1.exe -c license.txt license.lzw\n";
    std::cerr << "\tlzw_v1.exe -d license.lzw new_license.txt\n";
  }

  std::cerr << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    print_usage("Wrong number of arguments.");
    return EXIT_FAILURE;
  }

  enum class Mode { Compress, Decompress };

  Mode m;

  if (std::string(argv[1]) == "-c")
    m = Mode::Compress;
  else if (std::string(argv[1]) == "-d")
    m = Mode::Decompress;
  else {
    print_usage(std::string("flag `") + argv[1] + "' is not recognized.");
    return EXIT_FAILURE;
  }

  std::ifstream input_file(argv[2], std::ios_base::binary);

  if (!input_file.is_open()) {
    print_usage(std::string("input_file `") + argv[2] +
                "' could not be opened.");
    return EXIT_FAILURE;
  }

  std::ofstream output_file(argv[3], std::ios_base::binary);

  if (!output_file.is_open()) {
    print_usage(std::string("output_file `") + argv[3] +
                "' could not be opened.");
    return EXIT_FAILURE;
  }

  try {
    input_file.exceptions(std::ios_base::badbit);
    output_file.exceptions(std::ios_base::badbit | std::ios_base::failbit);

    if (m == Mode::Compress)
      compress(input_file, output_file);
    else if (m == Mode::Decompress)
      decompress(input_file, output_file);
  } catch (const std::ios_base::failure &f) {
    print_usage(std::string("File input/output failure: ") + f.what() + '.',
                false);
    return EXIT_FAILURE;
  } catch (const std::exception &e) {
    print_usage(std::string("Caught exception: ") + e.what() + '.', false);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
