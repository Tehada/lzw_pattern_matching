///
/// @file
/// @author Julius Pettersson
/// @copyright MIT/Expat License.
/// @brief LZW file compressor
/// @version 1
///
/// This is the C++11 implementation of a Lempel-Ziv-Welch single-file command-line compressor.
/// It uses the simpler fixed-width code compression method.
/// It was written with Doxygen comments.
///
/// @see http://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv%E2%80%93Welch
/// @see http://marknelson.us/2011/11/08/lzw-revisited/
/// @see http://www.cs.duke.edu/csed/curious/compression/lzw.html
/// @see http://warp.povusers.org/EfficientLZW/index.html
/// @see http://en.cppreference.com/
/// @see http://www.doxygen.org/
///

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

/// Type used to store and retrieve codes.
using CodeType = std::uint16_t;

namespace globals {

/// Dictionary Maximum Size (when reached, the dictionary will be reset)
const CodeType dms {std::numeric_limits<CodeType>::max()};

} // namespace globals

///
/// @brief Helper operator intended to simplify code.
/// @param vc   original vector
/// @param c    element to be appended
/// @returns vector resulting from appending `c` to `vc`
///
std::vector<char> operator + (std::vector<char> vc, char c)
{
    vc.push_back(c);
    return vc;
}

///
/// @brief Compresses the contents of `is` and writes the result to `os`.
/// @param [in] is      input stream
/// @param [out] os     output stream
///
void compress(std::istream &is, std::ostream &os)
{
    std::map<std::vector<char>, CodeType> dictionary;

    // "named" lambda function, used to reset the dictionary to its initial contents
    const auto reset_dictionary = [&dictionary] {
        dictionary.clear();

        const long int minc = std::numeric_limits<char>::min();
        const long int maxc = std::numeric_limits<char>::max();

        for (long int c = minc; c <= maxc; ++c)
        {
            // to prevent Undefined Behavior, resulting from reading and modifying
            // the dictionary object at the same time
            const CodeType dictionary_size = dictionary.size();

            dictionary[{static_cast<char> (c)}] = dictionary_size;
        }
    };

    reset_dictionary();

    std::vector<char> s; // String
    char c;

    while (is.get(c))
    {
        // dictionary's maximum size was reached
        if (dictionary.size() == globals::dms)
            reset_dictionary();

        s.push_back(c);

        if (dictionary.count(s) == 0)
        {
            // to prevent Undefined Behavior, resulting from reading and modifying
            // the dictionary object at the same time
            const CodeType dictionary_size = dictionary.size();

            dictionary[s] = dictionary_size;
            s.pop_back();
            os.write(reinterpret_cast<const char *> (&dictionary.at(s)), sizeof (CodeType));
            s = {c};
        }
    }

    if (!s.empty())
        os.write(reinterpret_cast<const char *> (&dictionary.at(s)), sizeof (CodeType));
}

///
/// @brief Decompresses the contents of `is` and writes the result to `os`.
/// @param [in] is      input stream
/// @param [out] os     output stream
///
void decompress(std::istream &is, std::ostream &os)
{
    std::vector<std::vector<char>> dictionary;

    // "named" lambda function, used to reset the dictionary to its initial contents
    const auto reset_dictionary = [&dictionary] {
        dictionary.clear();
        dictionary.reserve(globals::dms);

        const long int minc = std::numeric_limits<char>::min();
        const long int maxc = std::numeric_limits<char>::max();

        for (long int c = minc; c <= maxc; ++c)
            dictionary.push_back({static_cast<char> (c)});
    };

    reset_dictionary();

    std::vector<char> s; // String
    CodeType k; // Key

    while (is.read(reinterpret_cast<char *> (&k), sizeof (CodeType)))
    {
        // dictionary's maximum size was reached
        if (dictionary.size() == globals::dms)
            reset_dictionary();

        if (k > dictionary.size())
            throw std::runtime_error("invalid compressed code");

        if (k == dictionary.size())
            dictionary.push_back(s + s.front());
        else
        if (!s.empty())
            dictionary.push_back(s + dictionary.at(k).front());

        os.write(&dictionary.at(k).front(), dictionary.at(k).size());
        s = dictionary.at(k);
    }

    if (!is.eof() || is.gcount() != 0)
        throw std::runtime_error("corrupted compressed file");
}

///
/// @brief Prints usage information and a custom error message.
/// @param s    custom error message to be printed
/// @param su   Show Usage information
///
void print_usage(const std::string &s = "", bool su = true)
{
    if (!s.empty())
        std::cerr << "\nERROR: " << s << '\n';

    if (su)
    {
        std::cerr << "\nUsage:\n";
        std::cerr << "\tprogram -flag input_file output_file\n\n";
        std::cerr << "Where `flag' is either `c' for compressing, or `d' for decompressing, and\n";
        std::cerr << "`input_file' and `output_file' are distinct files.\n\n";
        std::cerr << "Examples:\n";
        std::cerr << "\tlzw_v1.exe -c license.txt license.lzw\n";
        std::cerr << "\tlzw_v1.exe -d license.lzw new_license.txt\n";
    }

    std::cerr << std::endl;
}

///
/// @brief Actual program entry point.
/// @param argc             number of command line arguments
/// @param [in] argv        array of command line arguments
/// @retval EXIT_FAILURE    for failed operation
/// @retval EXIT_SUCCESS    for successful operation
///
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        print_usage("Wrong number of arguments.");
        return EXIT_FAILURE;
    }

    enum class Mode {
        Compress,
        Decompress
    };

    Mode m;

    if (std::string(argv[1]) == "-c")
        m = Mode::Compress;
    else
    if (std::string(argv[1]) == "-d")
        m = Mode::Decompress;
    else
    {
        print_usage(std::string("flag `") + argv[1] + "' is not recognized.");
        return EXIT_FAILURE;
    }

    std::ifstream input_file(argv[2], std::ios_base::binary);

    if (!input_file.is_open())
    {
        print_usage(std::string("input_file `") + argv[2] + "' could not be opened.");
        return EXIT_FAILURE;
    }

    std::ofstream output_file(argv[3], std::ios_base::binary);

    if (!output_file.is_open())
    {
        print_usage(std::string("output_file `") + argv[3] + "' could not be opened.");
        return EXIT_FAILURE;
    }

    try
    {
        input_file.exceptions(std::ios_base::badbit);
        output_file.exceptions(std::ios_base::badbit | std::ios_base::failbit);

        if (m == Mode::Compress)
            compress(input_file, output_file);
        else
        if (m == Mode::Decompress)
            decompress(input_file, output_file);
    }
    catch (const std::ios_base::failure &f)
    {
        print_usage(std::string("File input/output failure: ") + f.what() + '.', false);
        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        print_usage(std::string("Caught exception: ") + e.what() + '.', false);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
