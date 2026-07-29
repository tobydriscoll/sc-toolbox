#pragma once
// Reads the plain-text golden files produced by
// tests/cpp/exportGoldensText.m (see that script for the format).
// Deliberately avoids any MATLAB-file-format library: the format is a
// whitespace-delimited text dump designed to be trivial to parse with
// std::istringstream / strtod.

#include <cstdlib>
#include <complex>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace golden {

using Matrix = Eigen::MatrixXcd;

struct Case {
    std::string desc;
    double tol = 0.0;
    std::map<std::string, Matrix> inputs;
    std::map<std::string, Matrix> outputs;
    std::map<std::string, std::string> inputStrings;
    std::map<std::string, std::string> outputStrings;

    // Convenience accessors. Throw std::out_of_range if missing (a
    // missing field means the golden generator didn't record it).
    const Matrix& in(const std::string& name) const { return inputs.at(name); }
    const Matrix& out(const std::string& name) const { return outputs.at(name); }
    const std::string& inStr(const std::string& name) const { return inputStrings.at(name); }
    const std::string& outStr(const std::string& name) const { return outputStrings.at(name); }

    double inReal(const std::string& name) const { return in(name)(0, 0).real(); }
    int inInt(const std::string& name) const { return static_cast<int>(std::lround(inReal(name))); }
};

using GroupMap = std::map<std::string, std::vector<Case>>;

// std::istream::operator>>(double&) does not reliably parse "NaN"/"Inf"
// tokens across standard library implementations (libc++ fails outright),
// even though the golden text format relies on exactly those tokens for
// non-finite values. Tokenize manually and parse each with strtod, which
// is required by the C standard to accept "nan"/"inf" case-insensitively.
inline double parseToken(const std::string& tok) {
    return std::strtod(tok.c_str(), nullptr);
}

inline Matrix readMatrix(std::istream& is, int rows, int cols, bool isComplex) {
    Matrix m(rows, cols);
    for (int r = 0; r < rows; ++r) {
        std::string line;
        if (!std::getline(is, line)) throw std::runtime_error("golden_reader: unexpected EOF reading matrix row");
        std::istringstream ls(line);
        for (int c = 0; c < cols; ++c) {
            std::string retok, imtok;
            ls >> retok;
            double im = 0.0;
            if (isComplex) {
                ls >> imtok;
                im = parseToken(imtok);
            }
            m(r, c) = std::complex<double>(parseToken(retok), im);
        }
    }
    return m;
}

inline GroupMap loadGoldens(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("golden_reader: could not open " + path);

    GroupMap groups;
    std::string line;
    std::string currentGroup;
    std::vector<Case>* currentVec = nullptr;
    Case current;
    bool inCase = false;

    while (std::getline(f, line)) {
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "GROUP") {
            std::string name;
            ls >> name;
            currentGroup = name;
            currentVec = &groups[name];
        } else if (tag == "CASE") {
            current = Case{};
            inCase = true;
        } else if (tag == "DESC") {
            std::string rest;
            std::getline(ls, rest);
            if (!rest.empty() && rest[0] == ' ') rest.erase(0, 1);
            current.desc = rest;
        } else if (tag == "TOL") {
            ls >> current.tol;
        } else if (tag == "INPUT" || tag == "OUTPUT") {
            std::string fname, kindOrRows;
            ls >> fname;
            // Either "<rows> <cols> <R|C>" or "STR" or "SKIP"
            std::string next;
            ls >> next;
            if (next == "STR") {
                std::string strval;
                std::getline(f, strval);
                if (tag == "INPUT") current.inputStrings[fname] = strval;
                else current.outputStrings[fname] = strval;
            } else if (next == "SKIP") {
                // field intentionally omitted (non-numeric MATLAB type)
            } else {
                int rows = std::stoi(next);
                int cols;
                std::string kind;
                ls >> cols >> kind;
                bool isComplex = (kind == "C");
                Matrix m = readMatrix(f, rows, cols, isComplex);
                if (tag == "INPUT") current.inputs[fname] = m;
                else current.outputs[fname] = m;
            }
        } else if (tag == "ENDCASE") {
            if (currentVec) currentVec->push_back(current);
            inCase = false;
        } else if (tag == "ENDGROUP") {
            currentVec = nullptr;
        }
    }
    (void)inCase;
    return groups;
}

}  // namespace golden
