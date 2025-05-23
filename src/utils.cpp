#include "../header/utils.h"
#include <string>
#include <algorithm>
#include <vector>

namespace utils {
    namespace strings {
        std::string toLower(const std::string& str) {
            std::string newStr = str;
            std::transform(newStr.begin(), newStr.end(), newStr.begin(), [](unsigned char c){ return std::tolower(c); });

            return newStr;
        }

        std::string toUpper(const std::string& str) {
            std::string newStr = str;
            std::transform(newStr.begin(), newStr.end(), newStr.begin(), [](unsigned char c) { return std::toupper(c); });

            return newStr;
        }

        std::string replace(const std::string& str, const char& target, const char& replaceWith) {
            std::string newStr = str;
            for (char& c : newStr) {
                if (target == c) {
                    c = replaceWith;
                }
            }
            return newStr;
        }

        std::string convertBytesToStr(const std::vector<unsigned char>& data) {
            return std::string(data.begin(), data.end());
        }
    }

    namespace math {
        int add(int a, int b) {
            return a + b;
        }
    }
}