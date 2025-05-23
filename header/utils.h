#pragma once

#include <string>
#include <algorithm>
#include <vector>

namespace utils
{
    namespace strings {
        std::string toLower(const std::string& str);
        std::string toUpper(const std::string& str);

        std::string replace(const std::string& str, const char& target, const char& replaceWith);

        std::string convertBytesToStr(const std::vector<unsigned char>& data);
    }

    namespace math {
        
    }
    
} // namespace utils
