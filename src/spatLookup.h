#ifndef SPATLOOKUP_GUARD
#define SPATLOOKUP_GUARD

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring>

struct SpatHashDouble {
    std::size_t operator()(double d) const {
        if (std::isnan(d)) return 0;
        return std::hash<double>{}(d);
    }
};

struct SpatEqDouble {
    bool operator()(double a, double b) const {
        if (std::isnan(a)) return std::isnan(b);
        return a == b;
    }
};

struct SpatHashVector {
    std::size_t operator()(const std::vector<double>& v) const {
        std::size_t seed = 0;
        for (double d : v) {
            std::size_t h = std::isnan(d) ? 0 : std::hash<double>{}(d);
            // Boost hash_combine logic
            seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct SpatEqVector {
    bool operator()(const std::vector<double>& a, const std::vector<double>& b) const {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::isnan(a[i])) {
                if (!std::isnan(b[i])) return false;
            } else if (a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }
};

template <typename K, typename V>
using SpatHashMap = std::unordered_map<K, V,
    typename std::conditional<std::is_same<K, double>::value, SpatHashDouble, SpatHashVector>::type,
    typename std::conditional<std::is_same<K, double>::value, SpatEqDouble, SpatEqVector>::type>;

template <typename K>
using SpatHashSet = std::unordered_set<K,
    typename std::conditional<std::is_same<K, double>::value, SpatHashDouble, SpatHashVector>::type,
    typename std::conditional<std::is_same<K, double>::value, SpatEqDouble, SpatEqVector>::type>;

template <typename K>
using SpatFrequencyTable = SpatHashMap<K, size_t>;

#endif // SPATLOOKUP_GUARD
