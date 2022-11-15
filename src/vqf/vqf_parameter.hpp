#pragma once

#include <string>

namespace filters::vqf {
    
    enum class FingerprintSize : size_t {
        _16bit = 16, _8bit = 8
    };

    /**
     * @brief Vector Quotient Filter parameter for the benchmark
     * note: the fingerprint size is currently set with the _k parameter in the Filter structure, not this parameter
     * 
     * @tparam fSize fingerprint size in bits
     */
    template <FingerprintSize fSize>
    struct VQFParameter {
        static constexpr size_t max_n_retries{16};
        static constexpr FingerprintSize fingerprintSize{fSize};

        static std::string to_string() {
            std::string s_fSize;
            switch (fSize) {
                case FingerprintSize::_16bit:
                    s_fSize = "16";
                    break;
                case FingerprintSize::_8bit:
                    s_fSize = "8";
                    break;
            }

            std::string s = "{";
            s += "\"fingerprint size\": \"" + s_fSize + "\", ";
            s += "\"max_n_retries\": " + std::to_string(max_n_retries) + "}";
            return s;
        }
    };

    template<size_t> using Standard = VQFParameter<FingerprintSize::_16bit>;
    template<size_t> using SmallFingerprints = VQFParameter<FingerprintSize::_8bit>;
}