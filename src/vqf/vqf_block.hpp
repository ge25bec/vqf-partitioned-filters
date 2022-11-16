#pragma once

#include <bitset>

using std::uint64_t;
using std::int64_t;
using std::uint32_t;
using std::uint16_t;
using std::uint8_t;
using std::int8_t;
using RegisterSize = filters::parameter::RegisterSize;
using SIMD = filters::parameter::SIMD;
using MT = filters::parameter::MultiThreading;

namespace filters::vqf {

/**
 * @brief mask for metadata interpretation in select(); 
 * 
 * taken from:
 * P. Pandey, A. Conway, and R. Johnson. vqf. Mar. 2021. https://github.com/splatlab/vqf (visited on 11/07/2022).
 * see LICENSE in src/vqf
 */
static uint64_t pdepMask[64 + 128] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   1ULL << 0, 1ULL << 1, 1ULL << 2, 1ULL << 3, 1ULL << 4, 1ULL << 5, 1ULL << 6, 1ULL << 7, 1ULL << 8, 1ULL << 9,
   1ULL << 10, 1ULL << 11, 1ULL << 12, 1ULL << 13, 1ULL << 14, 1ULL << 15, 1ULL << 16, 1ULL << 17, 1ULL << 18, 1ULL << 19, 
   1ULL << 20, 1ULL << 21, 1ULL << 22, 1ULL << 23, 1ULL << 24, 1ULL << 25, 1ULL << 26, 1ULL << 27, 1ULL << 28, 1ULL << 29, 
   1ULL << 30, 1ULL << 31, 1ULL << 32, 1ULL << 33, 1ULL << 34, 1ULL << 35, 1ULL << 36, 1ULL << 37, 1ULL << 38, 1ULL << 39, 
   1ULL << 40, 1ULL << 41, 1ULL << 42, 1ULL << 43, 1ULL << 44, 1ULL << 45, 1ULL << 46, 1ULL << 47, 1ULL << 48, 1ULL << 49, 
   1ULL << 50, 1ULL << 51, 1ULL << 52, 1ULL << 53, 1ULL << 54, 1ULL << 55, 1ULL << 56, 1ULL << 57, 1ULL << 58, 1ULL << 59, 
   1ULL << 60, 1ULL << 61, 1ULL << 62, 1ULL << 63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/**
 * @brief Structure of a single block inside the Vector Quotient Filter; each block contains a mini-filter 
 * 
 * @tparam SIMD determines the level of SIMD optimisation
 * @tparam MT determines if multi-threading is enabled (deprecated)
 * @tparam uint8_t fingerprint size in bits
 * @tparam uint8_t capacity of fingerprints in a block
 * @tparam uint8_t number of buckets in a block
 */
template <SIMD, MT, uint8_t, uint8_t, uint8_t>
struct VQFBlock { };

/**
 * @brief Structure of a single block inside the Vector Quotient Filter, scalar version
 * 
 * @tparam multiThreading determines if multi-threading is enabled (deprecated)
 * @tparam fingerprintSize fingerprint size in bits
 * @tparam blockCapacity capacity of fingerprints in a block
 * @tparam numBuckets number of buckets in a block
 */
template <MT multiThreading, uint8_t fingerprintSize, uint8_t blockCapacity, uint8_t numBuckets>
struct VQFBlock<SIMD::Scalar, multiThreading, fingerprintSize, blockCapacity, numBuckets> {
  using T = typename std::conditional<fingerprintSize == 16, uint16_t, uint8_t>::type;

  template <RegisterSize, SIMD, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
  friend struct VectorQuotientFilterContainer;

  private:
    std::bitset<numBuckets + blockCapacity> bucketSizeVector{}; //lists the number of fingerprints in each bucket in unary;
    std::array<T, blockCapacity> fingerprints{};

    /**
     * @brief Construct a new VQFBlock object by initialising the metadata
     * 
     */
    VQFBlock() {
      if constexpr (fingerprintSize == 16) {
        bucketSizeVector = 68719476735; //68719476735 is 2^36 - 1; for 36 buckets,
      } else {
        bucketSizeVector = std::move(std::bitset<numBuckets + blockCapacity>{"11111111111111111111111111111111111111111111111111111111111111111111111111111111"}); //80 1's for 80 buckets
      }
    }

    /**
     * @brief Finds the i'th 1 in bucketSizeVector, corresponding to the i'th bucket; 
     * 
     * taken from:
     * P. Pandey, A. Conway, and R. Johnson. vqf. Mar. 2021. https://github.com/splatlab/vqf (visited on 11/07/2022).
     * see LICENSE in src/vqf
     * 
     * @param i the required bucket
     * @return uint8_t index of the 1 if i is in [0, numBuckets), else numBuckets + blockCapacity
     */
    forceinline unroll_loops
    uint8_t select(uint8_t i) const {
      
      if constexpr (fingerprintSize == 16) {
        uint64_t val = _pdep_u64(pdepMask[i + 64], bucketSizeVector.to_ulong());
        return _tzcnt_u64(val);
      } else {
        std::bitset<numBuckets + blockCapacity> lowerMask{0xffffffffffffffff};

        //evaluates the number of 1's in the lower 64 bit of the metadata
        uint64_t lower_word = (bucketSizeVector & lowerMask).to_ulong();
        uint64_t lower_pdep = _pdep_u64(pdepMask[i + 64], lower_word);
        if (lower_pdep != 0) {
          return _tzcnt_u64(lower_pdep);
        }

        //evaluates the number of 1's in the higher 64 bit of the metadata
        i = i - _popcnt64(lower_word);
        uint64_t higher_word = (bucketSizeVector >> 64).to_ulong();
        uint64_t val = _pdep_u64(pdepMask[i + 64], higher_word);
        return _tzcnt_u64(val) + 64;
      }
    }

    /**
     * @brief looks up a fingerprint in the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of the input key
     * @return int8_t index of fingerprint in the fingerprint array, or -1 if not found 
     */
    forceinline unroll_loops
    int8_t lookup(uint8_t bucketIndex, uint16_t fingerprint) {
      uint8_t start{}, end{};

      //Computes start of run of elements in bucket at bucketIndex
      if (bucketIndex != 0) {
        start = select(bucketIndex - 1);

        if (start == numBuckets + blockCapacity)
          return -1;

        start = start - bucketIndex + 1;
      }

      //Computes end of elements in bucket
      end = select(bucketIndex);
      
      if (end == numBuckets + blockCapacity)
        return -1;

      end -= bucketIndex;

      //fingerprint comparison loop
      while(start < end) {
        if ((start >= fingerprints.size()) || start == numBuckets + blockCapacity) {
          std::cout << "Error in lookup!\n";
          return -1;
        }

        if (fingerprints[start] == fingerprint)
          return start;

        start++;
      }

      return -1;
    }

    /**
     * @brief inserts a fingerprint into the block
     * 
     * @param bucketIndex bucket index of fingerprint
     * @param fingerprint fingerprint of input key
     * @return bool true if insert was successful, false otherwise 
     */
    forceinline unroll_loops
    bool insert(uint8_t bucketIndex, uint16_t fingerprint) {
      const uint8_t metadataIndex{select(bucketIndex)}; //position for new 0 to be inserted into bucketSizeVector
      const uint8_t fingerprintSlot{(uint8_t) (metadataIndex - bucketIndex)}; //slot for new fingerprint to be inserted in
      uint8_t maxLength{blockCapacity - 1}; //loop bound
      uint8_t temp{}; //intermediate result
      std::bitset<numBuckets + blockCapacity> metadataMaskLo{0xFFFFFFFFFFFFFFFF}; //mask for low bits of metadata manipulation
      std::bitset<numBuckets + blockCapacity> metadataMaskHi{0xFFFFFFFFFFFFFFFF}; //mask for high bits of metadata manipulation
    
      if (metadataIndex == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty metadata index!\n";
        return false;
      }

      temp = select(numBuckets - 1);

      //filter is full
      if (temp == numBuckets + blockCapacity - 1) {
        //std::cout << "Error! The filter is full!\n";
        return false;

        //faulty metadata
      } else if (temp == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty index!\n";
        return false;
      }

      //adjusting the metadata masks for 128 bit metadata
      if constexpr (fingerprintSize != 16) {
        metadataMaskLo <<= 64;
        metadataMaskLo |= 0xFFFFFFFFFFFFFFFF;
        metadataMaskHi <<= 64;
        metadataMaskHi |= 0xFFFFFFFFFFFFFFFF;
      }

      //inserting a 0 at the metadata bit at metadataIndex
      metadataMaskLo <<= metadataIndex;
      metadataMaskLo = bucketSizeVector & metadataMaskLo.flip();
      bucketSizeVector <<= 1;
      bucketSizeVector &= metadataMaskHi << (metadataIndex + 1);
      bucketSizeVector |= metadataMaskLo;

      //shifting of the fingerprint array, to make space for new fingerprint
      while (maxLength > fingerprintSlot) {
        fingerprints[maxLength] = fingerprints[maxLength - 1];
        maxLength--;
      }
    
      //fingerprints[blockCapacity - 1] = fingerprint; (in original paper)
      fingerprints[fingerprintSlot] = fingerprint;

      return true;
    }

    /**
     * @brief removes a fingerprint from the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of input key
     * @return bool true if remove was successful, false otherwise 
     */
    forceinline unroll_loops
    bool remove(uint8_t bucketIndex, uint16_t fingerprint) {
      int8_t fingerprintPos{lookup(bucketIndex, fingerprint)}; //position of fingerprint to be removed
      const uint8_t metadataIndex{(uint8_t) (bucketIndex + fingerprintPos)};  //index of metadata bit to be deleted
      std::bitset<numBuckets + blockCapacity> metadataMaskLo{0xFFFFFFFFFFFFFFFF}; //mask for low bits of metadata manipulation
      std::bitset<numBuckets + blockCapacity> metadataMaskHi{0xFFFFFFFFFFFFFFFF}; //mask for high bits of metadata manipulation

      if (fingerprintPos < 0)
        return false;

      //adjusting the metadata masks for 128 bit metadata
      if constexpr (fingerprintSize != 16) {
        metadataMaskLo <<= 64;
        metadataMaskLo |= 0xFFFFFFFFFFFFFFFF;
        metadataMaskHi <<= 64;
        metadataMaskHi |= 0xFFFFFFFFFFFFFFFF;
      }

      //removing a 0 from the metadata bit at metadataIndex
      metadataMaskLo <<= metadataIndex;
      metadataMaskLo = bucketSizeVector & metadataMaskLo.flip();
      bucketSizeVector >>= 1;
      bucketSizeVector &= metadataMaskHi << metadataIndex;
      bucketSizeVector |= metadataMaskLo;

      //shifting the fingerprint array to effectively overwrite the fingerprint
      while (fingerprintPos < blockCapacity - 1) {
        fingerprints[fingerprintPos] = fingerprints[fingerprintPos + 1];
        fingerprintPos++;
      }

      return true;
    }
};

/**
 * @brief Structure of a single block inside the Vector Quotient Filter, AVX2 version
 * 
 * @tparam multiThreading determines if multi-threading is enabled (deprecated)
 * @tparam fingerprintSize fingerprint size in bits
 * @tparam blockCapacity capacity of fingerprints in a block
 * @tparam numBuckets number of buckets in a block
 */
template <MT multiThreading, uint8_t fingerprintSize, uint8_t blockCapacity, uint8_t numBuckets>
struct VQFBlock<SIMD::AVX2, multiThreading, fingerprintSize, blockCapacity, numBuckets> {
  using T = typename std::conditional<fingerprintSize == 16, uint16_t, uint8_t>::type;

  template <RegisterSize, SIMD, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
  friend struct VectorQuotientFilterContainer;

  private:
    std::array<T, blockCapacity> fingerprints{};
    std::bitset<numBuckets + blockCapacity> bucketSizeVector{}; //lists the number of fingerprints in each bucket in unary;

    VQFBlock() {
      if constexpr (fingerprintSize == 16) {
        bucketSizeVector = 68719476735; //68719476735 is 2^36 - 1 for 36 buckets,
      } else {
        bucketSizeVector = std::move(std::bitset<numBuckets + blockCapacity>{"11111111111111111111111111111111111111111111111111111111111111111111111111111111"}); //80 1's for 80 buckets
      }
    }

    /**
     * @brief Finds the i'th 1 in bucketSizeVector, corresponding to the i'th bucket; 
     * 
     * taken from:
     * P. Pandey, A. Conway, and R. Johnson. vqf. Mar. 2021. https://github.com/splatlab/vqf (visited on 11/07/2022).
     * see LICENSE in src/vqf
     * 
     * @param i the required bucket
     * @return uint8_t index of the 1 if i is in [0, numBuckets), else numBuckets + blockCapacity
     */
    forceinline unroll_loops
    uint8_t select(uint8_t i) const {
      
      if constexpr (fingerprintSize == 16) {
        uint64_t val = _pdep_u64(pdepMask[i + 64], bucketSizeVector.to_ullong());
        return _tzcnt_u64(val);
      } else {
        std::bitset<numBuckets + blockCapacity> lowerMask{0xffffffffffffffff};

        //evaluates the number of 1's in the low 64 bit of bucketSizeVector
        uint64_t lower_word = (bucketSizeVector & lowerMask).to_ulong();
        uint64_t lower_pdep = _pdep_u64(pdepMask[i + 64], lower_word);
        if (lower_pdep != 0) {
          return _tzcnt_u64(lower_pdep);
        }

        //evaluates the number of 1's in the high 64 bit of bucketSizeVector
        i = i - _popcnt64(lower_word);
        uint64_t higher_word = (bucketSizeVector >> 64).to_ullong();
        uint64_t val = _pdep_u64(pdepMask[i + 64], higher_word);
        return _tzcnt_u64(val) + 64;
      }
    }

    /**
     * @brief looks up a fingerprint in the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingeprint of the input key
     * @return int8_t index of fingerprint in fingerprint array, or -1 if not found 
     */
    forceinline unroll_loops
    int8_t lookup(uint8_t bucketIndex, uint16_t fingerprint) const {
      uint8_t start{}, end{};

      //Compute start of run of elements in bucket at bucketIndex
      if (bucketIndex != 0) {
        start = select(bucketIndex - 1);

        if (start == numBuckets + blockCapacity)
          return -1;

        start = start - bucketIndex + 1;
      }

      //Compute end of bucket
      end = select(bucketIndex);
      
      if (end == numBuckets + blockCapacity)
        return -1;

      end -= bucketIndex;

      //fingerprint comparisons
      while(start < end) {
        if ((start >= blockCapacity) || start == numBuckets + blockCapacity) {
          std::cout << "Error in lookup!\n";
          return -1;
        }

        if (fingerprints[start] == fingerprint)
          return start;

        start++;
      }

      return -1;
    }

    /**
     * @brief inserts a fingerprint into the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of the input key
     * @return bool true if insert was successful, false if not
     */
    forceinline unroll_loops
    bool insert(uint8_t bucketIndex, uint16_t fingerprint) {
      const uint8_t metadataIndex{select(bucketIndex)}; //position in bucketSizeVector where the new 0 will be inserted
      const uint8_t fingerprintSlot{(uint8_t) (metadataIndex - bucketIndex)}; //position in fingerprint array where the new fingerprint will be inserted
      uint8_t maxLength{blockCapacity - 1}; //loop bound
      uint8_t temp{}; //intermediate result
      std::bitset<numBuckets + blockCapacity> metadataMaskLo{0xFFFFFFFFFFFFFFFF}; //metadata manipulation mask for the low bit of bucketSizeVector
      std::bitset<numBuckets + blockCapacity> metadataMaskHi{0xFFFFFFFFFFFFFFFF}; //metadata manipulation mask for the high bit of bucketSizeVector
      __m256i _data{};
    
      if (metadataIndex == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty metadata index!\n";
        return false;
      }

      temp = select(numBuckets - 1);

      //filter is full
      if (temp == blockCapacity + numBuckets - 1) {
        //std::cout << "Error! The filter is full!\n";
        return false;

        //faulty metadata
      } else if (temp == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty index!\n";
        return false;
      }

      //adjusting the metadata masks for 128 bit metadata
      if constexpr (fingerprintSize != 16) {
        metadataMaskLo <<= 64;
        metadataMaskLo |= 0xFFFFFFFFFFFFFFFF;
        metadataMaskHi <<= 64;
        metadataMaskHi |= 0xFFFFFFFFFFFFFFFF;
      }

      //inserting a 0 at the metadata bit at metadataIndex
      metadataMaskLo <<= metadataIndex;
      metadataMaskLo = bucketSizeVector & metadataMaskLo.flip();
      bucketSizeVector <<= 1;
      bucketSizeVector &= metadataMaskHi << (metadataIndex + 1);
      bucketSizeVector |= metadataMaskLo;

      //shifting the fingerprint array by one to make space for the new fingerprint
      if constexpr (fingerprintSize == 16) {
        if (fingerprintSlot < 12) {
          _data = _mm256_loadu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + maxLength - 17));
      
          _mm256_storeu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + maxLength - 16), _data);

          maxLength -= 16;
        }
      } else {
        if (fingerprintSlot < 16) {
          _data = _mm256_loadu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + maxLength - 33));
      
          _mm256_storeu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + maxLength - 32), _data);

          maxLength -= 32;
        }
      }

      //shifting the rest of the fingerprint array scalarly
      while (maxLength > fingerprintSlot) {
        fingerprints[maxLength] = fingerprints[maxLength - 1];
        maxLength--;
      }
    
      //fingerprints[blockCapacity - 1] = fingerprint; (in orginal paper)
      fingerprints[fingerprintSlot] = fingerprint;

      return true;
    }

    /**
     * @brief removes a fingerprint from the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of the input key
     * @return bool true if remove was successful, false otherwise 
     */
    forceinline unroll_loops
    bool remove(uint8_t bucketIndex, uint16_t fingerprint) {
      int8_t fingerprintPos{lookup(bucketIndex, fingerprint)}; //position of fingerprint in the fingerprint array, that will be removed
      uint8_t metadataIndex{};  //index of metadata bit to be deleted
      std::bitset<numBuckets + blockCapacity> metadataMaskLo{0xFFFFFFFFFFFFFFFF}; //mask for low bit of the metedata manipulation
      std::bitset<numBuckets + blockCapacity> metadataMaskHi{0xFFFFFFFFFFFFFFFF}; //mask for high bit of the metadata manipulation
      __m256i _data{};

      if (fingerprintPos < 0)
        return false;

      metadataIndex = bucketIndex + (uint8_t) fingerprintPos;

      //adjusting the metadata masks for 128 bit metadata
      if constexpr (fingerprintSize != 16) {
        metadataMaskLo <<= 64;
        metadataMaskLo |= 0xFFFFFFFFFFFFFFFF;
        metadataMaskHi <<= 64;
        metadataMaskHi |= 0xFFFFFFFFFFFFFFFF;
      }

      //removing a 0 at the metadata bit at metadataIndex
      metadataMaskLo <<= metadataIndex;
      metadataMaskLo = bucketSizeVector & metadataMaskLo.flip();
      bucketSizeVector >>= 1;
      bucketSizeVector &= metadataMaskHi << metadataIndex;
      bucketSizeVector |= metadataMaskLo;

      //shifting the fingerprint array one to the right, to overwrite fingerprint
      if constexpr (fingerprintSize == 16) {
        if (fingerprintPos < 12) {
          _data = _mm256_loadu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + fingerprintPos + 1));
      
          _mm256_storeu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + fingerprintPos), _data);

          fingerprintPos += 16;
        }
      } else {
        if (fingerprintPos < 16) {
          _data = _mm256_loadu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + fingerprintPos + 1));
      
          _mm256_storeu_si256(reinterpret_cast<__m256i*>(fingerprints.data() + fingerprintPos), _data);

          fingerprintPos += 32;
        }
      }

      //shifting the rest of the fingerprint array scalarly
      while (fingerprintPos < blockCapacity - 1) {
        fingerprints[fingerprintPos] = fingerprints[fingerprintPos + 1];
        fingerprintPos++;
      }

      return true;
    }
};

/**
 * @brief Structure of a single block inside the Vector Quotient Filter, AVX512 version
 * 
 * @tparam multiThreading determines if multi-threading is enabled (deprecated)
 * @tparam fingerprintSize fingerprint size in bits
 * @tparam blockCapacity capacity of fingerprints in a block
 * @tparam numBuckets number of buckets in a block
 */
template <MT multiThreading, uint8_t fingerprintSize, uint8_t blockCapacity, uint8_t numBuckets>
struct alignas(64) VQFBlock<SIMD::AVX512, multiThreading, fingerprintSize, blockCapacity, numBuckets> {
  using T = typename std::conditional<fingerprintSize == 16, uint16_t, uint8_t>::type;

  template <RegisterSize, SIMD, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
  friend struct VectorQuotientFilterContainer;

  private:
    std::bitset<numBuckets + blockCapacity> bucketSizeVector{}; //lists the number of fingerprints in each bucket in unary;
    std::array<T, blockCapacity> fingerprints{};

    VQFBlock() {
      if constexpr (fingerprintSize == 16) {
        bucketSizeVector = 68719476735; //68719476735 is 2^36 - 1 for 36 buckets,
      } else {
        bucketSizeVector = std::move(std::bitset<numBuckets + blockCapacity>{"11111111111111111111111111111111111111111111111111111111111111111111111111111111"}); //80 1's for 80 buckets
      }
    }

    /**
     * @brief Finds the i'th 1 in bucketSizeVector, corresponding to the i'th bucket; 
     * 
     * taken from:
     * P. Pandey, A. Conway, and R. Johnson. vqf. Mar. 2021. https://github.com/splatlab/vqf (visited on 11/07/2022).
     * see LICENSE in src/vqf
     * 
     * @param i the required bucket
     * @return uint8_t index of the 1 if i is in [0, numBuckets), else numBuckets + blockCapacity
     */
    forceinline unroll_loops
    uint8_t select(uint8_t i) const {
      
      if constexpr (fingerprintSize == 16) {
        uint64_t val = _pdep_u64(pdepMask[i + 64], bucketSizeVector.to_ulong());
        return _tzcnt_u64(val);
      } else {
        std::bitset<numBuckets + blockCapacity> lowerMask{0xffffffffffffffff};

        //evaluates the number of 1's in the lower 64 bit of bucketSizeVector
        uint64_t lower_word = (bucketSizeVector & lowerMask).to_ulong();
        uint64_t lower_pdep = _pdep_u64(pdepMask[i + 64], lower_word);
        if (lower_pdep != 0) {
          return _tzcnt_u64(lower_pdep);
        }

        //evaluates the number of 1's in the high 64 bit of bucketSizeVector
        i = i - _popcnt64(lower_word);
        uint64_t higher_word = (bucketSizeVector >> 64).to_ulong();
        uint64_t val = _pdep_u64(pdepMask[i + 64], higher_word);
        return _tzcnt_u64(val) + 64;
      }
    }

    /**
     * @brief looks up a fingerprint in the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of the input key
     * @return int8_t index of the fingerprint, or -1 if not found 
     */
    forceinline
    int8_t lookup(uint8_t bucketIndex, uint16_t fingerprint) const {
      uint8_t start{}, end{};
      uint64_t lookupResult{}; 
      __mmask32 loadMask{}; //mask to load the fingerprints in bucket
      __mmask64 loadMask8bit{}; //mask to load the fingerprints in bucket for 8 bit fingerprints
      const __m512i _cmp{_mm512_set1_epi16(fingerprint)}; //fingerprint to compare other fingerprints to
      const __m512i _cmp8bit{_mm512_set1_epi8(fingerprint)}; //fingerprint to compare other fingerprints to, for 8 bit 
      __m512i _data{};

      //Compute start of run of elements in bucket at bucketIndex
      if (bucketIndex != 0) {
        start = select(bucketIndex - 1);

        if (start == numBuckets + blockCapacity)
          return -1;

        start = start - bucketIndex + 1;
      }

      //compute end of elements in bucket
      end = select(bucketIndex);
      
      if (end == numBuckets + blockCapacity)
        return -1;

      end -= bucketIndex;

      //comparing the fingerprints with avx512
      if constexpr (fingerprintSize == 16) {
        //shifting the mask so that only fingerprints from start to end are considered   
        loadMask = _cvtu32_mask32((268435455 >> (28 - end + start)) << (4 + start));

        //fingerprint comparison
        _data = _mm512_load_si512(reinterpret_cast<const __m512i *>(this));
        loadMask = _mm512_mask_cmpeq_epu16_mask(loadMask, _data, _cmp);
        lookupResult = _cvtmask32_u32(loadMask);
      } else {
        //shifting the mask so that only fingerprints from start to end are considered
        loadMask8bit= _cvtu64_mask64((281474976710655 >> (48 - end + start)) << (16 + start));
        
        //fingerprint comparisons
        _data = _mm512_load_si512(reinterpret_cast<const __m512i *>(this)); 
        loadMask8bit = _mm512_mask_cmpeq_epu8_mask(loadMask8bit, _data, _cmp8bit);
        lookupResult = _cvtmask64_u64(loadMask8bit);
      }

      //evaluation of comparison result
      if (lookupResult > 0) {
        return _bit_scan_forward(lookupResult) + start;
      }

      return -1;
    }

    /**
     * @brief inserts a fingerprint into the block
     * 
     * @param bucketIndex bucket index of the fingerprint
     * @param fingerprint fingerprint of the input key
     * @return bool true if insert was successful, otherwise false 
     */
    forceinline
    bool insert(uint8_t bucketIndex, uint16_t fingerprint) {
      const uint8_t metadataIndex{select(bucketIndex)}; //index in metadata where new 0 will be inserted
      const uint8_t fingerprintSlot{(uint8_t) (metadataIndex - bucketIndex)}; //index in fingerprint array where new fingerprint will be inserted
      uint8_t temp{}; //intermediate result
      std::bitset<numBuckets + blockCapacity> metadataMaskLo{0xFFFFFFFFFFFFFFFF}; //metadata mask for low 64 bit of metadata manipulation
      std::bitset<numBuckets + blockCapacity> metadataMaskHi{0xFFFFFFFFFFFFFFFF}; //metadata mask for high 64 bit of metadata manipulation
      __mmask32 loadMask{}; //mask to load the fingerprints in bucket
      __mmask64 loadMask8bit{}; //mask to load the fingerprints in bucket for 8 bit fingerprints
      __m512i _data{};
    
      if (metadataIndex == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty metadata index!\n";
        return false;
      }

      temp = select(numBuckets - 1);

      //filter is full
      if (temp == blockCapacity + numBuckets - 1) {
        //std::cout << "Error! The filter is full!\n";
        return false;

        //faulty metadata
      } else if (temp == numBuckets + blockCapacity) {
        std::cout << "Error in insert! Faulty index!\n";
        return false;
      }

      //adjusting the metadata masks for 128 bit metadata
      if constexpr (fingerprintSize != 16) {
        metadataMaskLo <<= 64;
        metadataMaskLo |= 0xFFFFFFFFFFFFFFFF;
        metadataMaskHi <<= 64;
        metadataMaskHi |= 0xFFFFFFFFFFFFFFFF;
      }

      //inserting a 0 at the metadata bit at metadataIndex
      metadataMaskLo <<= metadataIndex;
      metadataMaskLo = bucketSizeVector & metadataMaskLo.flip();
      bucketSizeVector <<= 1;
      bucketSizeVector &= metadataMaskHi << (metadataIndex + 1);
      bucketSizeVector |= metadataMaskLo;

      //shifting the fingerprint array to the right to make space for new fingerprint
      if constexpr (fingerprintSize == 16) {
        /** using permb (only draft, not tested because flag is not supported)
        #ifdef __AVX512VBMI__
        __m512i _idx{_mm512_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63)};
        __m512i _idxShifted{_mm512_set_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 0)};
        __mmask64 idxMask{_cvtu64_mask64(9223372036854775807UL << (8 + (fingerprintSlot << 1)))};

        __m512i _idxRes{_mm512_mask_blend_epi8(idxMask, _idx, _idxShifted)};

        _data = _mm512_load_si512(reinterpret_cast<const __m512i *>(this));
        _data = _mm512_permutexvar_epi8(_idxRes, _data);
        _data = _mm512_store_si512(reinterpret_cast<const __m512i *>(this));
        #endif
        **/

        loadMask = _cvtu32_mask32(134217727 >> fingerprintSlot);
        _data = _mm512_maskz_loadu_epi16(loadMask, fingerprints.data() + fingerprintSlot);
        _mm512_mask_storeu_epi16(fingerprints.data() + fingerprintSlot + 1, loadMask, _data);
      } else {
        loadMask8bit = _cvtu64_mask64(140737488355327 >> fingerprintSlot);
        _data = _mm512_maskz_loadu_epi8(loadMask8bit, fingerprints.data() + fingerprintSlot);
        _mm512_mask_storeu_epi8(fingerprints.data() + fingerprintSlot + 1, loadMask8bit, _data);
      }
    
      //fingerprints[blockCapacity - 1] = fingerprint; (in original paper)
      fingerprints[fingerprintSlot] = fingerprint;

      return true;
    }
};
}