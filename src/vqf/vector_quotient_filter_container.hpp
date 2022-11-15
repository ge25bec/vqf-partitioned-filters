#pragma once

#include <iostream>
#include <cstdint>
#include <vector>
#include <deque>
#include <immintrin.h>
#include <shared_mutex>
#include "vqf_block.hpp"

using std::uint64_t;
using std::uint32_t;
using std::uint16_t;
using std::uint8_t;

using RegisterSize = filters::parameter::RegisterSize;
using SIMD = filters::parameter::SIMD;

namespace filters::vqf {

/**
 * @brief struct for the entire Vector Quotient Filter container
 * 
 * @tparam RegisterSize determines whether the filter takes 64 bit or 32 bit input keys
 * @tparam SIMD determines the level of SIMD the filter uses
 * @tparam OptimizationParameter determines if the filter uses multi-threading
 * @tparam fSize sets the fingerprint size
 * @tparam Hasher determines the hashing procedure the filter uses for input keys
 * @tparam Vector enclosing type for SIMD registers with Hasher/Addresser
 * @tparam Addresser determines the addressing procedure the filter uses for first block index of an input key
 */
template <RegisterSize, SIMD, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
struct VectorQuotientFilterContainer { };

/**
 * @brief struct for the entire Vector Quotient Filter container without SIMD 
 * 
 * @tparam paramSize determines whether the filter takes 64 bit or 32 bit input keys
 * @tparam OptimizationParameter determines if the filter uses multi-threading
 * @tparam fSize sets the fingerprint size
 * @tparam Hasher determines the hashing procedure the filter uses for input keys
 * @tparam Vector enclosing type for SIMD registers with Hasher/Addresser
 * @tparam Addresser determines the addressing procedure the filter uses for first block index of an input key
 */
template <RegisterSize paramSize, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
struct VectorQuotientFilterContainer<paramSize, SIMD::Scalar, OptimizationParameter, fSize, Hasher, Vector, Addresser> {
  using T = typename Vector::T;
  using OP = OptimizationParameter;

  static constexpr uint8_t fingerprintSize{fSize}; //size of fingerprints in bit
  static constexpr uint8_t blockCapacity{(fSize == 16) ? 28 : 48}; //capacity of fingerprints in a block
  static constexpr uint8_t numBuckets{(fSize == 16) ? 36 : 80}; //number of buckets in a block
  static constexpr uint8_t shortcutOpt{(uint8_t) (0.75 * blockCapacity)}; //shortcut optimisation for insert: threshold for checking the occupancy of second block

  static_assert(sizeof(VQFBlock<SIMD::Scalar, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>) == 64); //static_assert to check that blocks always fit into a cache line
  std::vector<VQFBlock<SIMD::Scalar, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>> blocks{}; //vector of blocks for the filter

  std::deque<std::shared_mutex> mutexVector{}; //vector of shared mutexes, each one corresponding to a block
  Addresser addresser; //block addresser
  T numBlocks{2}; //number of blocks in the filter

  /**
   * @brief Construct a Vector Quotient Filter Container object with only two blocks
   * 
   */
  VectorQuotientFilterContainer() {
    //computes the block range in the filter
    addresser = std::move(Addresser(&numBlocks, 1));
    numBlocks = addresser.get_size(0);

    blocks.resize(numBlocks, VQFBlock<SIMD::Scalar, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
    mutexVector.resize(numBlocks);
  }

  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param histogram contains number of elements to store in the filter at index 0
   */
  VectorQuotientFilterContainer(size_t s, const T *histogram) {
    const T h{(T) ceil((size_t) histogram[0] * s / 100 / blockCapacity)}; //number of required blocks to store elements
  
    //computes the block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::Scalar, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
    mutexVector.resize(numBlocks);
  }
  
  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param length number of elements to store in the filter
   */
  VectorQuotientFilterContainer(size_t s, size_t length) {
    const T h{(T) ceil(length * s / 100 / blockCapacity)}; //number of required blocks to store elements

    //computes the block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::Scalar, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
    mutexVector.resize(numBlocks);
  }

  /**
   * @brief looks up an input key in the filter
   * 
   * @param x input key
   * @return bool true if lookup was positive, false otherwise 
   */
  forceinline
  bool lookup(uint64_t x) {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint64_t queryHash{(uint64_t) v.vector}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (queryHash & 65535) : (queryHash & 255))};
    uint64_t blockIndex1{(uint64_t) w.vector};
    uint64_t blockIndex2{blockIndex1 ^ queryHash};
    uint64_t bucketIndex{queryHash >> fingerprintSize};

    //bringing indexes into range
    blockIndex2 = ((__uint128_t) (T) blockIndex2 * numBlocks) >> (uint8_t) paramSize;  
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in lookup! Block Index invalid!\n";
      return false;
    }

    if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled) {
      //shared locking of the first block
      mutexVector.at(blockIndex1).lock_shared();
      bool result{blocks.at(blockIndex1).lookup(bucketIndex, fingerprint) != -1};
      mutexVector.at(blockIndex1).unlock_shared();

      //shared locking of the second block, if first lookup unsuccessful
      if (!result) {
        mutexVector.at(blockIndex2).lock_shared();
        result = result || (blocks.at(blockIndex2).lookup(bucketIndex, fingerprint) != -1);
        mutexVector.at(blockIndex2).unlock_shared();
      }

      return result;
    } else {
      //looking up fingerprint in both blocks
      return ((blocks.at(blockIndex1).lookup(bucketIndex, fingerprint) != -1) || (blocks.at(blockIndex2).lookup(bucketIndex, fingerprint) != -1));
    }
  }

  /**
   * @brief looks up an array of input keys in the filter
   * 
   * @param values array of input keys
   * @param length length of the array
   * @return size_t number of positive lookups 
   */
  forceinline unroll_loops
  size_t lookup(const T *values, size_t length) {
    size_t result{};

    for (size_t i{}; i < length; i++) {
      if (lookup(values[i]))
        result++;
    }

    return result;
  }

  /**
   * @brief inserts an input key into the filter
   * 
   * @param x input key
   * @return bool true if insert was successful, false otherwise 
   */
  forceinline
  bool insert(uint64_t x) {
    bool result{};
    uint64_t i{}; //index of the emptier block
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint64_t insertHash{(uint64_t) v.vector}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (insertHash & 65535) : (insertHash & 255))};
    uint64_t blockIndex1{(uint64_t) w.vector};
    uint64_t blockIndex2{insertHash ^ blockIndex1};
    uint64_t bucketIndex{insertHash >> fingerprintSize};

    //bringing indexes into range
    blockIndex2 = ((__uint128_t) (T) blockIndex2 * numBlocks) >> (uint8_t) paramSize;  
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;
    
    i = blockIndex1;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in Insert! Block Index invalid! \n";
      return false;
    }

    //locking the first block
    if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
      mutexVector.at(i).lock();

    //shortcut optimisation: block 2 is only checked if block1's occupancy is more than 75%
    if (blocks.at(i).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
      if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled) {
        //locking the second block, following lock order of increasing block index
        if (blockIndex2 < i) {
          mutexVector.at(i).unlock();
          mutexVector.at(blockIndex2).lock();
          mutexVector.at(i).lock();
        } else if (blockIndex2 > i) {
          mutexVector.at(blockIndex2).lock();
        }
      }
      
      //selecting index i of emptier block
      if (blocks.at(blockIndex2).select(numBuckets - 1) < blocks.at(blockIndex1).select(numBuckets - 1)) {
        i = blockIndex2;

        //unlocking irrelevant block
        if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled) {
          if (i != blockIndex1)
            mutexVector.at(blockIndex1).unlock();
        }
      } else {
        //unlocking irrelevant block
        if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled) {
          if (i != blockIndex2)
            mutexVector.at(blockIndex2).unlock();
        }
      }
    }

    //inserting fingerprint into emptier block
    result = blocks.at(i).insert(bucketIndex, fingerprint);
    
    //unlocking other block
    if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
      mutexVector.at(i).unlock();

    return result;
  }

  /**
   * @brief inserts an array of input keys into the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of successful inserts 
   */
  forceinline unroll_loops
  size_t insert(const T *values, size_t length) {
    size_t result{};

    for (size_t i{}; i < length; i++) {
      if (insert(values[i]))
        result++;
    }

    return result;
  }

  /**
   * @brief removes an input key from the filter
   * 
   * @param x input key
   * @return bool true if remove was successful, false otherwise 
   */
  forceinline
  bool remove(uint64_t x) {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint64_t removeHash{(uint64_t) v.vector}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (removeHash & 65535) : (removeHash & 255))};
    uint64_t blockIndex1{(uint64_t) w.vector};
    uint64_t blockIndex2{blockIndex1 ^ removeHash};
    uint64_t bucketIndex{removeHash >> fingerprintSize};

    //bringing indexes into range
    blockIndex2 = ((__uint128_t) (T) blockIndex2 * numBlocks) >> (uint8_t) paramSize;  
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in Remove! Block Index invalid!\n";
      return false;
    }

    //locking the first block
    if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
      mutexVector.at(blockIndex1).lock();

    //trying to remove fingerprint from first block
    if (blocks.at(blockIndex1).remove(bucketIndex, fingerprint)) {

      //unlocking first block
      if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
        mutexVector.at(blockIndex1).unlock();
      
      return true;
    } else {
      //locking second block
      if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
        mutexVector.at(blockIndex2).lock();

      //trying to remove fingerprint from second block, if first remove failed
      bool result{blocks.at(blockIndex2).remove(bucketIndex, fingerprint)};

      //unlocking second block
      if constexpr (OP::multiThreading == filters::parameter::MultiThreading::Enabled)
        mutexVector.at(blockIndex2).unlock();

      return result;
    }
  }

  /**
   * @brief returns total amount of space allocated for fingerprint arrays
   * 
   * @return size_t space allocated for the fingerprints in bytes
   */
  size_t size() const {
    return (size_t) numBlocks * blockCapacity * fingerprintSize/8;
  }
};

/**
 * @brief struct for the entire Vector Quotient Filter container using AVX2 (supports only 32 bit keys)
 * 
 * @tparam OptimizationParameter determines if the filter uses multi-threading
 * @tparam fSize sets the fingerprint size
 * @tparam Hasher determines the hashing procedure the filter uses for input keys
 * @tparam Vector enclosing type for SIMD registers with Hasher/Addresser
 * @tparam Addresser determines the addressing procedure the filter uses for first block index of an input key
 */
template <typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
struct VectorQuotientFilterContainer<RegisterSize::_32bit, SIMD::AVX2, OptimizationParameter, fSize, Hasher, Vector, Addresser> {
  using T = typename Vector::T;
  using OP = OptimizationParameter;

  static constexpr uint8_t fingerprintSize{fSize}; //size of fingerprints in bit
  static constexpr uint8_t blockCapacity{(fSize == 16) ? 28 : 48}; //capacity of fingerprints in a block
  static constexpr uint8_t numBuckets{(fSize == 16) ? 36 : 80}; //number of buckets in a block
  static constexpr uint8_t shortcutOpt{(uint8_t) (0.75 * blockCapacity)}; //shortcut optimisation for insert: threshold for checking the occupancy of second block

  static_assert(sizeof(VQFBlock<SIMD::AVX2, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>) == 64); //static_assert to check that blocks always fit into a cache line
  std::vector<VQFBlock<SIMD::AVX2, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>> blocks{}; //vector of blocks for the filter
  
  Addresser addresser; //block addresser
  T numBlocks{2}; //number of blocks in the filter

  /**
   * @brief Construct a new Vector Quotient Filter Container object with only 2 blocks
   * 
   */
  VectorQuotientFilterContainer() {
    //computes block range in the filter
    addresser = std::move(Addresser(&numBlocks, 1));
    numBlocks = addresser.get_size(0);

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX2, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param histogram contains number of elements to store in the filter at index 0
   */
  VectorQuotientFilterContainer(size_t s, const T *histogram) {
    const T h{(T) ceil((size_t) histogram[0] * s / 100 / blockCapacity)}; //number of blocks required to store elements
    
    //computes block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX2, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param length number of elements to store in the filter
   */
  VectorQuotientFilterContainer(size_t s, size_t length) {
    const T h{(T) ceil(length * s / 100 / blockCapacity)}; //number of blocks required to store elements

    //computes block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX2, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief returns total amount of space allocated for fingerprint arrays
   * 
   * @return size_t space allocated for fingerprint arrays in bytes
   */
  size_t size() const {
    return (size_t) numBlocks * blockCapacity * fingerprintSize/8;
  }

  /**
   * @brief looks up an input key in the filter
   * 
   * @param x input key
   * @return bool true if lookup was positive, false otherwise 
   */
  forceinline
  bool lookup(uint32_t x) const {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint32_t queryHash{(uint32_t) _mm256_extract_epi32(v.vector, 0)}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (queryHash & 65535) : (queryHash & 255))};
    uint32_t blockIndex1{(uint32_t) _mm256_extract_epi32(w.vector, 0)};
    uint32_t blockIndex2{(blockIndex1 ^ queryHash)};
    uint16_t bucketIndex{(uint16_t) (queryHash >> fingerprintSize)};

    //bringing indexes into range
    blockIndex2 = ((uint64_t) blockIndex2 * numBlocks) >> 32;
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in Lookup! Block Index invalid!\n";
      return false;
    }

    //looks up fingerprint in both blocks
    return ((blocks.at(blockIndex1).lookup(bucketIndex, fingerprint) != -1) || (blocks.at(blockIndex2).lookup(bucketIndex, fingerprint) != -1));
  }

  /**
   * @brief looks up an array of input keys in the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of positive lookups
   */
  forceinline unroll_loops
  size_t lookup(const uint32_t *values, size_t length) const {
    size_t result{};
    uint32_t *blockIndex1Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *blockIndex2Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *bucketIndexVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *fingerprintVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    Vector v, w;
    __m256i _data, _blockIndexes1{}, _blockIndexes2{}, _bucketIndexes{};
    const uint8_t loopRemainder{(uint8_t) (length & 7)}; //number of keys in values that need to be looked up scalarly
    const __m256i _mNumBlocks{_mm256_set1_epi32(numBlocks)}; //numBlocks in AVX2 register
    const __m256i _mNumBuckets{_mm256_set1_epi32(numBuckets)}; //numBuckets in AVX2 register
    const __m256i _fingerprintMask{(fingerprintSize == 16) ? _mm256_set1_epi32(65535) : _mm256_set1_epi32(255)}; //mask to compute fingerprint

    //vectorized loop; processes array in batches of eight
    for (size_t i{}; i < length - loopRemainder; i+=8) {
      v = Vector::load(values + i); //v <- input keys
      w = Hasher::hash(v);
      _data = w.vector;

      //computation of first block index
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      _blockIndexes1 = w.vector;

      //computation of second block index, bucket index and fingerprint
      _blockIndexes2 = _mm256_xor_si256(_blockIndexes1, _data);
      _bucketIndexes = _mm256_srli_epi32(_data, fingerprintSize);
      _data = _mm256_and_si256(_data, _fingerprintMask);

      //bringing indexes into range
      _blockIndexes2 = vectorizedMulhi_epu32(_blockIndexes2, _mNumBlocks);
      _bucketIndexes = _mm256_mulhi_epu16(_bucketIndexes, _mNumBuckets);
      _bucketIndexes = _mm256_slli_epi32(_bucketIndexes, 16);
      _bucketIndexes = _mm256_srli_epi32(_bucketIndexes, 16);
 
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex1Vector), _blockIndexes1);
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex2Vector), _blockIndexes2);
      _mm256_store_si256(reinterpret_cast<__m256i *>(bucketIndexVector), _bucketIndexes);
      _mm256_store_si256(reinterpret_cast<__m256i *>(fingerprintVector), _data);

      //looks up current content of parameter arrays
      for (size_t j{}; j < 8; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in Lookup! Block Index invalid!\n";
          return result;
        }

        //looks up current fingerprint in both blocks
        if ((blocks.at(blockIndex1Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1) || (blocks.at(blockIndex2Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1))
          result++;
      }
    }
    
    //looks up the last loopRemainder elements in values scalarly
    for (size_t i{length - loopRemainder}; i < length; i++) {
      if (lookup(values[i]))
        result++;
    }

    delete[] fingerprintVector;
    delete[] bucketIndexVector;
    delete[] blockIndex2Vector;
    delete[] blockIndex1Vector;

    return result;
  }

  /**
   * @brief inserts an input key into the filter
   * 
   * @param x input key
   * @return bool true if insert was successful, false otherwise 
   */
  forceinline
  bool insert(uint32_t x) {
    uint32_t i{}; //index of emptier block
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint32_t insertHash{(uint32_t) _mm256_extract_epi32(v.vector, 0)}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (insertHash & 65535) : (insertHash & 255))};
    uint32_t blockIndex1{(uint32_t) _mm256_extract_epi32(w.vector, 0)};
    uint32_t blockIndex2{(blockIndex1 ^ insertHash)};
    uint16_t bucketIndex{(uint16_t) (insertHash >> fingerprintSize)};

    //bringing indexes into range
    blockIndex2 = ((uint64_t) blockIndex2 * numBlocks) >> 32;
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;

    i = blockIndex1;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in Insert! Block Index invalid!\n";
      return false;
    }

    //shortcut optimization: block 2 is only checked if block1's occupancy is more than 75%
    if (blocks.at(blockIndex1).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
      //select index i of emptier block
      if (blocks.at(blockIndex2).select(numBuckets - 1) < blocks.at(blockIndex1).select(numBuckets - 1))
        i = blockIndex2;
    }

    //inserts fingerprint into emptier block
    return blocks.at(i).insert(bucketIndex, fingerprint);
  }

  /**
   * @brief inserts an array of input keys into the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of successful inserts
   */
  forceinline unroll_loops
  size_t insert(const uint32_t *values, size_t length) {
    uint32_t index{}; //index of emptier block
    size_t result{};
    uint32_t *blockIndex1Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *blockIndex2Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *bucketIndexVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *fingerprintVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    Vector v, w;
    __m256i _data, _blockIndexes1{}, _blockIndexes2{}, _bucketIndexes{};
    const uint8_t loopRemainder{(uint8_t) (length & 7)}; //remainder of values that needs to be inserted scalarly
    const __m256i _mNumBlocks{_mm256_set1_epi32(numBlocks)}; //numBlocks in AVX2 register
    const __m256i _mNumBuckets{_mm256_set1_epi32(numBuckets)}; //numBuckets in AVX2 register
    const __m256i _fingerprintMask{(fingerprintSize == 16) ? _mm256_set1_epi32(65535) : _mm256_set1_epi32(255)}; //mask to compute fingerprints

    //vectorized loop; processes array in batches of eight
    for (size_t i{}; i < length - loopRemainder; i+=8) {
      v = Vector::load(values + i);
      w = Hasher::hash(v);
      _data = w.vector;

      //computation of first block index
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      _blockIndexes1 = w.vector;

      //computation of second block index, bucket index and fingerprint
      _blockIndexes2 = _mm256_xor_si256(_blockIndexes1, _data);
      _bucketIndexes = _mm256_srli_epi32(_data, fingerprintSize);
      _data = _mm256_and_si256(_data, _fingerprintMask);

      //bringing indexes into range
      _blockIndexes2 = vectorizedMulhi_epu32(_blockIndexes2, _mNumBlocks);
      _bucketIndexes = _mm256_mulhi_epu16(_bucketIndexes, _mNumBuckets);
      _bucketIndexes = _mm256_slli_epi32(_bucketIndexes, 16);
      _bucketIndexes = _mm256_srli_epi32(_bucketIndexes, 16);
 
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex1Vector), _blockIndexes1);
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex2Vector), _blockIndexes2);
      _mm256_store_si256(reinterpret_cast<__m256i *>(bucketIndexVector), _bucketIndexes);
      _mm256_store_si256(reinterpret_cast<__m256i *>(fingerprintVector), _data);

      //inserts current contents of parameter arrays into filter
      for (size_t j{}; j < 8; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in Insert! Block Index invalid!\n";
          return result;
        }

        index = blockIndex1Vector[j];

        //shortcut optimization: block2 is only checked if block1's occupancy is more than 75%
        if (blocks.at(blockIndex1Vector[j]).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
          //select index of emptier block
          if (blocks.at(blockIndex2Vector[j]).select(numBuckets - 1) < blocks.at(index).select(numBuckets - 1))
            index = blockIndex2Vector[j];
        }

        //inserting current fingerprint into emptier block
        if (blocks.at(index).insert(bucketIndexVector[j], fingerprintVector[j]))
          result++;
      }
    }
    
    //inserts the last loopRemainder elements in values scalarly
    for (size_t i{length - loopRemainder}; i < length; i++) {
      if (insert(values[i]))
        result++;
    }

    delete[] fingerprintVector;
    delete[] bucketIndexVector;
    delete[] blockIndex2Vector;
    delete[] blockIndex1Vector;

    return result;
  }

  /**
   * @brief removes an input key from the filter
   * 
   * @param x input key
   * @return bool true if remove was successful, false otherwise 
   */
  forceinline
  bool remove(uint32_t x) {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const uint32_t removeHash{(uint32_t) _mm256_extract_epi32(v.vector, 0)}; //hash value of x
    const uint16_t fingerprint{(uint16_t) ((fingerprintSize == 16) ? (removeHash & 65535) : (removeHash & 255))};
    uint32_t blockIndex1{(uint32_t) _mm256_extract_epi32(w.vector, 0)};
    uint32_t blockIndex2{(blockIndex1 ^ removeHash)};
    uint16_t bucketIndex{(uint16_t) (removeHash >> fingerprintSize)};

    //bringing indexes into range
    blockIndex2 = ((uint64_t) blockIndex2 * numBlocks) >> 32;
    bucketIndex = ((uint32_t) (uint16_t) bucketIndex * numBuckets) >> 16;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in Remove! Block Index invalid!\n";
      return false;
    }

    //trying to remove fingerprint from first block
    if (blocks.at(blockIndex1).remove(bucketIndex, fingerprint)) {
      return true;
    } else {
      //trying to remove fingerprint from second block, if first remove was unsuccessful
      return blocks.at(blockIndex2).remove(bucketIndex, fingerprint);
    }
  }

  /**
   * @brief removes an array of input key from the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of successful removes 
   */
  forceinline unroll_loops
  size_t remove(const uint32_t *values, size_t length) {
    size_t result{};
    uint32_t *blockIndex1Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *blockIndex2Vector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *bucketIndexVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    uint32_t *fingerprintVector{reinterpret_cast<uint32_t *>(aligned_alloc(32, 32))};
    Vector v, w;
    __m256i _data, _blockIndexes1{}, _blockIndexes2{}, _bucketIndexes{};
    const uint8_t loopRemainder{(uint8_t) (length & 7)}; //remainder of values that need to be inserted scalarly
    const __m256i _mNumBlocks{_mm256_set1_epi32(numBlocks)}; //numBlocks in AVX2 register
    const __m256i _mNumBuckets{_mm256_set1_epi32(numBuckets)}; //numBuckets in AVX2 register
    const __m256i _fingerprintMask{(fingerprintSize == 16) ? _mm256_set1_epi32(65535) : _mm256_set1_epi32(255)}; //mask to compute fingerprints

    //vectorized loop; processes array in batches of eight
    for (size_t i{}; i < length - loopRemainder; i+=8) {
      v = Vector::load(values + i);
      w = Hasher::hash(v);
      _data = w.vector;

      //computation of first block index
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      _blockIndexes1 = w.vector;

      //computation of second block index, bucket index and fingerprint
      _blockIndexes2 = _mm256_xor_si256(_blockIndexes1, _data);
      _bucketIndexes = _mm256_srli_epi32(_data, fingerprintSize);
      _data = _mm256_and_si256(_data, _fingerprintMask);

      //bringing indexes into range
      _blockIndexes2 = vectorizedMulhi_epu32(_blockIndexes2, _mNumBlocks);
      _bucketIndexes = _mm256_mulhi_epu16(_bucketIndexes, _mNumBuckets);
      _bucketIndexes = _mm256_slli_epi32(_bucketIndexes, 16);
      _bucketIndexes = _mm256_srli_epi32(_bucketIndexes, 16);
 
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex1Vector), _blockIndexes1);
      _mm256_store_si256(reinterpret_cast<__m256i *>(blockIndex2Vector), _blockIndexes2);
      _mm256_store_si256(reinterpret_cast<__m256i *>(bucketIndexVector), _bucketIndexes);
      _mm256_store_si256(reinterpret_cast<__m256i *>(fingerprintVector), _data);

      //removes current contents of parameter arrays from the filter
      for (size_t j{}; j < 8; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in Remove! Block Index invalid!\n";
          return result;
        }

        //tries to remove current fingerprint from first block
        if (blocks.at(blockIndex1Vector[j]).remove(bucketIndexVector[j], fingerprintVector[j])) {
          result++;
          
          //tries to remove current fingerprint from second block
        } else if (blocks.at(blockIndex2Vector[j]).remove(bucketIndexVector[j], fingerprintVector[j])) {
          result++;
        }
      }
    }
    
    //removes the last loopRemainder elements in values from the filter scalarly
    for (size_t i{length - loopRemainder}; i < length; i++) {
      if (remove(values[i]))
        result++;
    }

    delete[] fingerprintVector;
    delete[] bucketIndexVector;
    delete[] blockIndex2Vector;
    delete[] blockIndex1Vector;

    return result;
  }
  
  private:
    /**
     * @brief computes mulhi (high 32 bit of multiplication) of packed 32 bit unsigned integers in two AVX2 registers
     * 
     * @param _x register of first multiplicands
     * @param _y register of second multiplicands
     * @return __m256i result of the mulhi operation
     */
    forceinline
    __m256i vectorizedMulhi_epu32(const __m256i _x, const __m256i _y) const {
      __m256i _xEven{_mm256_srli_epi64(_x, 32)}; //elements of _x at even indexes at required positions for multiplication
      __m256i _xOdd{_x}; //elements of _x at odd indexes at required positions for multiplication
      const __m256i _yEven{_mm256_srli_epi64(_y, 32)}; //elements of _y at even indexes at required positions for multiplication
      const __m256i _yOdd{_y}; //elements of _y at odd indexes at required positions for multiplication
      
      //multiplication
      _xEven = _mm256_mul_epu32(_xEven, _yEven);
      _xOdd = _mm256_mul_epu32(_xOdd, _yOdd);
      
      //shifts high 32 bit of odd indexes to correct positions
      _xOdd = _mm256_srli_epi64(_xOdd, 32);
      
      //blends result together
      return _mm256_blend_epi32(_xEven, _xOdd, 0b01010101);
    }
};

/**
 * @brief struct for the entire Vector Quotient Filter container using AVX512
 * 
 * @tparam paramSize determines whether the filter takes 64 bit or 32 bit input keys
 * @tparam OptimizationParameter determines if the filter uses multi-threading
 * @tparam fSize sets the fingerprint size
 * @tparam Hasher determines the hashing procedure the filter uses for input keys
 * @tparam Vector enclosing type for SIMD registers with Hasher/Addresser
 * @tparam Addresser determines the addressing procedure the filter uses for block indexes of an input key
 */
template <RegisterSize paramSize, typename OptimizationParameter, size_t fSize, typename Hasher, typename Vector, typename Addresser>
struct VectorQuotientFilterContainer<paramSize, SIMD::AVX512, OptimizationParameter, fSize, Hasher, Vector, Addresser> {
  using T = typename Vector::T;
  using OP = OptimizationParameter;

  static constexpr uint8_t fingerprintSize{fSize}; //fingerprint size in bits
  static constexpr uint8_t blockCapacity{(fSize == 16) ? 28 : 48}; //capacity of fingerprints in a block
  static constexpr uint8_t numBuckets{(fSize == 16) ? 36 : 80}; //number of buckets in a block
  static constexpr uint8_t shortcutOpt{(uint8_t) (0.75 * blockCapacity)}; //shortcut optimisation for insert: threshold for checking of occupancy of second block

  static_assert(sizeof(VQFBlock<SIMD::AVX512, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>) == 64); //static_assert to check that blocks always fit into a cache line
  std::vector<VQFBlock<SIMD::AVX512, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>> blocks{}; //vector of blocks in the filter
  Addresser addresser; //block addresser
  T numBlocks{2}; //number of blocks in the filter

  /**
   * @brief Construct a new Vector Quotient Filter Container object with only two blocks
   * 
   */
  VectorQuotientFilterContainer() {
    //computes block range in the filter
    addresser = std::move(Addresser(&numBlocks, 1));
    numBlocks = addresser.get_size(0);

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX512, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param histogram contains number of elements to store in the filter at index 0
   */
  VectorQuotientFilterContainer(size_t s, const T *histogram) {
    const T h{(T) ceil((size_t) histogram[0] * s / 100 / blockCapacity)}; //number of blocks required to store the elements in the filter

    //computes block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX512, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief Construct a new Vector Quotient Filter Container object
   * 
   * @param s scaling factor for space allocation to ensure build does not fail
   * @param length number of elements to store in the filter
   */
  VectorQuotientFilterContainer(size_t s, size_t length) {
    const T h{(T) ceil(length * s / 100 / blockCapacity)}; //number of blocks required to store the elements in the filter
  
    //computes block range in the filter
    if (h > 1) {
      addresser = std::move(Addresser(&h, 1));
      numBlocks = addresser.get_size(0);
    } else {
      addresser = std::move(Addresser(&numBlocks, 1));
      numBlocks = addresser.get_size(0);
    }

    blocks.resize(numBlocks, VQFBlock<SIMD::AVX512, OP::multiThreading, fingerprintSize, blockCapacity, numBuckets>());
  }

  /**
   * @brief returns total amount of space allocated for the fingerprint arrays
   * 
   * @return space allocated for the fingerprint arrays in bytes
   */
  size_t size() const {
    return (size_t) numBlocks * blockCapacity * fingerprintSize/8;
  }

  /**
   * @brief looks up an input key in the filter
   * 
   * @param x input key
   * @return bool true if lookup positive, false otherwise 
   */
  forceinline
  bool lookup(uint64 x) const {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const Vector y{addresser.compute_alternative_address_vertical(0, w, simd::extractAddressBits(Vector(x)))}; //vector containing second block index of x
    alignas(64) uint64_t queryHash{}; //hash value of x
    alignas(64) uint64_t blockIndex1{};
    alignas(64) uint64_t blockIndex2{};
    uint16_t fingerprint{};
    uint64_t bucketIndex{};

    //initialisation of variables
    if constexpr (paramSize == RegisterSize::_64bit) {
      _mm512_mask_store_epi64(&queryHash, 1, v.vector);
      _mm512_mask_store_epi64(&blockIndex1, 1, w.vector);
      _mm512_mask_store_epi64(&blockIndex2, 1, y.vector);
    } else {
      _mm512_mask_store_epi32(&queryHash, 1, v.vector);
      _mm512_mask_store_epi32(&blockIndex1, 1, w.vector);
      _mm512_mask_store_epi32(&blockIndex2, 1, y.vector);
    }

    fingerprint = (fingerprintSize == 16) ? (queryHash & 65535) : (queryHash & 255);
    
    //bringing bucket index into range
    bucketIndex = ((uint32_t) (uint16_t) (queryHash >> fingerprintSize) * numBuckets) >> 16;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in scalar lookup! Block Index invalid!\n";
      return false;
    }

    //looking up the fingerprint in both blocks
    return ((blocks.at(blockIndex1).lookup(bucketIndex, fingerprint) != -1) || (blocks.at(blockIndex2).lookup(bucketIndex, fingerprint) != -1));
  }

  /**
   * @brief looks up an array of input keys in the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of positive lookups 
   */
  forceinline unroll_loops
  size_t lookup(const T *values, size_t length) const {
    T index{}; //index of emptier block
    size_t i{}; //loop counter
    size_t result{};
    static constexpr uint8_t loopStep{((uint8_t) paramSize == 64) ? 8 : 16}; //number of keys that can be loaded into one 512 bit register
    static constexpr uint8_t arrLen = 64/sizeof(T); //
    const uint8_t loopRemainder{(uint8_t) (length & (loopStep - 1))}; //number of keys at the end of values that do not fill an entire 512 bit register
    T *blockIndex1Vector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *blockIndex2Vector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *bucketIndexVector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *fingerprintVector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    uint32_t loopRemainderMask{1}; //mask to load the loop remainder into a 512 bit register
    __m512i _data{}, _blockIndexes1{}, _blockIndexes2{}, _bucketIndexes{};
    __m512i _mNumBlocks{}; //numBlocks in AVX512 register
    __m512i _mNumBuckets{}; //numBuckets in AVX512 register
    __m512i _fingerprintMask{}; //mask to compute fingerprints
    Vector v, w, x;

    //initialisation of variables depending on paramSize
    if constexpr (paramSize == RegisterSize::_64bit) {
      _mNumBlocks = _mm512_set1_epi64(numBlocks);
      _mNumBuckets = _mm512_set1_epi64(numBuckets);
      _fingerprintMask = _mm512_set1_epi64((fingerprintSize == 16) ? 65535 : 255);
        
    } else {
      _mNumBlocks = _mm512_set1_epi32(numBlocks);
      _mNumBuckets = _mm512_set1_epi32(numBuckets);
      _fingerprintMask = _mm512_set1_epi32((fingerprintSize == 16) ? 65535 : 255);
    }

    //vectorized loop; processes values in batches of eight or 16, depending on paramSize
    for (i = 0; i < length - loopRemainder; i+=loopStep) {
      v = Vector::load(values + i);

      w = Hasher::hash(v);
      _data = w.vector;

      //computation of block indexes
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      x = addresser.compute_alternative_address_vertical(0, simd::extractAddressBits(w), v);
      _blockIndexes1 = w.vector;
      _blockIndexes2 = x.vector;
      
      //initialisation of bucket index
      if constexpr (paramSize == RegisterSize::_64bit) {
        _bucketIndexes = _mm512_srli_epi64(_data, fingerprintSize);
      } else {
        _bucketIndexes = _mm512_srli_epi32(_data, fingerprintSize);
      }

      //computation of fingerprint and bucket index
      _data = _mm512_and_si512(_data, _fingerprintMask);
      _bucketIndexes = vectorizedComputeBucketIndex(_bucketIndexes, _mNumBuckets);
 
      _mm512_store_si512(reinterpret_cast<__m512i *>(blockIndex1Vector), _blockIndexes1);
      _mm512_store_si512(reinterpret_cast<__m512i *>(blockIndex2Vector), _blockIndexes2);
      _mm512_store_si512(reinterpret_cast<__m512i *>(bucketIndexVector), _bucketIndexes);
      _mm512_store_si512(reinterpret_cast<__m512i *>(fingerprintVector), _data);

      //looks up the current contents of the parameter arrays
      for (size_t j{}; j < arrLen; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in vectorized lookup! Block Index invalid!\n";
          return result;
        }

        //looks up the current fingerprint in both blocks
        if ((blocks.at(blockIndex1Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1) || (blocks.at(blockIndex2Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1))
          result++;
      }

    }

    //masked computation for the loopRemainder
    if (loopRemainder > 0) {
      //initialisation of mask for the mask load
      loopRemainderMask <<= loopRemainder;
      loopRemainderMask -= 1;

      if constexpr (paramSize == RegisterSize::_64bit) {
        v = Vector(_mm512_maskz_load_epi64(loopRemainderMask, values + i));
      } else {
        v = Vector(_mm512_maskz_load_epi32(loopRemainderMask, values + i));
      }

      w = Hasher::hash(v);
      _data = w.vector;

      //computation of block indexes
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      x = addresser.compute_alternative_address_vertical(0, simd::extractAddressBits(w), v);
      _blockIndexes1 = w.vector;
      _blockIndexes2 = x.vector;
        
      //initialisation of bucket index
      if constexpr (paramSize == RegisterSize::_64bit) {
        _bucketIndexes = _mm512_srli_epi64(_data, fingerprintSize);
      } else {
        _bucketIndexes = _mm512_srli_epi32(_data, fingerprintSize);
      }

      //computation of fingerprint and bucket index
      _data = _mm512_and_si512(_data, _fingerprintMask);
      _bucketIndexes = vectorizedComputeBucketIndex(_bucketIndexes, _mNumBuckets);

      if constexpr (paramSize == RegisterSize::_64bit) {
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(blockIndex1Vector), loopRemainderMask, _blockIndexes1);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(blockIndex2Vector), loopRemainderMask, _blockIndexes2);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(bucketIndexVector), loopRemainderMask, _bucketIndexes);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(fingerprintVector), loopRemainderMask, _data);
      } else {
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(blockIndex1Vector), loopRemainderMask, _blockIndexes1);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(blockIndex2Vector), loopRemainderMask, _blockIndexes2);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(bucketIndexVector), loopRemainderMask, _bucketIndexes);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(fingerprintVector), loopRemainderMask, _data);
      }

      //looks up the current contents of the parameter arrays
      for (size_t j{}; j < loopRemainder; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in vectorized lookup! Block Index invalid!\n";
          return result;
        }

        //looks up the current fingerprint in both blocks
        if ((blocks.at(blockIndex1Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1) || (blocks.at(blockIndex2Vector[j]).lookup(bucketIndexVector[j], fingerprintVector[j]) != -1))
          result++;
      }
    }

    delete[] fingerprintVector;
    delete[] bucketIndexVector;
    delete[] blockIndex2Vector;
    delete[] blockIndex1Vector;

    return result;
  }

  /**
   * @brief inserts an input key into the filter
   * 
   * @param x input key
   * @return bool true if insert was successful, false otherwise
   */
  forceinline
  bool insert(uint64_t x) {
    const Vector v{Hasher::hash(Vector(x))}; //vector containing hash value of x
    const Vector w{addresser.compute_address_vertical(0, simd::extractAddressBits(Vector(x)))}; //vector containing first block index of x
    const Vector y{addresser.compute_alternative_address_vertical(0, w, simd::extractAddressBits(Vector(x)))}; //vector containing second block index of x
    alignas(64) uint64_t insertHash{}; //hash value of x
    alignas(64) uint64_t blockIndex1{};
    alignas(64) uint64_t blockIndex2{};
    uint64_t i{}; //index of emptier block
    uint16_t fingerprint{};
    uint64_t bucketIndex{};

    //initialisation of variables
    if constexpr (paramSize == RegisterSize::_64bit) {
      _mm512_mask_store_epi64(&insertHash, 1, v.vector);
      _mm512_mask_store_epi64(&blockIndex1, 1, w.vector);
      _mm512_mask_store_epi64(&blockIndex2, 1, y.vector);
    } else {
      _mm512_mask_store_epi32(&insertHash, 1, v.vector);
      _mm512_mask_store_epi32(&blockIndex1, 1, w.vector);
      _mm512_mask_store_epi32(&blockIndex2, 1, y.vector);
    }

    //bringing bucket index into range
    bucketIndex = ((uint32_t) (uint16_t) (insertHash >> fingerprintSize) * numBuckets) >> 16;
    
    fingerprint = (fingerprintSize == 16) ? (insertHash & 65535) : (insertHash & 255);
    i = blockIndex1;

    if ((blockIndex1 >= numBlocks) || (blockIndex2 >= numBlocks)) {
      std::cout << "Error in scalar insert! Block Index invalid!\n";
      return false;
    }

    //shortcut optimization: block 2 is only checked if block1's occupancy is more than 75%
    if (blocks.at(blockIndex1).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
      //select index i of emptier block
      if (blocks.at(blockIndex2).select(numBuckets - 1) < blocks.at(blockIndex1).select(numBuckets - 1))
        i = blockIndex2;
    }

    //inserts fingerprint into emptier block
    return blocks.at(i).insert(bucketIndex, fingerprint);
  }

  /**
   * @brief inserts an array of input keys into the filter
   * 
   * @param values array of input keys
   * @param length length of array
   * @return size_t number of successful inserts 
   */
  forceinline unroll_loops
  size_t insert(const T *values, size_t length) {
    T index{}; //index of emptier block
    size_t i{}; //loop counter
    size_t result{};
    static constexpr uint8_t loopStep{((uint8_t) paramSize == 64) ? 8 : 16}; //number of keys that fit into a 512 bit register
    static constexpr uint8_t arrLen = 64/sizeof(T);
    const uint8_t loopRemainder{(uint8_t) (length & (loopStep - 1))}; //remainder of values that needs to be accessed using a mask load
    T *blockIndex1Vector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *blockIndex2Vector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *bucketIndexVector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    T *fingerprintVector{reinterpret_cast<T *>(aligned_alloc(64, 64))};
    uint32_t loopRemainderMask{1}; //mask for the mask load of the loop remainder
    __m512i _data{}, _blockIndexes1{}, _blockIndexes2{}, _bucketIndexes{};
    __m512i _mNumBlocks{}; //numBlocks in AVX512 register
    __m512i _mNumBuckets{}; //numBuckets in AVX512 register
    __m512i _fingerprintMask{}; //mask to compute fingerprints
    Vector v, w, x;

    //initialisation of variables
    if constexpr (paramSize == RegisterSize::_64bit) {
      _mNumBlocks = _mm512_set1_epi64(numBlocks);
      _mNumBuckets = _mm512_set1_epi64(numBuckets);
      _fingerprintMask = _mm512_set1_epi64((fingerprintSize == 16) ? 65535 : 255);
        
    } else {
      _mNumBlocks = _mm512_set1_epi32(numBlocks);
      _mNumBuckets = _mm512_set1_epi32(numBuckets);
      _fingerprintMask = _mm512_set1_epi32((fingerprintSize == 16) ? 65535 : 255);
    }

    //vectorized loop; processes values in batches of 8 or 16, depending on paramSize
    for (i = 0; i < length - loopRemainder; i+=loopStep) {
      v = Vector::load(values + i);

      w = Hasher::hash(v);
      _data = w.vector;

      //computation of block indexes
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      x = addresser.compute_alternative_address_vertical(0, simd::extractAddressBits(w), v);
      _blockIndexes1 = w.vector;
      _blockIndexes2 = x.vector;
        
      //initialisation of bucket index
      if constexpr (paramSize == RegisterSize::_64bit) {
        _bucketIndexes = _mm512_srli_epi64(_data, fingerprintSize);
      } else {
        _bucketIndexes = _mm512_srli_epi32(_data, fingerprintSize);
      }

      //computation of fingerprint and bucket index
      _data = _mm512_and_si512(_data, _fingerprintMask);
      _bucketIndexes = vectorizedComputeBucketIndex(_bucketIndexes, _mNumBuckets);
 
      _mm512_store_si512(reinterpret_cast<__m512i *>(blockIndex1Vector), _blockIndexes1);
      _mm512_store_si512(reinterpret_cast<__m512i *>(blockIndex2Vector), _blockIndexes2);
      _mm512_store_si512(reinterpret_cast<__m512i *>(bucketIndexVector), _bucketIndexes);
      _mm512_store_si512(reinterpret_cast<__m512i *>(fingerprintVector), _data);

      //processes current contents of parameter arrays
      for (size_t j{}; j < arrLen; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in vectorized insert! Block Index invalid!\n";
          return result;
        }

        index = blockIndex1Vector[j];

        //shortcut optimization: block2 is only checked if block1's occupancy is more than 75%
        if (blocks.at(blockIndex1Vector[j]).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
          //select index of emptier block
          if (blocks.at(blockIndex2Vector[j]).select(numBuckets - 1) < blocks.at(index).select(numBuckets - 1))
            index = blockIndex2Vector[j];
        }

        //inserts current fingerprint into emptier block
        if (blocks.at(index).insert(bucketIndexVector[j], fingerprintVector[j]))
          result++;
      }
    }

    //masked computation for the rest of the values
    if (loopRemainder > 0) {
      //initialisation of mask for mask load
      loopRemainderMask <<= loopRemainder;
      loopRemainderMask -= 1;

      if constexpr (paramSize == RegisterSize::_64bit) {
        v = Vector(_mm512_maskz_load_epi64(loopRemainderMask, values + i));
      } else {
        v = Vector(_mm512_maskz_load_epi32(loopRemainderMask, values + i));
      }

      w = Hasher::hash(v);
      _data = w.vector;

      //computation of block indexes
      w = addresser.compute_address_vertical(0, simd::extractAddressBits(v));
      x = addresser.compute_alternative_address_vertical(0, simd::extractAddressBits(w), v);
      _blockIndexes1 = w.vector;
      _blockIndexes2 = x.vector;
      
      //initialisation of bucket index
      if constexpr (paramSize == RegisterSize::_64bit) {
        _bucketIndexes = _mm512_srli_epi64(_data, fingerprintSize);
      } else {
        _bucketIndexes = _mm512_srli_epi32(_data, fingerprintSize);
      }

      //computation of bucket index and fingerprint
      _data = _mm512_and_si512(_data, _fingerprintMask);
      _bucketIndexes = vectorizedComputeBucketIndex(_bucketIndexes, _mNumBuckets);

      if constexpr (paramSize == RegisterSize::_64bit) {
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(blockIndex1Vector), loopRemainderMask, _blockIndexes1);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(blockIndex2Vector), loopRemainderMask, _blockIndexes2);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(bucketIndexVector), loopRemainderMask, _bucketIndexes);
        _mm512_mask_store_epi64(reinterpret_cast<__m512i *>(fingerprintVector), loopRemainderMask, _data);
      } else {
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(blockIndex1Vector), loopRemainderMask, _blockIndexes1);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(blockIndex2Vector), loopRemainderMask, _blockIndexes2);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(bucketIndexVector), loopRemainderMask, _bucketIndexes);
        _mm512_mask_store_epi32(reinterpret_cast<__m512i *>(fingerprintVector), loopRemainderMask, _data);
      }

      //processes current contents of the parameter arrays
      for (size_t j{}; j < loopRemainder; j++) {
        if ((blockIndex1Vector[j] >= numBlocks) || (blockIndex2Vector[j] >= numBlocks)) {
          std::cout << "Error in vectorized insert! Block Index invalid!\n";
          return result;
        }

        index = blockIndex1Vector[j];

        //shortcut optimization: block2 is only checked if block1's occupancy is more than 75%
        if (blocks.at(blockIndex1Vector[j]).select(numBuckets - 1) - numBuckets + 1 >= shortcutOpt) {
          //select index of emptier block
          if (blocks.at(blockIndex2Vector[j]).select(numBuckets - 1) < blocks.at(index).select(numBuckets - 1))
            index = blockIndex2Vector[j];
        }

        //inserts current fingerprint into emptier block
        if (blocks.at(index).insert(bucketIndexVector[j], fingerprintVector[j]))
          result++;
      }
    }

    delete[] fingerprintVector;
    delete[] bucketIndexVector;
    delete[] blockIndex2Vector;
    delete[] blockIndex1Vector;

    return result;
  }
  
  private:

    /**
     * @brief brings AVX512 register of bucket indexes into range [0, numBuckets)
     * 
     * @param _index AVX512 register consisting of 64/32 bit out of range bucket indexes 
     * @param _mNumBuckets AVX512 register consisting of 64/32 entries of numBuckets
     * @return __m512i AVX512 register consisting of the 64/32 bit results of the operation 
     */
    forceinline
    __m512i vectorizedComputeBucketIndex(const __m512i _index, const __m512i _mNumBuckets) const {
      __m512i _temp{}; //intermediate result

      //16 bit mulhi operation, the unnecessary bits are shifted out
      if constexpr (paramSize == RegisterSize::_64bit) {
        _temp = _mm512_mulhi_epu16(_index, _mNumBuckets);
        _temp = _mm512_slli_epi64(_temp, 48);
        return _mm512_srli_epi64(_temp, 48);
      } else {
        _temp = _mm512_mulhi_epu16(_index, _mNumBuckets);
        _temp = _mm512_slli_epi32(_temp, 16);
        return _mm512_srli_epi32(_temp, 16);
      }
    }
  };
}
