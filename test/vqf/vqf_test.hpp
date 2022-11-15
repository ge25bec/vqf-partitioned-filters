#pragma once

#include <cstddef>
#include "vqf/vqf_parameter.hpp"
#include "../filter_test.hpp"

namespace test::vqf {

    namespace parameter = filters::parameter;
    namespace vqf = filters::vqf;

    static constexpr size_t n_s = 1000, n_l = 1000000;

    static constexpr size_t k_s = 8, k_l = 16;

    static constexpr size_t n_partition_s = 2, n_partition_l = 16;

    static constexpr filters::FilterType VectorQuotientFilter = filters::FilterType::VectorQuotientFilter;

    /*
     * Small Test Types
     */
    
    //power of two addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar641 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar642 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar643 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar644 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulScalar64, n_s, s, 0, 2, 2, expected_fp>;

    //magic addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar645 = FilterTestConfig<VectorQuotientFilter,
            FP, k_s, parameter::MagicMurmurScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar646 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar647 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicTwoIndependentMultiplyShiftScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar648 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulScalar64, n_s, s, 0, 2, 2, expected_fp>;

    //lemire addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar649 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar6410 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar6411 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireTwoIndependentMultiplyShiftScalar64, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar6412 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar64, n_s, s, 0, 2, 2, expected_fp>;


    //power of two addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar321 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar322 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar323 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar324 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulScalar32, n_s, s, 0, 2, 2, expected_fp>;
    
    //magic addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar325 = FilterTestConfig<VectorQuotientFilter,
            FP, k_s, parameter::MagicMurmurScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar326 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar327 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicTwoIndependentMultiplyShiftScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar328 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulScalar32, n_s, s, 0, 2, 2, expected_fp>;
        
    //lemire addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar329 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar3210 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar3211 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireTwoIndependentMultiplyShiftScalar32, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalar3212 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar32, n_s, s, 0, 2, 2, expected_fp>;


    //power of two addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX21 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX22 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX23 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX232, n_s, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX24 = FilterTestConfig<VectorQuotientFilter,
            FP, k_s, parameter::MagicMurmurAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX25 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX26 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX232, n_s, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX27 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX28 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX29 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX232, n_s, s, 0, 2, 2, expected_fp>;


    //power of two addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_s, parameter::MagicMurmurAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_5 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_6 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX51264, n_s, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_7 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_8 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX51264, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51264_9 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX51264, n_s, s, 0, 2, 2, expected_fp>;


    //power of two addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_5 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_6 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX51232, n_s, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_7 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_8 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX51232, n_s, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallAVX51232_9 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX51232, n_s, s, 0, 2, 2, expected_fp>;



    //power of two addressing scalar 64 bit multi-threaded
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalarMT1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar64MT, n_s, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalarMT2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar64MT, n_s, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalarMT3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar32MT, n_s, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFSmallScalarMT4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar32MT, n_s, s, 0, 4, 16, expected_fp>;

    /*
     * Large Test Types
     */

    //power of two addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar641 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar642 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar643 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar644 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulScalar64, n_l, s, 0, 2, 2, expected_fp>;

    //magic addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar645 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar646 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar647 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicTwoIndependentMultiplyShiftScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar648 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulScalar64, n_l, s, 0, 2, 2, expected_fp>;

    //lemire addressing scalar 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar649 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar6410 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar6411 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireTwoIndependentMultiplyShiftScalar64, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar6412 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar64, n_l, s, 0, 2, 2, expected_fp>;


    //power of two addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar321 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar322 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar323 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar324 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulScalar32, n_l, s, 0, 2, 2, expected_fp>;
    
    //magic addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar325 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar326 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar327 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicTwoIndependentMultiplyShiftScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar328 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulScalar32, n_l, s, 0, 2, 2, expected_fp>;
        
    //lemire addressing scalar 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar329 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar3210 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar3211 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireTwoIndependentMultiplyShiftScalar32, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalar3212 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar32, n_l, s, 0, 2, 2, expected_fp>;


    //power of two addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX21 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX22 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX23 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX232, n_l, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX24 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX25 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX26 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX232, n_l, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx2 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX27 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX28 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX29 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX232, n_l, s, 0, 2, 2, expected_fp>;

    
    //power of two addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_5 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_6 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX51264, n_l, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx512 64 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_7 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_8 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX51264, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51264_9 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX51264, n_l, s, 0, 2, 2, expected_fp>;


    //power of two addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoFasthashAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMulAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    
    //magic addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMurmurAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_5 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicFasthashAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_6 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::MagicMulAVX51232, n_l, s, 0, 2, 2, expected_fp>;

    //lemire addressing avx512 32 bit
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_7 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_8 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireFasthashAVX51232, n_l, s, 0, 2, 2, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeAVX51232_9 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulAVX51232, n_l, s, 0, 2, 2, expected_fp>;


    //power of two addressing scalar 64 bit multi-threaded
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalarMT1 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoMurmurScalar64MT, n_l, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalarMT2 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMurmurScalar64MT, n_l, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalarMT3 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::PowerOfTwoTwoIndependentMultiplyShiftScalar32MT, n_l, s, 0, 4, 16, expected_fp>;
    template<template<size_t> typename FP, size_t s, int64_t expected_fp> using VQFLargeScalarMT4 = FilterTestConfig<VectorQuotientFilter,
            FP, k_l, parameter::LemireMulScalar32MT, n_l, s, 0, 4, 16, expected_fp>;


    using VQFSmall64TestTypes = ::testing::Types<VQFSmallScalar641<vqf::Standard, 125, 1>, VQFSmallScalar642<vqf::Standard, 125, 0>, VQFSmallScalar643<vqf::Standard, 125, 0>, VQFSmallScalar644<vqf::Standard, 125, 0>, 
        VQFSmallScalar645<vqf::Standard, 125, 43>, VQFSmallScalar646<vqf::Standard, 125, 0>, VQFSmallScalar647<vqf::Standard, 125, 0>, VQFSmallScalar648<vqf::Standard, 125, 0>, VQFSmallScalar649<vqf::Standard, 125, 1>, 
        VQFSmallScalar6410<vqf::Standard, 125, 0>, VQFSmallScalar6411<vqf::Standard, 125, 0>, VQFSmallScalar6412<vqf::Standard, 125, 0>>;
    
    using VQFSmall32TestTypes = ::testing::Types<VQFSmallScalar321<vqf::Standard, 125, 0>, VQFSmallScalar322<vqf::Standard, 125, 0>, VQFSmallScalar323<vqf::Standard, 125, 0>, VQFSmallScalar324<vqf::Standard, 125, 0>, 
        VQFSmallScalar325<vqf::Standard, 125, 35>, VQFSmallScalar326<vqf::Standard, 125, 0>, VQFSmallScalar327<vqf::Standard, 125, 0>, VQFSmallScalar328<vqf::Standard, 125, 0>, VQFSmallScalar329<vqf::Standard, 125, 1>, 
        VQFSmallScalar3210<vqf::Standard, 125, 0>, VQFSmallScalar3211<vqf::Standard, 125, 0>, VQFSmallScalar3212<vqf::Standard, 125, 0>, VQFSmallAVX21<vqf::Standard, 125, 0>, VQFSmallAVX22<vqf::Standard, 125, 0>, 
        VQFSmallAVX23<vqf::Standard, 125, 0>, VQFSmallAVX24<vqf::Standard, 125, 35>, VQFSmallAVX25<vqf::Standard, 125, 0>, VQFSmallAVX26<vqf::Standard, 125, 0>, VQFSmallAVX27<vqf::Standard, 125, 1>, VQFSmallAVX28<vqf::Standard, 125, 0>, 
        VQFSmallAVX29<vqf::Standard, 125, 0>>;

    using VQFLarge64TestTypes = ::testing::Types<VQFLargeScalar641<vqf::Standard, 125, 132>, VQFLargeScalar642<vqf::Standard, 125, 120>, VQFLargeScalar643<vqf::Standard, 125, 57>, VQFLargeScalar644<vqf::Standard, 125, 3441770>, 
        VQFLargeScalar645<vqf::Standard, 125, 198>, VQFLargeScalar646<vqf::Standard, 125, 189>, VQFLargeScalar647<vqf::Standard, 125, 120>, VQFLargeScalar648<vqf::Standard, 125, 98>, VQFLargeScalar649<vqf::Standard, 125, 194>, 
        VQFLargeScalar6410<vqf::Standard, 125, 192>, VQFLargeScalar6411<vqf::Standard, 125, 87>, VQFLargeScalar6412<vqf::Standard, 125, 75>>;
    
    using VQFLarge32TestTypes = ::testing::Types<VQFLargeScalar321<vqf::Standard, 125, 148>, VQFLargeScalar322<vqf::Standard, 125, 67>, VQFLargeScalar323<vqf::Standard, 125, 60>, VQFLargeScalar324<vqf::Standard, 125, 3583382>, 
        VQFLargeScalar325<vqf::Standard, 125, 311>, VQFLargeScalar326<vqf::Standard, 125, 311>, VQFLargeScalar327<vqf::Standard, 125, 92>, VQFLargeScalar328<vqf::Standard, 125, 91>, VQFLargeScalar329<vqf::Standard, 125, 337>, 
        VQFLargeScalar3210<vqf::Standard, 125, 256>, VQFLargeScalar3211<vqf::Standard, 125, 94>, VQFLargeScalar3212<vqf::Standard, 125, 86>, VQFLargeAVX21<vqf::Standard, 125, 148>, VQFLargeAVX22<vqf::Standard, 125, 67>, 
        VQFLargeAVX23<vqf::Standard, 125, 3583382>, VQFLargeAVX24<vqf::Standard, 125, 311>, VQFLargeAVX25<vqf::Standard, 125, 311>, VQFLargeAVX26<vqf::Standard, 125, 91>, VQFLargeAVX27<vqf::Standard, 125, 337>, VQFLargeAVX28<vqf::Standard, 125, 256>, 
        VQFLargeAVX29<vqf::Standard, 125, 86>>;

   using VQFAVX512SmallTestTypes = ::testing::Types<VQFSmallAVX51264_1<vqf::Standard, 125, 0>, VQFSmallAVX51264_2<vqf::Standard, 125, 0>, VQFSmallAVX51264_3<vqf::Standard, 125, 0>, 
        VQFSmallAVX51264_4<vqf::Standard, 125, 47>, VQFSmallAVX51264_5<vqf::Standard, 125, 0>, VQFSmallAVX51264_6<vqf::Standard, 125, 0>, VQFSmallAVX51264_7<vqf::Standard, 125, 1>, VQFSmallAVX51264_8<vqf::Standard, 125, 0>, VQFSmallAVX51264_9<vqf::Standard, 125, 0>, 
        VQFSmallAVX51232_1<vqf::Standard, 125, 0>, VQFSmallAVX51232_2<vqf::Standard, 125, 0>, VQFSmallAVX51232_3<vqf::Standard, 125, 0>, VQFSmallAVX51232_7<vqf::Standard, 125, 0>, VQFSmallAVX51232_8<vqf::Standard, 125, 0>, VQFSmallAVX51232_9<vqf::Standard, 125, 0>>;

  using VQFAVX512LargeTestTypes = ::testing::Types<VQFLargeAVX51264_1<vqf::Standard, 125, 135>, VQFLargeAVX51264_2<vqf::Standard, 125, 138>, //VQFLargeAVX51264_3<vqf::Standard, 125, 0>,
        VQFLargeAVX51264_4<vqf::Standard, 125, 187>, VQFLargeAVX51264_5<vqf::Standard, 125, 200>, VQFLargeAVX51264_6<vqf::Standard, 125, 90>, VQFLargeAVX51264_7<vqf::Standard, 125, 174>, VQFLargeAVX51264_8<vqf::Standard, 125, 188>, VQFLargeAVX51264_9<vqf::Standard, 125, 91>, 
        VQFLargeAVX51232_1<vqf::Standard, 125, 152>, VQFLargeAVX51232_2<vqf::Standard, 125, 77>, VQFLargeAVX51232_3<vqf::Standard, 125, 3598917>, VQFLargeAVX51232_7<vqf::Standard, 125, 186>, VQFLargeAVX51232_8<vqf::Standard, 125, 110>, VQFLargeAVX51232_9<vqf::Standard, 125, 99>>;

  using VQFMTScalarTestTypes = ::testing::Types<VQFSmallScalarMT1<vqf::Standard, 125, 1>, VQFSmallScalarMT2<vqf::Standard, 125, 1>, VQFSmallScalarMT3<vqf::Standard, 125, 0>, VQFSmallScalarMT4<vqf::Standard, 125, 0>,
        VQFLargeScalarMT1<vqf::Standard, 125, 132>, VQFLargeScalarMT2<vqf::Standard, 125, 189>, VQFLargeScalarMT3<vqf::Standard, 125, 60>, VQFLargeScalarMT4<vqf::Standard, 125, 95>>;
}
