#include <gtest/gtest.h>
#include "vqf_test.hpp"

namespace test::vqf {

    INSTANTIATE_TYPED_TEST_CASE_P(VQFSmall64TestTypes, FilterTest, VQFSmall64TestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFSmall32TestTypes, FilterTest, VQFSmall32TestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFLarge64TestTypes, FilterTest, VQFLarge64TestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFLarge32TestTypes, FilterTest, VQFLarge32TestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFAVX512SmallTestTypes, FilterTest, VQFAVX512SmallTestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFAVX512LargeTestTypes, FilterTest, VQFAVX512LargeTestTypes);

    INSTANTIATE_TYPED_TEST_CASE_P(VQFMTScalarTestTypes, FilterTest, VQFMTScalarTestTypes);

}

MAIN();