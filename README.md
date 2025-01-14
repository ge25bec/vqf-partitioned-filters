# Partitioned Filters - Vector Quotient Filter
This project is a fork off the *partitioned-filters* repository [3] (see LICENSE in root directory), into which my implementation of the Vector Quotient Filter from the repository https://gitlab.db.in.tum.de/ge25bec/vector-quotient-filter was integrated. The Vector Quotient Filter was introduced by Pandey, et al. [1].
The VectorQuotientFilterContainer structure in src/vqf contains the whole filter. There are template-specialised versions for the structure that support AVX2 or AVX512 intrinsics. The scalar and AVX2 structure support lookups, insertions or removals of input keys, the AVX512 version only supports lookups and insertions.
The VQFBlock structure represents a block inside the filter. There are again template-specialised versions of it for the different SIMD levels. The "select(uint8_t i)" functions and "pdepMask" inside "vqf_block.hpp" are taken from the github repository of the original VQF prototype [2]. The License for the functions and mask can be found in src/vqf.
To circumvent most of the modulo instructions, the *Ross* [5] procedure is used.

The Vector Quotient Filter was integrated into the *partitioned-filters* repository [3], to test, benchmark and plot the filter. The *partitioned-filters* repository was the basis for *A four-dimensional Analysis of Partitioned Approximate Filters* by Schmidt, et al. [4], which sets the metrics of false positive rate, lookup/build performance and space usage to evaluate the filter performance.

## Related Work

[1] P. Pandey, A. Conway, J. Durie, M. A. Bender, M. Farach-Colton, and R. Johnson. “Vector Quotient Filters: Overcoming the Time/Space Trade-Off in Filter Design.” In: SIGMOD ’21: Proceedings of the 2021 International Conference on Management of Data, Virtual Event, China. New York, NY, United States: ACM SIGMOD, June 2021. url: https://doi.org/10.1145/3448016.3452841 (visited on 10/21/2022).

[2] P. Pandey, A. Conway, and R. Johnson. vqf. Mar. 2021. url: https://github.com/splatlab/vqf (visited on 11/07/2022).

[3] T. Schmidt. partitioned-filters. Aug. 2022. url: https://github.com/tum-db/partitioned-filters (visited on 08/21/2022).

[4] T. Schmidt, M. Bandle, and J. Giceva. “A four-dimensional Analysis of Partitioned Approximate Filters.” In: Proceedings of the VLDB Endowment, Volume 14, Issue 11. VLDB Endowment, ISSN:2150-8097, July 2021. url: https://doi.org/10.14778/3476249.3476286 (visited on 10/29/2022).

[5] K. A. Ross. “Efficient Hash Probes on Modern Processors.” In: IBM Research Report RC24100 (Nov. 2006). url: https://dominoweb.draco.res.ibm.com/reports/rc24100.pdf (visited on 11/13/2022).