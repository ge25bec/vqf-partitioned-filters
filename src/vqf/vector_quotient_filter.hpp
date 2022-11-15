#pragma once

#include <iostream>
#include <filter_base.hpp>
#include <vqf/vector_quotient_filter_container.hpp>
#include <vqf/vqf_parameter.hpp>

namespace filters {

    template<typename FilterParameter, size_t _k, typename OptimizationParameter>
    struct Filter<FilterType::VectorQuotientFilter, FilterParameter, _k, OptimizationParameter> {

        using FP = FilterParameter;
        using OP = OptimizationParameter;
        static constexpr size_t k{(_k == 8) ? 8 : 16};
        static constexpr size_t max_n_retries{FP::max_n_retries};
        static constexpr bool supports_add{true};
        static constexpr bool supports_add_partition{false};

        using Vector = simd::Vector<OP::registerSize, OP::simd>;
        using T = typename Vector::T;
        using M = typename Vector::M;
        using Hasher = hash::Hasher<OP::hashingMode, Vector, 3>;
        using Addresser = addresser::Addresser<OP::addressingMode, Vector>;
        using Container = vqf::VectorQuotientFilterContainer<OP::registerSize, OP::simd, OP, k, Hasher, Vector, Addresser>;

        // types used for fallback (always scalar)
        using FallbackVector = simd::Vector<OP::registerSize, parameter::SIMD::Scalar>;

        size_t s;
        size_t n_partitions;
        size_t n_retries;
        Container container;
        task::TaskQueue<OP::multiThreading> queue;

        Filter(size_t s, size_t n_partitions, size_t n_threads, size_t n_tasks_per_level) : s(s),
                                                                                            n_partitions(n_partitions),
                                                                                            n_retries(0),
                                                                                            queue(n_threads,
                                                                                                    n_tasks_per_level) {
        }

        forceinline
        void init(const T *histogram) {
            n_partitions = 1;
            container = std::move(Container(s, histogram));
        }

        forceinline
        void init(size_t length) {
            n_partitions = 1;
            container = std::move(Container(s, length));
        }

        forceinline
        bool contains(const T &value) {
            return contains(value, 0);
        }

        forceinline
        bool contains(const T &value, size_t index) {
            return container.lookup(value);
        }

        forceinline
        bool add(const T &value) {
            return container.insert(value);
        }

        forceinline
        bool add_all(const T *values, size_t length) {
            if constexpr ((OP::multiThreading == filters::parameter::MultiThreading::Enabled) && (OP::simd == filters::parameter::SIMD::Scalar)) {
                const size_t n_threads{queue.get_n_threads()};
                const size_t stepSize{length/n_threads};
                const size_t loopRemainder{length % n_threads};
                std::atomic<size_t> result{0};
                
                //dividing the input array into workloads for threads
                for (size_t i{}; i < n_threads; i++) {
                    queue.add_task([this, &result, values, i, stepSize](size_t) {
                        result += container.insert(values + i*stepSize, stepSize);
                    });
                }

                queue.add_task([this, &result, values, length, loopRemainder](size_t) {
                    result += container.insert(values + length - loopRemainder, loopRemainder);
                });

                queue.execute_tasks();

                return result;
            }
            return container.insert(values, length) == length;
        }

        forceinline
        bool add_partition(const T *values, size_t length, size_t index) {
            return false;
        }

        bool _construct_fallback(const FallbackVector &offset, const FallbackVector &histogram, const T *values,
                                 const size_t index) {
            return false;
        }

        bool _construct(const Vector &offset, const Vector &histogram, const T *values, const size_t index) {
            return false;
        }

        bool construct(T *values, size_t length) {
            bool success{};

            n_retries = 0;
            for (size_t i{}; (i < max_n_retries) && !success; i++) {
                init(length);
                success = add_all(values, length);

                if (!success)
                    n_retries++;
            }

            return success;
        }

        size_t count(T *values, size_t length) {
            if constexpr (OP::multiThreading == parameter::MultiThreading::Enabled) {
                const size_t n_threads{queue.get_n_threads()};
                const size_t stepSize{length/n_threads};
                const size_t loopRemainder{length % n_threads};
                std::atomic<size_t> result{0};

                //dividing the input array into workloads for threads
                for (size_t i{}; i < n_threads; i++) {
                    queue.add_task([this, &result, values, i, stepSize](size_t) {
                        result += container.lookup(values + i*stepSize, stepSize);
                    });
                }

                queue.add_task([this, &result, values, length, loopRemainder](size_t) {
                    result += container.lookup(values + length - loopRemainder, loopRemainder);
                });

                queue.execute_tasks();

                return result;
            } else {
                return container.lookup(values, length);
            }
        }

        size_t size() {
            return container.size();
        }

        size_t avg_size() {
            return size();
        }

        size_t retries() {
            return n_retries;
        }

        std::string to_string() {
            std::string s = "\n{\n";
            s += "\t\"k\": " + std::to_string(k) + ",\n";
            s += "\t\"size\": " + std::to_string(size() * 8) + " bits,\n";
            s += "\t\"n_partitions\": " + std::to_string(n_partitions) + ",\n";
            s += "\t\"n_retries\": " + std::to_string(n_retries) + ",\n";
            s += "\t\"filter_params\": " + FP::to_string() + ",\n";
            s += "\t\"optimization_params\": " + OP::to_string() + "\n";
            s += "}\n";

            return s;
        }
    };

} // filters