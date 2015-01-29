//
//  Random.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___Random_h
#define activemoleculesC___Random_h

namespace nologin {
    namespace utils {
        
        // Random number generator
        struct RNG {
            unsigned int MT[624];
            int index;
            
            RNG(int seed = 1) {
                init(seed);
            }
            
            void init(int seed = 1) {
                MT[0] = seed;
                for(int i = 1; i < 624; i++) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
                index = 0;
            }
            
            void generate() {
                const unsigned int MULT[] = {0, 2567483615UL};
                for(int i = 0; i < 227; i++) {
                    unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
                    MT[i] = MT[i+397] ^ (y >> 1);
                    MT[i] ^= MULT[y&1];
                }
                for(int i = 227; i < 623; i++) {
                    unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
                    MT[i] = MT[i-227] ^ (y >> 1);
                    MT[i] ^= MULT[y&1];
                }
                unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
                MT[623] = MT[623-227] ^ (y >> 1);
                MT[623] ^= MULT[y&1];
            }
            
            unsigned int rand() {
                if (index == 0) {
                    generate();
                }
                
                unsigned int y = MT[index];
                y ^= y >> 11;
                y ^= y << 7  & 2636928640UL;
                y ^= y << 15 & 4022730752UL;
                y ^= y >> 18;
                index = index == 623 ? 0 : index + 1;
                return y;
            }
            
            inline int next() {
                return rand();
            }
            
            inline int next(int x) {
                return rand() % x;
            }
            
            inline int next(int a, int b) {
                return a + (rand() % (b - a));
            }
            
            inline double nextDouble() {
                return (rand() + 0.5) * (1.0 / 4294967296.0);
            }
        };
        
        inline int uniform_random_index()
        {
#ifdef CUSTOM_UNIFORM_RANDOM_INDEX_FUNCTION
            return CUSTOM_UNIFORM_RANDOM_INDEX_FUNCTION % std::numeric_limits<IndexType>::max();
#else
            return std::rand();
#endif
        }
        
        inline int uniform_random_index_bounded(int upper)
        {
            return uniform_random_index() % upper;
        }
        
        /**
         * Returns uniformly distributed random value in range [0 : 1)
         */
        inline double uniform_random()
        {
#ifdef CUSTOM_UNIFORM_RANDOM_FUNCTION
            return CUSTOM_UNIFORM_RANDOM_FUNCTION;
#else
            return std::rand() / ((double)RAND_MAX + 1);
#endif
        }
        
        inline double gaussian_random()
        {
#ifdef CUSTOM_GAUSSIAN_RANDOM_FUNCTION
            return CUSTOM_GAUSSIAN_RANDOM_FUNCTION;
#else
            double x, y, radius;
            do {
                x = 2 * (std::rand() / ((double)RAND_MAX + 1)) - 1;
                y = 2 * (std::rand() / ((double)RAND_MAX + 1)) - 1;
                radius = (x * x) + (y * y);
            } while ((radius >= 1.0) || (radius == 0.0));
            radius = std::sqrt(-2 * std::log(radius) / radius);
            x *= radius;
            y *= radius;
            return x;
#endif
        }
        
        template <class RAI>
        inline void random_shuffle(RAI first, RAI last)
        {
            std::random_shuffle(first, last, uniform_random_index_bounded);
        }
        
        /**
         * Returns num_random samples of random values from range between begin and end.
         * I.e. it is effectively shuffle vector of data between begin and end iterators within num_random samples.
         */
        template<class bidiiter>
        inline bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
            size_t left = std::distance(begin, end);
            while (num_random--) {
                bidiiter r = begin;
                std::advance(r, rand()%left);
                std::swap(*begin, *r);
                ++begin;
                --left;
            }
            return begin;
        }
        
        class RandomSample {
            // generate "m_number" of data samples with the value within the range [0, m_max].
            int m_max;
            int m_number;
            
        public:
            RandomSample(int max, int number) : m_max(max), m_number(number) {}
            
            inline VI get_sample_index() {
                // fill vector with indices
                VI re_res(m_max);
                for (int i = 0; i < m_max; ++i)
                    re_res[i] = i;
                
                // suffle
                random_unique(re_res.begin(), re_res.end(), m_number);
                
                // resize vector
                re_res.resize(m_number);
                VI(re_res).swap(re_res);
                
                return re_res;
            }
        };
        
    }
}
#endif
