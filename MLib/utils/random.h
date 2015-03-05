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
        
        
        //
        // --------------------------------------------------------------------------------------------
        //
        
        // This is the Mersenne Twister random number generator MT19937, which
        // generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
        // starting from any odd seed in 0..(2^32 - 1).  This version is a recode
        // by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
        // Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
        // July-August 1997).
        //
        // Effectiveness of the recoding (on Goedel2.math.washington.edu, a DEC Alpha
        // running OSF/1) using GCC -O3 as a compiler: before recoding: 51.6 sec. to
        // generate 300 million random numbers; after recoding: 24.0 sec. for the same
        // (i.e., 46.5% of original time), so speed is now about 12.5 million random
        // number generations per second on this machine.
        //
        // According to the URL <http://www.math.keio.ac.jp/~matumoto/emt.html>
        // (and paraphrasing a bit in places), the Mersenne Twister is ``designed
        // with consideration of the flaws of various existing generators,'' has
        // a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
        // equidistributed, and ``has passed many stringent tests, including the
        // die-hard test of G. Marsaglia and the load test of P. Hellekalek and
        // S. Wegenkittl.''  It is efficient in memory usage (typically using 2506
        // to 5012 bytes of static data, depending on data type sizes, and the code
        // is quite short as well).  It generates random numbers in batches of 624
        // at a time, so the caching and pipelining of modern systems is exploited.
        // It is also divide- and mod-free.
        //
        // This library is free software; you can redistribute it and/or modify it
        // under the terms of the GNU Library General Public License as published by
        // the Free Software Foundation (either version 2 of the License or, at your
        // option, any later version).  This library is distributed in the hope that
        // it will be useful, but WITHOUT ANY WARRANTY, without even the implied
        // warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
        // the GNU Library General Public License for more details.  You should have
        // received a copy of the GNU Library General Public License along with this
        // library; if not, write to the Free Software Foundation, Inc., 59 Temple
        // Place, Suite 330, Boston, MA 02111-1307, USA.
        //
        // The code as Shawn received it included the following notice:
        //
        //   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When
        //   you use this, send an e-mail to <matumoto@math.keio.ac.jp> with
        //   an appropriate reference to your work.
        //
        // It would be nice to CC: <Cokus@math.washington.edu> when you write.
        //
        struct MT_RNG {
            //
            // uint32 must be an unsigned integer type capable of holding at least 32
            // bits; exactly 32 should be fastest, but 64 is better on an Alpha with
            // GCC at -O3 optimization so try your options and see whats best for you
            //
            typedef unsigned int uint32;
            
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v
            
            const uint32 K = 0x9908B0DFU; // a magic constant
            const uint32 N = 624; // length of state vector
            const uint32 M = 397; // a period parameter
            
            uint32   state[624 + 1];     // state vector + 1 extra to not violate ANSI C
            uint32   *next;          // next random value is computed from here
            int      left = -1;      // can *next++ this many times before reloading
            
            
            void seedMT(uint32 seed) {
                //
                // We initialize state[0..(N-1)] via the generator
                //
                //   x_new = (69069 * x_old) mod 2^32
                //
                // from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuths
                // _The Art of Computer Programming_, Volume 2, 3rd ed.
                //
                // Notes (SJC): I do not know what the initial state requirements
                // of the Mersenne Twister are, but it seems this seeding generator
                // could be better.  It achieves the maximum period for its modulus
                // (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
                // x_initial can be even, you have sequences like 0, 0, 0, ...;
                // 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
                // 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
                //
                // Even if x_initial is odd, if x_initial is 1 mod 4 then
                //
                //   the          lowest bit of x is always 1,
                //   the  next-to-lowest bit of x is always 0,
                //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
                //   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
                //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
                //    ...
                //
                // and if x_initial is 3 mod 4 then
                //
                //   the          lowest bit of x is always 1,
                //   the  next-to-lowest bit of x is always 1,
                //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
                //   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
                //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
                //    ...
                //
                // The generators potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
                // 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
                // also does well in the dimension 2..5 spectral tests, but it could be
                // better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
                //
                // Note that the random number user does not see the values generated
                // here directly since reloadMT() will always munge them first, so maybe
                // none of all of this matters.  In fact, the seed values made here could
                // even be extra-special desirable if the Mersenne Twister theory says
                // so-- thats why the only change I made is to restrict to odd seeds.
                //
                
                uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
                int    j;
                
                for (left = 0, *s++ = x, j = N; --j;
                     *s++ = (x*=69069U) & 0xFFFFFFFFU);
            }
            
            
            uint32 reloadMT(void) {
                uint32 *p0 = state, *p2 = state + 2, *pM = state + M, s0, s1;
                int    j;
                
                if (left < -1)
                    seedMT(4357U);
                
                left = N - 1, next = state + 1;
                
                for (s0 = state[0], s1 = state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++)
                    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
                
                for (pM = state, j = M; --j; s0 = s1, s1 = *p2++)
                    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
                
                s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
                s1 ^= (s1 >> 11);
                s1 ^= (s1 <<  7) & 0x9D2C5680U;
                s1 ^= (s1 << 15) & 0xEFC60000U;
                return(s1 ^ (s1 >> 18));
            }
            
            
            uint32 randomMT() {
                uint32 y;
                
                if(--left < 0)
                    return(reloadMT());
                
                y  = *next++;
                y ^= (y >> 11);
                y ^= (y <<  7) & 0x9D2C5680U;
                y ^= (y << 15) & 0xEFC60000U;
                y ^= (y >> 18);
                return(y);
            }
            
#define MAX_UINT_COKUS 4294967295  //basically 2^32-1
            double unif_rand(){
                return (((double)randomMT())/((double)MAX_UINT_COKUS));
            }
        };
        
        //
        // --------------------------------------------------------------------------------------------
        //

        
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
