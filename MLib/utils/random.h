//
//  Random.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___Random_h
#define activemoleculesC___Random_h

// Random number generator
struct RNG {
    unsigned int MT[624];
    int index;
    
    RNG(int seed = 1) {
        init(seed);
    }
    
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
    
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
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

inline IndexType uniform_random_index()
{
#ifdef CUSTOM_UNIFORM_RANDOM_INDEX_FUNCTION
    return CUSTOM_UNIFORM_RANDOM_INDEX_FUNCTION % std::numeric_limits<IndexType>::max();
#else
    return std::rand();
#endif
}

inline IndexType uniform_random_index_bounded(IndexType upper)
{
    return uniform_random_index() % upper;
}

inline ScalarType uniform_random()
{
#ifdef CUSTOM_UNIFORM_RANDOM_FUNCTION
    return CUSTOM_UNIFORM_RANDOM_FUNCTION;
#else
    return std::rand()/((double)RAND_MAX+1);
#endif
}

inline ScalarType gaussian_random()
{
#ifdef CUSTOM_GAUSSIAN_RANDOM_FUNCTION
    return CUSTOM_GAUSSIAN_RANDOM_FUNCTION;
#else
    ScalarType x, y, radius;
    do {
        x = 2*(std::rand()/((double)RAND_MAX+1)) - 1;
        y = 2*(std::rand()/((double)RAND_MAX+1)) - 1;
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
    std::random_shuffle(first,last,uniform_random_index_bounded);
}

class RandomSample {
    // class members
    // generate "m_number" of data with the value within the range [0, m_max].
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


#endif
