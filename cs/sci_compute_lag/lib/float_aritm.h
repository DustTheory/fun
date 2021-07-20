#include <cmath>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <exception>
#include <iostream>

#include "./essential.h"

#ifndef SCI_COMPUTE_LAG_LIB_FLOAT_ARITM_H_
#define SCI_COMPUTE_LAG_LIB_FLOAT_ARITM_H_

template <bool bSign, size_t bExp, size_t bMts>
struct Float;

typedef Float<1, 8, 23> Single;
typedef Float<1, 11, 52> Double;

template <bool bSign, size_t bExp, size_t bMts>
struct Float {
    static_assert(bExp +bMts+ bSign<= 64,
    "Invalid template arguments: total number of bits must be <= 64");
    uint64_t bin = 0;

 public:
    Float() {}
    explicit Float(double n) {
        uint64_t* p = reinterpret_cast<uint64_t*>(&n);
        bool sgn = n < 0;
        uint64_t exp, mts;
        uint64_t dExp = (*p & -0x8000000000000001) >> 52;
        uint64_t dMts = (*p << 12) >> 12;
        double M = static_cast<double>(dMts)/((uint64_t)1<<52);
        if (dExp == 0x7ff) {
            exp = ((uint64_t)1 << bExp) - 1;
            if (dMts == 0) {
                mts = 0;
            } else {
                mts = (uint64_t)1 | ((dMts>>22) << (bMts-1));
            }
        } else if (dExp == 0) {
            if (dMts != 0)
                throw std::runtime_error("Exponent overflow converting double to Float (subnormal doubles)");
            exp = mts = 0;
        } else {
            int dExpRel = (int)dExp - 1023;
            if (dExpRel < minExp())
                exp = 0, mts = 0;
            else if (dExpRel > maxExp())
                exp = ((uint64_t)1 << bExp) - 1, mts = 0;  // inf
            else
                exp = dExpRel + maxExp(),
                mts = (dMts >> (52-bMts)) + ((dMts >> (52-bMts-1)) & 1);
        }
        set(sgn, exp, mts);
    }

    Float(bool sign, uint64_t exp, uint64_t mts) {
        set(sign, exp, mts);
    }

    uint64_t getExp() const {
        return ((bin & (~((uint64_t)1 << (bMts + bExp)))) >> bMts);
    }

    uint64_t getMts() const {
        return (bin << (64-bMts)) >> (64-bMts);
    }

    bool getSign() const {
        return bin >> (bMts+bExp);
    }

    static std::string getBitStringStatic(uint64_t bin, uint8_t len = 64) {
        std::string bs("b" + std::string(len, '0'));
        int i = len;
        while (bin) {
            bs[i--] += (bin & 1);
            bin >>= 1;
        }
        return bs;
    }

    std::string getBitString() const {
        return getBitStringStatic(bin);
    }

    void set(bool sign, uint64_t exp, uint64_t mts){
        if (sign && !bSign)
                throw std::runtime_error("Float type is unsigned but negative value assigment was attempted");
        if (exp > (uint64_t)1 << bExp)
                throw std::runtime_error("Exponent out of range");
        if (mts > (uint64_t)1 << bMts)
                throw std::runtime_error("Mantissa out of range");
        if (bSign) {
            bin |= sign;
            bin <<= bExp;
        }
        bin |= exp;
        bin <<= bMts;
        bin |= mts;
    }

    uint64_t getRaw() const {
        return bin;
    }

    double getDouble() const {

        if (bSign == 1 && bExp == 11 && bMts == 52) {
            const double* p = reinterpret_cast<const double*>(&bin);
            return *p;
        }

        bool sign = getSign();
        int64_t exp = getExp();
        uint64_t mts = getMts();
        double M = 1 + static_cast<double>(mts) / ((uint64_t)1 << bMts);
        return (sign ? -1 : 1) * pow(2, exp-maxExp()) * M;
    }

    float getFloat() const {
         if (bSign == 1 && bExp == 8 && bMts == 23) {
            const float* p = reinterpret_cast<const float*>(&bin);
            return *p;
        }
        return getDouble();
    }

    uint64_t maxBiasedExp() const {
        return ((uint64_t)1 << bExp)-1;
    }

    int64_t maxExp() const {
        return ((uint64_t)1 << (bExp-1)) - 1;
    }

    int64_t minExp() const {
        return -(maxExp() - 1);
    }

    typedef Float<bSign, bExp, bMts> SameFloat;

    SameFloat operator+(const SameFloat &other) const {

        // a - pointer to number with larger exponent
        // b - pointer to number wiht smaller exponent
        const SameFloat* a = this;
        const SameFloat* b = &other;
        if (a->getExp() < b->getExp())
            std::swap(a, b);

        // get exponents of both number
        uint64_t exp1 = a->getExp();
        uint64_t exp2 = b->getExp();

        // get the mantissas too but shift left by one to leave room
        // for one more extra bit which will be used for rounding
        int64_t mts1 = a->getMts() << 1;
        int64_t mts2 = b->getMts() << 1;


        // Add on the implicit bits, taking care of subnormals
        mts1 |= (uint64_t)(exp1 != 0 || mts1 == 0) << (bMts+1);
        mts2 |= (uint64_t)(exp2 != 0 || mts2 == 0) << (bMts+1);

        // Make exponents the same
        mts2 >>= (exp1 - exp2);
        exp2 = exp1;

        int64_t implicitMask = (int64_t)1 << (bMts+1);
        int64_t carryMask = (int64_t)1 << (bMts+2);

        uint64_t newMts, newExp, sgn;

        // If the signs are the same, just do simple addition and copy the sign
        if (a->getSign() == b->getSign()) {
            
            sgn = a->getSign();  // just copy the sign

            // Sum the mantissas
            newMts = mts1+mts2;

            // If the sum overflows, shift the mantissa one bit
            // to the right and increment the exponent
            if (newMts & carryMask) {
                newMts >>= 1;
                exp1++;
            }

        } else {  // If signs differ it's more complicated
            int64_t sgn1 = a->getSign() ? -1 : 1;
            int64_t sgn2 = b->getSign() ? -1 : 1;

            if (mts1 < mts2) {
                std::swap(mts1, mts2);
                std::swap(sgn1, sgn2);
            }

            newMts = mts1-mts2;
            sgn = sgn1 == -1;   // Copy the sign of the larger number

            // Normalize the new mantissa by shifting left and
            // decrementing the exponent until a one is found at the
            // bMts+2 bit, and then cut the one off.
            // if the exponent becomes 0 before a one is found
            // idk i'm fucked i'll just save it as -inf idk what to do tbh

            while (!(newMts & implicitMask) && exp1 > 0) {
                newMts <<= 1;
                exp1--;
            }

            // If the exponent has been decremented to zero but
            // the mantissa is still not normalized, idk what to do rly
            // i'll just do nothing for now
            if (!(newMts & implicitMask)) 
                return SameFloat(0, 0, 0); // Guess I'll just return 0 :|
            
        } 
        // cut off the implicit bit, not handling subnormals
        newMts &= implicitMask-1;

        // Use the rounding bit to round;
        newMts = (newMts >> 1) + (newMts & 1);

        return SameFloat(sgn, exp1, newMts);

    }

    SameFloat operator*(const SameFloat &other) const {
        uint64_t exp1 = getExp();
        uint64_t exp2 = other.getExp();

        // get the mantissas too but shift left by one to leave room
        // for one more extra bit which will be used for rounding
        uint64_t mts1 = getMts();
        uint64_t mts2 = other.getMts();


        // Add on the implicit bits, taking care of subnormals
        mts1 |= (uint64_t)(exp1 != 0 || mts1 == 0) << bMts;
        mts2 |= (uint64_t)(exp2 != 0 || mts2 == 0) << bMts;

        int64_t implicitMask = (int64_t)1 << bMts;

        // Remove and count trailing zeros
        int tz1 = 0, tz2 = 0;
        while (!(mts1 & 1)) mts1 >>= 1, tz1++;
        while (!(mts2 & 1)) mts2 >>= 1, tz2++;

        uint64_t newMts = mts1*mts2;
        uint64_t newExp = exp1 + exp2 > maxExp() ? exp1 + exp2 - maxExp() : 0;
        if (newExp > maxBiasedExp() - 1)
            newExp = maxBiasedExp(), newMts = 0;

        if (((newMts >> (bMts * 2 - tz1 - tz2 + 1)) & 1))  // handle carry
            newExp++;

        while (newMts > (implicitMask << 1))
            newMts = (newMts>>1);

        if (newMts != 0)
            while (newMts < implicitMask)
                newMts <<= 1;

        // cut off the implicit bit, not handling subnormals
        newMts &= implicitMask-1;

        return SameFloat(getSign()^other.getSign(), newExp, newMts);
    }

    SameFloat operator/(const SameFloat& other) const {
        uint64_t exp1 = getExp();
        uint64_t exp2 = other.getExp();

        // get the mantissas too but shift left by one to leave room
        // for one more extra bit which will be used for rounding
        uint64_t mts1 = getMts();
        uint64_t mts2 = other.getMts();


        // Add on the implicit bits, taking care of subnormals
        mts1 |= (uint64_t)(exp1 != 0 || mts1 == 0) << bMts;
        mts2 |= (uint64_t)(exp2 != 0 || mts2 == 0) << bMts;

        int64_t implicitMask = (int64_t)1 << bMts;
        uint64_t lastBitMask = (int64_t)1 << 63;
        uint64_t newExp = maxExp() + exp1 - exp2 - (mts1 < mts2);

        // shift left until implicit bit is msb of uint64_t
        while (newExp > 1 && !(mts1 & lastBitMask))
            mts1 <<= 1;

        uint64_t newMts = mts1/mts2;
        while (newMts > (implicitMask << 1))
            newMts = (newMts>>1) + (newMts & 1);

        if (newMts != 0)
            while (newMts < implicitMask)
                newMts <<= 1;

        // cut off the implicit bit, not handling subnormals
        newMts &= implicitMask-1;

        return SameFloat(getSign()^other.getSign(), newExp, newMts);
    }

    SameFloat operator-() const {
        return SameFloat(!getSign(), getExp(), getMts());
    }

    SameFloat operator-(const SameFloat &other) const {
        return *this + -other;
    }

    // Approximate a/b using the Newton-Raphson method
    // Often doesn't work
    SameFloat approxDiv(const SameFloat &other, int iter = 10) const {
        return *this * approxReciprocal(other, iter);
    }

    static SameFloat approxReciprocal(SameFloat b, int iter = 10) {
        SameFloat bInv(b);
        if(b.getSign()) bInv = -bInv;
        SameFloat two(2.0);
        while (iter--){
            bInv = bInv * (two - b*bInv);
        }
        return bInv;
    }

    SameFloat cheatDiv(const SameFloat &other) const {
        return SameFloat(getDouble()/other.getDouble());
    }

    SameFloat timesPow2(int pow) const {
        return SameFloat(getSign(), getExp() + pow, getMts());
    }

    SameFloat sqrt(int iter = 10) const {
        SameFloat x(0.1);
        SameFloat half(0.5);
        while (iter--)
            x = (x + cheatDiv(x))*half;
        // Curretnly using cheating division,
        // I don't have a good division implememted yet
        // x = (x + *this/x)*half;
        // x = (x + approxDiv(x))*half;
        return SameFloat(x);
    }

    SameFloat invSqrt() const {
        SameFloat x(0.1);
        SameFloat three(3);
        SameFloat half(0.5);
        int iter(10);
        while (iter--)
            x = half*x*(three - (*this)*x*x);
        return x;
    }

    SameFloat sqrt2(int iter = 10) const {
        return approxReciprocal(invSqrt(), iter);
    } 

    bool operator<(const SameFloat &other) const {
        if (getSign() != other.getSign())
            return getSign();
        if (getExp() != other.getExp())
            return getExp() < other.getExp();
        return bin < other.bin;
    }

    bool operator<(double other) const {
        return *this < SameFloat(other);
    }

    SameFloat abs(){
        return SameFloat(0, getExp(), getMts());
    }

};

#endif  // SCI_COMPUTE_LAG_LIB_FLOAT_ARITM_H_
