#include <iostream>
#include<bitset>
#include <cmath>
#include <iomanip>

using namespace std;

template <bool bSign, size_t bExp, size_t bMts>
struct Float {
    static_assert(bExp +bMts+ bSign<= 64,
    "Invalid template arguments: total number of bits must be <= 64");
    uint64_t bin = 0;
 public:
    Float(){}
    explicit Float(double n) {
        uint64_t* p = (uint64_t*)(void*)&n;
        bool sgn = n < 0;
        uint64_t exp, mts;
        uint64_t dExp = (*p & -0x8000000000000001) >> 52;
        uint64_t dMts = (*p << 12) >> 12;
        double M = (double)dMts/((uint64_t)1<<52);
        if(dExp == 0x7ff){
            exp = ((uint64_t)1 << bExp) - 1;
            if(dMts == 0){
                mts = 0;
            }else{
                mts = (uint64_t)1 | ((dMts>>22) << (bMts-1));
            }
        }else if(dExp == 0){
            if(dMts != 0)
                throw runtime_error("Exponent overflow converting double to Float (subnormal doubles)");
            exp = mts = 0;
        }else{
            int maxExp = ((uint64_t)1 << (bExp-1)) - 1;
            int minExp = -(maxExp - 1);
            int dExpRel = dExp - 1023;
            if(dExpRel < minExp || dExpRel > maxExp)
                throw runtime_error("Exponent overflow converting double to Float");
            exp = dExpRel + maxExp;
            mts = (dMts >> (52-bMts)) + ((dMts >> (52-bMts-1)) & 1);
        }
        set(sgn, exp, mts);
    }

    Float(bool sign, uint64_t exp, uint64_t mts){
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
    static string getBitString(uint64_t bin, uint8_t len=64) {
        string bs("b" + string(len, '0'));
        int i = len;
        while (bin) {
            bs[i--] += (bin & 1);
            bin >>= 1;
        }
        return bs;
    }
    string getBitString() const {
        return getBitString(bin);
    }
    void set(bool sign, uint64_t exp, uint64_t mts){
        if(sign && !bSign)
                throw runtime_error("Float type is unsigned but negative value assigment was attempted");
        if(exp > (uint64_t)1<<bExp)
                throw runtime_error("Exponent out of range");
        if(mts > (uint64_t)1<<bMts)
                throw runtime_error("Mantissa out of range");
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
        bool sign = getSign();
        int64_t exp = getExp();
        uint64_t mts = getMts();
        double M = 1+(double)mts / ((uint64_t)1<<bMts);
        int64_t maxExp = ((int64_t)1 << (bExp-1))-1; 
        return (sign ? -1 : 1) * pow(2, exp-maxExp) * M;
    }

    typedef Float<bSign, bExp, bMts> SameFloat;

    SameFloat operator+(const SameFloat &other) const {
        // a - pointer to number with larger exponent
        // b - pointer to number wiht smaller exponent
        const SameFloat* a = this;
        const SameFloat* b = &other;
        if (a->getExp() < b->getExp())
            swap(a, b);

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
        if(a->getSign() == b->getSign()){
            
            sgn = a->getSign(); // just copy the sign

            // Sum the mantissas
            newMts = mts1+mts2;

            // If the sum overflows, shift the mantissa one bit 
            // to the right and increment the exponent
            if(newMts & carryMask){
                newMts >>= 1;
                exp1++;
            }

        } else { // If signs differ it's more complicated
            int64_t sgn1 = a->getSign() ? -1 : 1;
            int64_t sgn2 = b->getSign() ? -1 : 1;

            if(mts1 < mts2){
                swap(mts1, mts2);
                swap(sgn1, sgn2);
            }

            newMts = mts2-mts1;
            sgn = sgn1; // Copy the sign of the larger number

            if(sgn){ // If the difference is negative
                newMts = -newMts; // get the magnitute

                // Normalize the new mantissa by shifting left and
                // decrementing the exponent until a one is found at the
                // bMts+2 bit, and then cut the one off.
                // if the exponent becomes 0 before a one is found
                // idk i'm fucked i'll just save it as -inf idk what to do tbh

                while(!(newMts & implicitMask) && exp1 > 0){
                    newMts <<= 1;
                    exp1--;
                }

                // If the exponent has been decremented to zero but
                // the mantissa is still not normalized, idk what to do rly
                // i'll just do nothing for now
                if(!(newMts & implicitMask)) {
                    // Nothing
                    cout << "Undefined behaviour, search for this string in code to find error: 9832ihybgsag :)" << endl;
                }
            }
        } 
        // cut off the implicit bit, not handling subnormals
        newMts &= implicitMask-1;

        // Use the rounding bit to round;
        newMts = (newMts >> 1) + (newMts & 1);

        return SameFloat(sgn, exp1, newMts);

    }
};

int main(){


    Float<1, 8, 23> n1(0.6);
    Float<1, 8, 23> n2(0.1);
    auto sum = n1 + n2;

    cout << setprecision(10) << sum.getDouble() << endl;

/*     uint64_t raw = sum.getRaw();
    double *p = reinterpret_cast<double*>(&raw);
    cout << *p << endl; */
}