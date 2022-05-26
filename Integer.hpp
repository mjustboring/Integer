/*
    ################################################
    ##  THE BIG INTEGER LIBRARY | BY : SP MAURYA  ##
    ################################################

    DATE STARTED : 27th April 2021
    DATE FINISHED : 29th July 2021

    TWEAKS : 26th May 2022
*/

#ifndef BIG_INTEGER_LIBRARY_SP_MAURYA
#define BIG_INTEGER_LIBRARY_SP_MAURYA

#include <stack>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <unordered_map>

class Integer;
class Fin;
class Fout;
struct Div;

constexpr unsigned short digs = 8U;
constexpr unsigned short ksmt = 64U;
constexpr unsigned short bits = 0x10U;
constexpr unsigned short bmax = 0xFFFFU;
constexpr unsigned       base = 0x10000U;
constexpr unsigned       dmax = 100000000U;

constexpr unsigned digc(unsigned b) { return b * 0.30103 + 1; }
constexpr unsigned bitc(unsigned d) { return d * 3.3219281 + 1; }

constexpr unsigned p10[] = {
    1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000
};

std::unordered_map <unsigned, std::vector<unsigned>> dec_exp;
std::unordered_map <unsigned, std::vector<unsigned short>> bin_exp;

std::chrono::time_point<std::chrono::high_resolution_clock> start_tp;
std::chrono::time_point<std::chrono::high_resolution_clock> end_tp;

double get_time() noexcept {
    return std::chrono::duration<double>(end_tp - start_tp).count();
}

void start_timer() noexcept {
    start_tp = std::chrono::high_resolution_clock::now();
}

double end_timer() noexcept {
    end_tp = std::chrono::high_resolution_clock::now();
    return get_time();
}

constexpr char cset[] = {
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 65, 66, 67, 68, 69, 70, 71, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
    97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, -32, -31
};

constexpr int getv(const char c) noexcept {

    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
    if (c >= 'a' && c <= 'z') return c - 'a' + 36;

    if (c == -32) return 62;
    if (c == -31) return 63;

    return -1;
}

void seed_rand() noexcept {
    static bool seeded;
    if (!seeded) {
        srand((std::chrono::high_resolution_clock::now().time_since_epoch().count() >> 10) & bmax);
        seeded = true;
    }
}

long long unsigned gcd (
    long long unsigned a, long long unsigned b
) noexcept {
    return b == 0 ? a : gcd(b, a % b);
}

std::string ntos (long long unsigned n) {
    std::stringstream ss;
    std::string s;
    ss << n;
    ss >> s;
    return s;
}

long long unsigned ston (const std::string& s) {
    std::stringstream ss;
    long long unsigned n;
    ss << s;
    ss >> n;
    return n;
}

unsigned pop_count (unsigned long long n) noexcept {
    unsigned count = 0;
    while (n)
        count++, n >>= 1;
    return count;
}

std::pair<unsigned, unsigned> get_frac (
    const long double& num, const long double& eps = 1e-9
) {
    std::pair<unsigned, unsigned> lb {0, 1}, rb {1, 0}, nb;
    long double alpha, gamma;
    unsigned skip;

    alpha = gamma = std::abs(num);
    do {
        skip = gamma;
        nb.first = lb.first + skip * rb.first;
        nb.second = lb.second + skip * rb.second;
        lb = rb, rb = nb;
        gamma = 1 / (gamma - skip);
    } while (
        std::abs(alpha - (long double) nb.first / nb.second) > eps
    );
    return nb;
}

class Integer {

    // DATA
    std::vector<unsigned short> data;
    bool sign;

    // DATA INITIALIZERS
    Integer &init(unsigned long long, bool = false);
    Integer &init(const std::string &);

    // STATIC CLASS UTILITY FUNCTIONS
    static void trim (std::vector<unsigned short>&) noexcept;
    static void trim (std::vector<unsigned>&) noexcept;

    static void _1s_comp (std::vector<unsigned short>&) noexcept;
    static void _2s_comp (std::vector<unsigned short>&) noexcept;

    static bool vec_equal (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&
    ) noexcept;

    static bool vec_less (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&
    ) noexcept;

    static void left_shift (std::vector<unsigned short>&, unsigned);
    static void right_shift (std::vector<unsigned short>&, unsigned);

    static void add (
        std::vector<unsigned short>&, const std::vector<unsigned short>&,
        unsigned, unsigned, unsigned
    );

    static void add (
        std::vector<unsigned>&, const std::vector<unsigned>&,
        unsigned, unsigned, unsigned, const unsigned
    );

    static void sub (
        std::vector<unsigned short>&, const std::vector<unsigned short>&,
        unsigned, unsigned, unsigned
    );

    static void sub (
        std::vector<unsigned>&, const std::vector<unsigned>&,
        unsigned, unsigned, unsigned, const unsigned
    );

    static std::vector<unsigned short> mul (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&,
        unsigned, unsigned, unsigned, unsigned
    );

    static void mul (
        std::vector<unsigned short>&, unsigned
    );

    static std::vector<unsigned> mul (
        const std::vector<unsigned>&, const std::vector<unsigned>&,
        unsigned, unsigned, unsigned, unsigned, const unsigned
    );

    static std::vector<unsigned short> sqr (
        const std::vector<unsigned short>&, unsigned, unsigned
    );

    static std::vector<unsigned> sqr (
        const std::vector<unsigned>&, unsigned, unsigned, const unsigned
    );

    static std::vector<unsigned short> ksm_sqr (
        const std::vector<unsigned short>&, unsigned, unsigned
    );

    static std::vector<unsigned short> ksm_sqr (
        const std::vector<unsigned short>&
    );

    static std::vector<unsigned> ksm_sqr (
        const std::vector<unsigned>&, unsigned, unsigned,
        const unsigned
    );

    static std::vector<unsigned> ksm_sqr (
        const std::vector<unsigned>&, const unsigned
    );

    static std::vector<unsigned short> ksm (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&,
        unsigned, unsigned, unsigned
    );

    static std::vector<unsigned short> ksm (
        std::vector<unsigned short>, std::vector<unsigned short>,
        unsigned, unsigned
    );

    static std::vector<unsigned short> ksm (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&
    );

    static std::vector<unsigned> ksm (
        const std::vector<unsigned>&, const std::vector<unsigned>&,
        unsigned, unsigned, unsigned, const unsigned
    );

    static std::vector<unsigned> ksm (
        std::vector<unsigned>, std::vector<unsigned>,
        unsigned, unsigned, const unsigned
    );

    static std::vector<unsigned> ksm (
        const std::vector<unsigned>&, const std::vector<unsigned>&,
        const unsigned
    );

    static void div (
        const std::vector<unsigned short>&, unsigned,
        std::vector<unsigned short>&, unsigned&
    );

    static void div (
        const std::vector<unsigned short>&, const std::vector<unsigned short>&,
        std::vector<unsigned short>&, std::vector<unsigned short>&
    );

    static std::vector<unsigned short> pow (
        const std::vector<unsigned short>&, unsigned
    );

    static std::vector<unsigned> pow (
        const std::vector<unsigned>&, unsigned, const unsigned
    );

    static std::vector<unsigned> to (
        const std::vector<unsigned short>&, unsigned, unsigned,
        const unsigned, const std::vector<unsigned>&,
        std::unordered_map<unsigned, std::vector<unsigned>>&
    );

    static std::vector<unsigned> to (
        const std::vector<unsigned short>&, const unsigned = dmax
    );

    static std::vector<unsigned short> to (
        const std::vector<unsigned>&, unsigned, unsigned,
        const unsigned, const std::vector<unsigned short>&,
        std::unordered_map<unsigned, std::vector<unsigned short>>&
    );

    static std::vector<unsigned short> to (
        const std::vector<unsigned>&, const unsigned = dmax
    );

    static std::invalid_argument init_except (
        const std::string&, const unsigned, const std::string&
    );

    static std::vector<unsigned short> from_binp (
        const std::string&, const unsigned, const unsigned
    );

    static std::string to_binp (
        const std::vector<unsigned short>&, bool, const unsigned, bool
    );

    static std::vector<unsigned short> fact(const unsigned, const unsigned);

public:

    // CONSTRUCTORS
    Integer() noexcept;
    Integer(int) noexcept;
    Integer(long long) noexcept;
    Integer(unsigned) noexcept;
    Integer(unsigned long long) noexcept;
    Integer(const char *);
    Integer(const std::string &);

    // COPY CTOR
    Integer(const Integer &) noexcept;
    // MOVE CTOR
    Integer(Integer &&) noexcept;

    // COPY ASSIGNMENT OPERATOR
    Integer &operator=(const Integer &) noexcept;
    // MOVE ASSIGNMENT OPERATOR
    Integer &operator=(Integer &&) noexcept;

    // UTILITY
    unsigned bit_count() const noexcept;
    unsigned size() const noexcept;
    void swap(Integer &) noexcept;

    // ARITHMETIC OPERATORS
    Integer &operator+=(const Integer &);
    friend Integer operator+(const Integer &, const Integer &);

    Integer &operator-=(const Integer &);
    friend Integer operator-(const Integer &, const Integer &);

    Integer operator-() const noexcept;
    Integer &flip_sign() noexcept;

    Integer abs() const;
    Integer &make_abs();

    Integer &operator*=(const Integer &);
    Integer &operator*=(const unsigned &);
    Integer &operator*=(const int &);
    friend Integer operator*(const Integer &, const Integer &);
    friend Integer operator*(const Integer &, const unsigned &);
    friend Integer operator*(const unsigned &, const Integer &);
    friend Integer operator*(const Integer &, const int &);
    friend Integer operator*(const int &, const Integer &);

    Integer &operator/=(const Integer &);
    Integer &operator/=(const unsigned &);
    Integer &operator/=(const int &);
    Integer operator/(const unsigned &) const;
    Integer operator/(const int &) const;
    friend Integer operator/(const Integer &, const Integer &);

    Integer &operator%=(const Integer &);
    Integer &operator%=(const unsigned &);
    Integer &operator%=(const int &);
    Integer operator%(const unsigned &) const;
    Integer operator%(const int &) const;
    friend Integer operator%(const Integer &, const Integer &);

    friend Div div(const Integer &, const Integer &);

    // COMPARATORS
    friend bool operator<(const Integer &, const Integer &) noexcept;
    friend bool operator>(const Integer &, const Integer &) noexcept;
    friend bool operator<=(const Integer &, const Integer &) noexcept;
    friend bool operator>=(const Integer &, const Integer &) noexcept;
    friend bool operator==(const Integer &, const Integer &) noexcept;
    friend bool operator!=(const Integer &, const Integer &) noexcept;

    // INCREMENT & DECREMENT
    Integer &operator++() noexcept;
    Integer &operator--() noexcept;

    Integer operator++(int) noexcept;
    Integer operator--(int) noexcept;

    // BITWISE OPERATORS
    Integer &operator|=(const Integer &);
    Integer &operator&=(const Integer &);
    Integer &operator^=(const Integer &);

    friend Integer operator|(const Integer &, const Integer &);
    friend Integer operator&(const Integer &, const Integer &);
    friend Integer operator^(const Integer &, const Integer &);

    Integer &operator<<=(unsigned );
    Integer &operator>>=(unsigned );

    Integer operator<<(unsigned ) const;
    Integer operator>>(unsigned ) const;

    Integer operator~() const;
    Integer& compliment();

    // INPUT & OUTPUT
    friend std::ostream &operator<<(std::ostream &, const Integer &);
    friend std::istream &operator>>(std::istream &, Integer &);

    // MATH UTILITY
    Integer sqr() const;
    Integer nrt(uint8_t) const;
    Integer sqrt() const;
    Integer pow(uint8_t) const;
    Integer pow(uint8_t, uint8_t) const;
    Integer powf(float) const;
    Integer gcd(const Integer &) const;
    Integer lcm(const Integer &) const;
    bool odd() const noexcept;
    bool even() const noexcept;
    bool operator[](unsigned) const noexcept;
    Integer& set_bit(unsigned) noexcept;
    Integer& unset_bit(unsigned) noexcept;

    long double log(long double) const;
    long double log10() const;
    long double log2() const;
    long double loge() const;

    friend Integer fact(unsigned);
    friend Integer rand(unsigned);
    friend Integer rand_bit(unsigned);
    friend Integer rand_byte(unsigned);
    friend Integer rand_range(unsigned, unsigned);
    friend Integer rand_range(const Integer& , const Integer&);

    // CONVERSION OPERATORS
    explicit operator bool() const noexcept;
    operator std::string() const noexcept;

    // RADIX CONVERSIONS
    std::string str(const unsigned = 10, bool = true) const;
    std::string dec(bool = false) const;
    std::string bin(bool = false) const;
    std::string oct(bool = false) const;
    std::string hex(bool = false) const;

}; // END OF Integer

struct Div {

    Integer quot, rem;

    Div();
    Div(const Integer &, const Integer &) noexcept;
    Div(const Integer &, Integer &&) noexcept;
    Div(Integer &&, const Integer &) noexcept;
    Div(Integer &&, Integer &&) noexcept;

    Div(const Div &) noexcept;
    Div(Div &&) noexcept;

    Div &operator=(const Div &) noexcept;
    Div &operator=(Div &&) noexcept;

    friend std::ostream &operator<<(std::ostream &, const Div &);

}; // END OF Div

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Declarations to friend functions

Integer operator+(const Integer &, const Integer &);
Integer operator-(const Integer &, const Integer &);

Integer operator*(const Integer &, const Integer &);
Integer operator*(const Integer &, const unsigned &);
Integer operator*(const unsigned &, const Integer &);
Integer operator*(const Integer &, const int &);
Integer operator*(const int &, const Integer &);

Integer operator%(const Integer &, const Integer &);

Div div(const Integer &, const Integer &);

bool operator<(const Integer &, const Integer &) noexcept;
bool operator>(const Integer &, const Integer &) noexcept;
bool operator<=(const Integer &, const Integer &) noexcept;
bool operator>=(const Integer &, const Integer &) noexcept;
bool operator==(const Integer &, const Integer &) noexcept;
bool operator!=(const Integer &, const Integer &) noexcept;

Integer operator|(const Integer &, const Integer &);
Integer operator&(const Integer &, const Integer &);
Integer operator^(const Integer &, const Integer &);

std::ostream &operator<<(std::ostream &, const Integer &);
std::istream &operator>>(std::istream &, Integer &);

Integer fact(unsigned);
Integer rand(unsigned);
Integer rand_bit(unsigned);
Integer rand_byte(unsigned);
Integer rand_range(unsigned, unsigned);
Integer rand_range(const Integer& , const Integer&);

std::ostream &operator<<(std::ostream &, const Div &);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Integer DEFINITIONS

Integer &Integer::init(unsigned long long n, bool sign) {

    data.clear();
    this->sign = sign;

    while (n) {
        data.emplace_back(n);
        n >>= bits;
    }

    return *this;
}

Integer &Integer::init(const std::string &s) {

    data.clear(), sign = false;

    if (s.empty() || s == "0" || s == "-0")
        return *this;

    unsigned j = 0, i = s.size() - 1, base = 10, v;

    if (s[j] == '-')
        ++j, sign = true;

    if (j < i && s[j] == '0') {

        ++j, base = 0;

        while (j <= i && s[j] != '.') {

            v = getv(s[j]);

            if (v < 0 || v > 9)
                throw init_except("NON-INTEGER BASE", j, s);

            base = base * 10 + v;

            if (base > dmax)
                throw init_except("VERY LARGE BASE", j, s);

            ++j;
        }

        while (j <= i && s[j] == '.') ++j;
    }

    if (base < 2)
        throw std::invalid_argument("INVALID BASE");

    if (!(base & (base - 1))) {
        data = from_binp(s, j, base);
        return *this;
    }

    long long unsigned cur = 0;
    std::vector<long long unsigned> p {1};
    unsigned base_bits = -1, cur_bits = 0;

    while (p.back() <= dmax)
        p.push_back(p.back() * base),
        ++base_bits;

    std::vector<unsigned> new_data;

    if (base <= 64) {

        while (i - j < s.size()) {

            v = getv(s[i]);

            if (v == -1)
                throw init_except("INAPPROPRIATE VALUE", i, s);

            if (base <= 10 && (v < 0 || v > 9))
                throw init_except("FOR BASE " + ntos(base) + " NON-INTEGER VALUE", i, s);

            if (v >= base)
                throw init_except("FOR BASE " + ntos(base) + " VERY LARGE BIT", i, s);

            cur += v * p[cur_bits];
            --i, ++cur_bits;

            if (cur_bits >= base_bits) {
                cur_bits -= base_bits;
                new_data.push_back(cur % p[base_bits]);
                cur /= p[base_bits];
            }
        }
    } else {

        auto to = [&s, &j, &base] (unsigned& i) -> unsigned {

            unsigned n = 0, k = 0, v;

            while (i - j < s.size() && s[i] != '\'') {

                v = getv(s[i]);

                if (v == -1)
                    throw init_except("INAPPROPRIATE VALUE", i, s);

                if (v < 0 || v > 9)
                    throw init_except("FOR BASE " + ntos(base) + " (> 64) NON-INTEGER VALUE", i, s);

                n += v * p10[k];

                if (k > 7)
                    throw init_except("VERY LARGE BIT WIDTH (" + ntos(k + 1) + ")", i, s);

                if (n >= base)
                    throw init_except("BIT (" + ntos(n) + ") LARGER THAN BASE (" + ntos(base) + ")", i, s);

                ++k, --i;
            }
            return n;
        };
        while (i - j < s.size()) {

            cur += to(i) * p[cur_bits];
            --i, ++cur_bits;

            if (cur_bits >= base_bits) {
                cur_bits -= base_bits;
                new_data.push_back(cur % p[base_bits]);
                cur /= p[base_bits];
            }
        }
    }
    new_data.push_back(cur);

    trim(new_data);

    data = to(new_data, p[base_bits]);
    trim(data);

    return *this;
}

Integer::Integer() noexcept : sign(false) {}

Integer::Integer(unsigned n) noexcept { init(n); }
Integer::Integer(unsigned long long n) noexcept { init(n); }

Integer::Integer(int n) noexcept
    : sign(n < 0) {
    init(sign ? -n : n, sign);
}

Integer::Integer(long long n) noexcept
    : sign(n < 0) {
    init(sign ? -n : n, sign);
}

Integer::Integer(const char *s) { init(s); }
Integer::Integer(const std::string &s) { init(s); }

Integer::Integer(const Integer &n) noexcept
    : data(n.data), sign(n.sign) {
}
Integer::Integer(Integer &&n) noexcept
    : data(std::move(n.data)), sign(n.sign) {
}

Integer &Integer::operator=(const Integer &n) noexcept {

    if (&n != this) {
        data = n.data;
        sign = n.sign;
    }

    return *this;
}

Integer &Integer::operator=(Integer &&n) noexcept {

    if (&n != this) {
        data = std::move(n.data);
        sign = n.sign;
    }

    return *this;
}

unsigned Integer::bit_count() const noexcept {

    if (data.empty())
        return 0;

    return ((data.size() - 1) << 4) + pop_count(data.back());
}

unsigned Integer::size() const noexcept {

    if (data.empty())
        return 0;

    return data.size() << 1;
}

void Integer::swap(Integer &n) noexcept {

    Integer t = std::move(n);
    n = std::move(*this);
    *this = std::move(t);
}

Integer &Integer::operator+=(const Integer &n) {

    if (sign == n.sign) {

        data.resize(std::max(data.size(), n.data.size()) + 1);
        add(data, n.data, 0, 0, n.data.size());
    } else {

        if (vec_less(data, n.data)) {

            std::vector<unsigned short> temp = n.data;
            sub(temp, data, 0, 0, data.size());
            data = std::move(temp);
            sign = !sign;

        } else
            sub(data, n.data, 0, 0, n.data.size());
    }

    trim(data);

    if (data.empty())
        sign = false;

    return *this;
}

Integer operator+(const Integer &a, const Integer &b) {
    return std::move(Integer(a) += b);
}

Integer &Integer::operator-=(const Integer &n) {

    bool is_less = vec_less(data, n.data);

    if (sign == n.sign) {

        sign ^= is_less;

        if (is_less) {

            std::vector<unsigned short> temp = n.data;
            sub(temp, data, 0, 0, data.size());
            data = std::move(temp);
        } else
            sub(data, n.data, 0, 0, n.data.size());

    } else {

        data.resize(std::max(data.size(), n.data.size()) + 1);
        add(data, n.data, 0, 0, n.data.size());
    }
    trim(data);

    if (data.empty())
        sign = false;

    return *this;
}

Integer operator-(const Integer &a, const Integer &b) {
    return std::move(Integer(a) -= b);
}

Integer Integer::operator-() const noexcept {

    Integer temp = *this;
    temp.sign = !temp.sign;
    return temp;
}

Integer &Integer::flip_sign() noexcept {

    sign = !sign;
    return *this;
}

Integer Integer::abs() const {

    Integer temp = *this;
    temp.sign = false;
    return temp;
}

Integer &Integer::make_abs() {

    sign = false;
    return *this;
}

Integer &Integer::operator*=(const Integer &n) {

    data = ksm(data, n.data);
    sign ^= n.sign;
    return *this;
}

Integer &Integer::operator*=(const unsigned &n) {

    mul(data, n);
    return *this;
}

Integer &Integer::operator*=(const int &n) {

    sign ^= (n < 0);
    mul(data, n < 0 ? -n : n);
    return *this;
}

Integer operator*(const Integer &x, const Integer &y) {
    return std::move(Integer(x) *= y);
}
Integer operator*(const Integer &x, const unsigned &y) {
    return std::move(Integer(x) *= y);
}
Integer operator*(const unsigned &x, const Integer &y) {
    return std::move(Integer(y) *= x);
}
Integer operator*(const Integer &x, const int &y) {
    return std::move(Integer(x) *= y);
}
Integer operator*(const int &x, const Integer &y) {
    return std::move(Integer(y) *= x);
}

Integer &Integer::operator/=(const Integer &n) {

    sign ^= n.sign;
    std::vector<unsigned short> quot, rem;
    div(data, n.data, quot, rem);
    data = std::move(quot);
    return *this;
}

Integer &Integer::operator/=(const unsigned &n) {

    std::vector<unsigned short> quot;
    unsigned rem;
    div(data, n, quot, rem);
    data = std::move(quot);
    return *this;
}

Integer &Integer::operator/=(const int &n) {

    sign ^= (n < 0);
    std::vector<unsigned short> quot;
    unsigned rem;
    div(data, n < 0 ? -n : n, quot, rem);
    data = std::move(quot);
    return *this;
}

Integer Integer::operator/(const unsigned &n) const {
    return std::move(Integer(*this) /= n);
}
Integer Integer::operator/(const int &n) const {
    return std::move(Integer(*this) /= n);
}
Integer operator/(const Integer &a, const Integer &b) {
    return std::move(Integer(a) /= b);
}

Integer &Integer::operator%=(const Integer &n) {

    std::vector<unsigned short> quot, rem;
    div(data, n.data, quot, rem);
    data = std::move(rem);
    return *this;
}

Integer &Integer::operator%=(const unsigned &n) {

    std::vector<unsigned short> quot, rem;
    unsigned r;
    div(data, n, quot, r);

    if (r)
        rem.emplace_back(r), r >>= bits;
    if (r)
        rem.emplace_back(r);

    data = std::move(rem);
    return *this;
}

Integer &Integer::operator%=(const int &n) {

    std::vector<unsigned short> quot, rem;
    unsigned r;
    div(data, n < 0 ? -n : n, quot, r);

    if (r)
        rem.emplace_back(r), r >>= bits;
    if (r)
        rem.emplace_back(r);

    data = std::move(rem);
    return *this;
}

Integer Integer::operator%(const unsigned &n) const {
    return std::move(Integer(*this) %= n);
}
Integer Integer::operator%(const int &n) const {
    return std::move(Integer(*this) %= n);
}
Integer operator%(const Integer &a, const Integer &b) {
    return std::move(Integer(a) %= b);
}

Div div(const Integer &a, const Integer &b) {

    Div qr;
    qr.quot.sign = a.sign ^ b.sign;
    qr.rem.sign = a.sign;
    Integer::div(a.data, b.data, qr.quot.data, qr.rem.data);
    return qr;
}

bool operator==(const Integer &a, const Integer &b) noexcept {
    return (a.sign == b.sign) && Integer::vec_equal(a.data, b.data);
}

bool operator<(const Integer &a, const Integer &b) noexcept {
    if (a.sign != b.sign)
        return (a.sign < b.sign);
    return a.sign ^ Integer::vec_less(a.data, b.data);
}

bool operator>(const Integer &a, const Integer &b) noexcept { return b < a; }
bool operator!=(const Integer &a, const Integer &b) noexcept { return !(a == b); }
bool operator<=(const Integer &a, const Integer &b) noexcept { return !(a > b); }
bool operator>=(const Integer &a, const Integer &b) noexcept { return !(a < b); }

Integer &Integer::operator++() noexcept {

    if (sign) {

        for (auto &x : data) {
            if (x) {
                --x;
                break;
            }
            --x;
        }
    } else {

        data.emplace_back(0);

        for (auto &x : data) {
            ++x;
            if (x)
                break;
        }
    }

    trim(data);

    if (data.empty())
        sign = false;

    return *this;
}

Integer &Integer::operator--() noexcept {

    if (sign) {

        data.emplace_back(0);

        for (auto &x : data) {
            ++x;
            if (x)
                break;
        }
    } else {

        for (auto &x : data) {
            if (x) {
                --x;
                break;
            }
            --x;
        }
    }

    trim(data);

    if (data.empty())
        sign = false;

    return *this;
}

Integer Integer::operator++(int) noexcept {

    Integer temp = *this;
    ++(*this);

    return temp;
}

Integer Integer::operator--(int) noexcept {

    Integer temp = *this;
    --(*this);

    return temp;
}

Integer &Integer::operator|=(const Integer &n) {

    const unsigned m = std::max(data.size(), n.data.size());

    std::vector<unsigned short> x = data;
    std::vector<unsigned short> y = n.data;

    x.resize(m), y.resize(m);

    if (sign)
        _2s_comp(x);
    if (n.sign)
        _2s_comp(y);

    for (unsigned i = 0; i < m; ++i)
        x[i] |= y[i];

    if ((sign |= n.sign))
        _2s_comp(x);

    trim(x);
    data = std::move(x);

    return *this;
}

Integer operator|(const Integer &a, const Integer &b) {
    return std::move(Integer(a) |= b);
}

Integer &Integer::operator&=(const Integer &n) {

    const unsigned m = std::max(data.size(), n.data.size());

    std::vector<unsigned short> x = data;
    std::vector<unsigned short> y = n.data;

    x.resize(m), y.resize(m);

    if (sign)
        _2s_comp(x);
    if (n.sign)
        _2s_comp(y);

    for (unsigned i = 0; i < m; ++i)
        x[i] &= y[i];

    if (sign &= n.sign)
        _2s_comp(x);

    trim(x);
    data = std::move(x);

    return *this;
}

Integer operator&(const Integer &a, const Integer &b) {
    return std::move(Integer(a) &= b);
}

Integer &Integer::operator^=(const Integer &n) {

    const unsigned m = std::max(data.size(), n.data.size());

    std::vector<unsigned short> x = data;
    std::vector<unsigned short> y = n.data;

    x.resize(m), y.resize(m);

    if (sign)
        _2s_comp(x);
    if (n.sign)
        _2s_comp(y);

    for (unsigned i = 0; i < m; ++i)
        x[i] ^= y[i];

    if (sign ^= n.sign)
        _2s_comp(x);

    trim(x);
    data = std::move(x);

    return *this;
}

Integer operator^(const Integer &a, const Integer &b) {
    return std::move(Integer(a) ^= b);
}

Integer &Integer::operator<<=(unsigned shift) {

    if (data.empty() || shift <= 0)
        return *this;

    left_shift(data, shift);

    return *this;
}

Integer Integer::operator<<(unsigned shift) const {
    return std::move(Integer(*this) <<= shift);
}

Integer &Integer::operator>>=(unsigned shift) {

    if (data.empty() || shift <= 0)
        return *this;

    right_shift(data, shift);

    if (data.empty())
        sign = false;

    return *this;
}

Integer Integer::operator>>(unsigned shift) const {
    return std::move(Integer(*this) >>= shift);
}

Integer Integer::operator~() const {

    Integer temp = *this;
    temp.sign = !temp.sign;
    return --temp;
}

Integer& Integer::compliment() {

    sign = !sign;
    return --(*this);
}

std::ostream &operator<<(std::ostream &out, const Integer &n) {

    if (n.data.size()) {

        std::vector<unsigned> temp = Integer::to(n.data);
        unsigned i = temp.size() - 1;

        if (n.sign)
            out << '-';

        out << temp[i];

        for (--i; i < temp.size(); --i)
            out << std::setfill('0') << std::setw(digs) << temp[i];

    } else
        out << 0;

    return out;
}

std::istream &operator>>(std::istream &in, Integer &n) {

    std::string input;
    in >> input;
    n.init(input);
    return in;
}

Integer Integer::sqr() const {

    Integer temp;
    temp.data = ksm_sqr(data);
    temp.sign = false;
    return temp;
}

Integer Integer::nrt(uint8_t p) const {

    if (sign && ((p & 1) == 0))
        throw std::invalid_argument("NO REAL ROOT FOUND!");
    if (data.empty())
        return Integer();
    if (*this == 1 || p == 0)
        return 1;
    if (p == 1)
        return *this;

    Integer x, x0, a = abs();
    x.data.resize((data.size() / p) + 2);
    x.data.back() = 1;

    do {
        x0 = x;
        x = (p - 1) * x + (a / x.pow(p - 1));
        x /= p;
    } while (x < x0);

    x0.sign = sign;

    return x0;
}

Integer Integer::sqrt() const {

    if (sign)
        throw std::invalid_argument("NO REAL ROOT FOUND!");
    if (data.empty())
        return Integer();
    if (*this == 1)
        return 1;

    Integer x, x0;
    x.data.resize((data.size() >> 1) + 2);
    x.data.back() = 1;

    do {
        x0 = x;
        x += *this / x;
        x >>= 1;
    } while (x < x0);

    return x0;
}

Integer Integer::pow(uint8_t p) const {

    if (data.empty())
        return Integer();
    if (*this == 1)
        return 1;

    Integer ans = 1, temp = *this;
    while (p) {
        if (p & 1)
            ans *= temp;
        if (p >>= 1)
            temp = temp.sqr();
    }
    return ans;
}

Integer Integer::pow(uint8_t p, uint8_t q) const {

    if (q == 0)
        throw std::invalid_argument(ntos(p) + " / " + ntos(q) + " <--");
    if (p == 0)
        return 1;
    if (q == 1)
        return pow(p);
    if (p == 1)
        return nrt(q);

    uint8_t g = ::gcd(p, q);
    p /= g;
    q /= g;

    return pow(p).nrt(q);
}

Integer Integer::powf(float pf) const {
    std::pair<unsigned, unsigned> pq = get_frac(pf, 1e-2);
    return pow(pq.first).nrt(pq.second);
}

Integer Integer::gcd(const Integer &n) const {
    return n ? n.gcd(*this % n) : *this;
}

Integer Integer::lcm(const Integer &n) const {
    return (*this * n) / gcd(n);
}

bool Integer::odd() const noexcept {
    return (*this)[0];
}

bool Integer::even() const noexcept {
    return !(*this)[0];
}

bool Integer::operator[](unsigned index) const noexcept {

    if (index >= bit_count())
        return 0;
    return data[index >> 4] & (1 << (index & 15));
}

Integer& Integer::set_bit(unsigned index) noexcept {

    if (index >= bit_count())
        data.resize((index + bits) >> 4);
    data[index >> 4] |= (1 << (index & 15));

    return *this;
}

Integer& Integer::unset_bit(unsigned index) noexcept {

    if (index >= bit_count())
        return *this;
    data[index >> 4] &= ~(1 << (index & 15));
    trim(data);

    return *this;
}

std::vector<unsigned short> Integer::fact(const unsigned st, const unsigned en) {

    if(en - st < 2) {

        long long unsigned cur = st;

        if (st != en) cur *= en;

        std::vector<unsigned short> ans;

        while (cur)
            ans.push_back(cur),
            cur >>= bits;

        return ans;
    }

    const unsigned m = st + ((en - st) >> 1);

    std::vector<unsigned short> lhs = fact(st, m);
    std::vector<unsigned short> rhs = fact(m + 1, en);

    rhs = ksm(lhs, rhs);

    return rhs;
}

Integer fact(const unsigned n) {

    if (n == 0)
        return Integer(1);

    if (n <= 20) {
        long long unsigned cur = 1;
        for (unsigned i = 2; i <= n; ++i)
            cur *= i;
        return cur;
    }

    Integer f;
    f.data = Integer::fact(2, n);

    return f;
}

long double Integer::log2() const {

    if (data.empty())
        throw std::domain_error("LOG(0) is NOT DEFINED");

    if (size() == 1 && data.back() == 1)
        return 0;

    unsigned bitc = bit_count();
    long double exp = 1, temp = 0, ans = 0;

    for (unsigned i = 0; i < 4 && i < data.size(); ++i)
        temp = temp * base + *(data.end() - i - 1);

    if (bitc > 48) {
        ans = bitc - 48;
        temp /= (1 << pop_count(data.back()));
    }

    auto fp = [](long double n, unsigned p) -> long double {
        long double ans = 1, cur = n;
        while (p) {
            if (p & 1) ans *= cur;
            if (p >>= 1) cur *= cur;
        }
        return ans;
    };

    for (unsigned short i = 0; i <= 12; ++i) {

        while (temp >= 2)
            ans += exp, temp /= 2;

        if (temp <= 1) return ans;

        temp = fp(temp, 10);
        exp /= 10;
    }

    return ans;
}

long double Integer::log10() const {
    return log2() / 3.321928094887;
}

long double Integer::loge() const {
    return log2() / 1.442695040889;
}

long double Integer::log(long double base) const {
    return log2() / ::log2(base);
}

Integer rand(unsigned digit_count) {

    if (digit_count == 0)
        ++digit_count;

    seed_rand();

    const unsigned lb = bitc(digit_count-1), ub = bitc(digit_count);

    return rand_bit(1 + lb + rand() % (ub - lb - 1));
}

Integer rand_bit(unsigned bit_count) {

    if (bit_count == 0)
        return 0;

    seed_rand();

    Integer r;
    r.set_bit(bit_count - 1);

    const unsigned size = r.data.size();

    r.data.back() |= rand() & ((1 << (bit_count & 0xF)) - 1);

    for (unsigned i = size - 2; i < size; --i)
        r.data[i] = rand() & bmax;

    return r;
}

Integer rand_byte(unsigned byte_count) {

    if (byte_count == 0)
        return 0;

    seed_rand();

    Integer r = rand_bit(byte_count << 3);

    r.data.back() = rand();

    return r;
}

Integer rand_range(unsigned min_digit_count, unsigned max_digit_count) {

    if (min_digit_count > max_digit_count || max_digit_count == 0)
        return 0;

    seed_rand();

    const unsigned lb = bitc(min_digit_count-1), ub = bitc(max_digit_count);

    return rand_bit(1 + lb + rand() % (ub - lb - 1));
}

Integer rand_range(const Integer& lower, const Integer& upper) {

    if (lower > upper)
        return 0;

    seed_rand();

    Integer diff = upper - lower;

    if (diff < UINT64_MAX) {

        long long unsigned r = rand();
        r = (r << bits) | rand();
        r = (r << bits) | rand();
        r = (r << bits) | rand();

        return lower + r % diff;
    }

    unsigned dc = digc(diff.bit_count());

    if (rand() & 3)
        return lower + rand(dc) % diff;

    return lower + rand_range(rand() % dc, dc-1);
}

Integer::operator bool() const noexcept {
    return !data.empty();
}

Integer::operator std::string() const noexcept {
    return str(10, false);
}

std::string Integer::str(const unsigned base, bool prefix) const {
    if (base < 2 || base > dmax)
        throw std::invalid_argument("INVALID BASE");

    if (data.empty())
        return prefix ? ("0" + ntos(base) + "...0") : ("0");

    if (!(base & (base - 1)))
        return to_binp(data, sign, pop_count(base) - 1, prefix);

    std::string s;

    if (sign)
        s.push_back('-');

    if (prefix)
        s.push_back('0'),
        s.append(ntos(base)),
        s.append("...");

    unsigned base_bits = 0;
    long long unsigned mbase = 1;
    while (mbase * base <= dmax)
        mbase *= base, ++base_bits;

    std::vector<unsigned> converted_data = to(data, mbase);

    std::stack<unsigned> stk;

    for (auto cur : converted_data)
        for (unsigned i = 0; i < base_bits; ++i)
            stk.push(cur % base),
            cur /= base;

    while (stk.top() == 0) stk.pop();

    if (base < 64) {
        while (!stk.empty())
            s.push_back(cset[stk.top()]),
            stk.pop();
    } else {
        while (!stk.empty())
            s.append(ntos(stk.top())),
            s.push_back('\''),
            stk.pop();
        s.pop_back();
    }

    return s;
}

std::string Integer::dec(bool prefix) const {
    return str(10, prefix);
}

std::string Integer::bin(bool prefix) const {
    return to_binp(data, sign, 1, prefix);
}

std::string Integer::oct(bool prefix) const {
    return to_binp(data, sign, 3, prefix);
}

std::string Integer::hex(bool prefix) const {
    return to_binp(data, sign, 4, prefix);
}

    // Integer DEFINITIONS #END

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Div DEFINITION

Div::Div() {}
Div::Div(const Integer &q, const Integer &r) noexcept
    : quot(q), rem(r) {
}
Div::Div(const Integer &q, Integer &&r) noexcept
    : quot(q), rem(std::move(r)) {
}
Div::Div(Integer &&q, const Integer &r) noexcept
    : quot(std::move(q)), rem(r) {
}
Div::Div(Integer &&q, Integer &&r) noexcept
    : quot(std::move(q)), rem(std::move(r)) {
}

Div::Div(const Div &qr) noexcept
    : quot(qr.quot), rem(qr.rem) {
}
Div::Div(Div &&qr) noexcept
    : quot(std::move(qr.quot)), rem(std::move(qr.rem)) {
}

Div &Div::operator=(const Div &qr) noexcept {

    if (&qr != this) {
        quot = qr.quot;
        rem = qr.rem;
    }

    return *this;
}

Div &Div::operator=(Div &&qr) noexcept {

    if (&qr != this) {
        quot = std::move(qr.quot);
        rem = std::move(qr.rem);
    }

    return *this;
}

std::ostream &operator<<(std::ostream &out, const Div &qr) {
    out << "Quotient : " << qr.quot << '\n';
    out << "Remainder : " << qr.rem << '\n';
    return out;
}


void Integer::trim (std::vector<unsigned short> &data) noexcept {
    while (data.size() && data.back() == 0)
        data.pop_back();
}

void Integer::trim (std::vector<unsigned> &data) noexcept {
    while (data.size() && data.back() == 0)
        data.pop_back();
}

void Integer::_1s_comp (std::vector<unsigned short> &n) noexcept {
    for (unsigned i = 0; i < n.size(); ++i)
        n[i] = ~n[i];
}

void Integer::_2s_comp (std::vector<unsigned short> &n) noexcept {
    bool carry = true;
    for (unsigned i = 0; i < n.size(); ++i) {
        if (carry) {
            if (n[i])
                carry = false;
            --n[i];
        }
        n[i] = ~n[i];
    }
}

bool Integer::vec_equal (
    const std::vector<unsigned short> &a, const std::vector<unsigned short> &b
) noexcept {

    if (a.size() != b.size())
        return false;

    for (unsigned i = a.size() - 1; i < a.size(); --i)
        if (a[i] != b[i])
            return false;

    return true;
}

bool Integer::vec_less (
    const std::vector<unsigned short> &a, const std::vector<unsigned short> &b
) noexcept {

    if (a.size() != b.size())
        return (a.size() < b.size());

    for (unsigned i = a.size() - 1; i < a.size(); --i)
        if (a[i] != b[i])
            return (a[i] < b[i]);

    return false;
}

void Integer::left_shift (
    std::vector<unsigned short> &data, unsigned shift
) {
    unsigned elem_shift = (shift >> 4);
    unsigned short bits_shift = (shift & 15);
    unsigned short rem_bits_shift = bits - bits_shift;

    if (bits_shift == 0) {

        data.resize(data.size() + elem_shift);
        unsigned i = data.size() - 1, j = i - elem_shift;

        for (; j < data.size(); --i, --j)
            data[i] = data[j];
        for (; i < data.size(); --i)
            data[i] = 0;

        return trim(data);
    }

    data.resize(data.size() + elem_shift + 1);
    unsigned i = data.size() - 1, j = i - elem_shift - 1;

    if (j < data.size()) {

        data[i] = data[j] >> rem_bits_shift;

        for (--i; j; --i)
            data[i] = data[j] << bits_shift | data[--j] >> rem_bits_shift;

        data[i] = data[0] << bits_shift;

        while (--i < data.size())
            data[i] = 0;
    }

    return trim(data);
}

void Integer::right_shift (
    std::vector<unsigned short> &data, unsigned shift
) {
    unsigned elem_shift = (shift >> 4);
    unsigned short bits_shift = (shift & 15);
    unsigned short rem_bits_shift = bits - bits_shift;

    if (elem_shift >= data.size())
        return data.clear();

    if (bits_shift == 0) {

        unsigned i = 0, j = elem_shift;

        for (; j < data.size(); ++i, ++j)
            data[i] = data[j];
        data.resize(i);

        return trim(data);
    }

    unsigned i = 0, j = elem_shift;

    for (; j + 1 < data.size(); ++i)
        data[i] = data[j] >> bits_shift | data[++j] << rem_bits_shift;

    data[i] = data[j] >> bits_shift;
    data.resize(i + 1);

    return trim(data);
}

void Integer::add (
    std::vector<unsigned short> &x, const std::vector<unsigned short> &y,
    unsigned xst, unsigned yst, unsigned yn
) {
    bool carry = 0;
    unsigned cur = 0, i = 0;

    for (; i < x.size() && i < yn; ++i) {
        cur = (unsigned) x[i + xst] + y[i + yst] + carry;
        carry = (cur > bmax);
        x[i + xst] = cur;
    }

    for (; i < x.size() && carry; ++i) {
        if (x[i + xst] < bmax)
            carry = 0;
        ++x[i + xst];
    }
}

void Integer::add (
    std::vector<unsigned> &x, const std::vector<unsigned> &y,
    unsigned xst, unsigned yst, unsigned yn, const unsigned base
) {
    bool carry = 0;
    unsigned cur = 0, i = 0;

    for (; i < x.size() && i < yn; ++i) {
        cur = x[i + xst] + y[i + yst] + carry;
        carry = (cur >= base);
        x[i + xst] = cur - (carry ? base : 0);
    }
    for (; i < x.size() && carry; ++i) {
        if (x[i + xst] + 1 < base) {
            ++x[i + xst];
            break;
        }
        x[i + xst] = 0;
    }
}

void Integer::sub (
    std::vector<unsigned short> &x, const std::vector<unsigned short> &y,
    unsigned xst, unsigned yst, unsigned yn
) {
    int cur = 0;
    bool carry = 0;
    unsigned i = 0;

    for (; i < x.size() && i < yn; ++i) {
        cur = (int) x[i + xst] - y[i + yst] - carry;
        carry = (cur < 0);
        x[i + xst] = cur;
    }

    for (; i < x.size() && carry; ++i) {
        if (x[i + xst])
            carry = false;
        --x[i + xst];
    }
}

void Integer::sub (
    std::vector<unsigned> &x, const std::vector<unsigned> &y,
    unsigned xst, unsigned yst, unsigned yn, const unsigned base
) {
    int cur = 0;
    bool carry = 0;
    unsigned i = 0;

    for (; i < x.size() && i < yn; ++i) {
        cur = (int) x[i + xst] - y[i + yst] - carry;
        carry = (cur < 0);
        x[i + xst] = cur + (carry ? base : 0);
    }
    for (; i < x.size() && carry; ++i) {
        if (x[i + xst]) {
            --x[i + xst];
            break;
        }
        x[i + xst] = base - 1;
    }
}

std::vector<unsigned short> Integer::mul (
    const std::vector<unsigned short> &x, const std::vector<unsigned short> &y,
    unsigned xst, unsigned xn, unsigned yst, unsigned yn
) {
    unsigned long long cur = 0;
    const unsigned size = xn + yn - 1;
    std::vector<unsigned short> ans(size + 1);

    for (unsigned z = 0; z < size; ++z) {
        unsigned s = z >= yn ? z - yn + 1 : 0;
        unsigned e = std::min(xn - 1, z) + xst;
        unsigned j = z - s + yst;
        unsigned i = s + xst;

        for (; i <= e; ++i, --j)
            cur += (unsigned) x[i] * y[j];

        ans[z] = cur, cur >>= bits;
    }
    ans.back() = cur;

    return ans;
}

void Integer::mul (
    std::vector<unsigned short> &data, unsigned n
) {
    if (n == 1)
        return;
    if (data.empty())
        return;
    if (n == 0)
        return data.clear();

    unsigned long long cur = 0;
    for (unsigned i = 0; i < data.size(); ++i) {
        cur += (unsigned long long) data[i] * n;
        data[i] = cur;
        cur >>= bits;
    }

    if (cur)
        data.emplace_back(cur), cur >>= bits;
    if (cur)
        data.emplace_back(cur);
}

std::vector<unsigned> Integer::mul (
    const std::vector<unsigned> &x, const std::vector<unsigned> &y,
    unsigned xst, unsigned xn, unsigned yst, unsigned yn, const unsigned base
) {
    const unsigned size = xn + yn - 1;
    std::vector<unsigned> ans(size + 1);
    lldiv_t cur = {0, 0};

    for (unsigned z = 0; z < size; ++z) {
        unsigned s = z >= yn ? z - yn + 1 : 0;
        unsigned e = std::min(xn - 1, z) + xst;
        unsigned j = z - s + yst;
        unsigned i = s + xst;

        for (; i <= e; ++i, --j)
            cur.quot += (unsigned long long) x[i] * y[j];

        cur = lldiv(cur.quot, base);
        ans[z] = cur.rem;
    }
    ans.back() = cur.quot;

    return ans;
}

std::vector<unsigned short> Integer::sqr (
    const std::vector<unsigned short> &x, unsigned st, unsigned n
) {
    unsigned long long cur = 0;
    const unsigned size = (n << 1) - 1;
    std::vector<unsigned short> ans(size + 1);

    for (unsigned z = 0; z < size; ++z) {
        unsigned s = (z >= n ? z - n + 1 : 0) + st;
        unsigned e = std::min(n - 1, z) + st;

        for (; s < e; ++s, --e)
            cur += ((unsigned long long) x[s] * x[e]) << 1;

        if (s == e)
            cur += (unsigned) x[s] * x[e];

        ans[z] = cur, cur >>= bits;
    }
    ans.back() = cur;

    return ans;
}

std::vector<unsigned> Integer::sqr (
    const std::vector<unsigned> &x, unsigned st, unsigned n,
    const unsigned base
) {
    const unsigned size = (n << 1) - 1;
    std::vector<unsigned> ans(size + 1);
    lldiv_t cur = {0, 0};

    for (unsigned z = 0; z < size; ++z) {
        unsigned s = (z >= n ? z - n + 1 : 0) + st;
        unsigned e = std::min(n - 1, z) + st;

        for (; s < e; ++s, --e)
            cur.quot += ((unsigned long long) x[s] * x[e]) << 1;

        if (s == e)
            cur.quot += (unsigned long long) x[s] * x[e];

        cur = lldiv(cur.quot, base);
        ans[z] = cur.rem;
    }
    ans.back() = cur.quot;

    return ans;
}

std::vector<unsigned short> Integer::ksm_sqr (
    const std::vector<unsigned short> &x, unsigned st, unsigned n
) {
    if (n < ksmt)
        return sqr(x, st, n);

    unsigned m = n - (n >> 1);

    std::vector<unsigned short> xx(x.begin() + st, x.begin() + st + m);
    xx.emplace_back(0), add(xx, x, 0, st + m, n - m);

    std::vector<unsigned short> xx2 = ksm_sqr(xx, 0, xx.size());
    std::vector<unsigned short> xl2 = ksm_sqr(x, st, m);
    std::vector<unsigned short> xh2 = ksm_sqr(x, st + m, n - m);
    sub(xx2, xl2, 0, 0, xl2.size());
    sub(xx2, xh2, 0, 0, xh2.size());

    unsigned i = 0, j;
    std::vector<unsigned short> ans(n << 1);
    for (j = 0; j < xl2.size(); ++i, ++j)
        ans[i] = xl2[j];
    for (j = 0; j < xh2.size(); ++i, ++j)
        ans[i] = xh2[j];
    add(ans, xx2, m, 0, xx2.size());

    return ans;
}

std::vector<unsigned short> Integer::ksm_sqr (
    const std::vector<unsigned short> &x
) {
    const unsigned n = x.size();

    if (n == 0)
        return std::vector<unsigned short>{};

    if (n == 1 && x.back() == 1)
        return std::vector<unsigned short>{1};

    std::vector<unsigned short> ans;

    if (n < ksmt)
        ans = sqr(x, 0, n);
    else
        ans = ksm_sqr(x, 0, n);

    trim(ans);

    return ans;
}

std::vector<unsigned> Integer::ksm_sqr (
    const std::vector<unsigned> &x, unsigned st, unsigned n,
    const unsigned base
) {
    if (n < ksmt)
        return sqr(x, st, n, base);

    unsigned m = n - (n >> 1);

    std::vector<unsigned> xx(x.begin() + st, x.begin() + st + m);
    xx.emplace_back(0), add(xx, x, 0, st + m, n - m, base);

    std::vector<unsigned> xx2 = ksm_sqr(xx, 0, xx.size(), base);
    std::vector<unsigned> xl2 = ksm_sqr(x, st, m, base);
    std::vector<unsigned> xh2 = ksm_sqr(x, st + m, n - m, base);
    sub(xx2, xl2, 0, 0, xl2.size(), base);
    sub(xx2, xh2, 0, 0, xh2.size(), base);

    unsigned i = 0, j;
    std::vector<unsigned> ans(n << 1);
    for (j = 0; j < xl2.size(); ++i, ++j)
        ans[i] = xl2[j];
    for (j = 0; j < xh2.size(); ++i, ++j)
        ans[i] = xh2[j];
    add(ans, xx2, m, 0, xx2.size(), base);

    return ans;
}

std::vector<unsigned> Integer::ksm_sqr (
    const std::vector<unsigned> &x, const unsigned base
) {
    const unsigned n = x.size();

    if (n == 0)
        return std::vector<unsigned>{};
    if (n == 1 && x.back() == 1)
        return std::vector<unsigned>{1};

    std::vector<unsigned> ans;

    if (n < ksmt)
        ans = sqr(x, 0, n, base);
    else
        ans = ksm_sqr(x, 0, n, base);

    trim(ans);

    return ans;
}

std::vector<unsigned short> Integer::ksm (
    const std::vector<unsigned short> &x, const std::vector<unsigned short> &y,
    unsigned xst, unsigned yst, unsigned n
) {
    if (n < ksmt)
        return mul(x, y, xst, n, yst, n);

    unsigned m = n - (n >> 1);

    std::vector<unsigned short> xx(x.begin() + xst, x.begin() + xst + m);
    std::vector<unsigned short> yy(y.begin() + yst, y.begin() + yst + m);
    xx.emplace_back(0), add(xx, x, 0, xst + m, n - m);
    yy.emplace_back(0), add(yy, y, 0, yst + m, n - m);

    std::vector<unsigned short> xy = ksm(xx, yy, 0, 0, xx.size());
    std::vector<unsigned short> xyl = ksm(x, y, xst, yst, m);
    std::vector<unsigned short> xyh = ksm(x, y, xst + m, yst + m, n - m);
    sub(xy, xyl, 0, 0, xyl.size());
    sub(xy, xyh, 0, 0, xyh.size());

    unsigned i = 0, j;
    std::vector<unsigned short> ans(n << 1);
    for (j = 0; j < xyl.size(); ++i, ++j)
        ans[i] = xyl[j];
    for (j = 0; j < xyh.size(); ++i, ++j)
        ans[i] = xyh[j];
    add(ans, xy, m, 0, xy.size());

    return ans;
}

std::vector<unsigned short> Integer::ksm (
    std::vector<unsigned short> x, std::vector<unsigned short> y,
    unsigned xn, unsigned yn
) {
    if (xn < yn)
        return ksm(std::move(y), std::move(x), yn, xn);

    if (3 * yn > (xn << 1)) {
        y.resize(xn);
        return ksm(x, y, 0, 0, xn);
    }

    unsigned st = 0;
    x.resize(xn + (yn - xn % yn) % yn);
    std::vector<unsigned short> ans(x.size() + yn), carry(yn);

    while (st < xn) {

        std::vector<unsigned short> temp = ksm(x, y, st, 0, yn);
        add(temp, carry, 0, 0, yn);

        for (unsigned i = 0; i < yn; ++i, ++st)
            carry[i] = temp[i + yn], ans[st] = temp[i];
    }
    for (unsigned i = 0; i < yn; ++i, ++st)
        ans[st] = carry[i];

    return ans;
}

std::vector<unsigned short> Integer::ksm (
    const std::vector<unsigned short> &x, const std::vector<unsigned short> &y
) {
    const unsigned xn = x.size();
    const unsigned yn = y.size();

    if (xn == 0 || yn == 0)
        return {};
    if (yn == 1 && y.back() == 1)
        return x;
    if (xn == 1 && x.back() == 1)
        return y;

    std::vector<unsigned short> ans;

    if (std::min(xn, yn) < ksmt)
        ans = mul(x, y, 0, xn, 0, yn);
    else
        ans = ksm(x, y, xn, yn);

    trim(ans);

    return ans;
}

std::vector<unsigned> Integer::ksm (
    const std::vector<unsigned> &x, const std::vector<unsigned> &y,
    unsigned xst, unsigned yst, unsigned n, const unsigned base
) {
    if (n < ksmt)
        return mul(x, y, xst, n, yst, n, base);

    unsigned m = n - (n >> 1);

    std::vector<unsigned> xx(x.begin() + xst, x.begin() + xst + m);
    std::vector<unsigned> yy(y.begin() + yst, y.begin() + yst + m);
    xx.emplace_back(0), add(xx, x, 0, xst + m, n - m, base);
    yy.emplace_back(0), add(yy, y, 0, yst + m, n - m, base);

    std::vector<unsigned> xy = ksm(xx, yy, 0, 0, xx.size(), base);
    std::vector<unsigned> xyl = ksm(x, y, xst, yst, m, base);
    std::vector<unsigned> xyh = ksm(x, y, xst + m, yst + m, n - m, base);
    sub(xy, xyl, 0, 0, xyl.size(), base);
    sub(xy, xyh, 0, 0, xyh.size(), base);

    unsigned i = 0, j;
    std::vector<unsigned> ans(n << 1);
    for (j = 0; j < xyl.size(); ++i, ++j)
        ans[i] = xyl[j];
    for (j = 0; j < xyh.size(); ++i, ++j)
        ans[i] = xyh[j];
    add(ans, xy, m, 0, xy.size(), base);

    return ans;
}

std::vector<unsigned> Integer::ksm (
    std::vector<unsigned> x, std::vector<unsigned> y,
    unsigned xn, unsigned yn, const unsigned base
) {
    if (xn < yn)
        return ksm(std::move(y), std::move(x), yn, xn, base);

    if (3 * yn > (xn << 1)) {
        y.resize(xn);
        return ksm(x, y, 0, 0, xn, base);
    }

    unsigned st = 0;
    x.resize(xn + (yn - xn % yn) % yn);
    std::vector<unsigned> ans(x.size() + yn), carry(yn);

    while (st < xn) {

        std::vector<unsigned> temp = ksm(x, y, st, 0, yn, base);
        add(temp, carry, 0, 0, yn, base);

        for (unsigned i = 0; i < yn; ++i, ++st)
            carry[i] = temp[i + yn], ans[st] = temp[i];
    }
    for (unsigned i = 0; i < yn; ++i, ++st)
        ans[st] = carry[i];

    return ans;
}

std::vector<unsigned> Integer::ksm (
    const std::vector<unsigned> &x, const std::vector<unsigned> &y,
    const unsigned base
) {
    const unsigned xn = x.size();
    const unsigned yn = y.size();

    if (xn == 0 || yn == 0)
        return {};
    if (yn == 1 && y.back() == 1)
        return x;
    if (xn == 1 && x.back() == 1)
        return y;

    std::vector<unsigned> ans;

    if (std::min(xn, yn) < ksmt)
        ans = mul(x, y, 0, xn, 0, yn, base);
    else
        ans = ksm(x, y, xn, yn, base);

    trim(ans);

    return ans;
}

void Integer::div (
    const std::vector<unsigned short> &dividend, unsigned divisor,
    std::vector<unsigned short> &quot, unsigned &rem
) {
    if (divisor == 0)
        throw std::runtime_error("ZERO DIVISION ERROR");

    quot.clear(), rem = 0;

    if (dividend.empty() || divisor == 1) {
        quot = dividend;
        return;
    }

    unsigned long long cur = 0;
    quot.resize(dividend.size());

    for (unsigned i = dividend.size() - 1; i < dividend.size(); --i) {
        cur = cur << bits | dividend[i];
        quot[i] = cur / divisor;
        cur %= divisor;
    }
    trim(quot);
    rem = cur;
}

void Integer::div (
    const std::vector<unsigned short> &dividend, const std::vector<unsigned short> &divisor,
    std::vector<unsigned short> &quot, std::vector<unsigned short> &rem
) {
    if (divisor.empty())
        throw std::runtime_error("ZERO DIVISION ERROR");

    quot.clear(), rem.clear();

    if (divisor.size() <= 2) {

        unsigned dv = divisor.back(), r;
        if (divisor.size() == 2)
            dv = dv << bits | divisor[divisor.size() - 2];

        div(dividend, dv, quot, r);

        if (r)
            rem.emplace_back(r), r >>= bits;
        if (r)
            rem.emplace_back(r);

        return;
    }

    if (vec_less(dividend, divisor)) {
        rem = dividend;
        return;
    }

    rem = dividend;
    std::vector<unsigned short> dvr = divisor;

    unsigned short norm = bits, db = dvr.back();
    while (db)
        --norm, db >>= 1;

    if (norm) {
        left_shift(rem, norm);
        left_shift(dvr, norm);
    }

    const unsigned m = rem.size();
    const unsigned n = dvr.size();

    quot.resize(m - n + 1);
    rem.emplace_back(0);

    for (unsigned i = m - n; i < m; --i) {

        const unsigned cur = (unsigned) rem[i + n] << bits | rem[i + n - 1];

        unsigned q = cur / dvr[n - 1];
        unsigned r = cur % dvr[n - 1];

        if (q >= base || q * dvr[n - 2] > (r << bits | rem[i + n - 2])) {
            --q, r += dvr[n - 1];
            if (r < base && q * dvr[n - 2] > (r << bits | rem[i + n - 2]))
                --q;
        }

        if (q == 0)
            continue;

        long long x = 0;

        for (unsigned j = 0; j < n; ++j)
            x += rem[i + j], x -= q * dvr[j],
                rem[i + j] = x, x >>= bits;

        x += rem[i + n], rem[i + n] = x, x >>= bits;

        if (x < 0) {

            x = 0, --q;

            for (unsigned j = 0; j < n; ++j)
                x += rem[i + j], x += dvr[j],
                    rem[i + j] = x, x >>= bits;

            x += rem[i + n], rem[i + n] = x, x >>= bits;
        }

        quot[i] = q;
    }

    trim(quot);
    trim(rem);

    if (norm)
        right_shift(rem, norm);
}

std::vector<unsigned short> Integer::pow (
    const std::vector<unsigned short> &n, unsigned p
) {
    std::vector<unsigned short> ans{1}, cur = n;

    while (p) {
        if (p & 1)
            ans = ksm(ans, cur);

        if (p >>= 1)
            cur = ksm_sqr(cur);
    }

    return ans;
}

std::vector<unsigned> Integer::pow (
    const std::vector<unsigned> &n, unsigned p, const unsigned base
) {
    std::vector<unsigned> ans{1}, cur = n;

    while (p) {

        if (p & 1)
            ans = ksm(ans, cur, base);

        if (p >>= 1)
            cur = ksm_sqr(cur, base);
    }

    return ans;
}

std::vector<unsigned> Integer::to (
    const std::vector<unsigned short> &binary, unsigned st, unsigned n,
    const unsigned base, const std::vector<unsigned>& base_vec,
    std::unordered_map<unsigned, std::vector<unsigned>>& exp_map
) {
    if (n <= 4) {

        unsigned long long cur = 0;

        for (unsigned i = n - 1; i < n; --i)
            cur = (cur << bits) | binary[st + i];

        std::vector<unsigned> ans;

        while (cur) {
            ans.emplace_back(cur % base);
            cur /= base;
        }

        return ans;
    }

    unsigned m = (n >> 1);

    std::vector<unsigned> lhs = to(binary, st, m, base, base_vec, exp_map);
    std::vector<unsigned> rhs = to(binary, st + m, n - m, base, base_vec, exp_map);

    std::unordered_map<unsigned, std::vector<unsigned>>::iterator p = exp_map.find(m);

    if (p == exp_map.end())
        exp_map[m] = pow(base_vec, m, base);

    rhs = ksm(rhs, exp_map[m], base);

    rhs.resize(std::max(lhs.size(), rhs.size()) + 1);

    add(rhs, lhs, 0, 0, lhs.size(), base);

    return rhs;
}

std::vector<unsigned> Integer::to (
    const std::vector<unsigned short> &binary, const unsigned base
) {
    std::vector<unsigned> base_vec;
    for (unsigned b = ::base; b; b /= base)
        base_vec.push_back(b % base);

    std::unordered_map<unsigned, std::vector<unsigned>> exp_map;

    std::vector<unsigned> decimal = to(
        binary, 0, binary.size(), base, base_vec,
        base == dmax ? dec_exp : exp_map
    );
    trim(decimal);

    return decimal;
}

std::vector<unsigned short> Integer::to (
    const std::vector<unsigned> &decimal, unsigned st, unsigned n,
    const unsigned base, const std::vector<unsigned short>& base_vec,
    std::unordered_map<unsigned, std::vector<unsigned short>>& exp_map
) {
    if (n <= 2) {

        unsigned long long cur = 0;

        for (unsigned i = n - 1; i < n; --i)
            cur = (cur * base) + decimal[st + i];

        std::vector<unsigned short> ans;

        while (cur) {
            ans.emplace_back(cur);
            cur >>= bits;
        }

        return ans;
    }

    unsigned m = (n >> 1);

    std::vector<unsigned short> lhs = to(decimal, st, m, base, base_vec, exp_map);
    std::vector<unsigned short> rhs = to(decimal, st + m, n - m, base, base_vec, exp_map);

    std::unordered_map<unsigned, std::vector<unsigned short>>::iterator p = exp_map.find(m);

    if (p == exp_map.end())
        exp_map[m] = pow(base_vec, m);

    rhs = ksm(rhs, exp_map[m]);

    rhs.resize(std::max(lhs.size(), rhs.size()) + 1);

    add(rhs, lhs, 0, 0, lhs.size());

    return rhs;
}

std::vector<unsigned short> Integer::to (
    const std::vector<unsigned> &decimal, const unsigned base
) {
    std::vector<unsigned short> base_vec;
    for (unsigned b = base; b; b >>= bits)
        base_vec.push_back(b);

    std::unordered_map<unsigned, std::vector<unsigned short>> exp_map;

    std::vector<unsigned short> binary = to(
        decimal, 0, decimal.size(), base, base_vec,
        base == dmax ? bin_exp : exp_map
    );
    trim(binary);

    return binary;
}

std::invalid_argument Integer::init_except (
    const std::string& what, const unsigned where, const std::string& in
) {
    std::string under_line = "\n\t\t  ";
    std::string except = "\tAT THE TIME OF INITIALIZATION\n";

    except.append("\t\t");
    except.append(what);
    except.append(" \"");

    for (unsigned i = 0; i < what.size(); ++i)
        under_line.push_back(' ');

    unsigned w = where;

    for (unsigned i = 0; i < 8 && w; ++i)
        --w, under_line.push_back(' ');

    if (where > 8)
        except.append(". . . "), under_line.append("      ");

    under_line[under_line.size() - 1] = '~';
    under_line[under_line.size() - 2] = '~';
    under_line[under_line.size() - 3] = '~';
    under_line.append("^~~~");

    for (unsigned i = 0; i < 20 && w < in.size(); ++i, ++w)
        except.push_back(in[w]);

    w -= where, w += 23;

    if (in.size() - where > 12)
        except.append(" . . ."), under_line.append("      ");

    except.append("\" ");
    except.append("WAS ENCOUNTERED AT POSITION ");
    except.append(ntos(where + 1));

    for (unsigned i = 0; i < w; ++i)
        under_line.push_back(' ');

    under_line.append("~~~^~~~\n");

    except.append(under_line);

    return std::invalid_argument(except);
}

std::string Integer::to_binp (
    const std::vector<unsigned short>& data, bool sign,
    const unsigned base_bits, bool prefix
) {
    std::string s;

    if (sign)
        s.push_back('-');

    if (prefix)
        s.push_back('0'),
        s.append(ntos(1 << base_bits)),
        s.append("...");

    if (data.empty()) {
        s.push_back('0');
        return s;
    }

    const unsigned mask = (1 << base_bits) - 1;
    long long unsigned cur = 0;
    unsigned cur_bits = 0;

    std::stack<unsigned> stk;

    for (const auto& bit : data) {
        cur |= (long long unsigned) bit << cur_bits;
        cur_bits += bits;
        while (cur_bits >= base_bits) {
            stk.push(cur & mask);
            cur_bits -= base_bits;
            cur >>= base_bits;
        }
    }
    stk.push(cur);

    while (stk.top() == 0) stk.pop();

    if (base_bits <= 6) {
        while (!stk.empty())
            s.push_back(cset[stk.top()]),
            stk.pop();
    } else {
        while (!stk.empty())
            s.append(ntos(stk.top())),
            s.push_back('\''),
            stk.pop();
        s.pop_back();
    }

    return s;
}

std::vector<unsigned short> Integer::from_binp (
    const std::string& s, const unsigned j, const unsigned base
) {

    long long unsigned cur = 0;
    std::vector<unsigned short> data;
    unsigned i = s.size() - 1, cur_bits = 0, v;
    unsigned base_bits = pop_count(base) - 1;

    if (base_bits <= 6) {

        while (i - j < s.size()) {

            v = getv(s[i]);

            if (v == -1)
                throw init_except("INAPPROPRIATE VALUE", i, s);

            if (base <= 10 && (v < 0 || v > 9))
                throw init_except("FOR BASE " + ntos(base) + " NON-INTEGER VALUE", i, s);

            if (v >= base)
                throw init_except("FOR BASE " + ntos(base) + " VERY LARGE BIT", i, s);

            cur |= v << cur_bits;
            --i, cur_bits += base_bits;

            if (cur_bits >= bits) {
                cur_bits -= bits;
                data.push_back(cur);
                cur >>= bits;
            }
        }
    } else {

        auto to = [&s, &j, &base] (unsigned& i) -> unsigned {

            unsigned n = 0, k = 0, v;

            while (i - j < s.size() && s[i] != '\'') {

                v = getv(s[i]);

                if (v == -1)
                    throw init_except("INAPPROPRIATE VALUE", i, s);

                if (v < 0 || v > 9)
                    throw init_except("FOR BASE " + ntos(base) + " (> 64) NON-INTEGER VALUE", i, s);

                n += v * p10[k];

                if (k > 7)
                    throw init_except("VERY LARGE BIT WIDTH (" + ntos(k + 1) + ")", i, s);

                if (n >= base)
                    throw init_except("BIT (" + ntos(n) + ") LARGER THAN BASE (" + ntos(base) + ")", i, s);

                ++k, --i;
            }
            return n;
        };
        while (i - j < s.size()) {

            cur |= to(i) << cur_bits;
            --i, cur_bits += base_bits;

            if (cur_bits >= bits) {
                cur_bits -= bits;
                data.push_back(cur);
                cur >>= bits;
            }
        }
    }
    data.push_back(cur);
    trim(data);

    return data;
}

Integer abs(const Integer& n) noexcept {
    return n.abs();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // UTILITY CLASS FOR EASY OUTPUT TO A FILE

class Fout {

    std::string file_name;
    std::ofstream out;

public:

    Fout(const Fout&);
    Fout(const std::string&);
    Fout& operator= (const Fout&);
    ~Fout();

    void close() noexcept;
    void open(const std::string&);
    explicit operator bool() const noexcept;
    bool good() const noexcept;
    void reset();

    friend Fout& operator<< (Fout&, const Integer&);
    friend Fout& operator<< (Fout&, const char);
    friend Fout& operator<< (Fout&, const char*);
    friend Fout& operator<< (Fout&, const std::string&);
    friend Fout& operator<< (Fout&, const int);
    friend Fout& operator<< (Fout&, const long long);
    friend Fout& operator<< (Fout&, const unsigned);
    friend Fout& operator<< (Fout&, const long long unsigned);
    friend Fout& operator<< (Fout&, const long double);
    friend Fout& operator<< (Fout&, const Div&);
    friend Fout& fout(const char*);
};

Fout& operator<< (Fout&, const int);
Fout& operator<< (Fout&, const Integer&);
Fout& operator<< (Fout&, const char);
Fout& operator<< (Fout&, const char*);
Fout& operator<< (Fout&, const Div&);
Fout& operator<< (Fout&, const unsigned);
Fout& operator<< (Fout&, const long long);
Fout& operator<< (Fout&, const long double);
Fout& operator<< (Fout&, const std::string&);
Fout& operator<< (Fout&, const long long unsigned);
Fout& fout(const char*);

// ++++++++++++++++++++++++++++++++++++++++++++++++++

Fout::Fout(const std::string& file_name)
    : file_name(file_name), out(std::ofstream(this->file_name))
{}

Fout::Fout(const Fout& fout)
    : file_name(fout.file_name), out(std::ofstream(file_name))
{}

Fout& Fout::operator= (const Fout& fout) {
    if(&fout != this) {
        file_name = fout.file_name;
        close();
        out = std::ofstream(file_name);
    }
    return *this;
}

Fout::~Fout() {
    close();
}

void Fout::open(const std::string& s) {
    close();
    out.open(s);
    file_name = s;
}

void Fout::close() noexcept {
    if (*this)
        out.close();
}

Fout::operator bool() const noexcept {
    return out.is_open();
}

bool Fout::good() const noexcept {
    return bool(out);
}

void Fout::reset() {
    close();
    out.open(file_name);
}

Fout& operator<< (Fout& out, const Integer& n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const char c) {
    if (!out.good()) out.reset();
    out.out << c;
    return out;
}

Fout& operator<< (Fout& out, const char* s) {
    if (!out.good()) out.reset();
    out.out << s;
    return out;
}

Fout& operator<< (Fout& out, const std::string& s) {
    if (!out.good()) out.reset();
    out.out << s;
    return out;
}

Fout& operator<< (Fout& out, const int n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const long long n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const unsigned n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const long long unsigned n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const long double n) {
    if (!out.good()) out.reset();
    out.out << n;
    return out;
}

Fout& operator<< (Fout& out, const Div& qr) {
    if (!out.good()) out.reset();
    out.out << qr;
    return out;
}

Fout& fout(const char* s) {
    static Fout f = Fout(s);
    f.close();
    f.open(s);
    if (!f)
        throw std::runtime_error(
            "UNABLE TO WRITE IN '" + std::string(s) + "'..."
        );
    return f;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // UTILITY CLASS FOR EASY INPUT FROM A FILE

class Fin {

    std::string file_name;
    std::ifstream in;

public:

    Fin(const Fin&);
    Fin(const std::string&);
    Fin& operator= (const Fin&);
    ~Fin();

    void close() noexcept;
    void open(const std::string&);
    explicit operator bool() const noexcept;
    std::string read();
    std::string& read(std::string&);
    bool good() const noexcept;
    void reset();

    friend Fin& operator>> (Fin&, Integer&);
    friend Fin& operator>> (Fin&, char&);
    friend Fin& operator>> (Fin&, char*);
    friend Fin& operator>> (Fin&, std::string&);
    friend Fin& operator>> (Fin&, int&);
    friend Fin& operator>> (Fin&, long long&);
    friend Fin& operator>> (Fin&, unsigned&);
    friend Fin& operator>> (Fin&, long long unsigned&);
    friend Fin& operator>> (Fin&, long double&);
    friend Fin& fin(const char*);
};

Fin& operator>> (Fin&, int&);
Fin& operator>> (Fin&, Integer&);
Fin& operator>> (Fin&, char&);
Fin& operator>> (Fin&, char*);
Fin& operator>> (Fin&, unsigned&);
Fin& operator>> (Fin&, long long&);
Fin& operator>> (Fin&, std::string&);
Fin& operator>> (Fin&, long double&);
Fin& operator>> (Fin&, long long unsigned&);
Fin& fin(const char*);

// ++++++++++++++++++++++++++++++++++++++++++++++++++

Fin::Fin(const std::string& file_name)
    : file_name(file_name), in(std::ifstream(this->file_name))
{}

Fin::Fin(const Fin& fout)
    : file_name(fout.file_name), in(std::ifstream(file_name))
{}

Fin& Fin::operator= (const Fin& fin) {
    if(&fin != this) {
        file_name = fin.file_name;
        close();
        in = std::ifstream(file_name);
    }
    return *this;
}

Fin::~Fin() {
    close();
}

void Fin::open(const std::string& s) {
    close();
    in.open(s);
    file_name = s;
}

void Fin::close() noexcept {
    if (*this)
        in.close();
}

Fin::operator bool() const noexcept {
    return in.is_open();
}

bool Fin::good() const noexcept {
    return bool(in);
}

std::string Fin::read() {
    std::string s;
    return read(s);
}

std::string& Fin::read(std::string& s) {
    reset();
    while (in)
        s.push_back(in.get());
    return s;
}

void Fin::reset() {
    close();
    in.open(file_name);
}

Fin& operator>> (Fin& in, Integer& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& operator>> (Fin& in, char& c) {
    if (!in.good()) in.reset();
    in.in >> c;
    return in;
}

Fin& operator>> (Fin& in, char *s) {
    if (!in.good()) in.reset();

    int pos = in.in.tellg();

    in.in.seekg(0, std::ios::end);
    int size = in.in.tellg();

    in.in.seekg(pos);
    in.in.getline(s, size, ' ');

    return in;
}

Fin& operator>> (Fin& in, std::string& s) {
    if (!in.good()) in.reset();
    in.in >> s;
    return in;
}

Fin& operator>> (Fin& in, int& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& operator>> (Fin& in, long long& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& operator>> (Fin& in, unsigned& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& operator>> (Fin& in, long long unsigned& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& operator>> (Fin& in, long double& n) {
    if (!in.good()) in.reset();
    in.in >> n;
    return in;
}

Fin& fin(const char* s) {
    static Fin f = Fin(s);
    f.close();
    f.open(s);
    if (!f)
        throw std::runtime_error(
            "UNABLE TO READ FROM '" + std::string(s) + "'..."
        );
    return f;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif