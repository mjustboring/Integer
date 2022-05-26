#include <iostream>
#include "Integer.hpp"
using namespace std;

void multiply() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " * " << b << " = " << a * b << "\n\n";
}

void divide() {
    Integer a = rand(100);
    Integer b = rand(50);

    cout << "div(" << a << ", " << b << ") = " << "\n";
    cout << div(a, b) << "\n";

    cout << a << " / " << b << " = " << a / b << "\n\n";
}

void modulo() {
    Integer a = rand(100);
    Integer b = rand(50);
    cout << a << " % " << b << " = " << a % b << "\n\n";
}

void addition() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " + " << b << " = " << a + b << "\n\n";
}

void subtraction() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " - " << b << " = " << a - b << "\n\n";
}

void bitwise_and() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " & " << b << " = " << (a & b) << "\n\n";
}

void bitwise_or() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " | " << b << " = " << (a | b) << "\n\n";
}

void bitwise_xor() {
    Integer a = rand(100);
    Integer b = rand(100);
    cout << a << " ^ " << b << " = " << (a ^ b) << "\n\n";
}

void bitwise_compliment() {
    Integer a = rand(100);
    cout << "~" << a << " = " << ~a << "\n\n";
}

void math_util() {

    Integer n = rand(100);

    cout << "sqr(" << n << ") = " << n.sqr() << "\n\n";
    cout << "nrt(" << n << ", 5) = " << n.nrt(5) << "\n\n";
    cout << "sqrt(" << n << ") = " << n.sqrt() << "\n\n";
    cout << "pow(" << n << ", 3) = " << n.pow(3) << "\n\n";
    cout << "pow(" << n << ", 1/5) = " << n.pow(1, 5) << "\n\n";
    cout << "pow(" << n << ", 0.2) = " << n.powf(0.2) << "\n\n";

    Integer a = rand(100);
    Integer b = rand(100);

    cout << "gcd(" << a << ", " << b << ") = " << a.gcd(b) << "\n\n";
    cout << "lcm(" << a << ", " << b << ") = " << a.lcm(b) << "\n\n";

    cout << "odd(" << n << ") = " << n.odd() << "\n";
    cout << "even(" << n << ") = " << n.even() << "\n";

    cout << "set_bit(0, 64) = " << (n = Integer().set_bit(64)) << "\n";
    cout << "unset_bit(" << n << ", 64) = " << n.unset_bit(64) << "\n";

    n = a;

    cout << "log2(" << n << ") = " << n.log2() << "\n";
    cout << "log10(" << n << ") = " << n.log10() << "\n";
    cout << "loge(" << n << ") = " << n.loge() << "\n";
    cout << "log(" << n << ", 123) = " << n.log(123) << "\n\n";
    cout << "fact(100) = " << fact(100) << "\n\n";

    cout << "rand(100) = " << rand(100) << "\n";
    cout << "rand_bit(100) = " << rand_bit(100) << "\n";
    cout << "rand_byte(100) = " << rand_byte(100) << "\n";
    cout << "rand_range(70, 100) /* By number of digits */ = " << rand_range(70, 100) << "\n";
    cout << "rand_range(" << (a = rand(70)) << ", " << (b = rand(100)) << ") = " << rand_range(a, b) << "\n\n";
}

void utility() {

    Integer a = rand(100);
    Integer b = rand(100);

    cout << "bit_count(" << a << ") = " << a.bit_count() << "\n\n";
    cout << "size(" << a << ") = " << a.size() << "\n\n";

    cout << "Before swapping" << "\n";
    cout << "A = " << a << "\n";
    cout << "B = " << b << "\n\n";
    a.swap(b);
    cout << "After swapping" << "\n";
    cout << "A = " << a << "\n";
    cout << "B = " << b << "\n\n";
}

void radix() {

    Integer n = rand(100);

    cout << "dec(" << n << ") = " << n.dec(1) << "\n\n";
    cout << "bin(" << n << ") = " << n.bin(1) << "\n\n";
    cout << "oct(" << n << ") = " << n.oct(1) << "\n\n";
    cout << "hex(" << n << ") = " << n.hex(1) << "\n\n";
    cout << "str(" << n << ", 123) = " << n.str(123, 1) << "\n\n";
    cout << "str(" << n << ", 64) = " << n.str(64, 1) << "\n\n";
}

int main() {

    multiply();
    divide();
    modulo();
    addition();
    subtraction();
    bitwise_and();
    bitwise_or();
    bitwise_xor();
    bitwise_compliment();
    math_util();
    utility();
    radix();
}
