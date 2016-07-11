/*
 * Author: Becker Awqatty
 * 
 * Licensed under MIT Open-Source License
 */


///////////// Definitions ///////////////
#include "numeric_mod.h"
#include <limits>
#include <stdexcept>

//template <typename T>
//T gcd(T a, T b) {
//	T *p_a=&a, *p_A=&b;
//	while (*p_a != 0) {
//		(*p_A) %= (*p_a);
//		swap(p_A, p_a);
//	}
//	return *p_A;
//}

template <typename M_t, M_t MOD>
template <typename N_t>
inline M_t numeric_mod<M_t,MOD>::modulus(N_t x) {
	if (!std::numeric_limits<N_t>::is_signed) {
		return x%MOD;
	} else {
		N_t ret_val = x%MOD;
		if (ret_val < 0)
			ret_val+=MOD;
		return ret_val;
	}
}
	
/* Returns the integer d that satisfies
	xd = 1 (mod m), or d = (1/x mod m)
 -> xd = km+1 for some k
 - find a k such that x | km+1
	-> (gcd(x,m) must be 1)
	-> k = -1/m (mod x)
		= x-(1/m mod x)
	-> (recursively implements function)

 -> d = (km+1)/x

template <typename M_t, M_t MOD>
M_t numeric_mod<M_t,MOD>::inv_rec(M_t x, M_t m) {
	// (Inductively guarantees return <= m)
	if (x==1 || x==0)
		return x;
	const M_t q=(m/x),r=(m%x);
	const M_t k = x-inv_rec(r,x);
		// k != 0 (if original x has an inverse)
	return k*q + (k*r+1)/x;
}
 */

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>::numeric_mod() {}

template <typename M_t, M_t MOD>
template <typename N_t>
numeric_mod<M_t,MOD>::numeric_mod(N_t n) : number(modulus(n)) {}

template <typename M_t, M_t MOD, typename N_t>
static inline numeric_mod<M_t, MOD> make_w_o_value_check(N_t x) {
	numeric_mod<M_t, MOD> ret;
	ret.number = x;
	return ret;
}

template <typename M_t, M_t MOD>
template <typename N_t>
numeric_mod<M_t,MOD>::operator N_t() const {
	return N_t(number);
}
	
// Returns the value congruent to (*this) within the range
//	[offset, MOD+offset)
template <typename M_t, M_t MOD>
template <typename N_t>
N_t numeric_mod<M_t,MOD>::offset_by(N_t offset) const {
	return N_t(numeric_mod<M_t,MOD>(*this)-=offset) + offset;
}

// (In/De)crenemt
template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator++() {
	if (number == MOD-1)
		number = 0;
	else ++number;
	return *this;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD> numeric_mod<M_t,MOD>::operator++(int) {
	numeric_mod<M_t,MOD> tmp(*this);
	operator++();
	return tmp;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator--() {
	if (number == 0)
		number = MOD-1;
	else --number;
	return *this;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD> numeric_mod<M_t,MOD>::operator--(int) {
	numeric_mod<M_t,MOD> tmp(*this);
	operator--();
	return tmp;
}

// Arithmetic Assignment & Inverses
template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator+=(
		const numeric_mod<M_t,MOD>& arg) {
	if (number >= MOD-arg.number)
		number -= MOD-arg.number;
	else
		number += arg.number;
	return *this;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator-=(
		const numeric_mod<M_t,MOD>& arg) {
	if (number < arg.number)
		number += MOD-arg.number;
	else
		number -= arg.number;
	return *this;
}

template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> numeric_mod<M_t,MOD>::operator-() const {
	return MOD-number;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator*=(
		const numeric_mod<M_t,MOD>& arg) {
	// Shortcut situations
	if (number==0 || arg.number==1) return *this;
	if (arg.number==0 || number==1) {
		number = arg.number;
		return *this;
	}
	
	// Schrage's Method
	const M_t q = MOD/number, r = MOD%number;
	number *= arg.number%q;
	operator -= (((r<=q) ? r : numeric_mod<M_t,MOD>(r)) *= (arg.number/q));
	
	return *this;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD>& numeric_mod<M_t,MOD>::operator/=(
		const numeric_mod<M_t,MOD>& arg) {
	// Faster w/o shortcut check for arg.number==1
	const M_t q=MOD/arg.number, r=MOD%arg.number;
	std::pair<M_t,M_t> k;
	
	Alg_FactorDiff alg;
	alg.run(number, arg.number,r, k);
	
	number = k.first;
	operator+=(k.second*q);
		// kr = -1/m modx -> k<x -> kq < xq < m
	operator*=(alg.factor_is_pos()
		? alg.factor()
		: -numeric_mod<M_t,MOD>(alg.factor()) );
	
	return *this;
}

template <typename M_t, M_t MOD>
void numeric_mod<M_t,MOD>::Alg_FactorDiff::run(
		M_t d, M_t x,M_t y,
			// input
		std::pair<M_t,M_t> &k) {
			// output
	is_pos=true;
	diff=d;
	alg_rec(x,y, k.first,k.second);
}
template <typename M_t, M_t MOD>
void numeric_mod<M_t,MOD>::Alg_FactorDiff::alg_rec(
		M_t x,M_t y, M_t &k1,M_t &k2) {
	if (y==0) {
		if (diff%x != 0) {
			// Undefined if gcd(denominator,mod) does not divide numerator
			// Excludes case denom=0 (fails earlier)
			throw std::runtime_error("modulus quotient does not exist");
		}
		diff/=x;
			// factor is stored & applied once at very end
			// (avoids overflow for numeric_limits<M_t>::max==MOD)
		k1=1;
		k2=0;
			// k1*x - k2*y = diff
			// -> k1 = 1*(diff/x)
			// -> k2 undefined -> k2=0 for smallest result
	}
	else {
		is_pos = !is_pos;
		alg_rec(y,x%y, k2,k1);
		k2 += k1*(x/y);
	}
}
template <typename M_t, M_t MOD>
M_t numeric_mod<M_t,MOD>::Alg_FactorDiff::factor() const {
	return diff;
}
template <typename M_t, M_t MOD>
bool numeric_mod<M_t,MOD>::Alg_FactorDiff::factor_is_pos() const {
	return is_pos;
}

template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD> numeric_mod<M_t,MOD>::inv() const {
	return make_w_o_value_check<>(1) /= (*this);
}

// Powers
/*
template <typename M_t, M_t MOD, typename N_t>
inline numeric_mod<M_t,MOD> to_pow(const N_t& arg) {
	if (arg<0) {
		operator=(inv());
		return to_pow(numeric_mod<M_t,MOD-1>(-arg));
	}
	else return to_pow(numeric_mod<M_t,MOD-1>(arg));
}
template <typename M_t, M_t MOD>
numeric_mod<M_t,MOD> to_pow(const numeric_mod<M_t,MOD-1>& arg) {
	// TODO
}
*/
	

// Non-Member Operations

// Addition
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) += arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2) {
	return numeric_mod<M_t,MOD>(arg2) += arg1;
		//switched to avoid 1 implicit constructor call
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) += arg2;
}

// Subtraction
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) -= arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2) {
	return numeric_mod<M_t,MOD>(arg1) -= arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) -= arg2;
}

// Multiplication
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) *= arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2) {
	return numeric_mod<M_t,MOD>(arg2) *= arg1;
		//switched to avoid 1 implicit constructor call
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) *= arg2;
}

// Division
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) /= arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2) {
	return numeric_mod<M_t,MOD>(arg1) /= arg2;
}
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2) {
	return numeric_mod<M_t,MOD>(arg1) /= arg2;
}





// Testing Functions

#include <iostream>
#include <ctime>
using namespace std;

template <typename M>
void test_table(M (*func)(M a, M b), const M mod) {
	cout << "a  \\  b" << "\t";
	for (M i=0; i<mod; ++i) {
		cout << "|" << i << "\t";
	}
	cout << endl;
	for (M i=0; i<=mod; ++i) {
		cout << "--------";
	}
	cout << endl;
	
	for (M i=0; i<mod; ++i) {
		cout << i << "\t";
		for (M j=0; j<mod; ++j) {
			cout << "|" << (*func)(i,j) << "\t";
		}
		cout << endl;
	}
}

void test_time(
		void (*func)(unsigned long int i),
		const unsigned long int reps) {
	clock_t start, end;
	
	start = clock();
	for (unsigned int i=1; i<reps; ++i) {
		(*func)(i);
	}
	end = clock();
	
	cout << "clock cycles elapsed: " << end-start << endl;
}


const M1 mod=7;
const M2 mod2=1000000007;

const unsigned long int REPS=1000000;

// Test Addition
M1 mod_add(M1 arg1, M1 arg2) {
	return M1(numeric_mod<M1,mod>(arg1)+=arg2);
}
void test_add_table() {
	test_table(&mod_add, mod);
}

// Test Subtraction
M1 mod_subt(M1 arg1, M1 arg2) {
	return M1(numeric_mod<M1,mod>(arg1)-=arg2);
}
void test_subt_table() {
	test_table(&mod_subt, mod);
}

// Test Multiplication
M1 mod_mult(M1 arg1, M1 arg2) {
	return M1(numeric_mod<M1,mod>(arg1)*=arg2);
}
void test_mult_table() {
	test_table(&mod_mult, mod);
}

void mod_mult_time(unsigned long int i) {
	static numeric_mod<M2,mod2> m2=1;
	m2 *= (M2)i;
}
void test_mult_time() {
	test_time(&mod_mult_time, REPS);
}

// Test Division
M1 mod_div(M1 arg1, M1 arg2) {
	return ((arg2==0) ? mod : M1(numeric_mod<M1,mod>(arg1)/=arg2));
}
void test_div_table() {
	test_table(&mod_div, mod);
}

void mod_div_time(unsigned long int i) {
	static numeric_mod<M2,mod2> m2=1;
	m2 /= (M2)i;
}
void test_div_time() {
	test_time(&mod_div_time, REPS);
}
