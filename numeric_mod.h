/*
 * Author: Becker Awqatty
 * 
 * Licensed under MIT Open-Source License
 */


#ifndef NUMERIC_MOD_H
#define NUMERIC_MOD_H

///////////// Declarations ///////////////
#include <utility>

template <typename M_t, M_t MOD>
class numeric_mod {
private:
	M_t number;
	
	template <typename N_t>
	static inline M_t modulus(N_t x);

	template <typename N_t>
	static inline numeric_mod<M_t, MOD> make_w_o_value_check(N_t x);
	
	//static M_t inv_rec(M_t x, M_t m);
	// k1*x - k2*y = d, x>y
	class Alg_FactorDiff {
	private:
		M_t diff;
		bool is_pos;
		void alg_rec(M_t x,M_t y, M_t &k1,M_t &k2);
		
	public:
		void run(
				M_t d, M_t x,M_t y,
					// input
				std::pair<M_t,M_t> &k);
					// output
		M_t factor() const;
		bool factor_is_pos() const;
	};
	
public:
	numeric_mod();
		template <typename N_t>
	numeric_mod(N_t n);
		template <typename N_t>
	explicit operator N_t() const;
	
	// Returns the value congruent to (*this) within the range
	//	[offset, MOD+offset)
		template <typename N_t>
	N_t offset_by(N_t offset) const;
	
	// (In/De)crenemt
	numeric_mod<M_t,MOD>& operator++();
	numeric_mod<M_t,MOD> operator++(int);
	numeric_mod<M_t,MOD>& operator--();
	numeric_mod<M_t,MOD> operator--(int);
	
	// Arithmetic Assignment & Inverses
	numeric_mod<M_t,MOD>& operator+=(const numeric_mod<M_t,MOD>& arg);
	numeric_mod<M_t,MOD>& operator-=(const numeric_mod<M_t,MOD>& arg);
	numeric_mod<M_t,MOD> operator-() const;
	
	numeric_mod<M_t,MOD>& operator*=(const numeric_mod<M_t,MOD>& arg);
	numeric_mod<M_t,MOD>& operator/=(const numeric_mod<M_t,MOD>& arg);
	numeric_mod<M_t,MOD> inv() const;
	
	//	template <typename N_t>
	//numeric_mod<M_t,MOD>& to_pow(const N_t& arg);
	//	template <>
	//numeric_mod<M_t,MOD>& to_pow(const numeric_mod<M_t,MOD-1>& arg);
};

// Addition
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator+(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2);

// Subtraction
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator-(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2);

// Multiplication
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator*(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2);

// Division
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD>& arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		const numeric_mod<M_t,MOD>& arg1,
		M_t arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<M_t,MOD> operator/(
		M_t arg1,
		const numeric_mod<M_t,MOD>& arg2);

// Powers
/*
template <typename M_t, M_t MOD, typename N_t>
inline numeric_mod<>M_t,MOD> pow(
		const numeric_mod<M_t,MOD>& arg1,
		N_t arg2);
template <typename M_t, M_t MOD>
inline numeric_mod<>M_t,MOD> pow(
		const numeric_mod<M_t,MOD>& arg1,
		const numeric_mod<M_t,MOD-1>& arg2);		
*/


// Testing Functions
typedef unsigned short M1;
typedef unsigned long long int M2;

template <typename M>
void test_table(M (*func)(M a, M b), const M mod);
void test_time(
		void (*func)(unsigned long int i),
		const unsigned long int reps);

// Test Addition
M1 mod_add(M1 arg1, M1 arg2);
void test_add_table();

// Test Subtraction
M1 mod_subt(M1 arg1, M1 arg2);
void test_subt_table();

// Test Multiplication
M1 mod_mult(M1 arg1, M1 arg2);
void test_mult_table();
void mod_mult_assg_time(unsigned long int i);
void test_mult_assg_time();

// Test Division
M1 mod_div(M1 arg1, M1 arg2);
void test_div_table();
void mod_div_time(unsigned long int i);
void test_div_time();


#endif
