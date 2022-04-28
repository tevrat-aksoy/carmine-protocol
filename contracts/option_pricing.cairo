# Contract that calculates Black-Scholes with Choudhury's approximation to std normal CDF
# https://www.hrpub.org/download/20140305/MS7-13401470.pdf.

%lang starknet

from starkware.cairo.common.cairo_builtins import HashBuiltin
from starkware.cairo.common.math import sign

from contracts.Math64x61 import (
    Math64x61_fromFelt,
    Math64x61_toFelt,
    Math64x61_exp,
    Math64x61_sqrt,
    Math64x61_mul,
    Math64x61_div,
    Math64x61_add
)


func _const{syscall_ptr : felt*, pedersen_ptr : HashBuiltin*, range_check_ptr}(x: felt) -> (
    res : felt
):
    let (cons) = Math64x61_fromFelt(x)
    let (thousand) = Math64x61_fromFelt(1000)
    let (res) = Math64x61_div(cons, thousand)
    return (res)
end

# Calculates approximate value of standard normal CDF
@view
func std_normal_cdf{syscall_ptr : felt*, pedersen_ptr : HashBuiltin*, range_check_ptr}(x: felt) -> (
    res : felt
):
    alloc_locals

    # TODO: rename variables so that they are reasonably named
    # TODO: extract repeated code

    let (million) = Math64x61_fromFelt(10**6)
    let (MINUS_ONE) = Math64x61_fromFelt(-1)

    let (sign_value) = sign(x)
    if sign_value == -1:
        let (symmetric_value) = std_normal_cdf(-x)
        let (symmetric_value_) = Math64x61_fromFelt(symmetric_value)
        let (neg_symmetric_value_) = Math64x61_mul(symmetric_value_, MINUS_ONE)
        let (res) = Math64x61_add(million, neg_symmetric_value_)
        let (res_) = Math64x61_toFelt(res)
        return (res_)
    end

    # Can't have it outside, since it requires "range_check_ptr"
    # TODO: how to define the constants reasonably?... throws an error for "const PI_a", same for "const (PI_a)"
    # 3.14159265358979
    let (PI_a) = Math64x61_fromFelt(314159265358979)
    let (PI_b) = Math64x61_fromFelt(10**14)
    let (PI) = Math64x61_div(PI_a, PI_b)

    # TODO: check if the ONE, TWO,... have to be defined like this
    # TODO: take the constants out of this function... how? Math64x61_fromFelt requires range_check_ptr
    let (ONE) = Math64x61_fromFelt(1)
    let (TWO) = Math64x61_fromFelt(2)
    let (THREE) = Math64x61_fromFelt(3)
    let (two_pi) = Math64x61_mul(TWO, PI)
    let (root_of_two_pi) = Math64x61_sqrt(two_pi)
    let (inv_root_of_two_pi) = Math64x61_div(ONE, root_of_two_pi)
    let (MINUS_ONE) = Math64x61_fromFelt(-1)
    let (CDF_CONST) = Math64x61_mul(MINUS_ONE, inv_root_of_two_pi)
    let (const_a) = _const(226)
    let (const_b) = _const(640)
    let (const_c) = _const(330)

    let (x_) = Math64x61_fromFelt(x)
    let (x_squared) = Math64x61_mul(x_, x_)
    let (x_squared_half) = Math64x61_div(x_squared, TWO)
    let (minus_x_squared_half) = Math64x61_mul(MINUS_ONE, x_squared_half)
    let (numerator) = Math64x61_exp(minus_x_squared_half)

    let (denominator_b) = Math64x61_mul(const_b, x_)
    let (denominator_a) = Math64x61_add(const_a, denominator_b)
    let (sqrt_den_part) = Math64x61_sqrt(x_squared + THREE)
    let (denominator_c) = Math64x61_mul(const_c, sqrt_den_part)
    let (denominator) = Math64x61_add(denominator_a, denominator_c)

    let (res_a) = Math64x61_div(numerator, denominator)
    let (res_b) = Math64x61_mul(CDF_CONST, res_a)
    let (res) = Math64x61_add(ONE, res_b)

    let (res_a) = Math64x61_mul(res, million)
    let (res_ad) = Math64x61_toFelt(res_a)
    return (res_ad)
end

# Calculates approximate value for Black Scholes
# TBD
@view
func black_scholes{syscall_ptr : felt*, pedersen_ptr : HashBuiltin*, range_check_ptr}(x: felt) -> (
    res : felt
):
    let res = 4 * 4
    return (res)
end
