[![CI](https://github.com/vargaconsulting/tiny-crypto/actions/workflows/ci.yml/badge.svg)](https://github.com/vargaconsulting/tiny-crypto/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/vargaconsulting/tiny-crypto/branch/main/graph/badge.svg)](https://codecov.io/gh/vargaconsulting/tiny-crypto)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/950847209.svg)](https://doi.org/10.5281/zenodo.15492419)
[![GitHub release](https://img.shields.io/github/v/release/vargaconsulting/tiny-crypto.svg)](https://github.com/vargaconsulting/tiny-crypto/releases)
[![Documentation](https://img.shields.io/badge/docs-stable-blue)](https://vargaconsulting.github.io/tiny-crypto)

# Tiny Crypto – Exploring Cryptography with Small Prime Fields

**A small-scale cryptography playground using tiny prime fields for easy manual verification.**

## Overview

TinyCrypto is an educational project that implements fundamental cryptographic algorithms using **small prime fields**. The goal is to **simplify the math** and make manual verification feasible, allowing learners to focus on **concepts** without getting lost in large numbers.

The project includes:
- **Elliptic Curve Cryptography (ECC)** – Weierstrass curves, point arithmetic, and scalar multiplication.
- **Finite Field Arithmetic** – Type-safe `𝔽ₚ` implementation with modular ops and field-specific overloads.
- **Shamir's Secret Sharing** – Securely split and reconstruct secrets.
- **Distributed Key Management** – Collaborative cryptographic protocols.
- **Other Cryptographic Routines** – Prime testing, modular inverse, quadratic residue checks, and more.

## Features

- Strong **finite field typing** via `Fp{T,P}`, ensuring safe and consistent modular arithmetic.
- **Elliptic curve operations** over tiny prime fields:
  - Curve generation over `𝔽ₚ` for small primes
  - Point addition, doubling, scalar multiplication
  - Generator detection and validation
  - Point-at-infinity and AffinePoint construction
- Clean ANSI-formatted output for symbolic math (e.g., `𝔽₃₁`, `𝔾(6,2)`)
- Fully written in **Julia**, for clarity and performance.
- **Educational focus** – Ideal for teaching and cryptographic exploration.

## Example Usage

```julia
using TinyCrypto
import TinyCrypto: is_identity

# Find a suitable Weierstrass curve over small field primes and coefficient ranges with iterative method
# general syntax: curve-name(prime-field-range, parameter-range₁, ..., parameter-rangeₙ, max_cofactor=8)
curve = Weierstrass(97:103, 10:15, 2:7)  # (prime range, a range, b range)
# Output: Weierstrass{𝔽₉₇}: y² = x³ + 10x + 3 | 𝔾(0,10), q = 101, h = 1, #E = 10

E = curve_points(curve) # Get all curve points
# → 101-element Vector{AffinePoint₁₂₈}: (0𝔽₉₇,10𝔽₉₇), (0𝔽₉₇,87𝔽₉₇), ..., (96𝔽₉₇,63𝔽₉₇), (∞,∞)
S = subgroup_points(curve) # same as `curve_points` as #E(𝔽₉₇) ∈ primes

curve = Montgomery(30:40, 8:40, 3:40)
# Montgomery{𝔽₃₇}: 8y² = x³ + 3x² + x | 𝔾(15,2), q = 11, h = 4, #E = 44
S = subgroup_points(curve) # when cofactor greater than one, there are multiple sub groups
E = curve_points(curve)

curve = TwistedEdwards(50:100, 1:20, 1:20)
# TwistedEdwards{𝔽₉₇}: 3x² + y² = 1 + 10x²y² | 𝔾(48,93), q = 29, h = 4, #E = 116

curve = Edwards(50:100, 1:20, max_cofactor=8)
# Edwards{𝔽₇₉}: 1x² + y² = 1 + 6x²y² | 𝔾(7,58), q = 23, h = 4, #E = 92
is_singular(curve) # false

## Point arithmetic on curve
𝔾 = curve.G                          # (7𝔽₇₉,58𝔽₇₉) ∈ Edwards{Fp{UInt128, 79}}
𝔾 + 𝔾 + 𝔾 + 𝔾                        # (68𝔽₇₉,67𝔽₇₉)  point addition on a cyclic group
2𝔾                                   # (31𝔽₇₉,51𝔽₇₉) 
2𝔾 == 𝔾 + 𝔾                          # true

𝒪 = identity(curve)                   # identity point on Edwards curve variant: (0𝔽₇₉,1𝔽₇₉)
𝒪 + 𝔾 == 𝔾                            # true
𝔾 + 𝒪 == 𝔾                            # true  

# direct construction from paramters:
const 𝔽₃₁ = 𝔽ₚ{UInt8, 31}                   # define field prime 𝔽ₚ for the abelian group, or use non-unicode `Fp{base_type, prime_number}`
Weierstrass{𝔽₃₁}(6, 9, 37, 1, (0,3))        # Weierstrass{𝔽₃₁}: y² = x³ + 6x + 9 | 𝔾(0,3), q = 37, h = 1, #E = 37
# or pass the 𝔽ₚ field prime directly as in:
curve = Weierstrass(31, 6, 9, 37, 1, (0,3)) # Weierstrass{𝔽₃₁}: y² = x³ + 6x + 9 | 𝔾(0,3), q = 37, h = 1, #E = 37

𝒪 = identity(curve)                    # (∞,∞) Weierstrass curve identity is infinity
is_identity(𝒪)                         # true 
∞ = infinity(curve)                    # (∞,∞) as expected
is_infinity(∞)                         # true

E = curve_points(curve)                # 
is_point_on_curve(E[1], curve)         # true 
is_point_on_curve(curve.G, curve)      # true, G generator point defines the cyclic group, which in this case the entire abelien group 



## Tiny Hash (not for cryptographic use, for obvious reasons)
H("byte size hash of a string")        # hashes string to a byte, given the abelian group is byte size
H₁₆("16 bit hash")                     # in case you need more hash space
H₈("is same") == H("is same")          # true, it is just a alias, smae as H8 
```

## Installation
```
using Pkg
Pkg.add(url="https://github.com/vargaconsulting/tiny-crypto.git")
```
## Development
```
using Pkg
Pkg.develop(url="https://github.com/vargaconsulting/tiny-crypto.git")
```
Please follow these steps:
1. create new issue ie: `gh issue create `
2. branch on top of release tag: `gh issue develop your-new-isse`
3. for commits follow pattern: `[#14]:svarga:feature, short description of your change`
4. rebase to latest release tag and manage conflicts if any 
5. upload changes: `git push -u origin 14-svarga-feature-short-name`
6. submit pull/merge request `gh pr create --fill`
7. We'll review your pull request and coordinate merging. Thanks for contributing!

## Future Work
- **Typed finite field system** – `𝔽ₚ{T,P}` representation for safe modular arithmetic
- **Symbolic math output** – pretty-printing with `𝔽ₚ`, `𝔾(x,y)` and Unicode formatting
- **Robust EC point operations** – scalar multiplication, inversion, curve validation
- **Support for popular blockchain curves** – `secp256k1`, `ed25519`, `altbn128`, ...
- **Full signature schemes** – ECDSA, EdDSA, Schnorr signatures
- **Shamir’s Secret Sharing** – secret splitting & reconstruction
- **Proactive refresh & threshold MPC** – secure key management with distributed trust
- **Interactive playground** – simulate distributed signing protocols, visually

## Licence 
MIT License. Copyright © 2025 Steven Varga, Varga Consulting, Toronto, ON, Canada 🍁
For questions, reach out at info@vargaconsulting.ca.
