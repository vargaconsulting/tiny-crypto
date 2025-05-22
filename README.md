[![CI](https://github.com/vargaconsulting/tiny-crypto/actions/workflows/ci.yml/badge.svg)](https://github.com/vargaconsulting/tiny-crypto/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/vargaconsulting/tiny-crypto/branch/main/graph/badge.svg)](https://codecov.io/gh/vargaconsulting/tiny-crypto)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

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
  - Point-at-infinity and ECPoint construction
- Clean ANSI-formatted output for symbolic math (e.g., `𝔽₃₁`, `𝔾(6,2)`)
- Fully written in **Julia**, for clarity and performance.
- **Educational focus** – Ideal for teaching and cryptographic exploration.

## Example Usage

```julia
using TinyCrypto

# Find a suitable Weierstrass curve over small field primes and coefficient ranges
curve = Weierstrass(97:103, 10:15, 2:7)  # (prime range, a range, b range)
# Output: Weierstrass curve: y² = x³ + 10x + 3 |𝔽₉₇ with order: 101 and 𝔾(0,10)

# Get all curve points
E = curve_points(curve)
# → 101-element Vector{ECPoint₁₂₈}: (0𝔽₉₇,10𝔽₉₇), (0𝔽₉₇,87𝔽₉₇), ..., (96𝔽₉₇,63𝔽₉₇), (∞,∞)

# ECDSA example
private_key = 42
msg = 33

public_key = private_key2public(private_key, curve) # ECPoint on curve
r, s, v = ecdsa_sign(curve, private_key, msg)
ecdsa_verify(curve, public_key, r, s, v)  # true
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
