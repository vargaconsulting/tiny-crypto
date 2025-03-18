# Tiny Crypto – Exploring Cryptography with Small Prime Fields

**A small-scale cryptography playground using tiny prime fields for easy manual verification.**

## Overview

TinyCurve is an educational project that implements fundamental cryptographic algorithms using **small prime fields**. The goal is to **simplify the math** and make manual verification feasible, allowing learners to focus on **concepts** without getting lost in large numbers.

The project includes:
- **Elliptic Curve Cryptography (ECC)** – Weierstrass curves, point addition, and scalar multiplication.
- **Shamir's Secret Sharing** – Securely split and reconstruct secrets.
- **Distributed Key Management** – Collaborative cryptographic protocols.
- **Other Cryptographic Routines** – Modular arithmetic, prime field operations, and more.

## Features
- Uses **tiny prime fields** to make calculations manually verifiable.
- Implements **common elliptic curve operations**, including:
  - Curve generation over small primes
  - Point addition & scalar multiplication
  - Checking if a point is on a curve
  - Finding generators for prime order curves
- Written in **Julia**, a language known for its mathematical clarity.
- **Educational focus** – Ideal for teaching and understanding cryptography.

## Example Usage
```julia
using TinyCrypto

curve = Weierstrass(97:103, 0:2, 6:100)
# Weierstrass curve: y² = x³ + 7 mod 97 with order: 79 and generator point: (1,28)
E = curve_points(curve)
# 79-element Vector{Any}:  (1,28) (1,69) ...  (96,54) (∞,∞)

private_key, msg = 42, 33
public_key = private_key2public(private_key, curve) # (71,45)
r,s,v = ecdsa_sign(curve, private_key, msg)         # (9,57,0)
ecdsa_verify(curve, public_key,  r,s,v)             # true
```

## Installation
```bash
git clone https://github.com/vargaconsulting/tiny-crypto.git
cd tiny-crypto
```

## Future Work
- Implement ECDSA signing & verification.
- Add Shamir's Secret Sharing Scheme.
- Proactive Update 
- Explore distributed key management protocols.
- Add Playground for Distributed Signing Systems

## Licence 
MIT License. (c) 2010-2025 Varga Consulting, Toronto, ON.
For questions, reach out at info@vargaconsulting.ca.
