use crate::{Limb, Uint};

use super::reduction::montgomery_reduction;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_zkvm_platform::syscall::{bigint, sys_bigint};

// TODO(victor): Implement this method with sys_bigint.
pub(crate) const fn mul_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    b: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
) -> Uint<LIMBS> {
    let product = a.mul_wide(b);
    montgomery_reduction::<LIMBS>(&product, modulus, mod_neg_inv)
}

// TODO(victor): Implement this method with sys_bigint.
pub(crate) const fn square_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
) -> Uint<LIMBS> {
    let product = a.square_wide();
    montgomery_reduction::<LIMBS>(&product, modulus, mod_neg_inv)
}
