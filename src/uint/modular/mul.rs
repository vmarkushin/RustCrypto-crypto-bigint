use crate::{Limb, Uint};

use super::reduction::montgomery_reduction;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_zkvm_platform::syscall::{bigint, sys_bigint};

// TODO(victor): Are there any safety concerns with the fact that the output from sys_bigint may
// not be reduced. If so, we can introduce a check here to make sure it's less than the modulus.
pub(crate) fn mul_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    b: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
    r_inv: &Uint<LIMBS>,
) -> Uint<LIMBS> {
    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    if LIMBS == bigint::WIDTH_WORDS {
        let mut result = core::mem::MaybeUninit::<[u32; LIMBS]>::uninit();
        return Uint::<LIMBS>::from_words(unsafe {
            sys_bigint(
                result.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                a.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                b.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            // a and b are in montgomery form (a' * R) and (b' * R).
            // Getting the final result ((a' * b') * R) requires removing a multiple of R.
            sys_bigint(
                result.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                result.as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                r_inv.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            result.assume_init()
        });
    }

    let product = a.mul_wide(b);
    montgomery_reduction::<LIMBS>(&product, modulus, mod_neg_inv)
}

pub(crate) fn square_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
    r_inv: &Uint<LIMBS>,
) -> Uint<LIMBS> {
    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    if Uint::<LIMBS>::BITS == bigint::WIDTH_BITS && Uint::<LIMBS>::LIMBS == bigint::WIDTH_WORDS {
        let mut result = core::mem::MaybeUninit::<[u32; LIMBS]>::uninit();
        return Uint::<LIMBS>::from_words(unsafe {
            sys_bigint(
                result.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                a.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                a.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            // a and b are in montgomery form (a' * R) and (b' * R).
            // Getting the final result ((a' * b') * R) requires removing a multiple of R.
            sys_bigint(
                result.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                result.as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                r_inv.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            result.assume_init()
        });
    }

    let product = a.square_wide();
    montgomery_reduction::<LIMBS>(&product, modulus, mod_neg_inv)
}
