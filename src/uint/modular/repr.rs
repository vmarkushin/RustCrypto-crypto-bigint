use crate::{Limb, Uint};

use super::reduction::montgomery_reduction;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_zkvm_platform::syscall::{bigint, sys_bigint};

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use subtle::ConstantTimeLess;

pub(crate) fn into_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    r2: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
    _r: &Uint<LIMBS>,
) -> Uint<LIMBS> {
    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    if LIMBS == bigint::WIDTH_WORDS {
        let result = Uint::<LIMBS>::from_words(unsafe {
            let mut out = core::mem::MaybeUninit::<[u32; LIMBS]>::uninit();
            sys_bigint(
                out.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                a.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                _r.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            out.assume_init()
        });
        // Assert that the Prover returned the canonical representation of the result, i.e. that it
        // is fully reduced and has no multiples of the modulus included.
        assert!(bool::from(result.ct_lt(&modulus)));
        return result;
    }

    let product = a.mul_wide(r2);
    montgomery_reduction::<LIMBS>(&product, modulus, mod_neg_inv)
}

pub(crate) fn from_montgomery_form<const LIMBS: usize>(
    a: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
    mod_neg_inv: Limb,
    _r_inv: &Uint<LIMBS>,
) -> Uint<LIMBS> {
    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    if LIMBS == bigint::WIDTH_WORDS {
        let result = Uint::<LIMBS>::from_words(unsafe {
            let mut out = core::mem::MaybeUninit::<[u32; LIMBS]>::uninit();
            sys_bigint(
                out.as_mut_ptr() as *mut [u32; bigint::WIDTH_WORDS],
                bigint::OP_MULTIPLY,
                a.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                _r_inv.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
                modulus.as_words().as_ptr() as *const [u32; bigint::WIDTH_WORDS],
            );
            out.assume_init()
        });
        // Assert that the Prover returned the canonical representation of the result, i.e. that it
        // is fully reduced and has no multiples of the modulus included.
        assert!(bool::from(result.ct_lt(&modulus)));
        return result;
    }

    montgomery_reduction::<LIMBS>(&(*a, Uint::<LIMBS>::ZERO), modulus, mod_neg_inv)
}
