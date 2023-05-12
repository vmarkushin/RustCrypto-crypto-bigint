//! Test to ensure that `impl_modulus!` works from outside this crate.

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_zkvm as _;

use crypto_bigint::{impl_modulus, U64};

impl_modulus!(TestMod, U64, "30e4b8f030ab42f3");
