# Kira modular multiplication patch

`kira-3.1-mod-mul.patch` is the unmodified upstream commit diff for
[Kira MR !61](https://gitlab.com/kira-pyred/kira/-/merge_requests/61),
"Fix mod_mul silently miscomputing on AArch64 (Apple Silicon)", by
Paulo Cabral Sanz. Copyright and GPL-3.0-or-later terms of the affected Kira
source remain applicable; the upstream GPL text is included as `Kira-COPYING`.

- Kira base: `5760f42f8decd7eb5474c780db75b5e5e7b450a0` (3.1).
- Fix commit: `dad964cd79a39a9fca728bb2606f59537a95109c`.
- Immutable diff: <https://gitlab.com/kira-pyred/kira/-/commit/dad964cd79a39a9fca728bb2606f59537a95109c.diff>.
- Patch SHA-256: `ac74e82e9f1e5ebf3945d6738159d2dee251789f252900fafae3de2ec8ad03ec`.
- FireFly: `0e155a9d879362521e9157e8fdac10b9553fcacb` (2.0.3).
- FLINT: 2.8.4; the pinned upstream wrap contains SHA-256 hashes for the
  source archive and Meson adaptation.

The installer uses this local patch, checks its hash, verifies the Kira
checkout, and replaces FireFly's mutable `kira-2` wrap revision with the
commit above. System compiler and Homebrew/library versions can vary;
`install.log` records Meson's dependency discovery and `installation.json`
records source pins and the installed binary hash. This is reproducible
source selection, not a claim of byte-identical builds across machines.

The precision guard selects exact `__uint128_t` multiplication when
`LDBL_MANT_DIG < 64`. It is a precision condition, not a blanket ARM ban.
Keep the exact coefficient checks when moving to a future upstream release.
Review the release's source and rerun the platform tests before dropping
this patch. Never infer patch presence from `3.1-dirty` alone.
