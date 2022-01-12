# Clatter Rattles RS

![Build](https://github.com/bartOssh/clatter_rattles_rs/workflows/Rust/badge.svg?branch=master)

<p align="center">
    <img
    width="85%" height="85%"
    src="https://github.com/bartOssh/clatter_rattles_rs/blob/master/assets/logo.jpg"/>
</p>

## Description

:headphones: Sound pattern recognition library.

Features:

- Match a sample from a local file.
- Match a sample from a sound stream (Now only mp3 radio stream).

Strategies:

- Use a fingerprint hashing mechanism based on
  [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform)
  and [Accustic Fingerprint Hash](https://en.wikipedia.org/wiki/Acoustic_fingerprint).
- Focus on parallel stream chunk fingerprint hashing method.
- Try to offer a simple way to match a given sample against songs in database
  collection.
- Provide a built-in solution to handle fingerprint storage.
    This library offers [Redis](https://redis.io/) database support for storing
    hashes and finding matches. By using `trait Repository`, it is fairly easy
    to implement this feature in any database. Using a database engine that
    supports fast hash lookup is highly recommended. The best performance can be
    achieved on a fairly simple schema within database RAM.

Improvements:

- The weakest point so far is that only one thread can read and write to
  the database.
    The main problem is that the hash lookup wasn't showing any significant
    performance improvement when using parallel threads. Additionally, when the
    database size increases, parallel calls make hash lookups take longer.
    This factor made me stick to single-threaded database calls. (If you know
    how to improve this, please provide your solution or issue a proposal
    ticket).

### Testing

- Test against self cached file samples
- Run local redis database with `systemctl start redis-server.service`
- Test it with `cargo test -- --nocapture` or `cargo test` <- with less
  information

### Building

- Building. Please build with `RUSTFLAGS="--emit=asm"` flag enabled, which
  forces the compiler to use a single LLVM module for the entire crate, allowing
  better LLVM optimization.

### Test:

Hardware: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz, Intel Corporation Cannon Point-LP Shared SRAM Speed: 2400 MT/s,
Build with: rust stable-x86_64-unknown-linux-gnu.
Redis server: v=5.0.7 sha=00000000:0 malloc=jemalloc-5.2.1 bits=64 build=49eb440ab0786b6a

STATISTICS:

Total time to calculate and populate db with 31 songs took: 142626 ms
Average time to calculate and add to db 1 song 4600 ms
Average time to calculate and add to db 1000 fingerprints took: 807 ms
Calculating samples fingerprints

Track : Sample
Guns N' Roses - Welcome To The Jungle.mp3 : Guns N' Roses - Welcome To The Jungle.mp3
Match accuracy 58 %


Track : Sample
Maroon 5 - Memories.mp3 : Maroon 5 - Memories.mp3
Match accuracy 61 %


Track : Sample
Tove Lo - Habits .mp3 : Tove Lo - Habits .mp3
Match accuracy 70 %


Track : Sample
twenty one pilots - Ride.mp3 : twenty one pilots - Ride.mp3
Match accuracy 62 %

Total time to calculate and match 4 samples took: 2430 ms
Average time to calculate fingerprint for 1 sample and match it against fingerprints in db took 607 ms

All samples matched corresponding songs.


## Contributors

- [Claire Amalfitano](https://github.com/polypodioides)
- [Bartosz Lenart](https://github.com/bartOssh)
- [Oskar Piechowicz](https://github.com/opiechow)
- [Janusz Roll](https://github.com/janeek1995)
- [Marcin Szymczak](https://github.com/m-sz)

## License

[MIT](https://github.com/bartOssh/clatter_rattles_rs/blob/master/LICENSE)

## Dedication

:cherry_blossom: Dla Grzechotki :cherries:
