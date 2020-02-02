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

- A fingerprint hashing mechanism based on
  [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform)
  and [Fingerprint Hash](https://en.wikipedia.org/wiki/Fingerprint_(computing)).
- A focus on parallel stream chunk fingerprint hashing method.
- Try to offer a simple way to match a given sample against songs in database
  collection.
- A built-in solution to handle fingerprint storage.
    This library offers [Redis](https://redis.io/) database support for storing
    hashes and finding matches. By using `trait Repository`, it is fairly easy
    to implement this feature in any database. Using a database engine that
    supports fast hash lookup is highly recommended. The best performance can be
    achieved on a fairly simple schema within database RAM.
- Async stream request and parallel stream hashing.
    There is no need to build a parallel mechanism by yourself.
    The stream listener collects chunks asynchronously and decodes them in a
    separate thread, then pipes it to the receiver. The receiver is boxed by
    atomic type, so hashing and matching can be done in a separate thread than
    the listener. This design allows more processor resources to be used, and
    speeds up the matching algorithm significantly.

Improvements:

- The weakest point so far is that only one thread can read and write to
  database.
    The main problem is that the hash lookup wasn't showing any significant
    performance improvement when using parallel threads. Additionally, when the
    database size increases, parallel calls make hash lookups take longer.
    This factor made me stick to single-threaded database calls. (If you know
    how to improve this, please provide your solution or issue a proposal
    ticket).

### Testing

- Run local redis database with `systemctl start redis-server.service`
- Test it with `cargo test -- --nocapture` or `cargo test` <- with less
  information

### Building

- Building. Please build with `RUSTFLAGS="--emit=asm"` flag enabled, which
  forces the compiler to use a single LLVM module for the entire crate, allowing
  better LLVM optimization.

## License

[MIT](https://github.com/bartOssh/clatter_rattles_rs/blob/master/LICENSE)

## Dedication

:cherry_blossom: Dla Grzechotki :cherries:
