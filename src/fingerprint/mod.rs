use rustfft::algorithm::Radix4;
use rustfft::num_complex::Complex;
use rustfft::num_traits::Zero;
use rustfft::FFT;
use std::error::Error;
use std::fs::File;
use read_byte_slice::{ByteSliceIter, FallibleStreamingIterator};

const CHUNK_SIZE: usize = 1024; // chunk size from readable buffer (file or stream)
const CHANNELS: usize = 2; // number of channels to compose floting point number from 

const FREQ_BINS: &[usize] = &[32, 40, 80, 120, 180, 320]; // Each value in array is a top range frequency to calculate local maximum magnitude for
const FUZZ_FACTOR: usize = 2; // higher the value of this factor, lower the fingerprint entropy, and less bias the algorithm become to the sound noises

pub enum FingerprinterErr {
    WrongChunksSize(Box<dyn Error>),
    NotCalculated,
}

/// Helper struct for calculating acoustic fingerprint
///
#[allow(dead_code)]
pub struct FingerprintHandle {
    /// FFT algorithm
    fft: Radix4<f32>,
    fft_window_size: usize,
    input: Vec<Complex<f32>>
}

#[allow(dead_code)]
impl FingerprintHandle {
    pub fn new() -> FingerprintHandle {
        let fft_window_size = CHUNK_SIZE * 2;
        let input = Vec::new();
        FingerprintHandle {
            fft: Radix4::new(fft_window_size, false),
            fft_window_size,
            input
        }
    }

    /// Calculate discret representation of accustic file
    /// 
    /// # Arguments:
    /// * path - path to accustinc file
    /// 
    /// # Returns buffer of discretize accustic files if success or dynamic error otherwise
    /// 
    pub fn transform_file_to_discret(&mut self, path: &String) -> Result<Vec<usize>, Box<dyn Error>> {
        let f = File::open(path)?;
        let mut iter_chunk = ByteSliceIter::new(f, CHUNK_SIZE);
        let mut fingerprints = Vec::new();
        while let Some(chunk) = iter_chunk.next()? {
            match self.transform_chunk_to_fft_fingerprint(chunk) {
                Ok(fingerprint) => fingerprints.push(fingerprint), 
                Err(FingerprinterErr::NotCalculated) => (),
                Err(FingerprinterErr::WrongChunksSize(e)) => println!("{}", e),
            }
        }
        Ok(fingerprints)
    }

    /// Transforms buffer chunk to district fingerprint
    /// 
    /// #Arguments:
    /// * bites - chunk of size 1024 bites to be transformed
    /// 
    /// # Returns single fingerprint if success or error FingerprinterErr otherwie
    /// 
    pub fn transform_chunk_to_fft_fingerprint(&mut self, bites: &[u8]) -> Result<usize, FingerprinterErr> {
        if bites.len() == CHUNK_SIZE {
            for sample in bites.chunks_exact(CHANNELS) {
                self.input.push(Complex::from(f32::from(
                    sample.iter().fold(0, |sum, x| sum + *x / CHANNELS as u8)),
                ));
            }
            if self.input.len() == self.fft_window_size {
                let mut output: Vec<Complex<f32>> = vec![Complex::zero(); self.fft_window_size];
                self.fft.process(&mut self.input, &mut output);
                self.input.clear();
                return Ok(calculate_fingerprint(&output));
            }
            return Err(FingerprinterErr::NotCalculated);
        }
        Err(FingerprinterErr::WrongChunksSize(Box::from(format!("Wrong bites length, is: {}, should be: {}", bites.len() , CHUNK_SIZE))))
    }
}

/// Find points with max magnitude in each of the bins
/// 
#[allow(dead_code)]
fn calculate_fingerprint(arr: &[Complex<f32>]) -> usize {
    let mut high_scores: Vec<f32> = vec![0.0; FREQ_BINS.len()];
    let mut record_points: Vec<usize> = vec![0; FREQ_BINS.len()];

    for bin in FREQ_BINS[0]..=FREQ_BINS[FREQ_BINS.len() - 1] {
        let magnitude = arr[bin].re.hypot(arr[bin].im);

        let mut bin_idx = 0;
        while FREQ_BINS[bin_idx] < bin {
            bin_idx += 1;
        }

        if magnitude > high_scores[bin_idx] {
            high_scores[bin_idx] = magnitude;
            record_points[bin_idx] = bin;
        }
    }

    encode(&record_points)
}

/// Encodeing function with reverse order
/// 
fn encode(arr: &[usize]) -> usize {
    (arr[4] - (arr[4] % FUZZ_FACTOR)) * usize::pow(10, 10)
        + (arr[3] - (arr[3] % FUZZ_FACTOR)) * usize::pow(10, 8)
        + (arr[2] - (arr[2] % FUZZ_FACTOR)) * usize::pow(10, 5)
        + (arr[1] - (arr[1] % FUZZ_FACTOR)) * usize::pow(10, 2)
        + (arr[0] - (arr[0] % FUZZ_FACTOR))
}

#[cfg(test)]
mod tests {
    use rand::prelude::*;
    #[test]
    // #[ignore]
    fn test_hash() {
        let record_points_0 = vec![32, 45, 100, 140, 235, 300];
        let record_points_1 = vec![33, 45, 100, 145, 235, 300];
        assert_eq!(super::encode(&record_points_0), 2354010004432);
        assert_ne!(super::encode(&record_points_0), super::encode(&record_points_1));
    }
    #[test]
    fn test_calculate_fingerprint() {
        let mut rng = rand::thread_rng();
        let mut arr_f32: Vec<f32> = vec![0.0; super::CHUNK_SIZE * 2];
        arr_f32.iter_mut().for_each(|complex_num| {
            *complex_num = rng.gen::<f32>() * 10000_f32;
        });
        let arr: Vec<super::Complex<f32>> = arr_f32.iter().map(super::Complex::from).collect();
        let fingerprint: usize = super::calculate_fingerprint(&arr);
        let fingerprint_log10 = (fingerprint as f64).log10();
        assert_eq!(
            fingerprint_log10 > 12_f64 && fingerprint_log10 < 13_f64,
            true
        );
    }

    #[test]
    fn test_transform_file_to_discret() {
        use std::time::Instant;
        let now = Instant::now();
        let mut fingerprinter = super::FingerprintHandle::new();
        let fingerprints = fingerprinter.transform_file_to_discret(&format!("assets/space_cover.mp3")).unwrap();
        println!("Transfroming to discret first file took: {} ms", now.elapsed().as_millis());
        println!("Total number of fingerprints of first file is: {}", fingerprints.len());
        let fingerprints = fingerprinter.transform_file_to_discret(&format!("assets/The 8 second music video.mp3")).unwrap();
        println!("Transfroming to discret second file took: {} ms", now.elapsed().as_millis());
        println!("Total number of fingerprints of second file is: {}", fingerprints.len());
    }
}
