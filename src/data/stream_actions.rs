// use super:: PlaylistHelper;
// use crate::helpers::pick_most_likely;
// use crossbeam_channel::{unbounded, Receiver, Sender};
// use m3u8_rs::playlist::{MasterPlaylist, MediaPlaylist, Playlist, VariantStream};
// use reqwest::{get, Url};
// use std::collections::HashMap;
// use std::error::Error;
// use std::io::Cursor;
// use std::sync::{Arc, Mutex};
// use std::thread;
// use std::thread::JoinHandle;
// use tokio::runtime::Runtime;

// /// Helps to watch for matches of fingerprint findings from the stream
// /// 
// #[derive(Clone, Debug)]
// pub struct MatchesWatcher {
//     findings: HashMap<String, usize>,
//     fingerprints_count: usize,
//     fingerprints_threshold: usize,
//     chunks_count: usize,
//     chunks_threshold: usize,
// }

// /// Helps to reads stream and holds multi thread pipe to for sending stream between threads 
// /// 
// #[derive(Clone, Debug)]
// struct StreamListener {
//     uri: Url,
//     receiver: Receiver<Vec<f32>>,
//     sender: Sender<Vec<f32>>,
//     is_active: bool,
// }

// /// Casing StreamListener to allow threaded atomic and mutable access to its active state and receiver
// /// 
// /// Thanks to boxing by atomic type approach it is easy to leverage application performance
// /// by running stream listener that is a decoder - writer
// /// and stream reader that can be a fingerprint hasher - stream matcher in separate threads
// ///    
// #[derive(Clone, Debug)]
// pub struct ArcStreamListener(Arc<Mutex<StreamListener>>);

// #[allow(dead_code)]
// impl MatchesWatcher {
//     /// Create instance of MatchesWatcher
//     ///
//     /// # Arguments:
//     /// * fingerprints_threshold - threshold for fingerprints that should be totally matched for chunks number equal chunks_threshold
//     /// * chunks_threshold - threshold for max chunks to take to account for watching for matching fingerprints
//     ///
//     /// # Returns new instance of MatchesWatcher
//     /// 
//     pub fn new(fingerprints_threshold: usize, chunks_threshold: usize) -> Self {
//         let findings = HashMap::new();
//         let fingerprints_count = 0;
//         let chunks_count = 0;
//         Self {
//             findings,
//             fingerprints_count,
//             fingerprints_threshold,
//             chunks_count,
//             chunks_threshold,
//         }
//     }
//     /// Feeds matches watcher with fingerprint findings
//     ///
//     /// # Arguments:
//     /// * findings - collection of songs and value of matching fingerprints for one stream chunk
//     ///
//     /// # Returns Option with tuple of sum of findings and difference between most likely matching song for given findings
//     ///
//     pub fn feed(
//         &mut self,
//         findings: HashMap<String, usize>,
//     ) -> Option<(HashMap<String, usize>, usize)> {
//         for (song, value) in findings.iter() {
//             match self.findings.get_mut(song) {
//                 Some(count) => {
//                     *count += value;
//                 }
//                 None => {
//                     self.findings.insert(song.clone(), *value);
//                 }
//             };
//         }
//         self.fingerprints_count = pick_most_likely(&self.findings).1;
//         self.chunks_count += 1;
//         if self.fingerprints_count >= self.fingerprints_threshold
//             || self.chunks_count == self.chunks_threshold
//         {
//             let difference_from_threshold = self.fingerprints_threshold - self.fingerprints_count;
//             self.fingerprints_count = 0;
//             self.chunks_count = 0;
//             return Some((self.findings.clone(), difference_from_threshold));
//         }
//         None
//     }
// }

// #[allow(dead_code)]
// impl ArcStreamListener {
//     /// Create new instance of StreamListener
//     ///
//     /// # Arguments:
//     ///
//     /// * uri - address of external server emitting stream
//     ///
//     /// # Returns result of instance of StreamListener or dyn Error otherwise
//     ///
//     pub fn new(uri: String) -> Result<Self, Box<dyn Error>> {
//         let (sender, receiver) = unbounded();
//         let uri = Url::parse(&uri)?;
//         let is_active = false;
//         Ok(Self(Arc::new(Mutex::new(StreamListener {
//             uri,
//             receiver,
//             sender,
//             is_active,
//         }))))
//     }

//     /// Getter for receiver
//     ///
//     /// # Returns receiver pipe that listen for decoded stream chunk
//     ///
//     pub fn get_listener(&self) -> Receiver<Vec<f32>> {
//         self.0.lock().unwrap().receiver.clone()
//     }

//     /// Check if stream listener is active
//     ///
//     /// # Returns true if listener is in active state, false otherwise
//     /// 
//     pub fn is_active(&self) -> bool {
//         self.0.lock().unwrap().is_active
//     }

//     /// Deactivates stream listener
//     /// 
//     pub fn deactivate(&mut self) {
//         self.0.lock().unwrap().is_active = false;
//     }

//     /// Runs loop for fetching m3u8 stream
//     ///
//     /// # Returns Ok if endpoint responded with valid playlist and stream or dyn Error otherwise
//     ///
//     pub async fn run_m3u8(&mut self) -> Result<(), Box<dyn Error + 'static>> {
//         let master_playlist = fetch_master_playlist(&self.0.lock().unwrap().uri).await?;
//         if let Ok(uri) = master_playlist.find_uri() {
//             let uri = Url::parse(&uri)?;
//             let media_playlist = fetch_media_playlist(&uri).await?;
//         };
//         // TODO: need implementation of loop of stream read
//         Ok(())
//     }

//     /// Runs listener ins separate thread that collects stream and feeds pipe sender
//     ///
//     /// # Returns success of thread join handle if listener is activated and async block has no errors, dyn Error otherwise
//     ///
//     pub async fn run_mp3(&self) -> Result<JoinHandle<()>, Box<dyn Error>> {
//         if self.0.lock().unwrap().is_active == true {
//             return Err(Box::from(format!(
//                 "Listener is active and should be deactivated first"
//             )));
//         }
//         self.0.lock().unwrap().is_active = true;
//         let listener_clone = self.clone();
//         let stream_listener_proc = thread::spawn(move || {
//             if let Err(e) = Runtime::new()
//                 .unwrap()
//                 .block_on(listen_mp3_stream(listener_clone))
//             {
//                 panic!(format!("{:?}", e));
//             };
//         });
//         Ok(stream_listener_proc)
//     }
// }

// impl PlaylistHelper for MasterPlaylist {
//     fn find_uri(&self) -> Result<String, Box<dyn Error>> {
//         let variants: Vec<VariantStream> = self
//             .variants
//             .clone()
//             .into_iter()
//             .filter(|variant_stream| variant_stream.uri.len() > 0)
//             .collect();
//         if variants.len() > 0 {
//             return Ok(variants[0].uri.to_owned());
//         }
//         Err(Box::from(format!(
//             "No available uri in MasterPlaylist variants"
//         )))
//     }
// }

// async fn listen_mp3_stream(listener: ArcStreamListener) -> Result<(), Box<dyn Error + 'static>> {
//     let mut res = reqwest::get(listener.0.lock().unwrap().uri.clone()).await?;
//     while let Some(chunk) = res.chunk().await? {
//         if listener.0.lock().unwrap().is_active == false {
//             // stop loop and finish listening
//             break;
//         };
//         let readable_buffer = Cursor::new(chunk);
//         let decoded = decode_mp3_from_chunk(readable_buffer);
//         match decoded {
//             Ok(_result) => {
//                 if _result.len() > 0 {
//                     listener.0.lock().unwrap().sender.send(_result)?;
//                 }
//             }
//             Err(_) => (),
//         }
//     }
//     Ok(())
// }

// async fn fetch_master_playlist(uri: &Url) -> Result<MasterPlaylist, Box<dyn Error + 'static>> {
//     let text = get_from_as_string(uri).await?;
//     match m3u8_rs::parse_playlist_res(&text.as_bytes()) {
//         Ok(Playlist::MasterPlaylist(pl)) => Ok(pl),
//         Ok(Playlist::MediaPlaylist(_)) => Err(Box::from(format!("Wrong format of media playlist"))),
//         Err(e) => Err(Box::from(format!("{:?}", e))),
//     }
// }

// async fn fetch_media_playlist(uri: &Url) -> Result<MediaPlaylist, Box<dyn Error + 'static>> {
//     let text = get_from_as_string(uri).await?;
//     match m3u8_rs::parse_playlist_res(&text.as_bytes()) {
//         Ok(Playlist::MasterPlaylist(_)) => {
//             Err(Box::from(format!("Wrong format of media playlist")))
//         }
//         Ok(Playlist::MediaPlaylist(pl)) => Ok(pl),
//         Err(e) => Err(Box::from(format!("{:?}", e))),
//     }
// }

// async fn get_from_as_string(uri: &Url) -> Result<String, Box<dyn Error + 'static>> {
//     let text = get(uri.clone()).await?.text().await?;
//     Ok(text)
// }

// #[cfg(test)]
// mod test {
//     use super::{ArcStreamListener, Runtime};
//     use futures_await_test::async_test;
//     use std::thread;
//     use std::thread::sleep;
//     use std::time::Duration;

//     #[async_test]
//     #[ignore]
//     async fn test_get_m3u8_stream() {
//         dotenv::dotenv().ok();
//         if dotenv!("RADIO_STREAM_ENABLED") == "true" {
//             let mut listener = ArcStreamListener::new(
//                 format!("http://a.files.bbci.co.uk/media/live/manifesto/audio/simulcast/hls/uk/sbr_high/ak/bbc_radio_two.m3u8")
//             ).unwrap();
//             let a = Runtime::new().unwrap().block_on(listener.run_m3u8());
//             println!("\nGetting m3ue stream {:?}\n", &a);
//         } else {
//             println!("test_get_m3ue_stream does nothing");
//         }
//     }
//     #[async_test]
//     #[ignore]
//     async fn test_get_mp3_stream() {
//         dotenv::dotenv().ok();
//         if dotenv!("RADIO_STREAM_ENABLED") == "true" {
//             let mut listener = ArcStreamListener::new(
//                 format!("https://str2b.openstream.co/604?aw_0_1st.collectionid=3162&stationId=3162&publisherId=628&listenerid=1580311050432_0.47836979431904714&awparams=companionAds%3Atrue&aw_0_1st.version=1.1.4%3Ahtml5")
//             ).unwrap();
//             let receiver = listener.get_listener();
//             let reader = thread::spawn(move || {
//                 let mut tested = false;
//                 loop {
//                     let decoded = receiver.recv().unwrap();
//                     if !tested {
//                         tested = true;
//                         // testing is stream correct
//                         assert_eq!(decoded.len() > 0, true);
//                         println!(
//                             "\nReceived decoded stream of {:?} floats by crossbeam channel pipe",
//                             &decoded.len()
//                         );
//                     } else {
//                         break;
//                     }
//                 }
//             });
//             if let Ok(a) = Runtime::new().unwrap().block_on(listener.run_mp3()) {
//                 sleep(Duration::from_secs(3));
//                 listener.deactivate();
//                 a.join().unwrap();
//                 reader.join().unwrap();
//             };
//         } else {
//             println!("test_get_mp3_stream does nothing");
//             return ();
//         }
//     }
// }
