use std::collections::HashMap;

/// Chooses the song or acoustic sample that has highest value of matching hashes
///
/// # Arguments:
/// * findings - collection of all songs with matching hash count as a value
///
/// # Returns tuple of best matching song, the one with the highest value
///
#[allow(dead_code)]
pub fn pick_most_likely(findings: &HashMap<String, usize>) -> (String, usize) {
    // TODO: consider returning Option or Result
    
    let empty = "No match found".to_string();
    let best = findings.iter().fold(
        (&empty, 0),
        |mut current, (string, &count)| {
            if count > current.1 {
                current.1 = count;
                current.0 = string;
            }
            current
        }
    );

    (best.0.clone(), best.1)
}

#[cfg(test)]
mod test {
    use crate::helpers::pick_most_likely;
    #[test]
    fn test_pick_most_likely() {
        use std::collections::HashMap;

        let mut findings = HashMap::new();
        findings.insert("Mock song 1".to_string(), 1);
        findings.insert("Mock song 2".to_string(), 2);
        findings.insert("Mock song 3".to_string(), 0);
        findings.insert("Best mock song".to_string(), 4);

        let (best_fit, count) = pick_most_likely(&findings);
        assert_eq!(best_fit, "Best mock song");
        assert_eq!(count, 4);
    }
}
