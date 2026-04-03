#![allow(dead_code, unused_variables)]
use clap::Parser;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Result, Context};
use rayon::ThreadPoolBuilder;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input primer file (TSV or FASTA format)
    #[arg(short, long)]
    input: PathBuf,
    
    /// Output file path
    #[arg(short, long)]
    output: Option<PathBuf>,
    
    /// Maximum allowed mismatches in reverse complement alignment
    #[arg(short, long, default_value_t = 2)]
    mismatch: u32,
    
    /// Badness threshold (only output pairs with badness > threshold)
    #[arg(short, long, default_value_t = 10.0)]
    threshold: f64,
    
    /// Number of threads for parallel processing (must be ≥1)
    #[arg(short = 'T', long, default_value_t = 4)]
    threads: usize,
}

#[derive(Debug, Clone)]
struct Primer {
    id: String,
    sequence: String,
}

impl Primer {
    fn new(id: String, sequence: String) -> Self {
        Self {
            id,
            sequence: sequence.to_uppercase(),
        }
    }
}

fn read_primers(path: &PathBuf) -> Result<Vec<Primer>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open primer file: {:?}", path))?;
    let reader = BufReader::new(file);
    
    let extension = path.extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or("")
        .to_lowercase();
    
    match extension.as_str() {
        "tsv" => read_tsv_primers(reader),
        "fa" | "fasta" => read_fasta_primers(reader),
        _ => {
            // Try to auto-detect format by content
            let mut lines = reader.lines();
            let first_line = lines.next()
                .transpose()?
                .unwrap_or_default();
            
            if first_line.starts_with('>') {
                // Reopen file for FASTA reading
                let file = File::open(path)?;
                let reader = BufReader::new(file);
                read_fasta_primers(reader)
            } else if first_line.contains('\t') {
                // Reopen file for TSV reading
                let file = File::open(path)?;
                let reader = BufReader::new(file);
                read_tsv_primers(reader)
            } else {
                anyhow::bail!("Unknown primer file format. Expected TSV or FASTA")
            }
        }
    }
}

fn read_tsv_primers<R: BufRead>(reader: R) -> Result<Vec<Primer>> {
    let mut primers = Vec::new();
    
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 2 {
            anyhow::bail!("Line {}: expected at least 2 columns (ID and sequence), found {}", 
                         line_num + 1, parts.len());
        }
        
        let id = parts[0].trim().to_string();
        let sequence = parts[1].trim().to_string();
        
        if !sequence.chars().all(|c| "ATCGatcg".contains(c)) {
            anyhow::bail!("Line {}: sequence contains invalid characters: {}", 
                         line_num + 1, sequence);
        }
        
        primers.push(Primer::new(id, sequence));
    }
    
    Ok(primers)
}

fn read_fasta_primers<R: BufRead>(reader: R) -> Result<Vec<Primer>> {
    let mut primers = Vec::new();
    let mut current_id = None;
    let mut current_sequence = String::new();
    
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        
        if trimmed.is_empty() {
            continue;
        }
        
        if trimmed.starts_with('>') {
            // Save previous primer if exists
            if let Some(id) = current_id.take() {
                if !current_sequence.is_empty() {
                    primers.push(Primer::new(id, current_sequence));
                }
                current_sequence = String::new();
            }
            
            // Parse new ID (remove '>' and any whitespace)
            let id = trimmed[1..].trim().to_string();
            if id.is_empty() {
                anyhow::bail!("Line {}: empty primer ID", line_num + 1);
            }
            current_id = Some(id);
        } else {
            // Sequence line
            if current_id.is_none() {
                anyhow::bail!("Line {}: sequence line without preceding ID", line_num + 1);
            }
            
            let seq_part = trimmed.replace(' ', "").replace('\t', "");
            if !seq_part.chars().all(|c| "ATCGatcg".contains(c)) {
                anyhow::bail!("Line {}: sequence contains invalid characters: {}", 
                             line_num + 1, seq_part);
            }
            current_sequence.push_str(&seq_part);
        }
    }
    
    // Save the last primer
    if let Some(id) = current_id.take() {
        if !current_sequence.is_empty() {
            primers.push(Primer::new(id, current_sequence));
        }
    }
    
    Ok(primers)
}

// ----------------------------------------------------------------------------
// Reverse complement and alignment algorithms
// ----------------------------------------------------------------------------

fn reverse_complement(seq: &str) -> String {
    let mut result = String::with_capacity(seq.len());
    for c in seq.chars().rev() {
        match c {
            'A' | 'a' => result.push('T'),
            'T' | 't' => result.push('A'),
            'C' | 'c' => result.push('G'),
            'G' | 'g' => result.push('C'),
            _ => result.push(c), // Keep other characters as is
        }
    }
    result
}

/// Result of reverse complement alignment
#[derive(Debug, Clone)]
struct AlignmentResult {
    subsequence: String,      // The reverse complementary subsequence
    length: usize,            // Length of the subsequence
    gc_count: usize,          // Number of G/C nucleotides
    d1: usize,                // Distance to 3' end of primer1
    d2: usize,                // Distance to 3' end of primer2
    mismatches: usize,        // Number of mismatches in alignment
}

/// Find the longest common substring with allowed mismatches
/// Implements the same DP algorithm as sub.badness.reverse_complementary.mismatch2.py
fn find_longest_common_substring_with_mismatches(
    seq1: &str,
    seq2: &str,
    max_mismatches: u32,
) -> Option<AlignmentResult> {
    let len1 = seq1.len();
    let len2 = seq2.len();
    
    if len1 == 0 || len2 == 0 {
        return None;
    }
    
    let max_mismatches = max_mismatches as usize;
    
    // dp[i][j][k]: longest common suffix length ending at s1[i-1], s2[j-1] with k mismatches
    let mut dp = vec![vec![vec![0; max_mismatches + 1]; len2 + 1]; len1 + 1];
    
    let mut max_len = 0;
    let mut end_i = 0;
    let mut end_j = 0;
    let mut used_mismatches = 0;
    
    // Fill DP table
    for i in 1..=len1 {
        for j in 1..=len2 {
            for k in 0..=max_mismatches {
                if seq1.as_bytes()[i - 1] == seq2.as_bytes()[j - 1] {
                    // Match
                    dp[i][j][k] = dp[i - 1][j - 1][k] + 1;
                } else if k > 0 {
                    // Mismatch, consume one mismatch count
                    dp[i][j][k] = dp[i - 1][j - 1][k - 1] + 1;
                } else {
                    // No mismatches allowed and characters don't match
                    dp[i][j][k] = 0;
                }
                
                if dp[i][j][k] > max_len {
                    max_len = dp[i][j][k];
                    end_i = i;
                    end_j = j;
                    used_mismatches = k;
                }
            }
        }
    }
    
    if max_len == 0 {
        return None;
    }
    
    // Backtrack to extract the subsequence
    let mut subsequence_chars = Vec::with_capacity(max_len);
    let mut i = end_i;
    let mut j = end_j;
    let mut k = used_mismatches;
    let mut remaining_len = max_len;
    
    while i > 0 && j > 0 && remaining_len > 0 {
        subsequence_chars.push(seq1.as_bytes()[i - 1] as char);
        
        if seq1.as_bytes()[i - 1] == seq2.as_bytes()[j - 1] {
            i -= 1;
            j -= 1;
        } else if k > 0 {
            i -= 1;
            j -= 1;
            k -= 1;
        }
        remaining_len -= 1;
    }
    
    subsequence_chars.reverse();
    let subsequence: String = subsequence_chars.into_iter().collect();
    
    // Calculate distances (following Python script logic)
    let d1 = len1 - end_i;
    let d2 = end_j - max_len;
    
    // Count GC in subsequence
    let gc_count = subsequence.chars()
        .filter(|c| *c == 'G' || *c == 'C' || *c == 'g' || *c == 'c')
        .count();
    
    Some(AlignmentResult {
        subsequence,
        length: max_len,
        gc_count,
        d1,
        d2,
        mismatches: used_mismatches,
    })
}

/// Find the longest common substring (exact match, no mismatches)
/// This implements the same algorithm as the Python script
fn find_exact_longest_common_substring(
    seq1: &str,
    seq2: &str,
) -> Option<AlignmentResult> {
    let len1 = seq1.len();
    let len2 = seq2.len();
    
    if len1 == 0 || len2 == 0 {
        return None;
    }
    
    // Create DP table initialized to 0
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];
    let mut max_len = 0;
    let mut end_row = 0;
    let mut end_col = 0;
    
    // Fill DP table
    for i in 1..=len1 {
        for j in 1..=len2 {
            if seq1.as_bytes()[i - 1] == seq2.as_bytes()[j - 1] {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if dp[i][j] > max_len {
                    max_len = dp[i][j];
                    end_row = i;
                    end_col = j;
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }
    
    if max_len == 0 {
        return None;
    }
    
    // Extract the matching subsequence from seq1
    let start_pos = end_row - max_len;
    let subsequence: String = seq1.chars()
        .skip(start_pos)
        .take(max_len)
        .collect();
    
    // Calculate distances to 3' ends (following Python logic)
    // In Python: d1 = len1 - row, d2 = col - lcstr_len
    // where row = end_row, col = end_col, lcstr_len = max_len
    let d1 = len1 - end_row;
    let d2 = end_col - max_len;
    
    // Count GC in subsequence
    let gc_count = subsequence.chars()
        .filter(|c| *c == 'G' || *c == 'C' || *c == 'g' || *c == 'c')
        .count();
    
    Some(AlignmentResult {
        subsequence,
        length: max_len,
        gc_count,
        d1,
        d2,
        mismatches: 0,
    })
}

/// Calculate badness score based on alignment result
fn calculate_badness(alignment: &AlignmentResult) -> f64 {
    let len = alignment.length;
    let d1 = alignment.d1;
    let d2 = alignment.d2;
    
    // Calculate weight P
    let p = if d1 == 0 && d2 == 0 {
        100.0
    } else if d1 == 0 || d2 == 0 {
        10.0
    } else {
        1.0
    };
    
    p * (2.0_f64.powi(len as i32)) / ((d1 + 1) as f64 * (d2 + 1) as f64)
}

/// Determine risk level based on log10(badness)
fn risk_level(log_badness: f64) -> &'static str {
    if log_badness >= 3.5 {
        "高风险"
    } else if log_badness >= 3.0 {
        "中风险"
    } else {
        "低风险"
    }
}

/// Result for a primer pair evaluation
#[derive(Debug)]
struct PrimerPairResult {
    primer1_seq: String,
    primer2_seq: String,
    primer1_id: String,
    primer2_id: String,
    alignment: AlignmentResult,
    badness: f64,
}

impl PrimerPairResult {
    fn to_tsv_line(&self) -> String {
        let log_badness = self.badness.log10();
        let risk = risk_level(log_badness);
        format!("{}\t{}\t'{}'\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{:.4}\t{}",
                self.primer1_seq,
                self.primer2_seq,
                self.alignment.subsequence,
                self.alignment.length,
                self.alignment.d1,
                self.alignment.d2,
                self.badness,
                self.primer1_id,
                self.primer2_id,
                log_badness,
                risk)
    }
}

/// Evaluate all primer pairs
fn evaluate_primer_pairs(
    primers: &[Primer],
    max_mismatches: u32,
    badness_threshold: f64,
) -> Vec<PrimerPairResult> {
    let mut results = Vec::new();
    let num_primers = primers.len();
    
    for i in 0..num_primers {
        for j in (i + 1)..num_primers {
            let primer1 = &primers[i];
            let primer2 = &primers[j];
            
            // Get reverse complement of primer2
            let rc_primer2 = reverse_complement(&primer2.sequence);
            
            // Find longest common substring
            let alignment = if max_mismatches == 0 {
                // Exact match (same as Perl/Python version)
                find_exact_longest_common_substring(&primer1.sequence, &rc_primer2)
            } else {
                // Allow mismatches
                find_longest_common_substring_with_mismatches(
                    &primer1.sequence,
                    &rc_primer2,
                    max_mismatches,
                )
            };
            
            if let Some(alignment) = alignment {
                // Only consider alignments with length > 0
                if alignment.length > 0 {
                    let badness = calculate_badness(&alignment);
                    
                    // Filter by threshold
                    if badness > badness_threshold {
                        results.push(PrimerPairResult {
                            primer1_seq: primer1.sequence.clone(),
                            primer2_seq: primer2.sequence.clone(),
                            primer1_id: primer1.id.clone(),
                            primer2_id: primer2.id.clone(),
                            alignment,
                            badness,
                        });
                    }
                }
            }
        }
    }
    
    // Sort by badness (descending)
    results.sort_by(|a, b| b.badness.partial_cmp(&a.badness).unwrap());
    
    results
}

/// Evaluate all primer pairs in parallel using Rayon
fn evaluate_primer_pairs_parallel(
    primers: &[Primer],
    max_mismatches: u32,
    badness_threshold: f64,
) -> Vec<PrimerPairResult> {
    use rayon::prelude::*;
    
    let num_primers = primers.len();
    
    // Process combinations in parallel without storing all indices
    let results: Vec<Option<PrimerPairResult>> = (0..num_primers)
        .into_par_iter()
        .flat_map_iter(|i| {
            (i + 1..num_primers).map(move |j| (i, j))
        })
        .map(|(i, j)| {
            let primer1 = &primers[i];
            let primer2 = &primers[j];
            
            // Get reverse complement of primer2
            let rc_primer2 = reverse_complement(&primer2.sequence);
            
            // Find longest common substring
            let alignment = if max_mismatches == 0 {
                find_exact_longest_common_substring(&primer1.sequence, &rc_primer2)
            } else {
                find_longest_common_substring_with_mismatches(
                    &primer1.sequence,
                    &rc_primer2,
                    max_mismatches,
                )
            };
            
            alignment.and_then(|alignment| {
                // Only consider alignments with length > 0
                if alignment.length > 0 {
                    let badness = calculate_badness(&alignment);
                    
                    // Filter by threshold
                    if badness > badness_threshold {
                        Some(PrimerPairResult {
                            primer1_seq: primer1.sequence.clone(),
                            primer2_seq: primer2.sequence.clone(),
                            primer1_id: primer1.id.clone(),
                            primer2_id: primer2.id.clone(),
                            alignment,
                            badness,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
        })
        .collect();
    
    // Filter out None values and sort by badness (descending)
    let mut results: Vec<PrimerPairResult> = results.into_iter().flatten().collect();
    results.sort_by(|a, b| b.badness.partial_cmp(&a.badness).unwrap());
    
    results
}

/// Write results to output file
fn write_results(output_path: &PathBuf, results: &[PrimerPairResult]) -> Result<()> {
    use std::fs::File;
    use std::io::Write;
    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
    
    // Write header
    writeln!(file, "Primer_1\tPrimer_2\tsub_rev\trev.len\trev.d1\trev.d2\tbadness\tPrimer_1.id\tPrimer_2.id\tlog10(badness)\t风险等级")?;
    
    // Write results
    for result in results {
        writeln!(file, "{}", result.to_tsv_line())?;
    }
    
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Validate thread count
    if args.threads == 0 {
        anyhow::bail!("Thread count must be at least 1");
    }
    
    println!("Input file: {:?}", args.input);
    println!("Mismatch threshold: {}", args.mismatch);
    println!("Badness threshold: {}", args.threshold);
    
    // Set number of threads for Rayon
    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("Failed to initialize thread pool");
    println!("Using {} threads for parallel processing", args.threads);
    
    // Read primers
    let primers = read_primers(&args.input)
        .with_context(|| format!("Failed to read primers from {:?}", args.input))?;
    
    println!("Loaded {} primers", primers.len());
    
    // Determine output path
    let output_path = match &args.output {
        Some(path) => path.clone(),
        None => {
            let mut path = args.input.clone();
            path.set_file_name("oneTubePrimer.3EndBadness.txt");
            path
        }
    };
    println!("Output file: {:?}", output_path);
    
    // Evaluate primer pairs
    println!("Evaluating primer pairs...");
    let results = evaluate_primer_pairs_parallel(&primers, args.mismatch, args.threshold);
    println!("Found {} primer pairs with badness > {}", results.len(), args.threshold);
    
    // Write results
    write_results(&output_path, &results)
        .with_context(|| format!("Failed to write results to {:?}", output_path))?;
    
    println!("Results written to {:?}", output_path);
    
    // Print summary
    if !results.is_empty() {
        println!("\nTop results:");
        for (i, result) in results.iter().take(5).enumerate() {
            let log_badness = result.badness.log10();
            let risk = risk_level(log_badness);
            println!("{}. {} vs {}: badness={:.2}, log10={:.4}, {}",
                     i + 1,
                     result.primer1_id,
                     result.primer2_id,
                     result.badness,
                     log_badness,
                     risk);
        }
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("atcg"), "CGAT");
        assert_eq!(reverse_complement("AATTCCGG"), "CCGGAATT");
    }
    
    #[test]
    fn test_find_longest_common_substring_with_mismatches() {
        let seq1 = "ATCGATCG";
        let seq2 = "TAGCTAGC"; // Reverse complement of seq1
        
        // With 0 mismatches, should find full match
        let result = find_longest_common_substring_with_mismatches(seq1, seq2, 0);
        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.subsequence, "ATCGATCG");
        assert_eq!(result.length, 8);
        assert_eq!(result.mismatches, 0);
        
        // With 2 mismatches
        let seq3 = "ATCGXXXX";
        let result = find_longest_common_substring_with_mismatches(seq1, seq3, 2);
        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.subsequence, "ATCG");
        assert_eq!(result.length, 4);
        assert!(result.mismatches <= 2);
    }
    
    #[test]
    fn test_calculate_badness() {
        let alignment = AlignmentResult {
            subsequence: "ATCG".to_string(),
            length: 4,
            gc_count: 2,
            d1: 0,
            d2: 0,
            mismatches: 0,
        };
        
        // d1=0, d2=0, P=100
        let badness = calculate_badness(&alignment);
        assert!(badness > 0.0);
        
        let alignment2 = AlignmentResult {
            subsequence: "ATCG".to_string(),
            length: 4,
            gc_count: 2,
            d1: 2,
            d2: 3,
            mismatches: 0,
        };
        
        // d1=2, d2=3, P=1
        let badness2 = calculate_badness(&alignment2);
        assert!(badness2 > 0.0);
        assert!(badness2 < badness); // With distances > 0, badness should be smaller
    }
    
    #[test]
    fn test_risk_level() {
        assert_eq!(risk_level(3.6), "高风险");
        assert_eq!(risk_level(3.2), "中风险");
        assert_eq!(risk_level(2.5), "低风险");
    }
    
    #[test]
    fn test_exact_longest_common_substring() {
        // Test case from Perl output: AAATACGGTGGATTAAACAAAAGCAA vs CTCGGATGCTGTGGGTGTTT
        // Perl output: sub_rev='AAACA', len=5, d1=7, d2=0, badness=40.00
        
        let seq1 = "AAATACGGTGGATTAAACAAAAGCAA";
        let seq2_rev = "AACACCCACAGCATCCGAG"; // Reverse complement of CTCGGATGCTGTGGGTGTTT
        
        let result = find_exact_longest_common_substring(seq1, seq2_rev);
        assert!(result.is_some());
        let alignment = result.unwrap();
        
        // The longest common substring should be "AAACA" (5 bp)
        assert_eq!(alignment.subsequence, "AAACA");
        assert_eq!(alignment.length, 5);
        assert_eq!(alignment.d1, 7);
        assert_eq!(alignment.d2, 0);
        
        // Calculate badness
        let badness = calculate_badness(&alignment);
        assert!((badness - 40.0).abs() < 0.01); // Allow small floating point difference
    }
}
