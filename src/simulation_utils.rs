use fxhash::{FxHashMap, FxHashSet};
use rand::distributions::{Bernoulli, Distribution, Uniform};
use rand::{thread_rng, Rng};
use std::collections::HashSet;
use rand_distr::Geometric;
use rayon::prelude::*;

pub fn gen_rand_string(n: usize) -> Vec<u8> {
    const CHARSET: &[u8] = b"ATCG";
    let mut rng = thread_rng();
    let return_string: Vec<u8> = (0..n)
        .map(|_| {
            let idx = rng.gen_range(0..CHARSET.len());
            CHARSET[idx] as u8
        })
        .collect();
    return return_string;
}

fn dirichlet_111<R: Rng>(rng: &mut R) -> (f64, f64, f64) {
    let u1: f64 = rng.gen::<f64>().clamp(f64::MIN_POSITIVE, 1.0);
    let u2: f64 = rng.gen::<f64>().clamp(f64::MIN_POSITIVE, 1.0);
    let u3: f64 = rng.gen::<f64>().clamp(f64::MIN_POSITIVE, 1.0);
    let g1 = -u1.ln();
    let g2 = -u2.ln();
    let g3 = -u3.ln();
    let s = g1 + g2 + g3;
    (g1 / s, g2 / s, g3 / s)
}

pub fn gen_mutated_string(
    sequence: &[u8],
    offset: usize,
    theta: f64,
    gamma: f64,
) -> (Vec<u8>, HashSet<(usize, usize)>) {
    let mut rng = rand::thread_rng();
    const CHARSET: &[u8] = b"ATCG";
    let p = theta / 3.0;
    let (w_i, w_s, w_d) = dirichlet_111(&mut rng); 
    let theta_i = theta * w_i;
    let theta_s = theta * w_s;
    let theta_d = theta * w_d; 
    let mut sp: Vec<u8> = Vec::new();
    let mut path: HashSet<(usize, usize)> = HashSet::new();
    let mut x = offset;
    let mut y = 0;
    path.insert((x, y));

    for &c in sequence {
        // let ins = rng.gen::<f64>() < p;
        let ins = rng.gen::<f64>() < theta_i;
        let mut ins_len = 0;
        if ins {
            let geom = Geometric::new(gamma).unwrap();
            ins_len = geom.sample(&mut rng) as usize + 1;
            for _ in 0..ins_len {
                let ins_chr = CHARSET[rng.gen_range(0..CHARSET.len())];
                sp.push(ins_chr);
            }
        }
        sp.push(c);
        // if rng.gen::<f64>() < p
        if rng.gen::<f64>() < theta_s {
            loop {
                let cand = CHARSET[rng.gen_range(0..CHARSET.len())];
                if cand != c {
                    sp.pop();
                    sp.push(cand);
                    break;
                }
            }
        }
        // let del = if rng.gen::<f64>() < p 
        let del = if rng.gen::<f64>() < theta_d {
            sp.pop();
            true
        } else {
            false
        };
    
        if del && !ins {
            x += 1;
            path.insert((x, y));
        } else if !del && ins {
            for _ in 0..ins_len {
                y += 1;
                path.insert((x, y));
            }
            x += 1;
            y += 1;
            path.insert((x, y));
        } else if del && ins {
            for _ in 0..ins_len {
                y += 1;
                path.insert((x, y));
            }
            x += 1;
            path.insert((x, y));
        } else {
            x += 1;
            y += 1;
            path.insert((x, y));
        }
    }

    (sp, path)
}

pub fn extend_gap(anchor_l: (usize, usize), anchor_r: (usize, usize), k: usize) -> HashSet<(usize, usize)> {
    let (xl, yl) = anchor_l;
    let (xr, yr) = anchor_r;

    let mut extension = HashSet::with_capacity((xr - xl) * (yr - yl));
    for x in (xl + k)..xr {
        for y in (yl + k)..yr {
            extension.insert((x, y));
        }
    }
    extension
}

pub fn extend_gap_points(anchor_l: (usize, usize), anchor_r: (usize, usize), k: usize) -> impl Iterator<Item = (usize, usize)> {
    let (xl, yl) = anchor_l;
    let (xr, yr) = anchor_r;

    (xl + k..xr).flat_map(move |x| {
        (yl + k..yr).map(move |y| (x, y))
    })
}

pub fn get_align_c(chain: &[(usize, usize)], k: usize) -> HashSet<(usize, usize)> {
    chain
        .par_iter()
        .flat_map_iter(|&(x, y)| (0..k).map(move |l| (x + l, y + l))) // anchor spans
        .chain(
            (0..chain.len().saturating_sub(1))
                .into_par_iter()
                .flat_map_iter(|i| {
                    let anchor_l = chain[i + 1];
                    let anchor_r = chain[i];
                    extend_gap_points(anchor_l, anchor_r, k)
                }),
        )
        .collect()
}


pub fn overlap_ratio(align_c: &HashSet<(usize, usize)>, hom_path: &HashSet<(usize, usize)>) -> f64 {
    let intersection_size = align_c.intersection(&hom_path).count() as f64;
    let denom = hom_path.len() as f64 + 1e-3;
    intersection_size / denom
}

pub fn get_hom_count(
    hom_path: &HashSet<(usize, usize)>,
    s: &[u8],
    sp: &[u8],
    k: usize,
    offset: usize,
) -> usize {
    let mut hom_num = 0;
    let mut path_diff = 0.0f64;
    let mut xp = 0;
    let mut yp = 0;

    for &(x, y) in hom_path.iter() {
        let mut in_count = 0;
        let mut letters_match = true;

        for l in 0..k {
            if hom_path.contains(&(x + l, y + l)) {
                if x + l < s.len() && y + l < sp.len() {
                    if s[x + l] != sp[y + l] {
                        letters_match = false;
                        break;
                    }
                } else {
                    letters_match = false;
                    break;
                }
                in_count += 1;
            }
        }

        if in_count == k && letters_match {
            hom_num += 1;
        }

        path_diff = x as f64 - y as f64 - offset as f64;
        xp = x;
        yp = y;
    }
    hom_num
}


pub fn get_recoverability(
    best_chain: &[(usize, usize)],
    hom_path: &HashSet<(usize, usize)>,
    k: usize,
) -> (f64) {
    if hom_path.is_empty() {
        return 0.0;
    }
    let actual_first_anchor_start = best_chain.last().copied();
    let actual_last_anchor_start  = best_chain.first().copied();

    let pre: usize = match actual_first_anchor_start {
        Some((ax, ay)) => hom_path
            .iter()
            .filter(|&&(x, y)| {
                (x <= ax && y <= ay) && !(x == ax && y == ay)
            })
            .count(),
        None => 0,
    };
    let post: usize = match actual_last_anchor_start {
        Some((lx, ly)) => {
            let end_x = lx.saturating_add(k.saturating_sub(1));
            let end_y = ly.saturating_add(k.saturating_sub(1));
            hom_path
                .iter()
                .filter(|&&(x, y)| {
                    (x >= end_x && y >= end_y) && !(x == end_x && y == end_y)
                })
                .count()
        }
        None => 0,
    };

    let total_hp = hom_path.len();
    let inside   = total_hp.saturating_sub(pre + post);
    let recoverability = if total_hp > 0 {
        (inside as f64 / total_hp as f64).clamp(0.0, 1.0) // = 1 - (pre+post)/total
    } else {
        0.0
    };

    (
        recoverability
    )
}

pub fn check_context_dependent_mutation(
    seeds_orig: &FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: &FxHashMap<&[u8], FxHashSet<usize>>,
    string_mut: &[u8],
    string_orig: &[u8],
    k: usize,
) {
    let mut num_context_mutations = 0;
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let difference: FxHashSet<_> = set1.difference(&set2).collect();

    //We only take a subset of context mutated k-mers,
    //if somehow a mismatch k-mer was generated through a mutation,
    //this won't capture it.
    for kmer in difference {
        let positions_s = seeds_orig.get(kmer).unwrap();
        for pos in positions_s {
            //dbg!(String::from_utf8(string_mut[*pos..*pos+k].to_vec()),String::from_utf8(kmer.to_vec()));
            if string_mut[*pos..*pos + k] == **kmer {
                num_context_mutations += 1;
            }
            if string_orig[*pos..*pos + k] != **kmer {
                panic!("something went wrong");
            }
        }
    }

    dbg!(num_context_mutations);
}

pub fn get_conservation(
    seeds_orig: FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: FxHashMap<&[u8], FxHashSet<usize>>,
    k: usize,
    str_len: usize,
) -> (f64, f64) {
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let intersection: FxHashSet<_> = set1.intersection(&set2).collect();
    let mut homologous_positions: FxHashSet<usize> = FxHashSet::default();
    let mut spurious_positions: FxHashSet<usize> = FxHashSet::default();
    //If we get the same k-mer on both strings, make sure that
    //they are at the same positions.
    //dbg!(intersection.len());
    for kmer in intersection {
        let positions_s = seeds_orig.get(kmer).unwrap();
        let positions_sp = seeds_mut.get(kmer).unwrap();

        let int_pos: FxHashSet<_> = positions_s.intersection(&positions_sp).collect();
        for pos in int_pos {
            for i in 0..k {
                homologous_positions.insert(pos + i);
            }
        }

        let bad_matches = positions_sp.difference(&positions_s).collect::<Vec<_>>();
        for pos in bad_matches {
            for i in 0..k {
                spurious_positions.insert(pos + i);
            }
        }
    }

    //dbg!(num_spurious_matches);
    return (
        homologous_positions.len() as f64 / str_len as f64,
        spurious_positions.len() as f64 / str_len as f64,
    );
}

pub fn get_schem_prob_minimizers(
    seeds_orig: &FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: &FxHashMap<&[u8], FxHashSet<usize>>,
    positions_changed: &Vec<bool>,
    k: usize,
) -> Vec<f64> {
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let intersection: FxHashSet<_> = set1.intersection(&set2).collect();
    let mut union = FxHashSet::default();

    for key in intersection {
        let positions = seeds_orig.get(key).unwrap();
        for pos in positions {
            union.insert(pos);
        }
    }

    let mut successes = vec![0; k];
    let mut total_counts = vec![0; k];
    let mut running_length = 0;
    for (i, pos) in positions_changed.iter().enumerate() {
        if *pos == true && running_length < 2 * k - 1 {
            running_length += 1;
        } else {
            if running_length >= k {
                for j in i - running_length..i - k + 1 {
                    if union.contains(&j) {
                        successes[running_length - k] += 1;
                        break;
                    }
                }
                total_counts[running_length - k] += 1;
            }
            if *pos == false {
                running_length = 0;
            }
        }
    }

    let mut probabilities = vec![0.0; k];
    for i in 0..k {
        probabilities[i] = successes[i] as f64 / total_counts[i] as f64;
    }

    return probabilities;

    //dbg!(probabilities,total_counts);
}

pub fn get_kmers_from_string<'a>(
    string: &'a [u8],
    k: usize,
    pos_orig: &'a [usize],
) -> FxHashMap<Vec<u8>, FxHashSet<(usize, usize)>> {
    let mut kmers = FxHashMap::default();

    for i in 0..string.len() - k as usize + 1 {
        let pos_vec = kmers
            .entry(string[i..i + k as usize].to_vec())
            .or_insert(FxHashSet::default());
        pos_vec.insert((pos_orig[i], pos_orig[i + k - 1]));
    }

    let mut num_items = 0;
    let mut mult_items = 0;
    for (_key, value) in kmers.iter() {
        num_items += 1;
        mult_items += value.len();
    }

    println!("Avg dup kmers:{}", (mult_items as f64) / (num_items as f64));

    return kmers;
}

pub fn get_conservation_gap(
    seeds_orig: FxHashMap<Vec<u8>, FxHashSet<(usize, usize)>>,
    seeds_mut: FxHashMap<Vec<u8>, FxHashSet<(usize, usize)>>,
    _k: usize,
    str_len: usize,
) -> (f64, f64) {
    let mut num_spurious_matches = 0;
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let intersection: FxHashSet<_> = set1.intersection(&set2).collect();
    let mut homologous_positions: FxHashSet<usize> = FxHashSet::default();
    //If we get the same k-mer on both strings, make sure that
    //they are at the same positions.
    //dbg!(intersection.len());
    for kmer in intersection {
        let positions_s = seeds_orig.get(kmer).unwrap();
        let positions_sp = seeds_mut.get(kmer).unwrap();

        let int_pos: FxHashSet<_> = positions_s.union(&positions_sp).collect();
        for pos in int_pos {
            for i in pos.0..pos.1 + 1 {
                homologous_positions.insert(i);
            }
        }

        let bad_matches = positions_sp.difference(&positions_s).collect::<Vec<_>>();
        if bad_matches.len() > 0 {
            //dbg!(&positions_s, &positions_sp);
        }
        num_spurious_matches += bad_matches.len();
    }

    //dbg!(num_spurious_matches);
    return (
        homologous_positions.len() as f64 / str_len as f64,
        num_spurious_matches as f64 / seeds_orig.len() as f64,
    );
}
