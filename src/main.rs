use crate::avl_tree::SearchTree;
use basic_seed_chainer::avl_tree;
use basic_seed_chainer::seeding_methods;
use basic_seed_chainer::simulation_utils;
use bio::alignment::pairwise::*;
use clap::{App, Arg, SubCommand};
use rayon::prelude::*;
use statistical;
use std::cmp;
use std::str;
use std::sync::Mutex;
use std::time::Instant;
use std::collections::HashSet;

use std::fs::{OpenOptions, metadata};
use std::io::{BufWriter, Write};
use std::path::Path;

fn append_result_csv<P: AsRef<Path>>(
    path: P,
    n: usize,
    m: usize,
    k: usize,
    recoverability: f64,
    expected_runtime: f64,
) -> std::io::Result<()> {
    let new_file = !path.as_ref().exists() || metadata(&path)?.len() == 0;

    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)?;
    let mut w = BufWriter::new(file);

    if new_file {
        writeln!(w, "n,m,k,recoverability,expected_runtime")?;
    }

    writeln!(
        w,
        "{},{},{},{:.17},{:.17}",
        n, m, k, recoverability, expected_runtime
    )?;

    w.flush()?;
    Ok(())
}


fn main() {
    let matches = App::new("basic_seed_chainer")
        .version("0.1")
        .about("Basic seed chain extend aligner on simulated sequences.")
        .arg(
            Arg::with_name("sketch")
                .short("s")
                .help("Use sketching. Default density is 1/(k - 6) where k = C log n (default no sketching)"),
        )
        .arg(
            Arg::with_name("minimizer")
                .short("m")
                .help("Use minimizers instead of open syncmers"),
        )
        .arg(
            Arg::with_name("wfa")
                .long("wfa")
                .help("Use wavefront aligner instead of standard DP for extension"),
        )
        .arg(
            Arg::with_name("substring")
                .long("substring")
                .help("The mutated string is a substring instead of a full string (default align two full strings). m = sqrt(n) + 200 where n is the sequence length"),
        )
        .arg(
            Arg::with_name("num_iter")
                .help("Number of iterations per value of k")
                .required(true)
                .takes_value(true)
                .index(1),
        )
        .arg(
            Arg::with_name("k")
                .short("k")
                .takes_value(true)
                .help("Maximum value of k (default 20, minimum k is 10)"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .takes_value(true)
                .help("Number of threads (default 20)"),
        )
        .arg(
            Arg::with_name("theta")
                .required(true)
                .index(2)
                .takes_value(true)
                .help("Value of mutation rate theta"),
        )
        .get_matches();
    let print_all_debug = false;
    let use_minimizers = matches.is_present("minimizer");
    let use_wfa = matches.is_present("wfa");
    let sketch = matches.is_present("sketch");
    let m_is_substring = matches.is_present("substring");
    let num_iters = matches
        .value_of("num_iter")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let max_k = matches
        .value_of("k")
        .unwrap_or("20")
        .parse::<usize>()
        .unwrap();
    let threads = matches
        .value_of("threads")
        .unwrap_or("20")
        .parse::<usize>()
        .unwrap();
    let th = matches.value_of("theta").unwrap().parse::<f64>().unwrap();
    let thetas = vec![th];
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    for theta in thetas {
        let mut history_align = vec![];
        let mut history_chain = vec![];
        let alpha = -((1.0 - theta) as f64).log(4.0);
        let _C = 3. / (1. - 2. * alpha);
        let ks = 20..max_k; 
        println!("Theta = {}, alpha = {}, C = {}", theta, alpha, _C);
        for k in ks {
            let k = k as usize;
            let d;
            if sketch {
                d = k - 6;
            } else {
                d = 1;
            }
            let w = 10.;
            let window = 2 * d - 1;
            let coeff_scale = 1.0;
            let v = k - d + 1;
            let t = (k - v + 1) / 2 + 1;
            let n = 4.0_f64.powf(k as f64 * (1. - 2. * alpha) / (3. * coeff_scale)) + 10.; // adding an extra constant to n to avoid m > n
            let n = n as usize;
            let gamma = 0.50;
            let mut s = simulation_utils::gen_rand_string(n);
            let mut s_seeds = seeding_methods::open_sync_seeds(&s, k, v, t).0;
            let zeta_inter = 3. / (1. - 2. * alpha) * 50. / 8.
                * (n as f64).log(4.)
                * (n as f64).ln()
                * ((1. - theta) as f64).powi(-(k as i32));
            let zeta_sketch_inter = zeta_inter * 4. * 3. / 2.
                + 8. * zeta_inter / ((1. - theta) as f64).powi(-(k as i32));
            let c_0 = 3.;
            let zeta_inter = 1. / (n as f64);
            let zeta_sketch_inter = 1. / (6. * zeta_sketch_inter);
            let zeta;
            if sketch {
                zeta = zeta_sketch_inter;
            } else {
                zeta = zeta_inter;
            }
            let zeta = zeta * w;
            let align_times: Mutex<Vec<_>> = Mutex::new(vec![]);
            let chain_times: Mutex<Vec<_>> = Mutex::new(vec![]);
            let recov: Mutex<f64> = Mutex::new(0.);
            let s_p_avg: Mutex<f64> = Mutex::new(0.);
            let m;
            let start_ind;

            if m_is_substring {
                start_ind = 1;
                let avg = ((2.0 * _C * alpha) + 1.0) / 2.0;
                m = ((n as f64).powf(avg)) as usize;
            } else {
                start_ind = 0;
                m = n as usize;
            }
            let offset = 100;
            println!("n = {},  m = {}", n, m);
            
            (0..num_iters)
                .collect::<Vec<usize>>()
                .into_par_iter()
                .for_each(|_l| {
                    let s_repro;
                    if !m_is_substring {
                        s_repro = simulation_utils::gen_rand_string(n);
                    } else {
                        s_repro = simulation_utils::gen_rand_string(2 * k);
                    }
                    let mut s_p;
                    let mut hom_path;
                    if m_is_substring {
                        let result = simulation_utils::gen_mutated_string(
                            &s[start_ind..start_ind + m],
                            start_ind,
                            theta,
                            gamma
                        );
                        s_p = result.0;
                        hom_path = result.1;
                        let mut s_p_len = s_p_avg.lock().unwrap();
                        *s_p_len += (s_p.len() as f64)
                    } else {
                        let result = simulation_utils::gen_mutated_string(
                            &s_repro[start_ind..start_ind + m],
                            start_ind,
                            theta,
                            gamma
                        );
                        s_p = result.0;
                        hom_path = result.1;
                        let mut s_p_len = s_p_avg.lock().unwrap();
                        *s_p_len += (s_p.len() as f64)
                    }

                    let now = Instant::now();
                    let s_repro_seeds;
                    let sp_seeds;
                    if !use_minimizers {
                        s_repro_seeds = seeding_methods::open_sync_seeds(&s_repro, k, v, t).0;
                        sp_seeds = seeding_methods::open_sync_seeds(&s_p, k, v, t).0;
                    } else {
                        s_repro_seeds = seeding_methods::minimizer_seeds(&s_repro, window, k).0;
                        sp_seeds = seeding_methods::minimizer_seeds(&s_p, window, k).0;
                    }
                    if print_all_debug {
                        println!("Seeding time {}", now.elapsed().as_secs_f32());
                    }

                    let mut anchors = vec![];
                    anchors.reserve(n);

                    // FOR THE PURPOSES OF CALCULATING N_H VS EN_H -- WE DONT ACTUALLY NEED THE CHAIN

                    let now = Instant::now();
                    if m_is_substring {
                        for (kmer2, poses2) in sp_seeds.iter() {
                            if s_seeds.contains_key(kmer2) {
                                let poses1 = &s_seeds[kmer2];
                                for pos1 in poses1 {
                                    for pos2 in poses2 {
                                        anchors.push((*pos1, *pos2));
                                    }
                                }
                            }
                        }
                    } else {
                        for (kmer1, poses1) in s_repro_seeds.iter() {
                            if sp_seeds.contains_key(kmer1) {
                                let poses2 = &sp_seeds[kmer1];
                                for pos1 in poses1 {
                                    for pos2 in poses2 {
                                        anchors.push((*pos1, *pos2));
                                    }
                                }
                            }
                        }
                    }
                    anchors.sort();
                    if print_all_debug {
                        println!("Number of anchors is {}", anchors.len());
                        println!("Anchor finding time {}", now.elapsed().as_secs_f32());
                    }

                    let mut f = vec![w];
                    let mut pointer_array = vec![0; anchors.len()];
                    let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();

                    let now = Instant::now();
                    for (i, anchor) in anchors.iter().enumerate() {
                        avl_tree.insert([anchor.1, i]);
                    }
                    if (anchors.len() > 0) {
                    avl_tree.update_query_info(
                        [anchors[0].1, 0],
                        w + zeta * (anchors[0].1 + anchors[0].0) as f64,
                        0,
                        100 * anchors[0].0 as usize + offset,
                        100 * anchors[0].1 as usize + offset,
                    );
                }

                    let mut spur = false;
                    for i in 1..anchors.len() {
                        if anchors[i].1 > anchors[i].0 + 1000 || anchors[i].0 > anchors[i].1 + 1000
                        {
                            spur = true;
                        }
                        let best_f_i;
                        let best_j;

                        let (best_score, best_id) = avl_tree.mrq(
                            [0, 0],
                            [anchors[i].1, i],
                            100 * anchors[i].0 as usize + offset,
                            100 * anchors[i].1 as usize + offset,
                        );
                        if spur {
                        }
                        if best_score == i64::MIN {
                            best_f_i = w;
                            best_j = i;
                        } else {
                            best_j = best_id;
                            best_f_i =
                                best_score as f64 - zeta * (anchors[i].0 + anchors[i].1) as f64 + w;
                        }

                        f.push(best_f_i);
                        avl_tree.update_query_info(
                            [anchors[i].1, i],
                            best_f_i + zeta * (anchors[i].0 + anchors[i].1) as f64,
                            i,
                            100 * anchors[i].0 as usize + offset,
                            100 * anchors[i].1 as usize + offset,
                        );

                        if best_j != usize::MAX {
                            pointer_array[i] = best_j;
                        }
                    }
                    //            println!("Chaining time {}", now.elapsed().as_secs_f32());
                    chain_times
                        .lock()
                        .unwrap()
                        .push(now.elapsed().as_secs_f32());

                    let mut vec: Vec<_> = f.iter().enumerate().collect();
                    vec.sort_by(|(_, v0), (_, v1)| v1.partial_cmp(v0).unwrap());

                    let mut best_chain = Vec::new();

                    if let Some(&(mut curr_i, _)) = vec.get(0) {
                        while let Some(&prev_i) = pointer_array.get(curr_i) {
                            if prev_i == curr_i { break; }               
                            if let Some(a) = anchors.get(curr_i) {
                                best_chain.push(a.clone());               
                            } else { break; }                            
                            curr_i = prev_i;
                        }
                        if let Some(a) = anchors.get(curr_i) {          
                            best_chain.push(a.clone());
                        }
                    }

                    let mut break_length = 0;
                    let mut break_start = false;
                    let mut break_a = 0;
                    let mut break_b = 0;
                    if best_chain.len() > 1 {
                        let range = best_chain[0].1 - best_chain[best_chain.len() - 1].1;
                        for anchors in best_chain.iter() {
                            if anchors.0 != anchors.1 + start_ind {
                                if break_start == false {
                                    break_start = true;
                                    break_a = cmp::min(anchors.0, anchors.1);
                                    break_b = cmp::max(anchors.0, anchors.1);
                                } else {
                                    let min = cmp::min(anchors.0, anchors.1);
                                    let max = cmp::max(anchors.0, anchors.1);
                                    if min < break_a {
                                        break_a = min;
                                    }
                                    if max > break_b {
                                        break_b = max;
                                    }
                                }
                            } else {
                                if break_start {
                                    break_start = false;
                                    break_length += break_b - break_a;
                                }
                            }
                        }
                        
                        if break_length > 0 {
                            if print_all_debug {
                                dbg!(break_length);
                            }
                        }
                    };
                    let mut rec = recov.lock().unwrap();
                    // recoverability calculations
                    let recoverability = simulation_utils::get_recoverability(&best_chain, &hom_path, k);
                    *rec += recoverability as f64;
                    let mut gap_intervals = vec![];
                    let mut curr_r_end = 0;
                    let mut curr_q_end = 0;
                    for (iter, anchor) in best_chain.iter().rev().enumerate() {
                        if iter > 0 {
                            let rgap;
                            let qgap;
                            if anchor.0 > curr_r_end + 1 && anchor.1 > curr_q_end + 1 {
                                rgap = (curr_r_end, anchor.0 - 1);
                                qgap = (curr_q_end, anchor.1 - 1);
                                gap_intervals.push((rgap, qgap));
                            } else {
                            }
                        }
                        curr_r_end = anchor.0 + k - 1;
                        curr_q_end = anchor.1 + k - 1;
                    }

                    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                    let now = Instant::now();
                    let mut aligner = Aligner::new(-5, -1, &score);
                    for gap in gap_intervals {
                        let ref_slice = &s[gap.0 .0..gap.0 .1];
                        let query_slice = &s[gap.1 .0..gap.1 .1];
                        if !use_wfa {
                            let _alignment =
                                aligner.global(&s[gap.0 .0..gap.0 .1], &s_p[gap.1 .0..gap.1 .1]);
                        } else {
                            //
                        }
                    }

                    align_times
                        .lock()
                        .unwrap()
                        .push(now.elapsed().as_secs_f32());
                });
            let mut align_times = align_times.into_inner().unwrap();
            let mut chain_times = chain_times.into_inner().unwrap();
            let align_mean = align_times.iter().sum::<f32>();
            let chain_mean = chain_times.iter().sum::<f32>();
            history_align.push(align_times);
            history_chain.push(chain_times);

            let recov_val = *recov.lock().unwrap() / (num_iters as f64);
            let sp_av = ((*s_p_avg.lock().unwrap() / (num_iters as f64)).round()) as usize;;
            println!("Value for k {}", k);
            println!(
                "Mean lower bound for recoverability for k {} is {}",
                k, recov_val
            );
            let align_mean_per_iter = align_mean as f64 / num_iters as f64;
            let chain_mean_per_iter = chain_mean as f64 / num_iters as f64;
            let expected_runtime = align_mean_per_iter + chain_mean_per_iter;
            append_result_csv("recoverability.csv", n, m, k, recov_val, expected_runtime); 
            println!("Mean extend time {}", align_mean / num_iters as f32);
            println!("Mean chain time {}", chain_mean / num_iters as f32);
        }
    }
}
