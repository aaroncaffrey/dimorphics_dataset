﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class misc
    {
        public const string module_name = nameof(misc);

        internal static List<subseq_match> find_subsequences(string pdb_id, string subsequence, bool fix)
        {
            var pdb = atom.load_atoms_pdb(pdb_id, new load_atoms_pdb_options()
            {
                find_3d_intramolecular = false,

                load_1d_blast_pssms = false,
                load_1d_iup_data = false,
                load_1d_sable = false,
                load_1d_dna_binding = false,

                load_2d_mpsa_sec_struct_predictions = false,

                load_3d_rsa_data = false,
                load_3d_dssp_data = true,
                load_3d_stride_data = true,
                load_3d_ring_data = false,
                load_3d_foldx_ala_scan = false,
            });

            var chain_order = pdb.First().pdb_model_chain_atoms.First().chain_order_in_pdb_file;

            var alphabetical_chain_order = chain_order.OrderBy(a => a.chain_id).Select((a, i) => (array_index: i, chain_id: a.chain_id)).ToList();

            var result = new List<subseq_match>();

            for (var k = 0; k < pdb.Count; k++)
            {
                var master_atoms = atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, pdb[k].pdb_model_chain_atoms);

                pdb[k].pdb_model_chain_atoms.Clear();
                pdb[k].pdb_model_chain_atoms.AddRange(master_atoms);

            }

            if (fix)
            {

                var seqs = pdb.Select(ch =>
                {
                    var atoms = ch.pdb_model_chain_atoms;
                    var master_atoms = atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, ch.pdb_model_chain_atoms);
                    var seq = String.Concat(master_atoms.Select(a => a.amino_acid).ToList());
                    return seq;
                }).ToList();


                for (var i = 0; i < seqs.Count; i++)
                {
                    for (var j = 0; j < seqs.Count; j++)
                    {
                        if (i == j) continue;

                        var seqi = seqs[i];
                        var seqj = seqs[j];

                        var chain_id = pdb[i].chain_id;

                        var residsj = pdb[j].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();
                        var residsi = pdb[i].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();

                        // algorithm 1 to deal with missing terminals on one chain
                        if (seqj.Length > seqi.Length && seqj.Contains(seqi, StringComparison.Ordinal))
                        {
                            var pre_final_index = seqj.IndexOf(seqi, StringComparison.Ordinal);
                            var post_first_index = pre_final_index + seqi.Length;

                            var pre = seqj.Substring(0, pre_final_index);
                            var post = seqj[post_first_index..];

                            var pre_len = pre.Length;
                            var post_len = post.Length;



                            var resids_before = residsj.Take(pre_len).OrderBy(a => a).ToList();
                            residsj.Reverse();
                            var resids_after = residsj.Take(post_len).OrderBy(a => a).ToList();
                            residsj.Reverse();

                            var i_first_res_id = residsi.Min();
                            var i_last_res_id = residsi.Max();

                            if (resids_before.Count > 0)
                            {
                                var resids_before_last_id = resids_before.Max();

                            }

                            if (resids_after.Count > 0)
                            {
                                var resids_after_last_id = resids_after.Min();
                            }

                            pdb[i].pdb_model_chain_atoms.InsertRange(0, pdb[j].pdb_model_chain_atoms.Where(a => resids_before.Contains(a.residue_index)).ToList());
                            pdb[i].pdb_model_chain_atoms.AddRange(pdb[j].pdb_model_chain_atoms.Where(a => resids_after.Contains(a.residue_index)).ToList());

                            residsi = pdb[i].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();
                            // note: chain id will still be wrong


                            //io_proxy.WriteLine();
                        }

                        // algorithm 2 to deal with missing termins on one chain
                        var x = residsi.Intersect(residsj).ToList();

                        var pct_intersect = (double)x.Count / (double)((residsi.Count + residsj.Count) / 2);

                        if (pct_intersect >= 0.75)
                        {
                            var i_missing_resids = residsj.Except(residsi).ToList();

                            pdb[i].pdb_model_chain_atoms.AddRange(pdb[j].pdb_model_chain_atoms.Where(a => i_missing_resids.Contains(a.residue_index)).ToList());

                            var z = pdb[i].pdb_model_chain_atoms.OrderBy(a => a.residue_index).ToList();

                            pdb[i].pdb_model_chain_atoms.Clear();
                            pdb[i].pdb_model_chain_atoms.AddRange(z);

                        }
                    }
                }
            }

            for (var k = 0; k < pdb.Count; k++)
            {
                var ch = pdb[k];

                var atoms = ch.pdb_model_chain_atoms;
                var master_atoms = atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, ch.pdb_model_chain_atoms);

                var seq = String.Concat(master_atoms.Select(a => a.amino_acid).ToList());
                var dssp = String.Concat(master_atoms.Select(a => a.dssp3_multimer).ToList());
                var stride = String.Concat(master_atoms.Select(a => a.stride3_multimer).ToList());

                var match_indexes = feature_calcs.IndexOfAll(seq, subsequence);

                if (match_indexes == null || match_indexes.Count == 0) continue;

                foreach (var match_index in match_indexes)
                {
                    var match_atoms = master_atoms.Where((a, i) => i >= match_index && i < match_index + subsequence.Length).ToList();

                    //var match_seq = seq.Substring(match_index, subsequence.Length);
                    //var match_dssp = dssp.Substring(match_index, subsequence.Length);
                    //var match_stride = stride.Substring(match_index, subsequence.Length);

                    var match_seq = String.Concat(match_atoms.Select(a => a.amino_acid).ToList());
                    var match_dssp = String.Concat(match_atoms.Select(a => a.dssp3_multimer).ToList());
                    var match_stride = String.Concat(match_atoms.Select(a => a.stride3_multimer).ToList());

                    var res_id_first = match_atoms.Min(a => a.residue_index);


                    if (!string.Equals(match_seq, subsequence, StringComparison.Ordinal)) throw new Exception();

                    var match = new subseq_match()
                    {
                        pdb_id = ch.pdb_id,
                        chain_id = ch.chain_id,
                        pdb_model_index = ch.pdb_model_index,
                        chain_number = chain_order.First(a => a.chain_id == ch.chain_id).array_index,
                        alphabetical_chain_number = alphabetical_chain_order.First(a => a.chain_id == ch.chain_id).array_index,
                        seq = seq,
                        subsequence = match_seq,
                        dssp = match_dssp,
                        stride = match_stride,
                        residue_index = res_id_first,
                        array_index = match_index,
                        total_matches = match_indexes.Count

                    };

                    result.Add(match);
                }

            }


            return result;

        }



        internal static void fill_missing_chains(CancellationTokenSource cts = null)
        {
            const string method_name = nameof(fill_missing_chains);

            using var i_cts = new CancellationTokenSource();
            if (cts == null) cts = i_cts;

            var lines_all = io_proxy.ReadAllLines(/*program.string_debug*/(@"c:\phd\bioinf\fill_missing_chain.csv"), module_name, method_name);

            var lines_header = lines_all.First();

            var lines = lines_all.Skip(1).Select(b =>
            {

                var a = b.Split(',');

                var i = 0;
                return new class_file_line()
                {
                    pair_id = a[i++].Trim(),
                    pdb_id = a[i++].Trim(),
                    dimer_type = a[i++].Trim(),
                    class_name = a[i++].Trim(),
                    symmetry = a[i++].Trim(),
                    parallelism = a[i++].Trim(),
                    strand_seq_merged = a[i++].Trim(),
                    strand_seq_unmerged = a[i++].Trim(),
                    res_id = a[i++].Trim(),
                    chain_colour = a[i++].Trim()
                };

            }).Distinct().ToList();

            //0 PDB ID
            //1 Dimer Type
            //2 Class
            //3 Symmetry
            //4 Parallelism
            //5 Strand Sequence
            //6 Res ID

            // first, fill in the 'missing' values, since those are the 'first' 
            // that won't work, because it might be a heterodimer, with different sequences, but both containing the strand subsequence

            // 

            //lines = lines.Skip(lines.FindIndex(a=>a.pdb_id.Equals(/*program.string_debug*/($@"2NS9", StringComparison.CurrentCultureIgnoreCase))).ToList();
            //lines = lines.Where(a => a.res_id.Length > 0).ToList();


            var errors = 0;

            var tasks = new List<Task>();
            var tasks_start_time = DateTime.Now;
            
            for (var i = 0; i < lines.Count; i++)
            {
                var line = lines[i];


                io_proxy.WriteLine((string) line.pdb_id);

                var task = Task.Run(() =>
                {
                    if (cts != null && cts.IsCancellationRequested) return;

                    var fix = false;

                    do
                    {
                        if (cts != null && cts.IsCancellationRequested) return;

                        var lookup = find_subsequences(line.pdb_id, line.strand_seq_merged, fix);
                        fix = false;

                        // filter lookup results by res id, if specified
                        if (line.res_id.Length > 0)
                        {

                            var lookup2 = lookup.Where(a => a.residue_index == Int32.Parse((string) line.res_id, NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).ToList();

                            if (lookup2.Count == 0) throw new Exception();


                            lookup = lookup2;
                        }


                        // filter by pymol chain colour, 0 is red, 1 is pink

                        if (string.Equals(line.chain_colour, /*program.string_debug*/($@"Red"), StringComparison.OrdinalIgnoreCase) && lookup.Any(a => a.chain_number == 0))
                        {
                            var colour_lookup = lookup.Where(a => a.chain_number == 0).ToList();

                            if (colour_lookup.Count > 1)
                                io_proxy.WriteLine();

                            line.res_id = colour_lookup.First().residue_index.ToString(CultureInfo.InvariantCulture);
                            line.chain_id = colour_lookup.First().chain_id.ToString(CultureInfo.InvariantCulture);
                        }

                        else if (string.Equals(line.chain_colour, /*program.string_debug*/($@"Pink"), StringComparison.OrdinalIgnoreCase) && lookup.Any(a => a.chain_number == 1))
                        {
                            var colour_lookup = lookup.Where(a => a.chain_number == 1).ToList();

                            if (colour_lookup.Count > 1)
                                io_proxy.WriteLine();

                            line.res_id = colour_lookup.First().residue_index.ToString(CultureInfo.InvariantCulture);
                            line.chain_id = colour_lookup.First().chain_id.ToString(CultureInfo.InvariantCulture);
                        }

                        else
                        {
                            errors++;
                            io_proxy.WriteLine(/*program.string_debug*/($@"errors: {errors}, fix: {fix}"), module_name, method_name);
                            //throw new Exception();
                            fix = true;
                        }

                    } while (fix);

                }, cts.Token);

                tasks.Add(task);

                program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, -1, module_name, method_name, cts);

            }

            //Task.WaitAll(tasks.ToArray());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, 0, module_name, method_name, cts);

            var result = lines.Select(a => String.Join(/*program.string_debug*/($@","), new string[] { a.pair_id, a.pdb_id, a.dimer_type, a.class_name, a.symmetry, a.parallelism, a.strand_seq_merged, a.strand_seq_unmerged, a.res_id, a.chain_colour, a.chain_id })).ToList();
            result.Insert(0, /*program.string_debug*/($@"{lines_header},chain_id"));

            var fn = @"c:\phd\bioinf\fixed_dataset.csv";

            //Directory.CreateDirectory(Path.GetDirectoryName(fn));

            io_proxy.WriteAllLines(fn, result, nameof(program), nameof(fill_missing_chains));
            io_proxy.WriteLine(/*program.string_debug*/($@"finished."), module_name, method_name);
            //Console.ReadLine();

        }
        internal static void extract_all_pdbs()
        {
            var pdbs = Directory.GetFiles(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb")), /*program.string_debug*/($@"*.pdb")).ToList();

            foreach (var p in pdbs)
            {
                io_proxy.WriteLine(/*program.string_debug*/($@"Extracting {p}"));
                atom.extract_split_pdb_chains(p, null);
            }
        }

        internal static void repair_all_extracted_pdbs(bool run, CancellationTokenSource cts = null)
        {
            const string module_name = nameof(misc);
            const string method_name = nameof(repair_all_extracted_pdbs);

            using var i_cts = new CancellationTokenSource();
            if (cts == null) cts = i_cts;

            var split_pdbs = Directory.GetFiles(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb_split")), /*program.string_debug*/($@"*.pdb")).ToList();

            var tasks = new List<Task>();
            var tasks_start_time = DateTime.Now;

            foreach (var pp in split_pdbs)
            {
                var p = pp;

                var task = Task.Run(() =>
                {
                    if (cts != null && cts.IsCancellationRequested) return;
                    
                    io_proxy.WriteLine(/*program.string_debug*/($@"Repairing {p}"));
                    info_foldx.foldx_repair_pdb(Path.GetFileNameWithoutExtension(p), run, /*program.string_debug*/($@"{Path.GetDirectoryName(p)}\"));
                }, cts.Token);

                tasks.Add(task);

                program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, -1, module_name, method_name, cts);

                //program.print_eta(tasks.Count(a => a.IsCompleted), split_pdbs.Count, tasks_start_time, module_name, method_name);
            }

            //Task.WaitAll(tasks.ToArray());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, 0, module_name, method_name, cts);

        }

        internal static List<string> get_pdb_sequences(bool limited_to_dimorphics_and_standard = true)
        {

            var dimorphics_data = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"csv"), /*program.string_debug*/($@"distinct dimorphics list.csv")), nameof(misc), nameof(get_pdb_sequences))
                .Skip(1).Where(a => !String.IsNullOrWhiteSpace(a.Replace(/*program.string_debug*/($@","), /*program.string_debug*/($@""), StringComparison.Ordinal))).Select((a, i) =>
                {
                    var x = a.Split(',');
                    return (
                        pdb_id: x[0].ToUpperInvariant(),
                        dimer_type: x[1],
                        class_name: x[2],
                        symmetry_mode: x[3],
                        parallelism: x[4],
                        chain_number: Int32.Parse(x[5], NumberStyles.Integer, NumberFormatInfo.InvariantInfo) - 1,
                        strand_seq: x[6],
                        optional_res_index: x[7]
                    );
                }).ToList();

            if (limited_to_dimorphics_and_standard)
            {
                dimorphics_data = dimorphics_data.Where(a => (string.Equals(a.class_name, /*program.string_debug*/($@"Dimorphic"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.class_name, /*program.string_debug*/($@"Single"), StringComparison.OrdinalIgnoreCase)) || (string.Equals(a.class_name, /*program.string_debug*/($@"Standard"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.class_name, /*program.string_debug*/($@"Multiple"), StringComparison.OrdinalIgnoreCase))).ToList();
            }


            var pdb_id_list = dimorphics_data.Select(a => a.pdb_id.Trim().ToUpperInvariant()).Distinct().ToList();


            var pdb_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb"));


            var pdbs = Directory.GetFiles(pdb_folder, /*program.string_debug*/($@"*.pdb")).ToList();

            var seqs = new List<string>();

            var reverse_lookup = new List<(string pdb, string seq)>();
            foreach (var pdb_id in pdbs)
            {
                var pdb_id_code = Path.GetFileNameWithoutExtension(pdb_id).Trim().ToUpperInvariant().Substring(0, 4);

                //if (!pdb_id_list.Contains(pdb_id_code)) continue;

                io_proxy.WriteLine(/*program.string_debug*/($@"Extracting sequence from {pdb_id}"));


                var chains = atom.load_atoms_pdb(pdb_id,
                    new load_atoms_pdb_options()
                    {
                        find_3d_intramolecular = true,

                        load_1d_blast_pssms = true,
                        load_1d_iup_data = true,
                        load_1d_sable = true,
                        load_1d_dna_binding = true,

                        load_2d_mpsa_sec_struct_predictions = true,

                        load_3d_rsa_data = true,
                        load_3d_dssp_data = true,
                        load_3d_stride_data = true,
                        load_3d_ring_data = true,
                        load_3d_foldx_ala_scan = true,
                    }, pdb_folder);

                foreach (var chain in chains)
                {


                    var pdb_seqs = atom.amino_acid_sequence(pdb_id, chain.pdb_model_chain_atoms);

                    foreach (var s in pdb_seqs)
                    {
                        seqs.Add(s.aa_sequence);
                        reverse_lookup.Add((chain.pdb_id + chain.chain_id, s.aa_sequence));
                    }
                }
            }

            var output_filename = @"c:\phd\bioinf\betastrands_dataset_sequences.txt";
            var rename_filename = @"c:\phd\bioinf\rename_num_to_pdb.bat";


            if (limited_to_dimorphics_and_standard)
            {
                output_filename = @"c:\phd\bioinf\betastrands_dataset_sequences_limited.txt";
                rename_filename = @"c:\phd\bioinf\rename_num_to_pdb_limited.bat";

            }


            seqs = seqs.OrderBy(a => a).Distinct().ToList();

            //            var rename_list = seqs.SelectMany((a, i) => reverse_lookup.Where(b=>b.seq==a).Select(b=> /*program.string_debug*/($@"copy {i}.* {b.pdb}.*")).ToList();

            var rename_list = reverse_lookup.Select(a => (a.pdb, a.seq, index: seqs.IndexOf(a.seq), cmd: new[] {/*program.string_debug*/($@"del {a.pdb}.*"), /*program.string_debug*/($@"copy {seqs.IndexOf(a.seq)}.* {a.pdb}.*")})).ToList();

            //Directory.CreateDirectory(Path.GetDirectoryName(output_filename));
            //Directory.CreateDirectory(Path.GetDirectoryName(rename_filename));

            io_proxy.WriteAllLines(output_filename, seqs, nameof(program), nameof(get_pdb_sequences));
            io_proxy.WriteAllLines(rename_filename, rename_list.SelectMany(a => a.cmd).ToList(), nameof(program), nameof(get_pdb_sequences));

            //for (var i = 0; i < seqs.Count; i++)
            //{
            //    var file_seq = new List<string>();

            //    file_seq.Add(/*program.string_debug*/($@">") + i);
            //    file_seq.Add(seqs[i]);

            //    program.WriteAllLines(Path.Combine(program.data_root_folder,"iupred2a\" + i + /*program.string_debug*/($@".seq", file_seq);
            //}

            return seqs;

        }



        internal static void get_pssms(string blast_db = @"swissprot", string db_folder = @"c:\phd\blast\db\", bool remote = false, bool run = false)
        {
            int num_iterations = 3;

            // this method:
            // 0. loads 'txt_sequences' text file of dataset sequences (but only uses those sequences also found within 'txt_sequences_limited')
            // 1. if 'output_fasta_files' is true, outputs the fasta sequences in individual files
            // 2. makes a list of command lines to run blast, skipping any pssm files which already exist (& can also run blast if 'run' set to true)

            var txt_sequences = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder,/*program.string_debug*/($@"betastrands_dataset_sequences.txt")), nameof(misc), nameof(get_pssms));
            var txt_sequences_limited = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"betastrands_dataset_sequences_limited.txt")), nameof(misc), nameof(get_pssms));

            //var excluded = txt_sequences.Except(txt_sequences_limited).ToList();
            //var included = txt_sequences.Except(excluded).ToList();

            //var fasta_file = Path.Combine(program.data_root_folder,"fasta\{index1}.fasta";
            //var pssm_file = Path.Combine(program.data_root_folder,"blast_pssm_{blast_db}_{num_iterations}_{(remote ? /*program.string_debug*/($@"remote" : /*program.string_debug*/($@"local")}\{index1}.pssm";


            var fasta_blast_input_folder = /*program.string_debug*/($@"c:\phd\fasta\");
            var fasta_output_folder = /*program.string_debug*/($@"c:\phd\fasta\");

            //var pssm_blast_output_folder = /*program.string_debug*/($@"c:\phd\data\blast_pssm_{blast_db}_{(!remote ? num_iterations.ToString(CultureInfo.InvariantCulture) + /*program.string_debug*/($@"_") : /*program.string_debug*/($@""))}{(remote ? /*program.string_debug*/($@"new_remote" : /*program.string_debug*/($@"local")}";
            //var already_calculated_pssms = program.ReadAllLines(@"C:\phd\{}");





            var output_fasta_files = false;
            if (output_fasta_files)
            {
                for (var index = 0; index < txt_sequences.Length; index++)
                {
                    var fasta_file = Path.Combine(fasta_output_folder, /*program.string_debug*/($@"{index}.fasta"));

                    if (File.Exists(fasta_file) && new FileInfo(fasta_file).Length > 0) continue;

                    io_proxy.WriteLine(/*program.string_debug*/($@"Saving fasta {index}"));
                    var t_seq = txt_sequences[index];
                    var seq = new sequence(/*program.string_debug*/($@">seq_{index}"), t_seq);
                    sequence.Save(fasta_file, seq);
                }
            }

            //var tasks = new List<Task<List<string>>>();
            //var max_tasks = 1; // Environment.ProcessorCount
            var skipped_not_limited = 0;
            var skipped_pssm_exists = 0;

            var dna_pred = new List<string>();

            var emax = 100.0;
            var emin = 0.00001;
            var ediv = 10;

            for (var evalue = emax; evalue >= emin; evalue /= ediv)
            {
                for (var inclusion_ethresh = emax; inclusion_ethresh >= emin; inclusion_ethresh /= ediv)
                {
                    var cmd_lines = new List<string>() { };

                    var pssm_blast_input_folder = /*program.string_debug*/($@"c:\phd\pssm\blast_pssm_{(!remote ? num_iterations.ToString(CultureInfo.InvariantCulture) + /*program.string_debug*/($@"_") : /*program.string_debug*/($@""))}{(remote ? /*program.string_debug*/($@"remote") : /*program.string_debug*/($@"local"))}_{blast_db}_{evalue.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}");
                    var pssm_blast_output_folder = /*program.string_debug*/($@"c:\phd\pssm\blast_pssm_{(!remote ? num_iterations.ToString(CultureInfo.InvariantCulture) + /*program.string_debug*/($@"_") : /*program.string_debug*/($@""))}{(remote ? /*program.string_debug*/($@"remote") : /*program.string_debug*/($@"local"))}_{blast_db}_{evalue.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}");

                    Directory.CreateDirectory(pssm_blast_input_folder);
                    Directory.CreateDirectory(pssm_blast_output_folder);

                    for (var index = 0; index < txt_sequences.Length; index++)
                    {
                        var index1 = index;

                        if (!txt_sequences_limited.Contains(txt_sequences[index1]))
                        {
                            skipped_not_limited++;
                            continue;
                        }
                        //io_proxy.WriteLine(/*program.string_debug*/($@"Querying psi blast with sequence {index1}");

                        var dna_str = /*program.string_debug*/($@"python3 predict.py {index1} {txt_sequences[index1]}");
                        if (!dna_pred.Contains(dna_str)) dna_pred.Add(dna_str);

                        var fasta_file = Path.Combine(fasta_blast_input_folder, /*program.string_debug*/($@"{index1}.fasta"));

                        var pssm_local_file = Path.Combine(pssm_blast_input_folder, /*program.string_debug*/($@"{index1}.pssm"));

                        var pssm_copy_file = Path.Combine(pssm_blast_output_folder, /*program.string_debug*/($@"{index1}.pssm"));

                        if (File.Exists(pssm_local_file) && new FileInfo(pssm_local_file).Length > 0 && io_proxy.ReadAllLines(pssm_local_file, nameof(misc), nameof(get_pssms)).Where(a => !String.IsNullOrWhiteSpace(a)).ToList().Any(x => x.StartsWith(/*program.string_debug*/($@"PSI Gapped"), StringComparison.Ordinal)))
                        {
                            skipped_pssm_exists++;
                            continue;
                        }



                        var cmd_line = info_blast_pssm.run_psi_blast_get_pssm(fasta_file, pssm_local_file, blast_db, db_folder, num_iterations, evalue.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture), inclusion_ethresh.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture), remote, run);

                        cmd_lines.Add(cmd_line);

                        if (!String.Equals(pssm_local_file, pssm_copy_file, StringComparison.Ordinal))
                        {
                            cmd_lines.Add(/*program.string_debug*/($@"copy /y ""{pssm_local_file}"" ""{pssm_copy_file}"""));
                        }



                    }

                    var batch_pssm_local_file = Path.Combine(pssm_blast_input_folder, /*program.string_debug*/($@"blast_pssm_{(!remote ? num_iterations.ToString(CultureInfo.InvariantCulture) + /*program.string_debug*/($@"_") : /*program.string_debug*/($@""))}{(remote ? /*program.string_debug*/($@"remote") : /*program.string_debug*/($@"local"))}_{blast_db}_{evalue.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString(/*program.string_debug*/($@"000.00000"), CultureInfo.InvariantCulture)}.bat"));

                    //Directory.CreateDirectory(Path.GetDirectoryName(batch_pssm_local_file));
                    io_proxy.WriteAllLines(batch_pssm_local_file, cmd_lines, nameof(program), nameof(get_pssms));

                }
            }

            //Task.WaitAll(tasks.ToArray());

            //var all_cmd_lines = tasks.SelectMany(a => a.Result).ToList();

            var dna_pred_batch_filename = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"dna_pred.bat"));
            //Directory.CreateDirectory(Path.GetDirectoryName(dna_pred_batch_filename));
            io_proxy.WriteAllLines(dna_pred_batch_filename, dna_pred, nameof(program), nameof(get_pssms));

            io_proxy.WriteLine(/*program.string_debug*/($@"{blast_db}_{num_iterations}_{(remote ? /*program.string_debug*/($@"remote") : /*program.string_debug*/($@"local"))}: {nameof(skipped_not_limited)}:{skipped_not_limited} {nameof(skipped_pssm_exists)}:{skipped_pssm_exists}"));
        }

        internal static void load_pssms()
        {
            var txt_sequences = io_proxy.ReadAllLines(@"c:\phd\bioinf\betastrands_dataset_sequences.txt", nameof(misc), nameof(load_pssms));



            for (var index = 0; index < txt_sequences.Length; index++)
            {
                var index1 = index;

                io_proxy.WriteLine(/*program.string_debug*/($@"Loading pssm {index1}"));


                var pssm_matrix = info_blast_pssm.load_psi_blast_pssm(/*program.string_debug*/($@"c:\phd\bs_pssm\{index1}.pssm"));

                pssm_matrix = info_blast_pssm.normalise_pssm(pssm_matrix);

                io_proxy.WriteLine(/*program.string_debug*/($@"Loaded pssm {index1}"));

            }


        }

        internal static (string header, string sequence) get_uniprot_sequence(string uniprot_id)
        {
            var uniprot_url = /*program.string_debug*/($@"https://www.uniprot.org/uniprot/{uniprot_id}.fasta");


            using var wc = new WebClient();
            var fasta = wc.DownloadString(uniprot_url)
                .Split(new char[] {'\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);

            var header = String.Join(/*program.string_debug*/($@""), fasta.Where(b => b.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).ToList());
            var sequence = String.Join(/*program.string_debug*/($@""), fasta.Where(b => !b.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).ToList());

            //io_proxy.WriteLine(sequence);
            return (header, sequence);
        }

        internal static void get_uniprot_sequences(/*bool limited_to_dimorphics_and_standard = true*/)
        {
            //var dimorphics_data = program.ReadAllLines(Path.Combine(program.data_root_folder,"csv\distinct dimorphics list.csv")
            //    .Skip(1).Where(a => !string.IsNullOrWhiteSpace(a.Replace(/*program.string_debug*/($@","), /*program.string_debug*/($@"")))).Select((a, i) =>
            //    {
            //        var x = a.Split(',');
            //        return (
            //            pdb_id: x[0].ToUpperInvariant(),
            //            dimer_type: x[1],
            //            class_name: x[2],
            //            symmetry_mode: x[3],
            //            parallelism: x[4],
            //            chain_number: int.Parse(x[5]) - 1,
            //            strand_seq: x[6],
            //            optional_res_index: x[7]
            //        );
            //    }).ToList();

            //if (limited_to_dimorphics_and_standard)
            //{
            //    dimorphics_data = dimorphics_data.Where(a => a.class_name == /*program.string_debug*/($@"Single" || a.class_name == /*program.string_debug*/($@"Multiple").ToList();
            //}


            //var pdb_id_list = dimorphics_data.Select(a => a.pdb_id.Trim().ToUpperInvariant()).Distinct().ToList();


            var sifts_lookup = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"sifts"), /*program.string_debug*/($@"pdb_chain_uniprot.csv")), nameof(misc), nameof(get_uniprot_sequences)).Skip(2).Select(a =>
            {

                var x = a.Split(',');
                var i = 0;
                var y = (
                    PDB: x[i++].Substring(0, 4).ToUpperInvariant(),
                    CHAIN: x[i++][0],
                    SP_PRIMARY: x[i++],
                    RES_BEG: x[i++],
                    RES_END: x[i++],
                    PDB_BEG: x[i++],
                    PDB_END: x[i++],
                    SP_BEG: x[i++],
                    SP_END: x[i++]);

                return y;
            }).ToList();

            var pdb_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb_split"));
            var uniprot_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"uniprot"));


            var files = Directory.GetFiles(pdb_folder, /*program.string_debug*/($@"*.pdb"));

            foreach (var a in files)
            {
                io_proxy.WriteLine(a);

                var file = Path.GetFileNameWithoutExtension(a).ToUpperInvariant();
                var pdb = file.Substring(0, 4);
                var chain = file[4];
                var uniprot_id = sifts_lookup.FirstOrDefault(b => b.PDB.Equals(pdb, StringComparison.OrdinalIgnoreCase) && b.CHAIN == chain).SP_PRIMARY;

                if (String.IsNullOrWhiteSpace(uniprot_id)) continue;

                var uniprot_fasta = get_uniprot_sequence(uniprot_id);

                //var data = (pdb: pdb, chain: chain, uniprot_id: uniprot_id, uniprot_header:uniprot_fasta.header, uniprot_sequence:uniprot_fasta.sequence);

                var fn = Path.Combine(uniprot_folder, /*program.string_debug*/($@"{pdb}{chain}.fasta"));
                //Directory.CreateDirectory(Path.GetDirectoryName(fn));
                io_proxy.WriteAllLines(fn, new string[] { uniprot_fasta.header, uniprot_fasta.sequence }, nameof(program), nameof(get_uniprot_sequences));
            }






        }

        internal static void check_buildmodel_position_scan_files_are_complete()
        {
            var cmd_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"foldx_calc_buildmodel_position_scan.bat"));
            var cmd_list_file_missing = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"missing_foldx_calc_buildmodel_position_scan.bat"));
            var cmd_list_file_missing_del = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"del_missing_foldx_calc_buildmodel_position_scan.bat"));

            var cmd_list = io_proxy.ReadAllLines(cmd_list_file, nameof(misc), nameof(check_buildmodel_position_scan_files_are_complete)).Where(a => a.Contains(/*program.string_debug*/($@"--mutant-file="), StringComparison.Ordinal)).ToList();


            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();


                var mutant_file = s.First(b => b.StartsWith(/*program.string_debug*/($@"--mutant-file="), StringComparison.Ordinal)).Split('=').Last();
                if (!File.Exists(mutant_file) || new FileInfo(mutant_file).Length == 0)
                {
                    io_proxy.WriteLine(/*program.string_debug*/($@"missing mutant file: {mutant_file}"));
                }

                var dif_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal);
                var average_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal).Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Average"), StringComparison.Ordinal);
                var raw_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal).Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Raw"), StringComparison.Ordinal);



                var mutant_data = io_proxy.ReadAllLines(mutant_file, nameof(misc), nameof(check_buildmodel_position_scan_files_are_complete)).Where(b => !String.IsNullOrWhiteSpace(b)).Select(b => b.Trim()).ToList();


                var dif_fxout_data = (!File.Exists(dif_fxout_file) || new FileInfo(dif_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(dif_fxout_file, nameof(misc), nameof(check_buildmodel_position_scan_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();
                var average_fxout_data = (!File.Exists(average_fxout_file) || new FileInfo(average_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(average_fxout_file, nameof(misc), nameof(check_buildmodel_position_scan_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();
                var raw_fxout_data = (!File.Exists(raw_fxout_file) || new FileInfo(raw_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(raw_fxout_file, nameof(misc), nameof(check_buildmodel_position_scan_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();


                return (dif_fxout_file, dif_fxout_data, average_fxout_file, average_fxout_data, raw_fxout_file, raw_fxout_data, mutant_file, mutant_data, a);

            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.mutant_data.Count != a.dif_fxout_data.Count).ToList();
            var not_missing = cmd_list_parsed.Except(missing).ToList();

            foreach (var m in missing)
            {
                var total_missing = m.mutant_data.Count - m.dif_fxout_data.Count;
                io_proxy.WriteLine(/*program.string_debug*/($@"missing {total_missing} / {m.mutant_data.Count}: {m.mutant_file}"));
            }



            io_proxy.WriteLine(/*program.string_debug*/($@"total missing: {missing.Count}"));
            io_proxy.WriteLine(/*program.string_debug*/($@"total not missing: {not_missing.Count}"));






            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith(/*program.string_debug*/($@"if not exist"), StringComparison.OrdinalIgnoreCase)).Select(b => /*program.string_debug*/($@"del {b.Split()[3]}")).ToList();
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Average"), StringComparison.Ordinal)).ToList());
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Raw"), StringComparison.Ordinal)).ToList());

            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing));
            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing_del));

            io_proxy.WriteAllLines(cmd_list_file_missing, missing_list, nameof(program), nameof(check_buildmodel_position_scan_files_are_complete));
            io_proxy.WriteAllLines(cmd_list_file_missing_del, missing_list_del, nameof(program), nameof(check_buildmodel_position_scan_files_are_complete));


            io_proxy.WriteLine();

        }


        internal static void check_buildmodel_subseq_mutant_files_are_complete()
        {
            var cmd_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"foldx_calc_buildmodel_subsequence_mutant.bat"));
            var cmd_list_file_missing = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"missing_foldx_calc_buildmodel_subsequence_mutant.bat"));
            var cmd_list_file_missing_del = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"del_missing_foldx_calc_buildmodel_subsequence_mutant.bat"));

            var cmd_list = io_proxy.ReadAllLines(cmd_list_file, nameof(misc), nameof(check_buildmodel_subseq_mutant_files_are_complete)).Where(a => a.Contains(/*program.string_debug*/($@"--mutant-file="), StringComparison.Ordinal)).ToList();



            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();


                var mutant_file = s.First(b => b.StartsWith(/*program.string_debug*/($@"--mutant-file="), StringComparison.Ordinal)).Split('=').Last();
                if (!File.Exists(mutant_file) || new FileInfo(mutant_file).Length == 0)
                {
                    io_proxy.WriteLine(/*program.string_debug*/($@"missing mutant file: {mutant_file}"));
                }

                var dif_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal);
                var average_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal).Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Average"), StringComparison.Ordinal);
                var raw_fxout_file = s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal).Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Raw"), StringComparison.Ordinal);



                var mutant_data = io_proxy.ReadAllLines(mutant_file, nameof(misc), nameof(check_buildmodel_subseq_mutant_files_are_complete)).Where(b => !String.IsNullOrWhiteSpace(b)).Select(b => b.Trim()).ToList();


                var dif_fxout_data = (!File.Exists(dif_fxout_file) || new FileInfo(dif_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(dif_fxout_file, nameof(misc), nameof(check_buildmodel_subseq_mutant_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();
                var average_fxout_data = (!File.Exists(average_fxout_file) || new FileInfo(average_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(average_fxout_file, nameof(misc), nameof(check_buildmodel_subseq_mutant_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();
                var raw_fxout_data = (!File.Exists(raw_fxout_file) || new FileInfo(raw_fxout_file).Length == 0) ? new List<string>() : io_proxy.ReadAllLines(raw_fxout_file, nameof(misc), nameof(check_buildmodel_subseq_mutant_files_are_complete)).Skip(9).Where(b => !String.IsNullOrWhiteSpace(b)).ToList();


                return (dif_fxout_file, dif_fxout_data, average_fxout_file, average_fxout_data, raw_fxout_file, raw_fxout_data, mutant_file, mutant_data, a);

            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.mutant_data.Count != a.dif_fxout_data.Count).ToList();
            var not_missing = cmd_list_parsed.Except(missing).ToList();

            foreach (var m in missing)
            {
                var total_missing = m.mutant_data.Count - m.dif_fxout_data.Count;
                io_proxy.WriteLine(/*program.string_debug*/($@"missing {total_missing} / {m.mutant_data.Count}: {m.mutant_file}"));
            }



            io_proxy.WriteLine(/*program.string_debug*/($@"total missing: {missing.Count}"));
            io_proxy.WriteLine(/*program.string_debug*/($@"total not missing: {not_missing.Count}"));


            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith(/*program.string_debug*/($@"if not exist"), StringComparison.OrdinalIgnoreCase)).Select(b => /*program.string_debug*/($@"del {b.Split()[3]}")).ToList();
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Average"), StringComparison.Ordinal)).ToList());
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace(/*program.string_debug*/($@"Dif"), /*program.string_debug*/($@"Raw"), StringComparison.Ordinal)).ToList());

            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing));
            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing_del));

            io_proxy.WriteAllLines(cmd_list_file_missing, missing_list, nameof(program), nameof(check_buildmodel_subseq_mutant_files_are_complete));
            io_proxy.WriteAllLines(cmd_list_file_missing_del, missing_list_del, nameof(program), nameof(check_buildmodel_subseq_mutant_files_are_complete));


            io_proxy.WriteLine();

        }

        internal static void check_ala_scanning_files_are_complete()
        {
            var ala_scan_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx_ala_scan"));

            var cmd_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"foldx_calc_ala_scanning.bat.skip"));
            var cmd_list_file_missing = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"missing_foldx_calc_ala_scanning.bat.skip"));
            var cmd_list_file_missing_del = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"del_missing_foldx_calc_ala_scanning.bat.skip"));

            var cmd_list = io_proxy.ReadAllLines(cmd_list_file, nameof(misc), nameof(check_ala_scanning_files_are_complete)).Where(a => a.Contains(/*program.string_debug*/($@"--pdb="), StringComparison.Ordinal)).ToList();

            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();

                var id = Path.GetFileNameWithoutExtension(s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal).Replace(/*program.string_debug*/($@"--pdb="), /*program.string_debug*/($@""), StringComparison.Ordinal));

                var pdb_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb_split_repair"), /*program.string_debug*/($@"{id}.pdb"));

                var ala_file = Path.Combine(ala_scan_folder, /*program.string_debug*/($@"{id}_AS.fxout"));

                var pdb_res = io_proxy
                    .ReadAllLines(pdb_file, nameof(misc), nameof(check_ala_scanning_files_are_complete))
                    .Where(b => b.StartsWith(/*program.string_debug*/($@"ATOM"), StringComparison.Ordinal))
                    .Select(c => /*program.string_debug*/($@"{c.Substring(17, 3)} {Int32.Parse(c.Substring(22, 4).Trim(), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)}"))
                    .Distinct()
                    .ToList();

                var ala_res = io_proxy
                    .ReadAllLines(ala_file, nameof(misc), nameof(check_ala_scanning_files_are_complete))
                    .Where(c => !String.IsNullOrWhiteSpace(c))
                    .Select(c => c.Substring(0, c.IndexOf(/*program.string_debug*/($@" to "), StringComparison.Ordinal)).Trim())
                    .ToList();

                ala_res = ala_res.Select(b => b.Replace(/*program.string_debug*/($@"H1S"), /*program.string_debug*/($@"HIS"), StringComparison.Ordinal)).ToList();
                ala_res = ala_res.Select(b => b.Replace(/*program.string_debug*/($@"H2S"), /*program.string_debug*/($@"HIS"), StringComparison.Ordinal)).ToList();

                var present = pdb_res.Intersect(ala_res).ToList();
                var missing1 = pdb_res.Except(ala_res).ToList();
                var extra = ala_res.Except(pdb_res).ToList();

                return (id: id, pdb_res: pdb_res, ala_res: ala_res, present: present, missing1: missing1, extra: extra, a);



                //positions: s[5].Split('=').Last().Split(','));
            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.missing1.Count > 0 || a.extra.Count > 0).ToList();

            io_proxy.WriteLine(/*program.string_debug*/($@"missed: {missing.Count}"));
            foreach (var m in missing)
            {
                io_proxy.WriteLine(/*program.string_debug*/($@"{m.id}, extra: {String.Join(/*program.string_debug*/($@"|"), m.extra)}, missing: {String.Join(/*program.string_debug*/($@"|"), m.missing1)}"));
            }


            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith(/*program.string_debug*/($@"if not exist"), StringComparison.OrdinalIgnoreCase)).Select(b => /*program.string_debug*/($@"del {b.Split()[3]}")).ToList();

            ///Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing));
            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing_del));

            io_proxy.WriteAllLines(cmd_list_file_missing, missing_list, nameof(program), nameof(check_ala_scanning_files_are_complete));
            io_proxy.WriteAllLines(cmd_list_file_missing_del, missing_list_del, nameof(program), nameof(check_ala_scanning_files_are_complete));

            io_proxy.WriteLine();
        }

        internal static void check_foldx_position_scan_files_are_completed()
        {

            var missing_res = compare_repaired_pdb_to_pdb();

            var ps_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"ps"));
            var file_mask = /*program.string_debug*/($@"PS_*.txt");

            var ps_files = Directory.GetFiles(ps_folder, file_mask, SearchOption.TopDirectoryOnly);

            var cmd_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"foldx_calc_position_scanning.bat"));
            var cmd_list_file_missing = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"missing_foldx_calc_position_scanning.bat"));
            var cmd_list_file_missing_del = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"bat"), /*program.string_debug*/($@"del_missing_foldx_calc_position_scanning.bat"));

            var cmd_list = io_proxy.ReadAllLines(cmd_list_file, nameof(misc), nameof(check_foldx_position_scan_files_are_completed)).Where(a => a.Contains(/*program.string_debug*/($@"--positions="), StringComparison.Ordinal)).ToList();

            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();

                return (file: s[3].Replace(/*program.string_debug*/($@""""), /*program.string_debug*/($@""), StringComparison.Ordinal), positions: s[5].Split('=').Last().Split(','), a);
            }).ToList();

            var incomplete_files = new List<string>();
            var complete_files = new List<string>();

            var missing = new List<string>();

            foreach (var c in cmd_list_parsed)
            {
                var file = c.file;

                var pdb = Path.GetFileNameWithoutExtension(file).Substring(3, 5).Trim();

                var pdb_missing_res = missing_res.Where(a => String.Equals(a.id, pdb, StringComparison.OrdinalIgnoreCase)).ToList();

                var actual_miss = pdb_missing_res.SelectMany(a => a.original_except_repair).ToList();

                var x = c.positions.Select(a => a[0] + a[1..]).ToList();

                // convert 'd' to the 25 codes

                x = x.SelectMany(a =>
                {
                    var rid = Int32.Parse(string.Join(/*program.string_debug*/($@""), a.Where(b => Char.IsDigit(b)).ToList()), NumberStyles.Integer, NumberFormatInfo.InvariantInfo);

                    if (actual_miss.Contains(rid))
                    {

                        return new List<string>() { };
                    }

                    if (a.Last() == 'd')
                    {                //123456789012345678901234567890
                        var d_codes = /*program.string_debug*/($@"GALVIPRTSCMKEQDNWYFypsh"); //e

                        var output = d_codes.Select(d => a[0..^1] + d).ToList();
                        return output;
                    }
                    else
                    {
                        return new List<string>() { a };

                    }
                }).Where(a => !String.IsNullOrWhiteSpace(a)).Distinct().ToList();

                var x_res = x.Select(a => Int32.Parse(string.Join(/*program.string_debug*/($@""), a.Where(b => Char.IsDigit(b)).ToList()), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).Distinct().ToList();

                var scanning_data = File.Exists(c.file) && new FileInfo(c.file).Length > 0 ? io_proxy.ReadAllLines(c.file, nameof(misc), nameof(check_foldx_position_scan_files_are_completed)).Select(a =>
                {
                    var z = a.Split();
                    z[0] = info_foldx.foldx_residues_aa_mutable.First(d => d.foldx_aa_code3 == z[0].Substring(0, 3)).standard_aa_code1 + z[0][3..];
                    return z;
                }).Distinct().ToList() : new List<string[]>();


                var data_codes = scanning_data.Select(a => a.First()).Distinct().ToList();
                var data_codes_res = data_codes.Select(a => Int32.Parse(string.Join(/*program.string_debug*/($@""), a.Where(b => Char.IsDigit(b)).ToList()), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).Distinct().ToList();

                var missing1 = x.Except(data_codes).Distinct().ToList();

                if (missing1.Count > 0)
                {
                    missing.Add(c.a);

                    incomplete_files.Add(file);

                    io_proxy.WriteLine();
                    io_proxy.WriteLine(/*program.string_debug*/($@"File: {file}"));
                    io_proxy.WriteLine();
                    io_proxy.WriteLine(/*program.string_debug*/($@"Missing {missing1.Count}: $@{string.Join(/*program.string_debug*/($@", "), missing1)}"));
                    io_proxy.WriteLine();
                }
                else
                {
                    complete_files.Add(file);

                }

            }

            io_proxy.WriteLine(/*program.string_debug*/($@"incomplete_file_count: {incomplete_files.Count}"));
            io_proxy.WriteLine(/*program.string_debug*/($@"complete_file_count: {complete_files.Count}"));

            //program.WriteAllLines(Path.Combine(program.data_root_folder,"foldx\ps_incomplete_files.txt", incomplete_files);



            var missing_list = missing.Select(a => a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith(/*program.string_debug*/($@"if not exist"), StringComparison.OrdinalIgnoreCase)).Select(b => /*program.string_debug*/($@"del {b.Split()[3]}")).ToList();

            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing));
            //Directory.CreateDirectory(Path.GetDirectoryName(cmd_list_file_missing_del));

            io_proxy.WriteAllLines(cmd_list_file_missing, missing_list, nameof(program), nameof(check_foldx_position_scan_files_are_completed));
            io_proxy.WriteAllLines(cmd_list_file_missing_del, missing_list_del, nameof(program), nameof(check_foldx_position_scan_files_are_completed));

            io_proxy.WriteLine();
        }

        internal static List<(string id, string id_repair, List<int> res_id_original, List<int> res_id_repair, List<int> intersect, List<int> original_except_repair, List<int> repair_except_original)> compare_repaired_pdb_to_pdb()
        {
            var split_pdbs = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb_split"));
            var split_pdbs_repair = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb_split_repair"));

            var pdbs_original = Directory.GetFiles(split_pdbs, /*program.string_debug*/($@"*.pdb"));
            var pdbs_repair = Directory.GetFiles(split_pdbs_repair, /*program.string_debug*/($@"*.pdb"));

            var ids = pdbs_original.Select(a => Path.GetFileNameWithoutExtension(a) ?? /*program.string_debug*/($@"")).ToList();

            var pdbs = ids.Select((a, i) => (id: Path.GetFileNameWithoutExtension(pdbs_original[i]), id_repair: Path.GetFileNameWithoutExtension(pdbs_repair[i]), pdb_original: pdbs_original[i], pdb_repair: pdbs_repair[i])).ToList();

            List<(string id, string id_repair, List<int> res_id_original, List<int> res_id_repair, List<int> intersect, List<int> original_except_repair, List<int> repair_except_original)> diff = pdbs.Select(a =>
            {
                var o = io_proxy.ReadAllLines(a.pdb_original, nameof(misc), nameof(compare_repaired_pdb_to_pdb)).Where(b => b.StartsWith(/*program.string_debug*/($@"ATOM"), StringComparison.Ordinal)).Select(b => Int32.Parse(b.Substring(22, 4).Trim(), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).Distinct().ToList();
                var r = io_proxy.ReadAllLines(a.pdb_repair, nameof(misc), nameof(compare_repaired_pdb_to_pdb)).Where(b => b.StartsWith(/*program.string_debug*/($@"ATOM"), StringComparison.Ordinal)).Select(b => Int32.Parse(b.Substring(22, 4).Trim(), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).Distinct().ToList();
                var i = o.Intersect(r).ToList();
                var e1 = o.Except(r).ToList();
                var e2 = r.Except(o).ToList();

                return (id: a.id, id_repair: a.id_repair, res_id_original: o, res_id_repair: r, intersect: i, original_except_repair: e1,
                    repair_except_original: e2);

            }).ToList();

            diff = diff.Where(a => a.intersect.Count != a.res_id_original.Count || a.intersect.Count != a.res_id_repair.Count).ToList();

            return diff;
        }

        internal static void get_limited_ids()
        {

            var seqs_all = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder,/*program.string_debug*/($@"betastrands_dataset_sequences.txt")), nameof(misc), nameof(get_limited_ids));
            var seqs_limited = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder,/*program.string_debug*/($@"betastrands_dataset_sequences_limited.txt")), nameof(misc), nameof(get_limited_ids));

            var ids = seqs_all.Select((a, i) => seqs_limited.Contains(a) ? i : -1).Distinct().Where(a => a != -1).ToList();

            var fn = Path.Combine(program.data_root_folder,/*program.string_debug*/($@"copy_limited.bat"));

            //Directory.CreateDirectory(Path.GetDirectoryName(fn));

            var ids_lines = ids.Select(a => /*program.string_debug*/($@"copy {a}.* limited\")).ToList();

            io_proxy.WriteAllLines(fn, ids_lines, nameof(program), nameof(get_limited_ids));
            //Console.ReadKey();

        }



        //internal static void get_missing_data()
        //{
        //    var lines = program.ReadAllLines(@"c:\phd\bioinf\dimorphics_data_set.csv").Skip(1).Select(a => a.Split(',').ToList()).ToList();

        //    var pdb_id_col = 1;
        //    var seq_col = 5;

        //    var list = lines.Select(a => (a[pdb_id_col], a[seq_col], a)).ToList();

        //    //var results = new List<string>();

        //    var tasks = new List<Task<List<string>>>();

        //    var max_tasks = 1000;

        //    foreach (var seq in list)
        //    {
        //        var z = seq;

        //        var task = Task.Run(() =>
        //        {
        //            var x = find_subsequences(z.Item1, z.Item2);

        //            for (var i = 0; i < x.Count; i++)
        //            {
        //               // x[i] = string.Join(/*program.string_debug*/($@","), z.a) + /*program.string_debug*/($@",") + x[i];
        //            }
        //            //results.AddRange(x);

        //            return x;

        //        });

        //        tasks.Add(task);

        //        while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
        //        {
        //            Task.WaitAny(tasks.ToArray<Task>());
        //        }
        //    }

        //    Task.WaitAll(tasks.ToArray<Task>());


        //    var results = tasks.SelectMany(a => a.Result).ToList();

        //    program.WriteAllLines(@"c:\phd\bioinf\res_ids_list.csv", results);

        //    io_proxy.WriteLine(/*program.string_debug*/($@"Finished.");
        //    Console.ReadLine();
        //}

        // xtodo: encode the foldx scores
        //internal static object _console_lock = new object();




        //subsequence_classification_data.encode_sequence(/*program.string_debug*/($@"x", /*program.string_debug*/($@"ALGY");

        //xtodo: create dhc dataset with dhc found by dssp3, dssp7, stride3, stride7, for comparison
        //xtodo: to see which secondary structure definition gives the best performance in the classifier


        //xtodo: ensure that data is placed in the correct column, in case of any features being empty/null/not included in some rows.

        //var pdb_files = Directory.GetFiles(Path.Combine(program.data_root_folder,"pdb\", /*program.string_debug*/($@"*.pdb");
        //pdb_files.ToList().ForEach(a=>Atom.run_dssp(Path.GetFileNameWithoutExtension(a)));
        //pdb_files.ToList().ForEach(a=>Atom.split_pdb_chains(Path.GetFileNameWithoutExtension(a)));
        //pdb_files.ToList().ForEach(a=>Program.WriteLine(/*program.string_debug*/($@"dssp -i {Path.GetFileNameWithoutExtension(a)}.pdb -o {Path.GetFileNameWithoutExtension(a)}.dssp"));
        //return;

        //if (run_larks_unknown)
        //{
        //    class_list.Add(
        //        (
        //        nameof(run_larks_unknown).EndsWith(/*program.string_debug*/($@"unknown") ? +1 : -1,
        //        nameof(run_larks_unknown).Replace(/*program.string_debug*/($@"run_", /*program.string_debug*/($@"")),
        //        /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(run_larks_unknown).Replace(/*program.string_debug*/($@"run_", /*program.string_debug*/($@""))}{(sequence_only_features ? /*program.string_debug*/($@"_seq_only" : /*program.string_debug*/($@""))}_{uid}.txt",
        //        /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(run_larks_unknown).Replace(/*program.string_debug*/($@"run_", /*program.string_debug*/($@""))}{(sequence_only_features ? /*program.string_debug*/($@"_seq_only" : /*program.string_debug*/($@""))}_{uid}.txt")

        //    );
        //}

        //else if (class_name == nameof(run_larks_unknown).Replace(/*program.string_debug*/($@"run_", /*program.string_debug*/($@"")))
        //{
        //    var larks = new List<string>() { /*program.string_debug*/($@"STGGYG", /*program.string_debug*/($@"GYNGFG", /*program.string_debug*/($@"SYSGYS", /*program.string_debug*/($@"SYSSYGQS", /*program.string_debug*/($@"GFGNFGTS" };
        //    sequence_list = larks;
        //}

        internal static void test2()
        {

            /*for (var i = 10; i >= 0; i--)
            {
                var aa_seq_pse_aac_options = new pse_aac_options()
                {
                    oaac = true,
                    oaac_binary = true,
                    motifs = true,
                    motifs_binary = true,
                    dipeptides = true,
                    dipeptides_binary = true,
                    //saac = true,
                    //saac_binary = true,
                    average_seq_position = true,
                    average_dipeptide_distance = true,
                };

                var ps = subsequence_classification_data.enum_protein_data_source.neighbourhood_1d;
                var category_prefix = new string('S', i);
                var group_prefix = new string('S', i);
                var aa_subsequence = new string('S', i);
                //var pse_ssc = subsequence_classification_data.calculate_aa_or_ss_sequence_classification_data(ps,category_prefix, group_prefix, ss_seq, feature_calcs.seq_type.secondary_structure_sequence, mpsa_pse_aac_options);
                //var pse_ssc_dssp_classification_data = subsequence_classification_data.calculate_aa_or_ss_sequence_classification_data(ps, category_prefix, /*program.string_debug* /($@"dssp_monomer", scd.dssp_monomer_subsequence, feature_calcs.seq_type.secondary_structure_sequence, aa_seq_pse_aac_options);
                var pse_aac_sequence_classification_data = subsequence_classification_data.calculate_aa_or_ss_sequence_classification_data(ps, category_prefix, group_prefix, aa_subsequence, feature_calcs.seq_type.secondary_structure_sequence, aa_seq_pse_aac_options, 100);

                io_proxy.WriteLine(pse_aac_sequence_classification_data.Count);
            }

            return;*/

            //for (var b = 0; b <= 1; b++)
            //{
            //    for (var i = 0; i < 100; i++)
            //    {
            //        var inp = new string('A', i);
            //        var x = feature_calcs.split_sequence(inp,3,0,b==1);
            //        io_proxy.WriteLine(inp + /*program.string_debug*/($@": " + string.Join(/*program.string_debug*/($@", "), x));
            //    }
            //}
            //return;

            //var aaa= subsequence_classification_data.aaindex_subset_templates_search();


            //return;
            //r_peptides.init_r();
            //io_proxy.WriteLine(r_peptides.get_values(/*program.string_debug*/($@"AAAAAAAAAALLLLLLL"));
            //for (var i = 0; i < 100; i++) Console.ReadLine();

            //return;

            //get_limited_ids();
            //return;

            //var xxxxx = new mpsa_reader(Path.Combine(program.data_root_folder,"ss_meta\0.psipred-swissprot");

            //return;
            ////pssm test
            //var blast_db_folder = @"c:\phd\blast\db\";
            //var blast_dbs =new[]{"nr", /*program.string_debug*/($@"swissprot", /*program.string_debug*/($@"uniref90"};
            //var remotes = new[]{true, false};
            //var run = false;
            //foreach (var blast_db in blast_dbs)
            //{
            //    foreach (var remote in remotes)
            //    {
            //        get_pssms(blast_db, blast_db_folder, remote, run);
            //    }
            //}
            //io_proxy.WriteLine(/*program.string_debug*/($@"Finished");
            //Console.ReadLine();
            //return;



            //dhc_dataset_maker.run_dhc_dataset_maker(substructure_type.dimorphic, 1, 2, true, false);
            //return;


            //for (var i = 3; i <= 5; i++)
            //{
            //    var a = new string(/*program.string_debug*/($@"ALGYALGY".ToCharArray(),0, i);
            //    var pepValues = r_peptides.get_values(a);
            //    io_proxy.WriteLine(/*program.string_debug*/($@"!!! length {i}: sequence {a}: features: {pepValues.Count}");
            //}

            //Console.ReadKey();
            //return;

            //REngine engine = REngine.GetInstance();
            ////engine.Evaluate(/*program.string_debug*/($@"install.packages(""Peptides"")");

            //var r_init_cmds = /*program.string_debug*/($@"
            //#install.packages(""devtools"", repos=""http://cran.us.r-project.org"")
            //#library(devtools)
            //#install_github(""dosorio/Peptides"")
            //library(Peptides)
            //";
            //r_init_cmds.Split(new char[] { '\r', '\n'}).Where(a=>!string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith(/*program.string_debug*/($@"#")).ToList().ForEach(a=> engine.Evaluate(a));

            //var xyz_e = engine.Evaluate(/*program.string_debug*/($@"xyz <- lengthpep(seq = ""GLPRKILCAIAKKKGKCKGPLKLVCKC"")").AsDataFrame();
            //NumericVector xyz = engine.GetSymbol(/*program.string_debug*/($@"xyz").AsNumeric();



            //var edge_files = Directory.GetFiles(@"C:\phd\ring_pdb\pdb\", /*program.string_debug*/($@"*.edges");
            ////var node_files = Directory.GetFiles(@"C:\phd\ring_pdb\pdb\","*.nodes");


            //var ring_edges = edge_files.SelectMany(a=> ring.ring_edge.load(a));

            //var i_types0 = ring_edges.SelectMany(a => new[]
            //{
            //    a.Interaction.text
            //}).ToList();

            //var ix = i_types0.Distinct().OrderBy(a => a).Select(a=>(a,i_types0.Count(b=>b==a))).ToList();

            //////var i_types1 = ring_edges.Select(a => a.Interaction.subtype_node1).Distinct().ToList();
            //////var i_types2 = ring_edges.Select(a => a.Interaction.subtype_node2).Distinct().ToList();

            //io_proxy.WriteLine(string.Join($"\r\n", ix));

            //SC MC





            //Console.ReadKey();
            //return;
            //compare_repaired_pdb_to_pdb();
            //return;

            //check_ala_scanning_files_are_complete();
            //check_foldx_position_scan_files_are_completed();
            //check_buildmodel_position_scan_files_are_complete();
            //check_buildmodel_subseq_mutant_files_are_complete();
            //Console.ReadLine();
            //return;


            // sable test
            //sable.load(@"c:\phd\sable\test.txt");
            // return;

            // uniprot test
            //get_uniprot_sequences();
            //return;

            // pdb sequences test
            //get_pdb_sequences(limited_to_dimorphics_and_standard: false);
            //get_pdb_sequences(limited_to_dimorphics_and_standard: true);
            //io_proxy.WriteLine(/*program.string_debug*/($@"Finished");
            //Console.ReadLine();

            //var fghx = feature_calcs.feature_pse_aac(/*program.string_debug*/($@"CCCCCCCCCCCCCEEEEEEEEHHHHHHT", feature_calcs.seq_type.secondary_structure_sequence, new feature_calcs.pse_aac_options(), false, true);
            //io_proxy.WriteLine();
            //get_pdb_sequences();
            //return;
            //var xxx = feature_calcs.feature_pse_aac(/*program.string_debug*/($@"AAA", false, false);

            //return;


            //fill_missing_chains();
            //get_missing_data();
            //return;

            //extract_all_pdbs();
            //repair_all_extracted_pdbs();

            //get_pdb_sequences();

            //get_pssms();
            //load_pssms();

            //return;
            // unknown classes
            //var run_larks_unknown = true;
            //if (run_dimorphics_pos && run_dhc_pos) throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: Is that right?");
            //if (run_coils_neg && run_strands_neg) throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: Is that right?");
            //var uid = 2;
            //class_list.Add((-1, substructure_type.coil, 3, nameof(substructure_type.coil), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.coil)}_{uid}", /*program.string_debug*/($@"c:\bioinf\d_{nameof(substructure_type.coil)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand, 2, nameof(substructure_type.standard_beta_strand), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.standard_beta_strand)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.standard_beta_strand)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand_with_host_coil, 3, nameof(substructure_type.standard_beta_strand_with_host_coil), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.standard_beta_strand_with_host_coil)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.standard_beta_strand_with_host_coil)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand_full_protein_sequence, 3, nameof(substructure_type.standard_beta_strand_full_protein_sequence), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.standard_beta_strand_full_protein_sequence)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.standard_beta_strand_full_protein_sequence)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics, 2, nameof(substructure_type.dimorphics), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.dimorphics)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.dimorphics)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics_with_host_coil, 3, nameof(substructure_type.dimorphics_with_host_coil), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.dimorphics_with_host_coil)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.dimorphics_with_host_coil)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics_full_protein_sequence, 3, nameof(substructure_type.dimorphics_full_protein_sequence), /*program.string_debug*/($@"c:\phd\bioinf\e_{nameof(substructure_type.dimorphics_full_protein_sequence)}_{uid}", /*program.string_debug*/($@"c:\phd\bioinf\d_{nameof(substructure_type.dimorphics_full_protein_sequence)}_{uid}"));


            //var make_distance_matrices = false;// true;
            //var encode_features = true;

            //var sequence_only_features = false;
            //var calculate_1d_features = true;
            //var calculate_3d_features = false;
            //var calculate_1d_nh_features = false;
            //var calculate_3d_nh_features = false;
        }

        //if (run_strands_neg)
        //{
        //    class_info_list.Add(new class_info()
        //    {
        //        class_id = -1,
        //        substructure_type = substructure_type.standard_strand,
        //        min_subsequence_length = 2,
        //        class_name = nameof(substructure_type.standard_strand),
        //    });
        //}

        //if (run_dimorphics_pos)
        //{
        //    class_info_list.Add(new class_info()
        //    {
        //        class_id = +1,
        //        substructure_type = substructure_type.dimorphic,
        //        min_subsequence_length = 2,
        //        class_name = nameof(substructure_type.dimorphic),
        //    });
        //}

        //if (run_dhc_pos)
        //{
        //    class_info_list.Add(new class_info()
        //    {
        //        class_id = +1,
        //        substructure_type = substructure_type.dimorphic_coil,
        //        min_subsequence_length = 3,
        //        class_name = nameof(substructure_type.dimorphic_coil),
        //    });
        //}

        //if (run_coils_neg)
        //{
        //    class_info_list.Add(new class_info()
        //    {
        //        class_id = -1,
        //        substructure_type = substructure_type.standard_coil,
        //        min_subsequence_length = 3,
        //        class_name = nameof(substructure_type.standard_coil),
        //    });
        //}


        //dimorphics_data.ForEach(a=>Atom.run_dssp(a));

        //return null;

        //dimorphics_data = dimorphics_data.Take(10).ToList();
        //dimorphics_data = dimorphics_data.Where(a => a[6].Length >= min_dimorphic_length_required).ToList();
        //var dimorphics_ids = dimorphics_data.Select(a => a.First().ToUpperInvariant()).Distinct().OrderBy(a=>a).ToList();

        //0   1   2   3   4   5   6   7   8   9
        //Pdb Id  Dimer Type  Class Name  Symmetry Mode   Parallelism Chain   Sequence Res Index Pdb Seq Pdb Seq Len

        //var dhc_features = new List<List<(string id, object value)>>();


        //int cl;
        //int ct;

        //lock (program._console_lock)
        //{
        //    io_proxy.WriteLine();
        //    cl = Console.CursorLeft + 1;
        //    ct = Console.CursorTop;
        //    io_proxy.WriteLine(/*program.string_debug*/($@"[{new string('.', dimorphics_data.Count)}]");
        //    io_proxy.WriteLine();
        //}
        //var console_lock = new object();

        //Console.SetCursorPosition(0,ct+1);
    }
}
