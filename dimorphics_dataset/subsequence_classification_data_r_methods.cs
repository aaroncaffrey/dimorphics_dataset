﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;

namespace dimorphics_dataset
{
    internal static class subsequence_classification_data_r_methods
    {
        public const string module_name = nameof(subsequence_classification_data_r_methods);

        public const string cache_r_protr_filename = @"c:\phd\dataset\r_protr_cache.csv";
        public const string cache_r_peptides_filename = @"c:\phd\dataset\r_peptides_cache.csv";

        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_protr;// = cache_r_load(cache_r_protr_filename);
        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_peptides;// = cache_r_load(cache_r_peptides_filename);

        internal static void cache_r_servers()
        {
            const string method_name = nameof(cache_r_servers);
            /*
            var files = new string[] { /*program.string_debug* /($@"c:\phd\dataset\l__(standard_coil).csv", /*program.string_debug* /($@"c:\phd\dataset\l__(dimorphic_coil).csv" };

            var psi_list = files.SelectMany(a => protein_subsequence_info.load(a)).ToList();

            var pdb_opt = new load_atoms_pdb_options()
            {
                first_model_only = true,
                first_icode_only = true,
                load_uniprot = false,

                load_2d_mpsa_sec_struct_predictions = false,
                load_2d_blast_pssms = false,
                load_2d_iup_data = false,
                load_2d_sable = false,
                load_2d_dna_binding = false,

                find_3d_intramolecular = true,
                load_3d_rsa_data = false,
                load_3d_dssp_data = false,
                load_3d_stride_data = false,
                load_3d_ring_data = false,
                load_3d_foldx_ala_scan = false,
            };

            var aa_list = psi_list.AsParallel().AsOrdered().SelectMany(a =>
            {
                var b = get_subsequence_classificiation_data(a, new feature_types(), null, pdb_opt, false, false);
                
                var c = b.get_regions().Select(c => c.region.aa_sequence).Distinct().ToList();

                return c;
            }).Distinct().ToList();

            aa_list = aa_list.Distinct().OrderBy(a => a.Length).ThenBy(a => a).ToList();
            */
            var aa_list = cache_r_get_seqs();
            aa_list = aa_list.OrderBy(a => a).ToList();

            {
                var r_protr_results = ParallelEnumerable.Select<string, (string aa_seq, List<feature_info> features)>(aa_list.AsParallel().AsOrdered(), a =>
                    (aa_seq: a, features: call_r_protr(a))).ToList();
                var r_protr_csv = r_protr_results.SelectMany(a =>
                        a.features.Select(b =>
                            /*program.string_debug*/($@"{a.aa_seq},{String.Join(/*program.string_debug*/($@","), b.key_value_list(true).Select(c => c.value).ToList())}")).ToList())
                    .ToList();
                io_proxy.WriteAllLines(cache_r_protr_filename, r_protr_csv);
            }

            {
                var r_peptides_results = ParallelEnumerable.Select<string, (string aa_seq, List<feature_info> features)>(aa_list.AsParallel().AsOrdered(), a =>
                    (aa_seq: a, features: call_r_peptides(a))).ToList();
                var r_peptides_csv = r_peptides_results.SelectMany(a =>
                        a.features.Select(b =>
                            /*program.string_debug*/($@"{a.aa_seq},{String.Join(/*program.string_debug*/($@","), b.key_value_list(true).Select(c => c.value).ToList())}")).ToList())
                    .ToList();
                io_proxy.WriteAllLines(cache_r_peptides_filename, r_peptides_csv);
            }

            //Console.WriteLine();
        }

        internal static List<string> cache_r_get_seqs()
        {
            const string method_name = nameof(cache_r_get_seqs);

            var files = new string[] { cache_r_protr_filename, cache_r_peptides_filename };
            
            var aa_seqs = files.Select(a => File.ReadAllLines(a).Select(b => b.Substring(0, b.IndexOf(','))).Distinct().OrderBy(c => c.Length).ThenBy(c => c).ToList()).ToList();

            if (!aa_seqs[0].SequenceEqual(aa_seqs[1])) throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: cache source sequences differ"));

            return aa_seqs.FirstOrDefault();
        }

        internal static void cache_r_servers_check()
        {
            const string method_name = nameof(cache_r_servers_check);

            var files = new string[] { cache_r_protr_filename, cache_r_peptides_filename };

            var aa_seqs = files.Select(a => File.ReadAllLines(a).Select(b => b.Substring(0, b.IndexOf(','))).Distinct().OrderBy(c => c.Length).ThenBy(c => c).ToList()).ToList();

            var equal = aa_seqs[0].SequenceEqual(aa_seqs[1]);

            io_proxy.WriteLine(/*program.string_debug*/($@"SequenceEqual: {equal}"));

            //Console.WriteLine();

        }

        internal static void cache_r_load()
        {
            const string method_name = nameof(cache_r_load);

            if (File.Exists(subsequence_classification_data_r_methods.cache_r_protr_filename) && new FileInfo(subsequence_classification_data_r_methods.cache_r_protr_filename).Length > 0)
            {
                subsequence_classification_data_r_methods.cache_r_protr = subsequence_classification_data_r_methods.cache_r_load(subsequence_classification_data_r_methods.cache_r_protr_filename);
            }

            if (File.Exists(subsequence_classification_data_r_methods.cache_r_peptides_filename) && new FileInfo(subsequence_classification_data_r_methods.cache_r_peptides_filename).Length > 0)
            {
                subsequence_classification_data_r_methods.cache_r_peptides = subsequence_classification_data_r_methods.cache_r_load(subsequence_classification_data_r_methods.cache_r_peptides_filename);
            }
        }

        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_load(string cache_file)
        {
            const string method_name = nameof(cache_r_load);

            // note: features must be in original order.  more important that it is correct than fast.  parallel was not any faster relatively.

            var sw1 = new Stopwatch();
            sw1.Start();
            var lines = io_proxy.ReadAllLines(cache_file);

            var lines_parsed = lines.AsParallel().AsOrdered().Select(line =>
            {
                var line_split = line.Split(',');

                var col_index = 1;

                var fi = new feature_info()
                {
                    alphabet = line_split[col_index++],
                    stats = line_split[col_index++],
                    dimension = int.Parse(line_split[col_index++], NumberStyles.Integer, NumberFormatInfo.InvariantInfo),
                    category = line_split[col_index++],
                    source = line_split[col_index++],
                    @group = line_split[col_index++],
                    member = line_split[col_index++],
                    perspective = line_split[col_index++],
                    feature_value = double.Parse(line_split[col_index++], NumberStyles.Float, NumberFormatInfo.InvariantInfo)
                };

                return (aa_sequence: line_split[0], feature: fi);
            }).ToList();

            var groupings = lines_parsed.AsParallel().AsOrdered().GroupBy(a => a.aa_sequence).ToList();
            var features_grouped = groupings.AsParallel().AsOrdered().Select(a => (aa_sequence: a.Key, features: a.Select(b => b.feature).ToList())).ToList();

            var counts = features_grouped.Select(a => a.features.Count).Distinct().ToList();
            sw1.Stop();

            io_proxy.WriteLine(/*program.string_debug*/($@"Loaded: {features_grouped.Count} sequences ({string.Join(/*program.string_debug*/($@", "), counts)} features) in {sw1.Elapsed:dd\:hh\:mm\:ss\.fff}"), module_name, method_name);
            return features_grouped;
        }

        internal static List<feature_info> call_r_peptides(string sequence/*, string alphabet_name, enum_protein_data_source source*/, bool priority_boost = false, ProcessPriorityClass priority = ProcessPriorityClass.Idle, CancellationTokenSource cts = null)
        {
            const string method_name = nameof(call_r_peptides);

            using var i_cts = new CancellationTokenSource();
            if (cts == null) cts = i_cts;

            if (cache_r_peptides != null && cache_r_peptides.Count > 0)
            {
                var cache_item = cache_r_peptides.FirstOrDefault(a => string.Equals(a.aa_sequence, sequence, StringComparison.Ordinal));
                if (cache_item.features != null && cache_item.features.Count > 0)
                {
                    return cache_item.features;
                }
            }

            const int min_sequence_length = 1;
            var enable_cache = false;

            if (String.IsNullOrWhiteSpace(sequence) || sequence.Length < min_sequence_length)
            {
                if (subsequence_classification_data_templates._peptides_data_template == null)
                {
                    subsequence_classification_data_templates._peptides_data_template = call_r_peptides(/*program.string_debug*/($@"ALG")/*, alphabet_name, source*/);
                    subsequence_classification_data_templates._peptides_data_template.ForEach(a => { /*a.source = /*program.string_debug* /($@"");*/ a.feature_value = 0; });
                }

                if (subsequence_classification_data_templates._peptides_data_template == null)
                {
                    throw new Exception();
                }


                var template = subsequence_classification_data_templates._peptides_data_template.Select(a => new feature_info(a) {/* alphabet = alphabet_name, source = /*program.string_debug* /($@"{source}"),*/ feature_value = 0 }).ToList();

                return template;
            }

            if (enable_cache)
            {
                lock (subsequence_classification_data_totals._peptides_cache_lock)
                {
                    var cache_index =
                        subsequence_classification_data_totals._peptides_cache.FindIndex(a => string.Equals(a.sequence, sequence, StringComparison.Ordinal));
                    if (cache_index > -1)
                    {
                        var cached_item = subsequence_classification_data_totals._peptides_cache[cache_index].features
                            .Select(a => new feature_info(a)
                            {
                                /*alphabet = alphabet_name, source = source.ToString()*/
                            }).ToList();


                        if (cached_item != null && cached_item.Count > 0)
                        {
                            return cached_item;
                        }
                    }
                }
            }

            var call_count = 0;
            lock (subsequence_classification_data_totals._peptides_call_count_lock)
            {
                call_count = subsequence_classification_data_totals._peptides_call_count++;
            }

            var this_exe = System.Diagnostics.Process.GetCurrentProcess().MainModule.FileName;
            var exe = Path.Combine(Path.GetDirectoryName(this_exe.Replace(/*program.string_debug*/($@"dimorphics_dataset\dimorphics_dataset"), /*program.string_debug*/($@"dimorphics_dataset\peptides_server"), StringComparison.Ordinal)), Path.GetFileName(this_exe).Replace(/*program.string_debug*/($@"dimorphics_dataset"), /*program.string_debug*/($@"peptides_server"), StringComparison.Ordinal));

            var start = new ProcessStartInfo
            {
                FileName = exe,
                //Arguments = /*program.string_debug*/($@"{call_count} {alphabet_name} {source.ToString()} {sequence}",
                Arguments = /*program.string_debug*/($@"{call_count} {sequence}"),
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                //WindowStyle = ProcessWindowStyle.Hidden,

            };

            var result = new List<feature_info>();

            using (var process = Process.Start(start))
            {

                if (process != null)
                {
                    try
                    {
                        process.PriorityBoostEnabled = priority_boost;
                    }
                    catch (Exception e)
                    {
                        io_proxy.log_exception(e, $@"", module_name, method_name);
                    }

                    try
                    {
                        process.PriorityClass = priority;
                    }
                    catch (Exception e)
                    {
                        io_proxy.log_exception(e, $@"", module_name, method_name);
                    }


                    var stdout = process.StandardOutput.ReadToEnd();
                    var stderr = process.StandardError.ReadToEnd();

                    stdout = stdout[(stdout.IndexOf('\n', StringComparison.Ordinal) + 1)..];

                    process.WaitForExit();

                    //io_proxy.WriteLine(/*program.string_debug*/($@"Data: " + data);

                    var r = feature_info_container.deserialise(stdout);

                    result = r.feautre_info_list;
                }
            }

            if (enable_cache)
            {
                if (result != null && result.Count > 0)
                {
                    lock (subsequence_classification_data_totals._peptides_cache_lock)
                    {
                        subsequence_classification_data_totals._peptides_cache.Add((sequence, result.Select(a => new feature_info(a)).ToList()));
                    }
                }
            }

            if (subsequence_classification_data_templates._peptides_data_template == null)
            {
                subsequence_classification_data_templates._peptides_data_template = result.Select(a => new feature_info(a) { /*source = /*program.string_debug* /($@""),*/ feature_value = 0 }).ToList();
            }

            //result.ForEach(a =>
            //{
            //    a.source = source.ToString();
            //    a.alphabet = alphabet_name;
            //});


            return !cts.IsCancellationRequested ? result : default;
        }

        internal static List<feature_info> call_r_protr(string sequence/*, string alphabet_name, enum_protein_data_source source*/, bool priority_boost = false, ProcessPriorityClass priority = ProcessPriorityClass.Idle, CancellationTokenSource cts = null)
        {
            const string method_name = nameof(call_r_protr);

            using var i_cts = new CancellationTokenSource();
            if (cts == null) cts = i_cts;

            if (cache_r_protr != null && cache_r_protr.Count > 0)
            {
                var cache_item = cache_r_protr.FirstOrDefault(a => string.Equals(a.aa_sequence, sequence, StringComparison.Ordinal));
                if (cache_item.features != null && cache_item.features.Count > 0)
                {
                    return cache_item.features;
                }
            }

            //Console.WriteLine(sequence);
            const int min_sequence_length = 1;
            var enable_cache = false;

            if (String.IsNullOrWhiteSpace(sequence) || sequence.Length < min_sequence_length)
            {
                if (subsequence_classification_data_templates._protr_data_template == null)
                {
                    subsequence_classification_data_templates._protr_data_template = call_r_protr(/*program.string_debug*/($@"ALG")/*, alphabet_name, source*/);
                    subsequence_classification_data_templates._protr_data_template.ForEach(a => {/* a.source = /*program.string_debug* /($@"");*/ a.feature_value = 0; });
                }

                if (subsequence_classification_data_templates._protr_data_template == null)
                {
                    throw new Exception();
                }

                var template = subsequence_classification_data_templates._protr_data_template.Select(a => new feature_info(a) {/* alphabet = alphabet_name, source = /*program.string_debug* /($@"{source}"),*/ feature_value = 0 }).ToList();

                return template;
            }

            if (enable_cache)
            {
                lock (subsequence_classification_data_totals._protr_cache_lock)
                {
                    var cache_index =
                        subsequence_classification_data_totals._protr_cache.FindIndex(a => string.Equals(a.sequence, sequence, StringComparison.Ordinal));
                    if (cache_index > -1)
                    {
                        var cached_item = subsequence_classification_data_totals._protr_cache[cache_index].features
                            .Select(a => new feature_info(a)
                            {
                                /* alphabet = alphabet_name, source = source.ToString()*/
                            }).ToList();

                        if (cached_item != null && cached_item.Count > 0)
                        {
                            return cached_item;
                        }
                    }
                }
            }

            var call_count = 0;
            lock (subsequence_classification_data_totals._protr_call_count_lock)
            {
                call_count = subsequence_classification_data_totals._protr_call_count++;
            }

            var this_exe = System.Diagnostics.Process.GetCurrentProcess().MainModule.FileName;
            var exe = Path.Combine(Path.GetDirectoryName(this_exe.Replace(/*program.string_debug*/($@"dimorphics_dataset\dimorphics_dataset"), /*program.string_debug*/($@"dimorphics_dataset\protr_server"), StringComparison.Ordinal)), Path.GetFileName(this_exe).Replace(/*program.string_debug*/($@"dimorphics_dataset"), /*program.string_debug*/($@"protr_server"), StringComparison.Ordinal));

            var start = new ProcessStartInfo
            {
                FileName = exe,
                //Arguments = /*program.string_debug*/($@"{call_count} {alphabet_name} {source.ToString()} {sequence}",
                Arguments = /*program.string_debug*/($@"{call_count} {sequence}"),
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                //WindowStyle = ProcessWindowStyle.Hidden,

            };


            var result = new List<feature_info>();

            using (var process = Process.Start(start))
            {

                if (process != null)
                {
                    try
                    {
                        process.PriorityBoostEnabled = priority_boost;
                    }
                    catch (Exception e)
                    {
                        io_proxy.log_exception(e, $@"", module_name, method_name);
                    }

                    try
                    {
                        process.PriorityClass = priority;
                    }
                    catch (Exception e)
                    {
                        io_proxy.log_exception(e, $@"", module_name, method_name);
                    }


                    var stdout = process.StandardOutput.ReadToEnd();
                    var stderr = process.StandardError.ReadToEnd();

                    stdout = stdout[(stdout.IndexOf('\n', StringComparison.Ordinal) + 1)..];



                    process.WaitForExit();

                    //io_proxy.WriteLine(/*program.string_debug*/($@"Data: " + data);

                    var r = feature_info_container.deserialise(stdout);

                    result = r.feautre_info_list;
                }
            }

            if (enable_cache)
            {
                if (result != null && result.Count > 0)
                {
                    lock (subsequence_classification_data_totals._protr_cache_lock)
                    {
                        subsequence_classification_data_totals._protr_cache.Add((sequence,
                            result.Select(a => new feature_info(a)).ToList()));
                    }
                }
            }

            if (subsequence_classification_data_templates._protr_data_template == null)
            {
                subsequence_classification_data_templates._protr_data_template = result.Select(a => new feature_info(a) { source = /*program.string_debug*/($@""), feature_value = 0 }).ToList();
            }

            //result.ForEach(a =>
            //{
            //    a.source = source.ToString();
            //    a.alphabet = alphabet_name;
            //});

            return !cts.IsCancellationRequested ? result : default;
        }
    }
}
