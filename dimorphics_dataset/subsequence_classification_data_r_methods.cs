using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace dimorphics_dataset
{
    internal static class subsequence_classification_data_r_methods
    {
        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_protr;// cache_r_load(@"e:\dataset\r_protr_cache.csv");
        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_peptides;//cache_r_load(@"e:\dataset\r_peptides_cache.csv");

        internal static void cache_r_servers()
        {
            /*
            var files = new string[] { $@"E:\dataset\l__(standard_coil).csv", $@"E:\dataset\l__(dimorphic_coil).csv" };

            var psi_list = files.SelectMany(a => protein_subsequence_info.load(a)).ToList();

            var pdb_opt = new atom.load_atoms_pdb_options()
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
                            $"{a.aa_seq},{String.Join(",", b.AsArray(true).Select(c => c.value).ToList())}").ToList())
                    .ToList();
                io_proxy.WriteAllLines($@"e:\dataset\r_protr_cache.csv", r_protr_csv);
            }

            {
                var r_peptides_results = ParallelEnumerable.Select<string, (string aa_seq, List<feature_info> features)>(aa_list.AsParallel().AsOrdered(), a =>
                    (aa_seq: a, features: call_r_peptides(a))).ToList();
                var r_peptides_csv = r_peptides_results.SelectMany(a =>
                        a.features.Select(b =>
                            $"{a.aa_seq},{String.Join(",", b.AsArray(true).Select(c => c.value).ToList())}").ToList())
                    .ToList();
                io_proxy.WriteAllLines($@"e:\dataset\r_peptides_cache.csv", r_peptides_csv);
            }

            //Console.WriteLine();
        }

        internal static List<string> cache_r_get_seqs()
        {
            var files = new string[] { $@"e:\dataset\r_protr_cache.csv" };

            var aa_seqs = files.SelectMany(a => File.ReadAllLines(a).Select(b => b.Substring(0, b.IndexOf(','))).Distinct().OrderBy(a => a.Length).ThenBy(a => a).ToList()).ToList();

            return aa_seqs;
        }

        internal static void cache_r_servers_check()
        {
            var files = new string[] { $@"e:\dataset\r_protr_cache.csv", $@"e:\dataset\r_peptides_cache.csv" };

            var aa_seqs = files.Select(a => File.ReadAllLines(a).Select(b => b.Substring(0, b.IndexOf(','))).Distinct().OrderBy(a => a.Length).ThenBy(a => a).ToList()).ToList();

            var equal = aa_seqs[0].SequenceEqual(aa_seqs[1]);

            io_proxy.WriteLine("SequenceEqual: " + equal);

            Console.WriteLine();

        }




        internal static List<(string aa_sequence, List<feature_info> features)> cache_r_load(string cache_file)
        {
            // note: features must be in original order.  more important that it is correct than fast.  parallel was not any faster relatively.

            var sw1 = new Stopwatch();
            sw1.Start();
            var a = io_proxy.ReadAllLines(cache_file);

            var b = a/*AsParallel().AsOrdered()*/.Select(y =>
            {
                var x = y.Split(',');
                var z = new feature_info()
                {
                    alphabet = x[1],
                    dimension = Int32.Parse(x[2], NumberStyles.Integer, NumberFormatInfo.InvariantInfo),
                    category = x[3],
                    source = x[4],
                    @group = x[5],
                    member = x[6],
                    perspective = x[7],
                    feature_value = Double.Parse(x[8], NumberStyles.Float, NumberFormatInfo.InvariantInfo)
                };
                return (aa_sequence: x[0], feature: z);
            }
            ).ToList();
            var c = b/*.AsParallel()./*AsOrdered()*/.GroupBy(a => a.aa_sequence).ToList();
            var d = c/*.AsParallel()./*AsOrdered().*/.Select(a => (aa_sequence: a.Key, features: a.Select(b => b.feature).ToList())).ToList();

            var counts = d.Select(a => a.features.Count).Distinct().ToList();
            sw1.Stop();

            io_proxy.WriteLine($@"Loaded: {d.Count} sequences ({String.Join(", ", counts)} features) in {sw1.Elapsed:dd\:hh\:mm\:ss\.fff}", nameof(program), nameof(cache_r_load));
            return d;
        }


        internal static List<feature_info> call_r_peptides(string sequence/*, string alphabet_name, enum_protein_data_source source*/, bool priority_boost = false, ProcessPriorityClass priority = ProcessPriorityClass.Idle)
        {
            if (cache_r_peptides != null && cache_r_peptides.Count > 0)
            {
                var cache_item = cache_r_peptides.FirstOrDefault(a => a.aa_sequence == sequence);
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
                    subsequence_classification_data_templates._peptides_data_template = call_r_peptides("ALG"/*, alphabet_name, source*/);
                    subsequence_classification_data_templates._peptides_data_template.ForEach(a => { /*a.source = "";*/ a.feature_value = 0; });
                }

                if (subsequence_classification_data_templates._peptides_data_template == null)
                {
                    throw new Exception();
                }


                var template = subsequence_classification_data_templates._peptides_data_template.Select(a => new feature_info(a) {/* alphabet = alphabet_name, source = source.ToString(),*/ feature_value = 0 }).ToList();

                return template;
            }

            if (enable_cache)
            {
                lock (subsequence_classification_data_totals._peptides_cache_lock)
                {
                    var cache_index =
                        subsequence_classification_data_totals._peptides_cache.FindIndex(a => a.sequence == sequence);
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
            var exe = Path.Combine(Path.GetDirectoryName(this_exe.Replace($@"dimorphics_dataset\dimorphics_dataset", $@"dimorphics_dataset\peptides_server", StringComparison.InvariantCulture)), Path.GetFileName(this_exe).Replace("dimorphics_dataset", "peptides_server", StringComparison.InvariantCulture));

            var start = new ProcessStartInfo
            {
                FileName = exe,
                //Arguments = $"{call_count} {alphabet_name} {source.ToString()} {sequence}",
                Arguments = $"{call_count} {sequence}",
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
                    catch (Exception)
                    {
                    }

                    try
                    {
                        process.PriorityClass = priority;
                    }
                    catch (Exception)
                    {
                    }


                    var stdout = process.StandardOutput.ReadToEnd();
                    var stderr = process.StandardError.ReadToEnd();

                    stdout = stdout[(stdout.IndexOf('\n', StringComparison.InvariantCulture) + 1)..];



                    process.WaitForExit();

                    //io_proxy.WriteLine("Data: " + data);

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
                subsequence_classification_data_templates._peptides_data_template = result.Select(a => new feature_info(a) { /*source = "",*/ feature_value = 0 }).ToList();
            }

            //result.ForEach(a =>
            //{
            //    a.source = source.ToString();
            //    a.alphabet = alphabet_name;
            //});


            return result;
        }

        internal static List<feature_info> call_r_protr(string sequence/*, string alphabet_name, enum_protein_data_source source*/, bool priority_boost = false, ProcessPriorityClass priority = ProcessPriorityClass.Idle)
        {
            if (cache_r_protr != null && cache_r_protr.Count > 0)
            {
                var cache_item = cache_r_protr.FirstOrDefault(a => a.aa_sequence == sequence);
                if (cache_item.features != null && cache_item.features.Count > 0)
                {
                    return cache_item.features;
                }
            }

            Console.WriteLine(sequence);
            const int min_sequence_length = 1;
            var enable_cache = false;

            if (String.IsNullOrWhiteSpace(sequence) || sequence.Length < min_sequence_length)
            {
                if (subsequence_classification_data_templates._protr_data_template == null)
                {
                    subsequence_classification_data_templates._protr_data_template = call_r_protr("ALG"/*, alphabet_name, source*/);
                    subsequence_classification_data_templates._protr_data_template.ForEach(a => {/* a.source = "";*/ a.feature_value = 0; });
                }

                if (subsequence_classification_data_templates._protr_data_template == null)
                {
                    throw new Exception();
                }

                var template = subsequence_classification_data_templates._protr_data_template.Select(a => new feature_info(a) {/* alphabet = alphabet_name, source = source.ToString(),*/ feature_value = 0 }).ToList();

                return template;
            }

            if (enable_cache)
            {
                lock (subsequence_classification_data_totals._protr_cache_lock)
                {
                    var cache_index =
                        subsequence_classification_data_totals._protr_cache.FindIndex(a => a.sequence == sequence);
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
            var exe = Path.Combine(Path.GetDirectoryName(this_exe.Replace($@"dimorphics_dataset\dimorphics_dataset", $@"dimorphics_dataset\protr_server", StringComparison.InvariantCulture)), Path.GetFileName(this_exe).Replace("dimorphics_dataset", "protr_server", StringComparison.InvariantCulture));

            var start = new ProcessStartInfo
            {
                FileName = exe,
                //Arguments = $"{call_count} {alphabet_name} {source.ToString()} {sequence}",
                Arguments = $"{call_count} {sequence}",
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
                    catch (Exception)
                    {

                    }

                    try
                    {
                        process.PriorityClass = priority;
                    }
                    catch (Exception)
                    {

                    }


                    var stdout = process.StandardOutput.ReadToEnd();
                    var stderr = process.StandardError.ReadToEnd();

                    stdout = stdout[(stdout.IndexOf('\n', StringComparison.InvariantCulture) + 1)..];



                    process.WaitForExit();

                    //io_proxy.WriteLine("Data: " + data);

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
                subsequence_classification_data_templates._protr_data_template = result.Select(a => new feature_info(a) { source = "", feature_value = 0 }).ToList();
            }

            //result.ForEach(a =>
            //{
            //    a.source = source.ToString();
            //    a.alphabet = alphabet_name;
            //});

            return result;
        }
    }
}
