using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Loader;
using System.Security.Cryptography;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class program
    {
        // todo: re-check tortusoity code, intramolecular distance code, and averageatom-atom distance cpde

        // todo: check which blast databases, output options are disabled.
        // todo: check blast distance option.

        internal static bool verbose = true;

        internal static string data_root_folder = $@"C:\betastrands_dataset\";

        internal static List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> get_dataset_pdb_id_list(string class_name_in_file = @"Single")
        {
            const string module_name = nameof(program);
            const string method_name = nameof(get_dataset_pdb_id_list);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            //List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)>

            var dimorphics_data1 = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, $@"csv", $@"distinct dimorphics list.csv"), module_name, method_name)
                .Skip(1)
                .Where(a => !string.IsNullOrWhiteSpace(a.Replace($@",", $@"", StringComparison.Ordinal)))
                .Select((a, i) =>
                {
                    var k = 0;
                    var x = a.Split(',');
                    return (
                        pdb_id: x[k++].ToUpperInvariant(),
                        dimer_type: x[k++],
                        class_name: x[k++],
                        symmetry_mode: x[k++],
                        parallelism: x[k++],
                        chain_number: int.Parse(x[k++], NumberStyles.Integer, NumberFormatInfo.InvariantInfo) - 1,
                        strand_seq: x[k++],
                        optional_res_index: x[k++]
                    );
                }).ToList();

            if (!string.IsNullOrWhiteSpace(class_name_in_file))
            {
                dimorphics_data1 = dimorphics_data1.Where(a => string.Equals(a.class_name, class_name_in_file, StringComparison.OrdinalIgnoreCase)).ToList();
            }

            //var pdb_id_list = dimorphics_data1.Select(a => (dimer_type: a.dimer_type, pdb_id: a.pdb_id)).Distinct().ToList();
            //
            //return pdb_id_list;

            return dimorphics_data1;
        }

        internal static subsequence_classification_data get_subsequence_classificiation_data(cmd_params p, protein_subsequence_info psi, protein_subsequence_info _template_protein = null, load_atoms_pdb_options pdb_options = null, bool load_ss_predictions = true, bool load_foldx_energy = true)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(get_subsequence_classificiation_data);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            // todo: check if pdb_folder should be pdb or pdb_repair?

            if (psi == null)
            {
                throw new ArgumentNullException(nameof(psi));
            }

            //if (feature_types == null)
            //{
            //    throw new ArgumentNullException(nameof(feature_types));
            //}

            //var features_1d = (feature_types?.feature_types_2d_interface?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_neighbourhood_2d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_chain_2d?.AsArray()?.Any(a => a.value) ?? false);

            //var features_3d = (feature_types?.feature_types_3d_interface?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_chain_3d?.AsArray()?.Any(a => a.value) ?? false);

            var pdb_atoms = atom.load_atoms_pdb
                (
                    pdb_id: psi.pdb_id,
                    pdb_options ??
                    new load_atoms_pdb_options()
                    {
                        //load_2d_mpsa_sec_struct_predictions = features_1d || features_3d,
                        //load_2d_blast_pssms = features_1d || features_3d,
                        //load_2d_iup_data = features_1d || features_3d,
                        //load_2d_sable = features_1d || features_3d,
                        //load_2d_dna_binding = features_1d || features_3d,

                        //find_3d_intramolecular = features_1d || features_3d,
                        //load_3d_dssp_data = features_1d || features_3d,
                        //load_3d_stride_data = features_1d || features_3d,
                        //load_3d_ring_data = features_1d || features_3d,
                        //load_3d_foldx_ala_scan = features_1d || features_3d,
                        //load_3d_rsa_data = features_1d || features_3d,


                        load_1d_blast_pssms = true,
                        load_1d_iup_data = true,
                        load_1d_sable = true,
                        load_1d_dna_binding = true,

                        load_2d_mpsa_sec_struct_predictions = true,
                        
                        find_3d_intramolecular = true,
                        load_3d_dssp_data = true,
                        load_3d_stride_data = true,
                        load_3d_ring_data = true,
                        load_3d_foldx_ala_scan = true,
                        load_3d_rsa_data = true,
                    }
                    )
                .Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();

            var chain_atoms = pdb_atoms.Where(a => a.chain_id == psi.chain_id).ToList();
            var interface_atoms = psi.res_ids.SelectMany(a => chain_atoms.Where(b => a.res_id == b.residue_index && a.amino_acid == b.amino_acid && a.i_code == b.i_code).ToList()).ToList();

            var scd = new subsequence_classification_data
            {
                dimer_type = psi.dimer_type ?? $@"",
                symmetry_mode = psi.symmetry_mode ?? $@"",
                parallelism = psi.parallelism ?? $@"",

                class_id = psi.class_id,
                class_name = psi.class_name,

                pdb_id = psi.pdb_id,
                chain_id = psi.chain_id,
            };

            scd.interface_region = new subsequence_classification_data_region(scd, interface_atoms, load_ss_predictions, load_foldx_energy);
            scd.chain_region = new subsequence_classification_data_region(scd, chain_atoms, load_ss_predictions, load_foldx_energy);

            // after the interface and chain properties are set, can then find neighborhood...
            scd.init_nh_flanking(6, load_ss_predictions, load_foldx_energy);
            scd.init_nh_contacts(5.0, false, load_ss_predictions, load_foldx_energy);

            var need_template = (
                                    string.IsNullOrEmpty(scd.interface_region.aa_sequence) ||
                                    string.IsNullOrEmpty(scd.chain_region.aa_sequence) ||
                                    string.IsNullOrEmpty(scd.nh_flank_region.aa_sequence) ||
                                    string.IsNullOrEmpty(scd.nh_contact_region.aa_sequence)
                                    );



            //if (feature_types?.feature_types_neighbourhood_1d?.AsArray()?.Any(a => a.value)??false)
            //{
            //    scd.neighbourhood_1d = atom.get_intramolecular_neighbourhood_1d(scd);
            //    if (string.IsNullOrEmpty(scd.neighbourhood_1d.aa_subsequence) || scd.neighbourhood_1d.interface_atoms?.Count == 0) {need_template = true;}

            //}

            //if (feature_types?.feature_types_protein_1d?.AsArray()?.Any(a => a.value)??false)
            //{
            //    scd.protein_1d = atom.get_intramolecular_protein_1d(scd);
            //    if (string.IsNullOrEmpty(scd.protein_1d.aa_subsequence) || scd.protein_1d.interface_atoms?.Count == 0) {need_template = true;}

            //}

            //if (feature_types?.feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value)??false)
            //{
            //    scd.neighbourhood_3d = atom.get_intramolecular_neighbourhood_3d(scd);
            //    if (string.IsNullOrEmpty(scd.neighbourhood_3d.aa_subsequence) || scd.neighbourhood_3d.interface_atoms?.Count == 0) {need_template = true;}

            //}

            //if (feature_types?.feature_types_protein_3d?.AsArray()?.Any(a => a.value)??false)
            //{
            //    scd.protein_3d = atom.get_intramolecular_protein_3d(scd);
            //    if (string.IsNullOrEmpty(scd.protein_3d.aa_subsequence) || scd.protein_3d.interface_atoms?.Count == 0) {need_template = true;}
            //}

            if (need_template && _template_protein != null && subsequence_classification_data_templates._template_scd == null)
            {
                subsequence_classification_data_templates._template_scd = get_subsequence_classificiation_data(p, _template_protein, null, pdb_options, load_ss_predictions, load_foldx_energy);
            }

            return scd;
        }

        internal static string join(IList<string> items)
        {
            if (items == null || items.Count == 0) return $@"";

            return string.Join($@"_", items.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()).Replace($@".", $@"_", StringComparison.Ordinal);
        }

        internal static void print_eta(int num_complete, int num_total, DateTime start_time, string module_name = @"", string method_name = @"")
        {
            var num_incomplete = num_total - num_complete;

            var num_complete_pct = num_total > 0 ? ((double)num_complete / (double)num_total) * 100 : 0;

            var time_taken_ticks = DateTime.Now.Subtract(start_time).Ticks;
            var time_taken = TimeSpan.FromTicks(time_taken_ticks);

            var est_time_each = num_complete > 0 ? (double)time_taken_ticks / (double)num_complete : 0;
            var est_time_remain = (double)est_time_each * (double)num_incomplete;

            var time_remaining = TimeSpan.FromTicks((long)est_time_remain);

            io_proxy.WriteLine($@"ETA: tasks complete: {num_complete}/{num_total} ({num_complete_pct:0.00}%) [time taken: {time_taken:dd\:hh\:mm\:ss\.fff}, time remaining: {time_remaining:dd\:hh\:mm\:ss\.fff}]", module_name, method_name);
        }

        internal static void wait_tasks(Task[] tasks, DateTime start_time, string module_name = @"", string method_name = @"")
        {
            if (tasks == null || tasks.Length == 0) return;

            do
            {
                var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToList();

                if (incomplete_tasks.Count > 0)
                {
                    var completed_task_index = Task.WaitAny(incomplete_tasks.ToArray<Task>());

                    if (completed_task_index > -1)
                    {
                        var num_total = tasks.Length;
                        var num_incomplete = tasks.Count(a => !a.IsCompleted);
                        var num_complete = num_total - num_incomplete;

                        print_eta(num_complete, num_total, start_time, module_name, method_name);
                    }
                }

            } while (tasks.Any(a => !a.IsCompleted));

            Task.WaitAll(tasks.ToArray<Task>());
        }

        internal static void close_notifications()
        {
            Console.CancelKeyPress += (sender, eventArgs) =>
            {
                io_proxy.WriteLine($@"Console.CancelKeyPress", nameof(program), nameof(close_notifications));
            };
            AssemblyLoadContext.Default.Unloading += context =>
            {
                io_proxy.WriteLine($@"AssemblyLoadContext.Default.Unloading", nameof(program), nameof(close_notifications));
            };
            AppDomain.CurrentDomain.ProcessExit += (sender, eventArgs) =>
            {
                io_proxy.WriteLine($@"AppDomain.CurrentDomain.ProcessExit", nameof(program), nameof(close_notifications));
            };
        }




        //var pdbs = get_dataset_pdb_id_list().Select(a => a.pdb_id).Distinct().ToList();
        //var at = pdbs.SelectMany(a => atom.load_atoms_pdb(a, new load_atoms_pdb_options()
        //{
        //    load_2d_mpsa_sec_struct_predictions = false,
        //    load_3d_dssp_data = false,
        //    load_2d_iup_data = false,
        //    load_2d_blast_pssms = false,

        //    first_model_only = true,
        //    first_icode_only = true,

        //    find_3d_intramolecular = false,
        //    find_3d_intermolecular = false,


        //    load_3d_foldx_ala_scan = false,

        //    load_3d_stride_data = false,
        //    load_sable = false,
        //    load_3d_ring_data = false,

        //    load_dna_binding_vars = false,
        //    load_2d_rsa_data = false,
        //})).SelectMany(a => a.pdb_model_chain_atoms.Select(b => (b.pdb_id,b.atom_type)).ToList()).ToList();

        //var atd = at.Select(a=>a.atom_type).Distinct().ToList();
        //var atc = atd.Select(a => (atom_type:a, count:at.Count(b => b.atom_type == a), files:at.Where(b=>b.atom_type==a).Select(b=>b.pdb_id).Distinct().Count())).OrderByDescending(a=>a.count).ToList();
        //atc.ForEach(a=>io_proxy.WriteLine($@"{a.atom_type},{a.count},{((double)a.count / (double)atd.Count)},{a.files},{((double)a.files / (double)pdbs.Count)}"));
        //return;
        //foreach (var pdb in pdbs)
        //{
        //    var ppp = atom.load_atoms_pdb(pdb, new load_atoms_pdb_options()
        //    {
        //        load_2d_mpsa_sec_struct_predictions = false,
        //        load_3d_dssp_data = false,
        //        load_2d_iup_data = false,
        //        load_2d_blast_pssms = false,

        //        first_model_only = true,
        //        first_icode_only = true,

        //        find_3d_intramolecular = false,
        //        find_3d_intermolecular = false,


        //        load_3d_foldx_ala_scan = false,

        //        load_3d_stride_data = false,
        //        load_sable = false,
        //        load_3d_ring_data = false,

        //        load_dna_binding_vars = false,
        //        load_2d_rsa_data = false,

        //    });

        //foreach (var pp in ppp)
        //{
        //    var p = pp.pdb_model_chain_atoms;

        //    var d = new double[p.Count, p.Count];
        //    for (var i = 0; i < p.Count; i++)
        //    {
        //        for (var j = 0; j < p.Count; j++)
        //        {
        //            if (i < j) continue;

        //            d[i, j] = atom.Distance3D(p[i], p[j]);
        //            d[j, i] = d[i, j];
        //        }
        //    }

        //    io_proxy.WriteLine();
        //}
        //}

        //return;





        internal static void Main(string[] args)
        {
            //var x2 = descriptive_stats_encoding_options.options_all_plus;
            //var x1 = descriptive_stats_encoding_options.options_all;
            //var x3 = descriptive_stats_encoding_options.options_default;

            //return;
            //var nums = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 10 };//{1.34, 353253.3, 4343.2344, 12.11111199, 9.1, 9.99, 0.01, 0.0001, 54654, 745637, 234234, 745634, 24, 54563476, 8090, 6453};

            //nums = new double[] { 222, 1122, 45444 };

            ////nums = Enumerable.Range(0,10).SelectMany(j=> nums.Concat(nums.Select((a, i) => a + ((i+1) * (j+1))).ToArray()).ToArray()).ToArray();
            ////nums = nums.Where(a => a != 0).ToArray();

            //var nums_ds = descriptive_stats.get_stat_values(data: nums, group_id_name: $@"", member_id_name: $@"", presorted: false, interval: true, distance: true, interquartile: true, abs: true, rescale:true);
            //var nums_ds_e = nums_ds.encode(new descriptive_stats_encoding_options($@"", true),  presorted: false, interval: true, distance: true, interquartile: true, abs: true, rescale:true);
            //nums_ds_e.ForEach(a=> Console.WriteLine($"{a.group_id} {a.member_id} {a.perspective_id} = {a.perspective_value}"));
            //return;


            const string module_name = nameof(program);
            const string method_name = nameof(Main);

            //var swx = new Stopwatch();
            //swx.Start();
            //var x = cache_r_load();
            //swx.Stop();
            //Console.WriteLine(swx.Elapsed.ToString($@"g"));

            //cache_r_servers();
            //cache_r_servers_check();
            //return;

            close_notifications();
            var child_priority = ProcessPriorityClass.Idle;
            var child_priority_boost = false;
            var p = new cmd_params(args);

            //var tag = p.get_output_file_tag();

            if (!p.parse_ok)
            {
                return;
            }

            program.verbose = p?.verbose.Value ?? true;

            //var feature_types2 = feature_types.feature_types_params(p?.area ?? null);

            //return;
            var sw1 = new Stopwatch();
            sw1.Start();

            var pdb_id_list = get_dataset_pdb_id_list();

            // 1. find subsequence details from analysis of pdb files
            var psi_list = part1_find_samples(p, pdb_id_list);

            // load r cache, if not using children (takes too many resources to load cache for individual items)
            if (!p.use_children.Value)
            {
                subsequence_classification_data_r_methods.cache_r_protr = subsequence_classification_data_r_methods.cache_r_load(subsequence_classification_data_r_methods.cache_r_protr_filename);
                subsequence_classification_data_r_methods.cache_r_peptides = subsequence_classification_data_r_methods.cache_r_load(subsequence_classification_data_r_methods.cache_r_peptides_filename);
            }

            var _template_protein = psi_list.First();
            //Parallel.For(0, psi_list.Count, index =>

            var tasks = new List<Task>();

            var tasks_start_time = DateTime.Now;

            for (var index = 0; index < psi_list.Count; index++)
            {
                var _index = index;

                if (p.first_index != null && p.first_index > -1)
                {
                    if (_index < p.first_index) continue;
                }

                if (p.last_index != null && p.last_index > -1)
                {
                    if (_index > p.last_index) continue;
                }

                var fns = get_output_filenames(p, _index);
                if (File.Exists(fns.fn_headers) && File.Exists(fns.fn_features) && File.Exists(fns.fn_comments)) continue;


                if (p.first_index == null && p.last_index == null)
                {
                    if (p.use_children.Value)
                    {
                        var task = Task.Run(() =>
                        {
                            while (true)
                            {
                                var cmd_line = $@"{Process.GetCurrentProcess().MainModule.FileName} {string.Join($@" ", args)} -first_index={_index} -last_index={_index}";

                                io_proxy.WriteLine(cmd_line, module_name, method_name);

                                using var process = Process.Start(Process.GetCurrentProcess().MainModule.FileName, $@"{string.Join($@" ", args)} -first_index={_index} -last_index={_index}");

                                if (process != null)
                                {
                                    try { process.PriorityClass = child_priority; } catch (Exception) { }
                                    try { process.PriorityBoostEnabled = child_priority_boost; } catch (Exception) { }

                                    process.WaitForExit();

                                    return;
                                }
                                else
                                {
                                    Task.Delay(new TimeSpan(0, 0, 5)).Wait();
                                }
                            }
                        });
                        tasks.Add(task);

                        var num_total = psi_list.Count;

                        var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToArray();
                        if (incomplete_tasks.Length >= Environment.ProcessorCount)
                        {
                            print_eta(num_total - incomplete_tasks.Length, num_total, tasks_start_time, module_name, method_name);
                            Task.WaitAny(incomplete_tasks);
                        }

                        var completed_tasks = tasks.Where(a => a.IsCompleted).ToArray();
                        print_eta(completed_tasks.Length, num_total, tasks_start_time, module_name, method_name);

                        continue;
                    }
                    /*else
                    {
                        var task = Task.Run(() =>
                        {
                            var args2 = args.Union(new string[] { $@"-first_index={_index}", $@"-last_index={_index}" }).ToArray();
                            Main(args2);
                        });
                        tasks.Add(task);

                        var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToArray();
                        if (incomplete_tasks.Length >= Environment.ProcessorCount)
                        {
                            Task.WaitAny(incomplete_tasks);
                        }

                        continue;
                    }*/
                    
                   
                }




                var psi_list2 = new List<protein_subsequence_info>() { psi_list[_index] };

                // 2. load available data sources
                var class_data_list = part2_load_data(p, psi_list2, _template_protein);

                // 3. encode the data as svm classification features
                var data_encoded_list = part3_encode_features(p, class_data_list);

                // 4. Save to file
                part4_save_outputs(p, data_encoded_list, _index);
            }

            wait_tasks(tasks.ToArray<Task>(), tasks_start_time, module_name, method_name);

            if (p.first_index == null && p.last_index == null)
            {
                part5_check_output_hashes(psi_list, p);

                part6_check_and_merge_outputs(psi_list, p);
            }

            sw1.Stop();

            if (p.first_index == null && p.last_index == null)
            {

                io_proxy.WriteLine($@"{p.class_id} {p.class_name}: Finished: ({sw1.Elapsed:dd\:hh\:mm\:ss\.fff})", module_name, method_name);
            }
        }
         
        private static void part6_check_and_merge_outputs(List<protein_subsequence_info> psi_list, cmd_params p)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part6_check_and_merge_outputs);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);

            var files = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(p, i);
                return (fn_headers: b.fn_headers, fn_features: b.fn_features, fn_comments: b.fn_comments);
            }).ToList();

            var text_header = io_proxy.ReadAllLines(files.First().fn_headers, module_name, method_name);
            var text_features = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_features, module_name, method_name).Skip(i == 0 ? 0 : 1).ToList()).ToList();
            var text_comments = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_comments, module_name, method_name).Skip(i == 0 ? 0 : 1).ToList()).ToList();

            // check for any lines with different total number of features
            var num_features = text_features.Select((a, i) => (line: i, num_features: a.Count(b => b == ','))).OrderByDescending(a => a.num_features).ToList();

            if (num_features.Select(a => a.num_features).Distinct().Count() > 1)
            {
                io_proxy.WriteLine($@"Feature counts don't match.", module_name, method_name);
            }
            else
            {
                io_proxy.WriteLine($@"Feature counts match.", module_name, method_name);
            }

            // check for any non-double type features
            var invalid_features = text_features
                .SelectMany(a =>
                    a.Split(',').Select((b, i) => (string.IsNullOrEmpty(b) || string.Equals(b, $@"∞", StringComparison.Ordinal) || string.Equals(b, $@"+∞", StringComparison.Ordinal) || string.Equals(b, $@"-∞", StringComparison.Ordinal) || string.Equals(b, $@"NaN", StringComparison.Ordinal)) ? i : -1)
                        .Where(a => a != -1).ToList()).Distinct().OrderBy(a => a).ToList();

            if (invalid_features.Count > 0)
            {
                io_proxy.WriteLine($@"Invalid feature ids: {string.Join($@", ", invalid_features)}", module_name, method_name);

                foreach (var i in invalid_features)
                {
                    var invalid_feature_header =
                        text_header.FirstOrDefault(a => a.StartsWith($@"{i},", StringComparison.Ordinal));
                    io_proxy.WriteLine($@"{i}: {invalid_feature_header}");
                }
            }
            else
            {
                io_proxy.WriteLine($@"No invalid feature ids.", module_name, method_name);
            }

            var fns = get_output_filenames(p, null);
            io_proxy.WriteAllLines(fns.fn_headers, text_header);
            io_proxy.WriteAllLines(fns.fn_features, text_features);
            io_proxy.WriteAllLines(fns.fn_comments, text_comments);
        }

        private static void part5_check_output_hashes(List<protein_subsequence_info> psi_list, cmd_params p)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part5_check_output_hashes);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            var file_hashes = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(p, i);
                return (fn_headers: file_sha_256(b.fn_headers), fn_features: file_sha_256(b.fn_features), fn_comments: file_sha_256(b.fn_comments));
            }).ToList();

            var file_hashes_headers = file_hashes.Select(a => a.fn_headers).ToList();
            var file_hashes_features = file_hashes.Select(a => a.fn_features).ToList();
            var file_hashes_comments = file_hashes.Select(a => a.fn_comments).ToList();

            var file_hashes_headers_cnt = file_hashes_headers.Distinct()
                .Select(a => (hash: a, count: file_hashes_headers.Count(b => string.Equals(a, b, StringComparison.Ordinal))))
                .OrderByDescending(a => a.count).ToList();

            var file_hashes_features_cnt = file_hashes_features.Distinct()
                .Select(a => (hash: a, count: file_hashes_features.Count(b => string.Equals(a, b, StringComparison.Ordinal))))
                .OrderByDescending(a => a.count).ToList();

            var file_hashes_comments_cnt = file_hashes_comments.Distinct()
                .Select(a => (hash: a, count: file_hashes_comments.Count(b => string.Equals(a, b, StringComparison.Ordinal))))
                .OrderByDescending(a => a.count).ToList();

            if (file_hashes_headers_cnt.Count != 1)
            {
                io_proxy.WriteLine($@"Hashes for headers:");
                file_hashes_headers_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }

            if (file_hashes_features_cnt.Count != psi_list.Count)
            {
                io_proxy.WriteLine($@"Hashes for features:");
                file_hashes_features_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }

            if (file_hashes_comments_cnt.Count != psi_list.Count)
            {
                io_proxy.WriteLine($@"Hashes for comments:");
                file_hashes_comments_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }
        }

        private static string file_sha_256(string filename)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(file_sha_256);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            using var stream = new BufferedStream(File.OpenRead(filename), 1024 * 1024 * 4);
            using var sha256 = new SHA256Managed();
            var hash_bytes = sha256.ComputeHash(stream);
            var hash_str = BitConverter.ToString(hash_bytes);//.Replace($@"-", String.Empty);
            return hash_str;
        }

       

        private static string get_input_filenames(cmd_params p, int? index = null)
        {
            var fn_input = Path.Combine(p.output_folder, $@"l__({(p.class_id > 0 ? $@"+" : $@"")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            return fn_input;
        }

        private static (string fn_headers, string fn_comments, string fn_features) get_output_filenames(cmd_params p, int? index)
        {
            if (p==null) throw new ArgumentNullException(nameof(p));
            //if (feature_types == null) throw new ArgumentNullException(nameof(feature_types));

            var tag = p.get_output_file_tag();

            //var fn_input = Path.Combine(cmd_params.output_folder, $@"l__({cmd_params.class_name}){(tag != null ? $@"_{tag}" : $@"")}.csv");
            var fn_headers =  Path.Combine(p.output_folder, $@"h_({(tag ?? $@"")})_({(p.class_id > 0 ? $@"+" : $@"")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            var fn_comments = Path.Combine(p.output_folder, $@"c_({(tag ?? $@"")})_({(p.class_id > 0 ? $@"+" : $@"")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            var fn_features = Path.Combine(p.output_folder, $@"f_({(tag ?? $@"")})_({(p.class_id > 0 ? $@"+" : $@"")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");

            return (fn_headers, fn_comments, fn_features);
        }

        private static void part4_save_outputs(cmd_params p, List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> data_encoded_list, int item_index)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part4_save_outputs);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);



            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine($@"{p.class_id} {p.class_name}: 4. Saving encoded data to file for {data_encoded_list.Count} items...", module_name, method_name);
            }

            var fns = get_output_filenames(p, item_index);

            // get header row indexes in csv format
            var row_feature_header_csv = string.Join($@",", Enumerable.Range(0, data_encoded_list.First().feature_info.Count));

            // get comments file header in csv format
            var row_comments_header_csv = string.Join($@",", subsequence_classification_data.get_row_comments_headers(data_encoded_list.First().instance_meta_data));

            // get list of the feature headers in csv format
            // save headers
            var feature_headers = feature_info.get_feature_headers_lines_csv(data_encoded_list.First().feature_info);
            io_proxy.WriteAllLines(fns.fn_headers, feature_headers, module_name, method_name);


            var tasks4 = new List<Task<(string row_feature_values_csv, string row_comments_csv)>>();
            var tasks4_start_time = DateTime.Now;

            for (var i = 0; i < data_encoded_list.Count; i++)
            {
                //var row_index = i;
                var row_index = i + item_index;

                var data_encoded_item = data_encoded_list[i];

                data_encoded_list[i] = default;

                var task = Task.Run(() =>
                {
                    // get feature values
                    var row_feature_values = data_encoded_item.feature_info.Select((a, fid) => $"{a.feature_value:G17}").ToList();

                    // convert feature values to csv format
                    var row_feature_values_csv = string.Join($@",", row_feature_values);

                    // get meta data about the example instance
                    var row_comments = subsequence_classification_data.get_row_comments(row_index, data_encoded_item.instance_meta_data);

                    var row_comments_csv = string.Join($@",", row_comments);

                    return (row_feature_values_csv, row_comments_csv);
                });

                tasks4.Add(task);

                wait_tasks(tasks4.ToArray<Task>(), tasks4_start_time, module_name, method_name);
            }

            data_encoded_list.Clear();
            data_encoded_list = null;
            wait_tasks(tasks4.ToArray<Task>(), tasks4_start_time, module_name, method_name);

            // save comments
            var comments_lines = tasks4.Select(a => a.Result.row_comments_csv).ToList();
            comments_lines.Insert(0, row_comments_header_csv);
            io_proxy.WriteAllLines(fns.fn_comments, comments_lines, module_name, method_name);

            // save features
            var features_lines = tasks4.Select(a => a.Result.row_feature_values_csv).ToList();
            features_lines.Insert(0, row_feature_header_csv);
            io_proxy.WriteAllLines(fns.fn_features, features_lines, module_name, method_name);
        }

        private static List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> part3_encode_features(cmd_params p, List<subsequence_classification_data> class_data_list)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part3_encode_features);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine($@"{p.class_id} {p.class_name}: 3. Encoding data for {class_data_list.Count} items...", module_name, method_name);
            }

            var tasks3 =
                new List<Task<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)>>();
            var tasks3_start_time = DateTime.Now;

            for (var i = 0; i < class_data_list.Count; i++)
            {
                var class_data_item = class_data_list[i];

                class_data_list[i] = null;

                var task = Task.Run(() =>
                {
                    var encoded_class_data_item = (class_data_item, subsequence_classification_data_methods.encode_subsequence_classification_data_row(p, class_data_item, p.max_features));

                    return encoded_class_data_item;
                });

                tasks3.Add(task);

                wait_tasks(tasks3.ToArray<Task>(), tasks3_start_time, module_name, method_name);
            }

            class_data_list.Clear();
            class_data_list = null;
            wait_tasks(tasks3.ToArray<Task>(), tasks3_start_time, module_name, method_name);
            var data_encoded_list = tasks3.Select(a => a.Result).ToList();
            return data_encoded_list;
        }

        private static List<subsequence_classification_data> part2_load_data(cmd_params p, List<protein_subsequence_info> psi_list, protein_subsequence_info _template_protein)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part2_load_data);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);



            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: 2. Loading available data for {psi_list.Count} items...",
                    nameof(program), nameof(part2_load_data));
            }

            var tasks2 = new List<Task<subsequence_classification_data>>();
            var tasks2_start_time = DateTime.Now;

            for (var i = 0; i < psi_list.Count; i++)
            {
                var psi = psi_list[i];

                psi_list[i] = null;

                var task = Task.Run(() =>
                {
                    var classification_data = get_subsequence_classificiation_data(p, psi, _template_protein);

                    return classification_data;
                });

                tasks2.Add(task);

                wait_tasks(tasks2.ToArray<Task>(), tasks2_start_time, module_name, method_name);
            }

            psi_list.Clear();
            psi_list = null;
            wait_tasks(tasks2.ToArray<Task>(), tasks2_start_time, module_name, method_name);
            var class_data_list = tasks2.Select(a => a.Result).ToList();
            return class_data_list;
        }

        private static List<protein_subsequence_info> part1_find_samples(cmd_params p, List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> pdb_id_list)
        {
            const string module_name = nameof(program);
            const string method_name = nameof(part1_find_samples);

            io_proxy.WriteLine($@"running {method_name}()...", module_name, method_name);


            const string standard_coil = @"standard_coil";
            const string dimorphic_coil = @"dimorphic_coil";

            var pdb_id_list2 = pdb_id_list;

            if (string.Equals(p.class_name, standard_coil, StringComparison.OrdinalIgnoreCase))
            {
                pdb_id_list2 = pdb_id_list.GroupBy(a => (a.pdb_id, a.chain_number)).Select(a => a.First()).ToList();
            }

            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine($@"{p.class_id} {p.class_name}: 1. Loading pdb info from {pdb_id_list2.Count} files...", module_name, method_name);
            }

            //var cache_filename = Path.Combine(cmd_params.output_folder, $@"l__({cmd_params.class_name}).csv");

            var cache_filename = get_input_filenames(p, null);

            var cache = protein_subsequence_info.load(cache_filename);

            if (cache != null && cache.Count > 0)
            {
                return cache;
            }

            List<(string pdb_id, char chain_id, List<int> res_ids)> invalid_res_ids = null;

            if (string.Equals(p.class_name, standard_coil, StringComparison.OrdinalIgnoreCase))
            {
                var p2 = new cmd_params(p)
                {
                    class_id = -p.class_id,
                    class_name = dimorphic_coil,
                    
                };

                var pdb_id_list3 = pdb_id_list.Select(a => (a.pdb_id, a.dimer_type, a.class_name, a.symmetry_mode, a.parallelism, a.chain_number, a.strand_seq, a.optional_res_index)).ToList();
                var dimorphic_coils_psi_list = part1_find_samples(p2, pdb_id_list3);

                invalid_res_ids = dimorphic_coils_psi_list.Select(a => (a.pdb_id, a.chain_id, a.res_ids.Select(b => b.res_id).ToList())).ToList();
                invalid_res_ids = invalid_res_ids.GroupBy(a => (a.pdb_id, a.chain_id)).Select(a =>
                    (a.Key.pdb_id,
                        a.Key.chain_id,
                        a.SelectMany(b => b.res_ids).Distinct().ToList())).ToList();
            }

            var tasks1 = new List<Task<List<protein_subsequence_info>>>();
            var tasks1_start_time = DateTime.Now;
            for (var i = 0; i < pdb_id_list2.Count; i++)
            {
                var pdb_id_item = pdb_id_list2[i];

                pdb_id_list2[i] = default;

                var task = Task.Run(() =>
                {
                    var psi_list = new List<protein_subsequence_info>();

                    if (string.Equals(p.class_name, standard_coil, StringComparison.OrdinalIgnoreCase))
                    {
                        psi_list = dataset_gen_coils.find_coils(pdb_id_item.dimer_type, pdb_id_item.pdb_id, pdb_id_item.chain_number, p.class_id, p.class_name, p.use_dssp3);
                    }
                    else if (string.Equals(p.class_name, dimorphic_coil, StringComparison.OrdinalIgnoreCase))
                    {
                        psi_list.Add(dataset_gen_dimorphic.get_dhc_item(p.class_id, p.class_name, p.use_dssp3, true, false, pdb_id_item));
                    }
                    else
                    {
                        throw new Exception();
                    }

                    if (invalid_res_ids != null && invalid_res_ids.Count > 0)
                    {
                        psi_list = psi_list.Where(pl => !invalid_res_ids.Any(iv =>
                            iv.pdb_id == pl.pdb_id &&
                            iv.chain_id == pl.chain_id &&
                            pl.res_ids.Any(c => iv.res_ids.Contains(c.res_id))
                        )).ToList();
                    }

                    return psi_list;
                });

                tasks1.Add(task);

                //wait_tasks(tasks1.ToArray<Task>(), nameof(program), nameof(part1));
            }

            pdb_id_list2?.Clear();
            pdb_id_list2 = null;

            pdb_id_list?.Clear();
            pdb_id_list = null;

            wait_tasks(tasks1.ToArray<Task>(), tasks1_start_time, module_name, method_name);
            var psi_list2 = tasks1.SelectMany(a => a.Result).ToList();

            psi_list2 = psi_list2.Where(a => a.aa_subsequence.Length >= p.min_sequence_length).ToList();

            protein_subsequence_info.save(cache_filename, psi_list2);

            return psi_list2;
        }
    }
}

