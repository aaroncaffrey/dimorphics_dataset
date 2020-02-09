using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;
using System.Runtime.Loader;
using System.Security.Cryptography;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class program
    {
        // todo: re-check tortusoity code, intramolecular distance code, and averageatom-atom distance cpde

        // todo: check which blast databases, output options are disabled.
        // todo: check blast distance option.

        internal static bool verbose = true;

        internal static string data_root_folder = $@"C:\betastrands_dataset\";

        public static List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> get_dataset_pdb_id_list(string class_name_in_file = "Single")
        {
            //List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)>

            var dimorphics_data1 = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, $@"csv", $@"distinct dimorphics list.csv"), nameof(program), nameof(get_dataset_pdb_id_list))
                .Skip(1)
                .Where(a => !string.IsNullOrWhiteSpace(a.Replace(",", "", StringComparison.InvariantCulture)))
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
                        chain_number: int.Parse(x[k++], NumberStyles.Integer, CultureInfo.InvariantCulture) - 1,
                        strand_seq: x[k++],
                        optional_res_index: x[k++]
                    );
                }).ToList();

            if (!string.IsNullOrWhiteSpace(class_name_in_file))
            {
                dimorphics_data1 = dimorphics_data1.Where(a => a.class_name == class_name_in_file).ToList();
            }

            //var pdb_id_list = dimorphics_data1.Select(a => (dimer_type: a.dimer_type, pdb_id: a.pdb_id)).Distinct().ToList();
            //
            //return pdb_id_list;

            return dimorphics_data1;
        }

        public static subsequence_classification_data get_subsequence_classificiation_data(protein_subsequence_info psi, feature_types feature_types, protein_subsequence_info _template_protein = null, atom.load_atoms_pdb_options pdb_options = null, bool load_ss_predictions = true, bool load_foldx_energy = true)
        {
            // todo: check if pdb_folder should be pdb or pdb_repair?

            if (psi == null)
            {
                throw new ArgumentNullException(nameof(psi));
            }

            if (feature_types == null)
            {
                throw new ArgumentNullException(nameof(feature_types));
            }

            //var features_1d = (feature_types?.feature_types_interface_2d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_neighbourhood_2d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_chain_2d?.AsArray()?.Any(a => a.value) ?? false);

            //var features_3d = (feature_types?.feature_types_interface_3d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value) ?? false) ||
            //                  (feature_types?.feature_types_chain_3d?.AsArray()?.Any(a => a.value) ?? false);

            var pdb_atoms = atom.load_atoms_pdb
                (
                    pdb_id: psi.pdb_id,
                    pdb_options ??
                    new atom.load_atoms_pdb_options()
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


                        load_2d_mpsa_sec_struct_predictions = true,
                        load_2d_blast_pssms = true,
                        load_2d_iup_data = true,
                        load_2d_sable = true,
                        load_2d_dna_binding = true,

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
                dimer_type = psi.dimer_type ?? "",
                symmetry_mode = psi.symmetry_mode ?? "",
                parallelism = psi.parallelism ?? "",

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
                subsequence_classification_data_templates._template_scd = get_subsequence_classificiation_data(_template_protein, feature_types, null, pdb_options, load_ss_predictions, load_foldx_energy);
            }

            return scd;
        }

        public static void wait_tasks(Task[] tasks, string module_name = "", string function_name = "")
        {
            if (tasks == null || tasks.Length == 0) return;

            //var sw1 = new Stopwatch();
            //sw1.Start();

            var start_time = DateTime.Now;

            do
            {
                var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToList();

                if (incomplete_tasks.Count > 0)
                {
                    var completed_task_index = Task.WaitAny(incomplete_tasks.ToArray<Task>());

                    if (completed_task_index > -1)
                    {
                        var incomplete = tasks.Count(a => !a.IsCompleted);
                        var complete = tasks.Length - incomplete;
                        var pct = ((double)complete / (double)tasks.Length) * 100;

                        var time_remaining = TimeSpan.FromTicks((long)(DateTime.Now.Subtract(start_time).Ticks * ((double)incomplete / (double)(complete == 0 ? 1 : complete))));

                        io_proxy.WriteLine($@"{module_name}.{function_name} -> {complete} / {tasks.Length} ( {pct:0.00} % ) [ {time_remaining:dd\:hh\:mm\:ss\.fff} ]",
                            nameof(program), nameof(Main));
                        //io_proxy.WriteLine($@"{module_name}.{function_name} -> {complete} / {tasks.Length} ( {pct:0.00} % )  Elapsed: {el_ts:dd\:hh\:mm\:ss\.fff}  ETA {togo_ts:dd\:hh\:mm\:ss\.fff}", nameof(program), nameof(Main));
                    }
                }

                //GC.Collect();
            } while (tasks.Any(a => !a.IsCompleted));

            Task.WaitAll(tasks.ToArray<Task>());
        }

        public static void close_notifications()
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
        //var at = pdbs.SelectMany(a => atom.load_atoms_pdb(a, new atom.load_atoms_pdb_options()
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
        //    var ppp = atom.load_atoms_pdb(pdb, new atom.load_atoms_pdb_options()
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





        public static void Main(string[] args)
        {

            //var swx = new Stopwatch();
            //swx.Start();
            //var x = cache_r_load();
            //swx.Stop();
            //Console.WriteLine(swx.Elapsed.ToString("g"));

            //cache_r_servers();
            //cache_r_servers_check();
            //return;

            close_notifications();
            var child_priority = ProcessPriorityClass.Idle;
            var child_priority_boost = false;
            var p = cmd_params.get_params(args);
            
            if (p == null)
            {
                return;
            }

            program.verbose = p?.verbose.Value ?? true;

            var feature_types2 = feature_types.feature_types_params(p?.area ?? null);

            //return;
            var sw1 = new Stopwatch();
            sw1.Start();

            var pdb_id_list = get_dataset_pdb_id_list();

            // 1. find subsequence details from analysis of pdb files
            var psi_list = part1_find_samples(p, pdb_id_list);

            // load r cache, if not using children (takes too many resources to load cache for individual items)
            if (!p.use_children.Value)
            {
                subsequence_classification_data_r_methods.cache_r_protr = subsequence_classification_data_r_methods.cache_r_load(@"e:\dataset\r_protr_cache.csv");
                subsequence_classification_data_r_methods.cache_r_peptides = subsequence_classification_data_r_methods.cache_r_load(@"e:\dataset\r_peptides_cache.csv");
            }

            var _template_protein = psi_list.First();
            //Parallel.For(0, psi_list.Count, index =>

            var tasks = new List<Task>();

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

                var fns = get_output_filenames(p, feature_types2, _index);
                if (File.Exists(fns.fn_headers) && File.Exists(fns.fn_features) && File.Exists(fns.fn_comments)) continue;


                if (p.first_index == null && p.last_index == null)
                {
                    if (p.use_children.Value)
                    {
                        var task = Task.Run(() =>
                        {
                            while (true)
                            {
                                using var process = Process.Start(Process.GetCurrentProcess().MainModule.FileName, $"{string.Join(" ", args)} -first_index={_index} -last_index={_index}");

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

                        var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToArray();
                        if (incomplete_tasks.Length >= Environment.ProcessorCount)
                        {
                            Task.WaitAny(incomplete_tasks);
                        }

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
                var class_data_list = part2_load_data(p, psi_list2, feature_types2, _template_protein);

                // 3. encode the data as svm classification features
                var data_encoded_list = part3_encode_features(p, class_data_list, feature_types2);

                // 4. Save to file
                part4_save_outputs(feature_types2, p, data_encoded_list, _index);
            }

            wait_tasks(tasks.ToArray<Task>(), nameof(program), nameof(Main));

            if (p.first_index == null && p.last_index == null)
            {
                part5_check_output_hashes(psi_list, feature_types2, p);

                part6_check_and_merge_outputs(psi_list, feature_types2, p);
            }

            sw1.Stop();

            if (p.first_index == null && p.last_index == null)
            {

                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: Finished: ({sw1.Elapsed:dd\:hh\:mm\:ss\.fff})",
                    nameof(program), nameof(Main));
            }
        }
         
        private static void part6_check_and_merge_outputs(List<protein_subsequence_info> psi_list, feature_types feature_types, cmd_params p)
        {
            var files = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(p, feature_types, i);
                return (fn_headers: b.fn_headers, fn_features: b.fn_features, fn_comments: b.fn_comments);
            }).ToList();

            var text_header = io_proxy.ReadAllLines(files.First().fn_headers, nameof(program), nameof(part6_check_and_merge_outputs));
            var text_features = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_features, nameof(program), nameof(part6_check_and_merge_outputs)).Skip(i == 0 ? 0 : 1).ToList()).ToList();
            var text_comments = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_comments, nameof(program), nameof(part6_check_and_merge_outputs)).Skip(i == 0 ? 0 : 1).ToList()).ToList();

            // check for any lines with different total number of features
            var num_features = text_features.Select((a, i) => (line: i, num_features: a.Count(b => b == ','))).OrderByDescending(a => a.num_features).ToList();

            if (num_features.Select(a => a.num_features).Distinct().Count() > 1)
            {
                io_proxy.WriteLine($@"Feature counts don't match.", nameof(program), nameof(part6_check_and_merge_outputs));
            }
            else
            {
                io_proxy.WriteLine($@"Feature counts match.", nameof(program), nameof(part6_check_and_merge_outputs));
            }

            // check for any non-double type features
            var invalid_features = text_features
                .SelectMany(a =>
                    a.Split(',').Select((b, i) => (string.IsNullOrEmpty(b) || b == "∞" || b == "+∞" || b == "-∞" || b == "NaN") ? i : -1)
                        .Where(a => a != -1).ToList()).Distinct().OrderBy(a => a).ToList();

            if (invalid_features.Count > 0)
            {
                io_proxy.WriteLine($@"Invalid feature ids: " + string.Join(", ", invalid_features), nameof(program), nameof(part6_check_and_merge_outputs));

                foreach (var i in invalid_features)
                {
                    var invalid_feature_header =
                        text_header.FirstOrDefault(a => a.StartsWith($"{i},", StringComparison.InvariantCulture));
                    io_proxy.WriteLine($@"{i}: {invalid_feature_header}");
                }
            }
            else
            {
                io_proxy.WriteLine($@"No invalid feature ids.", nameof(program), nameof(part6_check_and_merge_outputs));
            }

            var fns = get_output_filenames(p, feature_types, null);
            io_proxy.WriteAllLines(fns.fn_headers, text_header);
            io_proxy.WriteAllLines(fns.fn_features, text_features);
            io_proxy.WriteAllLines(fns.fn_comments, text_comments);
        }

        private static void part5_check_output_hashes(List<protein_subsequence_info> psi_list, feature_types feature_types,cmd_params p)
        {
            var file_hashes = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(p, feature_types, i);
                return (fn_headers: file_sha_256(b.fn_headers), fn_features: file_sha_256(b.fn_features), fn_comments: file_sha_256(b.fn_comments));
            }).ToList();

            var file_hashes_headers = file_hashes.Select(a => a.fn_headers).ToList();
            var file_hashes_features = file_hashes.Select(a => a.fn_features).ToList();
            var file_hashes_comments = file_hashes.Select(a => a.fn_comments).ToList();

            var file_hashes_headers_cnt = file_hashes_headers.Distinct()
                .Select(a => (hash: a, count: file_hashes_headers.Count(b => a == b)))
                .OrderByDescending(a => a.count).ToList();

            var file_hashes_features_cnt = file_hashes_features.Distinct()
                .Select(a => (hash: a, count: file_hashes_features.Count(b => a == b)))
                .OrderByDescending(a => a.count).ToList();

            var file_hashes_comments_cnt = file_hashes_comments.Distinct()
                .Select(a => (hash: a, count: file_hashes_comments.Count(b => a == b)))
                .OrderByDescending(a => a.count).ToList();

            if (file_hashes_headers_cnt.Count != 1)
            {
                io_proxy.WriteLine("Hashes for headers:");
                file_hashes_headers_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }

            if (file_hashes_features_cnt.Count != psi_list.Count)
            {
                io_proxy.WriteLine("Hashes for features:");
                file_hashes_features_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }

            if (file_hashes_comments_cnt.Count != psi_list.Count)
            {
                io_proxy.WriteLine("Hashes for comments:");
                file_hashes_comments_cnt.Where(a => a.count > 1).ToList().ForEach(a => io_proxy.WriteLine(a.ToString()));
            }
        }

        private static string file_sha_256(string filename)
        {
            using var stream = new BufferedStream(File.OpenRead(filename), 1024 * 1024 * 4);
            using var sha256 = new SHA256Managed();
            var hash_bytes = sha256.ComputeHash(stream);
            var hash_str = BitConverter.ToString(hash_bytes);//.Replace("-", String.Empty);
            return hash_str;
        }

       

        private static string get_input_filenames(cmd_params p, int? index = null)
        {
            var fn_input = Path.Combine(p.output_folder, $"l__({(p.class_id > 0 ? "+" : "")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            return fn_input;
        }

        private static (string fn_headers, string fn_comments, string fn_features) get_output_filenames(cmd_params p, feature_types feature_types, int? index)
        {
            if (p==null) throw new ArgumentNullException(nameof(p));
            if (feature_types == null) throw new ArgumentNullException(nameof(feature_types));

            var tag = feature_types.get_output_file_tag();

            //var fn_input = Path.Combine(cmd_params.output_folder, $"l__({cmd_params.class_name}){(tag != null ? $@"_{tag}" : $@"")}.csv");
            var fn_headers = Path.Combine(p.output_folder, $"h_{(tag ?? "")}_({(p.class_id > 0 ? "+" : "")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            var fn_comments = Path.Combine(p.output_folder, $"c_{(tag ?? "")}_({(p.class_id > 0 ? "+" : "")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");
            var fn_features = Path.Combine(p.output_folder, $"f_{(tag ?? "")}_({(p.class_id > 0 ? "+" : "")}{p.class_id})_({p.class_name}){(index != null ? $@"_{index}" : $@"")}.csv");

            return (fn_headers, fn_comments, fn_features);
        }

        private static void part4_save_outputs(
            feature_types feature_types, cmd_params p, List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> data_encoded_list, int tag)
        {
            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: 4. Saving encoded data to file for {data_encoded_list.Count} items...",
                    nameof(program), nameof(part4_save_outputs));
            }

            var fns = get_output_filenames(p, feature_types, tag);

            // get header row indexes in csv format
            var row_feature_header_csv = string.Join(",", Enumerable.Range(0, data_encoded_list.First().feature_info.Count));

            // get comments file header in csv format
            var row_comments_header_csv = string.Join(",", subsequence_classification_data.get_row_comments_headers(data_encoded_list.First().instance_meta_data));

            // get list of the feature headers in csv format
            // save headers
            var feature_headers = feature_info.get_feature_headers_lines_csv(data_encoded_list.First().feature_info);
            io_proxy.WriteAllLines(fns.fn_headers, feature_headers, nameof(program), nameof(part4_save_outputs));


            var tasks4 = new List<Task<(string row_feature_values_csv, string row_comments_csv)>>();

            for (var i = 0; i < data_encoded_list.Count; i++)
            {
                var row_index = i;
                var data_encoded_item = data_encoded_list[i];

                data_encoded_list[i] = default;

                var task = Task.Run(() =>
                {
                    // get feature values
                    var row_feature_values = data_encoded_item.feature_info.Select((a, fid) =>
                        a.feature_value.ToString("G17", CultureInfo.InvariantCulture)).ToList();

                    // convert feature values to csv format
                    var row_feature_values_csv = string.Join(",", row_feature_values);

                    // get meta data about the example instance
                    var row_comments = subsequence_classification_data.get_row_comments(row_index, data_encoded_item.instance_meta_data);

                    var row_comments_csv = string.Join(",", row_comments);

                    return (row_feature_values_csv, row_comments_csv);
                });

                tasks4.Add(task);

                wait_tasks(tasks4.ToArray<Task>(), nameof(program), nameof(part4_save_outputs));
            }

            data_encoded_list.Clear();
            data_encoded_list = null;
            wait_tasks(tasks4.ToArray<Task>(), nameof(program), nameof(part4_save_outputs));

            // save comments
            var comments_lines = tasks4.Select(a => a.Result.row_comments_csv).ToList();
            comments_lines.Insert(0, row_comments_header_csv);
            io_proxy.WriteAllLines(fns.fn_comments, comments_lines, nameof(program), nameof(part4_save_outputs));

            // save features
            var features_lines = tasks4.Select(a => a.Result.row_feature_values_csv).ToList();
            features_lines.Insert(0, row_feature_header_csv);
            io_proxy.WriteAllLines(fns.fn_features, features_lines, nameof(program), nameof(part4_save_outputs));
        }

        private static List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> part3_encode_features(
            cmd_params p, List<subsequence_classification_data> class_data_list, feature_types feature_types)
        {
            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: 3. Encoding data for {class_data_list.Count} items...",
                    nameof(program), nameof(part3_encode_features));
            }

            var tasks3 =
                new List<Task<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)>>();

            for (var i = 0; i < class_data_list.Count; i++)
            {
                var class_data_item = class_data_list[i];

                class_data_list[i] = null;

                var task = Task.Run(() =>
                {
                    var encoded_class_data_item = (class_data_item, subsequence_classification_data_methods.encode_subsequence_classification_data_row(class_data_item, p.max_features, feature_types));

                    return encoded_class_data_item;
                });

                tasks3.Add(task);

                wait_tasks(tasks3.ToArray<Task>(), nameof(program), nameof(part3_encode_features));
            }

            class_data_list.Clear();
            class_data_list = null;
            wait_tasks(tasks3.ToArray<Task>(), nameof(program), nameof(part3_encode_features));
            var data_encoded_list = tasks3.Select(a => a.Result).ToList();
            return data_encoded_list;
        }

        private static List<subsequence_classification_data> part2_load_data(cmd_params p, List<protein_subsequence_info> psi_list, feature_types feature_types, protein_subsequence_info _template_protein)
        {
            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: 2. Loading available data for {psi_list.Count} items...",
                    nameof(program), nameof(part2_load_data));
            }

            var tasks2 = new List<Task<subsequence_classification_data>>();

            for (var i = 0; i < psi_list.Count; i++)
            {
                var psi = psi_list[i];

                psi_list[i] = null;

                var task = Task.Run(() =>
                {
                    var classification_data = get_subsequence_classificiation_data(psi, feature_types, _template_protein);

                    return classification_data;
                });

                tasks2.Add(task);

                wait_tasks(tasks2.ToArray<Task>(), nameof(program), nameof(part2_load_data));
            }

            psi_list.Clear();
            psi_list = null;
            wait_tasks(tasks2.ToArray<Task>(), nameof(program), nameof(part2_load_data));
            var class_data_list = tasks2.Select(a => a.Result).ToList();
            return class_data_list;
        }

        private static List<protein_subsequence_info> part1_find_samples(cmd_params p, List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> pdb_id_list)
        {
            const string standard_coil = "standard_coil";
            const string dimorphic_coil = "dimorphic_coil";

            var pdb_id_list2 = pdb_id_list;

            if (string.Equals(p.class_name, standard_coil, StringComparison.InvariantCultureIgnoreCase))
            {
                pdb_id_list2 = pdb_id_list.GroupBy(a => (a.pdb_id, a.chain_number)).Select(a => a.First()).ToList();
            }

            if (p.first_index == null && p.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{p.class_id} {p.class_name}: 1. Loading pdb info from {pdb_id_list2.Count} files...",
                    nameof(program), nameof(part1_find_samples));
            }

            //var cache_filename = Path.Combine(cmd_params.output_folder, $@"l__({cmd_params.class_name}).csv");

            var cache_filename = get_input_filenames(p, null);

            var cache = protein_subsequence_info.load(cache_filename);

            if (cache != null && cache.Count > 0)
            {
                return cache;
            }

            List<(string pdb_id, char chain_id, List<int> res_ids)> invalid_res_ids = null;

            if (string.Equals(p.class_name, standard_coil, StringComparison.InvariantCultureIgnoreCase))
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

            for (var i = 0; i < pdb_id_list2.Count; i++)
            {
                var pdb_id_item = pdb_id_list2[i];

                pdb_id_list2[i] = default;

                var task = Task.Run(() =>
                {
                    var psi_list = new List<protein_subsequence_info>();

                    if (string.Equals(p.class_name, standard_coil, StringComparison.InvariantCultureIgnoreCase))
                    {
                        psi_list = dataset_gen_coils.find_coils(pdb_id_item.dimer_type, pdb_id_item.pdb_id, pdb_id_item.chain_number, p.class_id, p.class_name, p.use_dssp3);
                    }
                    else if (string.Equals(p.class_name, dimorphic_coil, StringComparison.InvariantCultureIgnoreCase))
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

            wait_tasks(tasks1.ToArray<Task>(), nameof(program), nameof(part1_find_samples));
            var psi_list2 = tasks1.SelectMany(a => a.Result).ToList();

            psi_list2 = psi_list2.Where(a => a.aa_subsequence.Length >= p.min_sequence_length).ToList();

            protein_subsequence_info.save(cache_filename, psi_list2);

            return psi_list2;
        }
    }
}

