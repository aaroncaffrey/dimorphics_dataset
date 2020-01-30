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

        internal static bool verbose;

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

        public static subsequence_classification_data get_subsequence_classificiation_data(protein_subsequence_info psi, feature_types feature_types, protein_subsequence_info _template_protein)
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

            var features_1d = (feature_types?.feature_types_interface_2d?.AsArray()?.Any(a => a.value) ?? false) ||
                              (feature_types?.feature_types_neighbourhood_2d?.AsArray()?.Any(a => a.value) ?? false) ||
                              (feature_types?.feature_types_chain_2d?.AsArray()?.Any(a => a.value) ?? false);

            var features_3d = (feature_types?.feature_types_interface_3d?.AsArray()?.Any(a => a.value) ?? false) ||
                              (feature_types?.feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value) ?? false) ||
                              (feature_types?.feature_types_chain_3d?.AsArray()?.Any(a => a.value) ?? false);

            var pdb_atoms = atom.load_atoms_pdb
                (
                    pdb_id: psi.pdb_id,
                    new atom.load_atoms_pdb_options()
                    {
                        load_2d_mpsa_sec_struct_predictions = features_1d || features_3d,
                        load_2d_blast_pssms = features_1d || features_3d,
                        load_2d_iup_data = features_1d || features_3d,
                        load_2d_sable = features_1d || features_3d,
                        load_2d_dna_binding = features_1d || features_3d,

                        find_3d_intramolecular = features_1d || features_3d,
                        load_3d_dssp_data = features_1d || features_3d,
                        load_3d_stride_data = features_1d || features_3d,
                        load_3d_ring_data = features_1d || features_3d,
                        load_3d_foldx_ala_scan = features_1d || features_3d,
                        load_3d_rsa_data = features_1d || features_3d,
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
            
            scd.interface_region = new subsequence_classification_data_region(scd, interface_atoms);
            scd.chain_region = new subsequence_classification_data_region(scd, chain_atoms);
            
            // after the interface and chain properties are set, can then find neighborhood...
            scd.init_nh_flanking();
            scd.init_nh_contacts();

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
                subsequence_classification_data_templates._template_scd = get_subsequence_classificiation_data(_template_protein, feature_types, null);
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

        public static string get_param(string name, List<string> args)
        {
            if (args == null)
            {
                return null;
            }

            var ix = args.FindIndex(a => string.Equals(a, $"-{name}", StringComparison.InvariantCultureIgnoreCase));
            
            if (ix == -1) return null;

            var value = (ix == args.Count - 1) ? "" : args[ix + 1];

            if (ix < args.Count - 1) args.RemoveAt(ix + 1);
            
            args.RemoveAt(ix);
            
            return value;
        }

        
        public static (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string output_folder, int? first_index, int? last_index) get_params(string[] args)
        {
            if (args == null || args.Length == 0)
            {
                var exe = Path.GetFileName(Process.GetCurrentProcess().MainModule.FileName);
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} usage:");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"  {exe} -area=[2i,2n,2p,3i,3n,3p] -use_dssp3=[true|false] -class_id=[-1|+1] -class_name=[class_name] -min_sequence_length=[3] -max_features=[100] -output_folder=[path] [-first_index=[i] -last_index=[j]]");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} examples:");
                io_proxy.WriteLine($@"");

                foreach (var a in new [] {"2i", "2n", "2p", "3i", "3n", "3p"})
                {
                    io_proxy.WriteLine($@"{a[0]}d {(a[1]=='i'?"interface subsequence":"")}{(a[1]=='n'?"neighbourhood":"")}{(a[1]=='p'?"protein":"")} area:");
                    io_proxy.WriteLine($@"  {exe} -area={a} -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{a}\dimorphic_coil\ 1> e:\dataset\{a}\dimorphic_coil\stdout.txt 2> e:\dataset\{a}\dimorphic_coil\stderr.txt");
                    io_proxy.WriteLine($@"  {exe} -area={a} -use_dssp3=true -class_id=-1 -class_name=standard_coil  -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{a}\standard_coil\ 1> e:\dataset\{a}\standard_coil\stdout.txt 2> e:\dataset\{a}\standard_coil\stderr.txt");
                    io_proxy.WriteLine($@"");
                }

                Environment.Exit(0);
                //return default;
                //throw new Exception($@"No {nameof(args)} specified.");
            }

            io_proxy.WriteLine(string.Join(" ", args), nameof(program), nameof(get_params));

            var args_list = args.SelectMany(a => a.Split(new char[] {'='}, StringSplitOptions.RemoveEmptyEntries)).ToList();

            (string[] area, bool? use_dssp3, int? class_id, string class_name, int? min_sequence_length, int? max_features, string output_folder, int? first_index, int? last_index) x;

            x = (
                get_param("area", args_list)?.ToLowerInvariant().Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries).OrderBy(a=>a).ToArray() ?? null,
                bool.Parse(get_param("use_dssp3", args_list)),
                int.Parse(get_param("class_id", args_list), NumberStyles.Integer, CultureInfo.InvariantCulture),
                get_param("class_name", args_list)?.ToLowerInvariant(),
                int.Parse(get_param("min_sequence_length", args_list), NumberStyles.Integer, CultureInfo.InvariantCulture),
                int.Parse(get_param("max_features", args_list), NumberStyles.Integer, CultureInfo.InvariantCulture),
                get_param("output_folder", args_list),
                int.TryParse(get_param("first_index", args_list), NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_first_index) ? tp_first_index : (int?)null,
                int.TryParse(get_param("last_index", args_list), NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_last_index) ? tp_last_index : (int?)null
            );
            
            if (args_list.Any())
            {
                io_proxy.WriteLine($@"Args unknown: {string.Join(" ", args_list)}", nameof(program), nameof(get_params));

                Environment.Exit(0);
            }

            if (x.area == null || x.area.Any(a => string.IsNullOrWhiteSpace(a))) throw new ArgumentNullException(nameof(x.area));
            if (x.area.Except(new string[] { "2i", "2n", "2p", "3i", "3n", "3p" }).Any()) throw new ArgumentOutOfRangeException(nameof(x.area));
            if (x.use_dssp3 == null) throw new ArgumentNullException(nameof(x.use_dssp3));
            if (x.class_id == null) throw new ArgumentNullException(nameof(x.class_id));
            if (string.IsNullOrWhiteSpace(x.class_name)) throw new ArgumentNullException(nameof(x.class_name));
            if (x.min_sequence_length == null) throw new ArgumentNullException(nameof(x.min_sequence_length));
            if (x.max_features == null) throw new ArgumentNullException(nameof(x.max_features));
            if (string.IsNullOrWhiteSpace(x.output_folder)) throw new ArgumentNullException(nameof(x.output_folder));

            var ret = (
                area:x.area, 
                use_dssp3:x.use_dssp3.Value, 
                class_id:x.class_id.Value,
                class_name:x.class_name,
                min_sequence_length:x.min_sequence_length.Value,
                max_features:x.max_features.Value, 
                output_folder:x.output_folder,
                first_index:x.first_index, 
                last_index:x.last_index ?? x.first_index
                           );

            if (x.first_index == null && x.last_index == null)
            {
                io_proxy.WriteLine($@"{nameof(ret.area)} = ""{string.Join(", ", ret.area)}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.use_dssp3)} = ""{ret.use_dssp3}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.class_id)} = ""{ret.class_id}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.class_name)} = ""{ret.class_name}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.min_sequence_length)} = ""{ret.min_sequence_length}""",
                    nameof(program), nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.max_features)} = ""{ret.max_features}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.output_folder)} = ""{ret.output_folder}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.first_index)} = ""{ret.first_index}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.last_index)} = ""{ret.last_index}""", nameof(program),
                    nameof(get_params));
            }

            return ret;
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
            close_notifications();

            var child_priority = ProcessPriorityClass.Idle;

            var cmd_params = get_params(args);

            program.verbose = (cmd_params.first_index == null && cmd_params.last_index == null);
            
            var feature_types = feature_types_params(cmd_params);

            //return;
            var sw1 = new Stopwatch();
            sw1.Start();

            var pdb_id_list = get_dataset_pdb_id_list();
            
            // 1. find subsequence details from analysis of pdb files
            var psi_list = part1(cmd_params, pdb_id_list);

            var _template_protein = psi_list.First();
            //Parallel.For(0, psi_list.Count, index =>

            var tasks = new List<Task>();

            for (var index = 0; index < psi_list.Count; index++)
            {
                var _index = index;

                if (cmd_params.first_index != null && cmd_params.first_index > -1)
                {
                    if (_index < cmd_params.first_index) continue;
                }

                if (cmd_params.last_index != null && cmd_params.last_index > -1)
                {
                    if (_index > cmd_params.last_index) continue;
                }

                var fns = get_output_filenames(cmd_params, _index);
                if (File.Exists(fns.fn_headers) && File.Exists(fns.fn_features) && File.Exists(fns.fn_comments)) continue;

                if (cmd_params.first_index == null && cmd_params.last_index == null)
                {
                    var task = Task.Run(() =>
                    {
                        while (true)
                        {
                            using var process = Process.Start(Process.GetCurrentProcess().MainModule.FileName, $"{string.Join(" ", args)} -first_index={_index} -last_index={_index}");

                            if (process != null)
                            {
                                try
                                {
                                    process.PriorityClass = child_priority;
                                }
                                catch (Exception)
                                {

                                }
                                
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

                
                

               
                var psi_list2 = new List<protein_subsequence_info>() { psi_list[_index] };

                // 2. load available data sources
                var class_data_list = part2(cmd_params, psi_list2, feature_types, _template_protein);

                // 3. encode the data as svm classification features
                var data_encoded_list = part3(cmd_params, class_data_list, feature_types);

                // 4. Save to file
                part4(cmd_params, data_encoded_list, _index);
            }

            wait_tasks(tasks.ToArray<Task>(), nameof(program), nameof(Main));

            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {
                part5(psi_list, cmd_params);

                part6(psi_list, cmd_params);
            }

            sw1.Stop();
            
            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {

                io_proxy.WriteLine(
                    $@"{cmd_params.class_id} {cmd_params.class_name}: Finished: ({sw1.Elapsed:dd\:hh\:mm\:ss\.fff})",
                    nameof(program), nameof(Main));
            }
        }

        private static void part6(List<protein_subsequence_info> psi_list,
          (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
              output_folder, int? first_index, int? last_index) cmd_params)
        {
            var files = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(cmd_params, i);
                return (fn_headers: b.fn_headers, fn_features: b.fn_features, fn_comments: b.fn_comments);
            }).ToList();

            var text_header = io_proxy.ReadAllLines(files.First().fn_headers, nameof(program), nameof(part6));
            var text_features = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_features, nameof(program), nameof(part6)).Skip(i == 0 ? 0 : 1).ToList()).ToList();
            var text_comments = files.SelectMany((a, i) => io_proxy.ReadAllLines(a.fn_comments, nameof(program), nameof(part6)).Skip(i == 0 ? 0 : 1).ToList()).ToList();

            // check for any lines with different total number of features
            var num_features = text_features.Select((a, i) => (line: i, num_features: a.Count(b => b == ','))).OrderByDescending(a => a.num_features).ToList();

            if (num_features.Select(a => a.num_features).Distinct().Count() > 1)
            {
                io_proxy.WriteLine($@"Feature counts don't match.", nameof(program), nameof(part6));
            }
            else
            {
                io_proxy.WriteLine($@"Feature counts match.", nameof(program), nameof(part6));
            }

            // check for any non-double type features
            var invalid_features = text_features
                .SelectMany(a =>
                    a.Split(',').Select((b, i) => (string.IsNullOrEmpty(b) || b == "∞" || b == "+∞" || b == "-∞" || b == "NaN") ? i : -1)
                        .Where(a => a != -1).ToList()).Distinct().OrderBy(a => a).ToList();

            if (invalid_features.Count > 0)
            {
                io_proxy.WriteLine($@"Invalid feature ids: " + string.Join(", ", invalid_features), nameof(program), nameof(part6));

                foreach (var i in invalid_features)
                {
                    var invalid_feature_header =
                        text_header.FirstOrDefault(a => a.StartsWith($"{i},", StringComparison.InvariantCulture));
                    io_proxy.WriteLine($@"{i}: {invalid_feature_header}");
                }
            }
            else
            {
                io_proxy.WriteLine($@"No invalid feature ids.", nameof(program), nameof(part6));
            }

            var fns = get_output_filenames(cmd_params, null);
            io_proxy.WriteAllLines(fns.fn_headers, text_header);
            io_proxy.WriteAllLines(fns.fn_features, text_features);
            io_proxy.WriteAllLines(fns.fn_comments, text_comments);
        }

        private static void part5(List<protein_subsequence_info> psi_list,
            (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
                output_folder, int? first_index, int? last_index) cmd_params)
        {
            var file_hashes = psi_list.AsParallel().AsOrdered().Select((a, i) =>
            {
                var b = get_output_filenames(cmd_params, i);
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

        private static feature_types feature_types_params(
            (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
                output_folder, int? first_index, int? last_index) cmd_params)
        {
            var do_2d_interface = cmd_params.area.Any(a => a == "2i");
            var do_2d_nh = cmd_params.area.Any(a => a == "2n");
            var do_2d_protein = cmd_params.area.Any(a => a == "2p");
            var do_3d_interface = cmd_params.area.Any(a => a == "3i");
            var do_3d_nh = cmd_params.area.Any(a => a == "3n");
            var do_3d_protein = cmd_params.area.Any(a => a == "3p");


            var feature_types = new feature_types()
            {
                feature_types_interface_2d = !do_2d_interface
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_2d_interface,
                        sequence_geometry_classification_data = do_2d_interface,
                        mpsa_classification_data_subsequence = do_2d_interface,
                        intrinsically_unordered_data = do_2d_interface,
                        aa_index_classification_data = do_2d_interface,
                        sable_classification_data = do_2d_interface,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_2d_interface,
                        r_peptides = do_2d_interface,
                        r_protr = do_2d_interface,
                    },
                feature_types_neighbourhood_2d = !do_2d_nh
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_2d_nh,
                        sequence_geometry_classification_data = do_2d_nh,
                        mpsa_classification_data_subsequence = do_2d_nh,
                        intrinsically_unordered_data = do_2d_nh,
                        aa_index_classification_data = do_2d_nh,
                        sable_classification_data = do_2d_nh,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_2d_nh,
                        r_peptides = do_2d_nh,
                        r_protr = do_2d_nh,
                    },
                feature_types_chain_2d = !do_2d_protein
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_2d_protein,
                        sequence_geometry_classification_data = do_2d_protein,
                        mpsa_classification_data_subsequence = do_2d_protein,
                        intrinsically_unordered_data = do_2d_protein,
                        aa_index_classification_data = do_2d_protein,
                        sable_classification_data = do_2d_protein,
                        dna_binding_prediction_data = do_2d_protein,
                        blast_pssm_subsequence_classification_data = do_2d_protein,
                        r_peptides = do_2d_protein,
                        r_protr = do_2d_protein,
                    },
                feature_types_interface_3d = !do_3d_interface
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_interface,
                        intramolecular_classification_data = do_3d_interface,

                        foldx_classification_data = do_3d_interface,
                        tortuosity_classification_data = do_3d_interface,
                        ring_classification_data = do_3d_interface,
                        pse_ssc_dssp_classification_data = do_3d_interface,
                    },
                feature_types_neighbourhood_3d = !do_3d_nh
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_nh,
                        intramolecular_classification_data = do_3d_nh,

                        foldx_classification_data = do_3d_nh,
                        tortuosity_classification_data = do_3d_nh,
                        ring_classification_data = do_3d_nh,
                        pse_ssc_dssp_classification_data = do_3d_nh,
                    },
                feature_types_chain_3d = !do_3d_protein
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_protein,
                        intramolecular_classification_data = do_3d_protein,

                        foldx_classification_data = do_3d_protein,
                        tortuosity_classification_data = do_3d_protein,
                        ring_classification_data = do_3d_protein,
                        pse_ssc_dssp_classification_data = do_3d_protein,
                    }
            };

            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_interface_2d)} = {feature_types.feature_types_interface_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_neighbourhood_2d)} = {feature_types.feature_types_neighbourhood_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_2d)} = {feature_types.feature_types_chain_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_interface_3d)} = {feature_types.feature_types_interface_3d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_3d)} = {feature_types.feature_types_chain_3d}",
                nameof(program), nameof(feature_types_params));
            return feature_types;
        }

        
        private static string get_input_filenames((string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string output_folder, int? first_index, int? last_index) cmd_params)
        {
            var fn_input = Path.Combine(cmd_params.output_folder, $"l__({cmd_params.class_name}).csv");

            return fn_input;
        }

        private static (string fn_headers, string fn_comments, string fn_features) get_output_filenames((string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string output_folder, int? first_index, int? last_index) cmd_params, int? tag)
        {
            
            var fn_headers = Path.Combine(cmd_params.output_folder, $"h__({cmd_params.class_name}){(tag != null ? $@"_{tag}" : $@"")}.csv");
            var fn_comments = Path.Combine(cmd_params.output_folder, $"c__({cmd_params.class_name}){(tag != null ? $@"_{tag}" : $@"")}.csv");
            var fn_features = Path.Combine(cmd_params.output_folder, $"f__({cmd_params.class_name}){(tag != null ? $@"_{tag}" : $@"")}.csv");

            return (fn_headers, fn_comments, fn_features);
        }

        private static void part4(
            (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
                output_folder, int? first_index, int? last_index) cmd_params, List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> data_encoded_list, int tag)
        {
            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{cmd_params.class_id} {cmd_params.class_name}: 4. Saving encoded data to file for {data_encoded_list.Count} items...",
                    nameof(program), nameof(part4));
            }

            var fns = get_output_filenames(cmd_params, tag);

            // get header row indexes in csv format
            var row_feature_header_csv = string.Join(",", Enumerable.Range(0, data_encoded_list.First().feature_info.Count));

            // get comments file header in csv format
            var row_comments_header_csv = string.Join(",", data_controller.get_row_comments_headers(data_encoded_list.First()));

            // get list of the feature headers in csv format
            var feature_headers = data_controller.get_feature_headers_lines_csv(data_encoded_list.First());
            
            io_proxy.WriteAllLines(fns.fn_headers, feature_headers, nameof(program), nameof(part4));

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
                    var row_comments = data_controller.get_row_comments(row_index, data_encoded_item);

                    var row_comments_csv = string.Join(",", row_comments);

                    return (row_feature_values_csv, row_comments_csv);
                });

                tasks4.Add(task);

                wait_tasks(tasks4.ToArray<Task>(), nameof(program), nameof(part4));
            }

            data_encoded_list.Clear();
            data_encoded_list = null;
            wait_tasks(tasks4.ToArray<Task>(), nameof(program), nameof(part4));


            var comments_lines = tasks4.Select(a => a.Result.row_comments_csv).ToList();
            comments_lines.Insert(0, row_comments_header_csv);
            
            io_proxy.WriteAllLines(fns.fn_comments, comments_lines, nameof(program), nameof(part4));

            var features_lines = tasks4.Select(a => a.Result.row_feature_values_csv).ToList();
            features_lines.Insert(0, row_feature_header_csv);
            
            io_proxy.WriteAllLines(fns.fn_features, features_lines, nameof(program), nameof(part4));
        }

        private static List<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)> part3(
            (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
                output_folder, int? first_index, int? last_index) cmd_params, List<subsequence_classification_data> class_data_list, feature_types feature_types)
        {
            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{cmd_params.class_id} {cmd_params.class_name}: 3. Encoding data for {class_data_list.Count} items...",
                    nameof(program), nameof(part3));
            }

            var tasks3 =
                new List<Task<(subsequence_classification_data instance_meta_data, List<feature_info> feature_info)>>();

            for (var i = 0; i < class_data_list.Count; i++)
            {
                var class_data_item = class_data_list[i];

                class_data_list[i] = null;

                var task = Task.Run(() =>
                {
                    var encoded_class_data_item = (class_data_item, subsequence_classification_data_methods.encode_subsequence_classification_data_row(class_data_item, cmd_params.max_features, feature_types));

                    return encoded_class_data_item;
                });

                tasks3.Add(task);

                wait_tasks(tasks3.ToArray<Task>(), nameof(program), nameof(part3));
            }
            
            class_data_list.Clear();
            class_data_list = null;
            wait_tasks(tasks3.ToArray<Task>(), nameof(program), nameof(part3));
            var data_encoded_list = tasks3.Select(a => a.Result).ToList();
            return data_encoded_list;
        }

        private static List<subsequence_classification_data> part2(
            (string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string
                output_folder, int? first_index, int? last_index) cmd_params, List<protein_subsequence_info> psi_list, feature_types feature_types, protein_subsequence_info _template_protein)
        {
            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{cmd_params.class_id} {cmd_params.class_name}: 2. Loading available data for {psi_list.Count} items...",
                    nameof(program), nameof(part2));
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

                wait_tasks(tasks2.ToArray<Task>(), nameof(program), nameof(part2));
            }
            
            psi_list.Clear();
            psi_list = null;
            wait_tasks(tasks2.ToArray<Task>(), nameof(program), nameof(part2));
            var class_data_list = tasks2.Select(a => a.Result).ToList();
            return class_data_list;
        }

        private static List<protein_subsequence_info> part1((string[] area, bool use_dssp3, int class_id, string class_name, int min_sequence_length, int max_features, string output_folder, int? first_index, int? last_index) cmd_params, List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> pdb_id_list)
        {
            if (cmd_params.class_name == "standard_coil")
            {
                pdb_id_list = pdb_id_list.GroupBy(a => a.pdb_id).Select(a => a.First()).ToList();
            }

            if (cmd_params.first_index == null && cmd_params.last_index == null)
            {
                io_proxy.WriteLine(
                    $@"{cmd_params.class_id} {cmd_params.class_name}: 1. Loading pdb info from {pdb_id_list.Count} files...",
                    nameof(program), nameof(part1));
            }

            var cache_filename = Path.Combine(cmd_params.output_folder, $@"l__({cmd_params.class_name}).csv");

            var cache = protein_subsequence_info.load(cache_filename);

            if (cache != null && cache.Count > 0)
            {
                return cache;
            }
            
            var tasks1 = new List<Task<List<protein_subsequence_info>>>();

            for (var i = 0; i < pdb_id_list.Count; i++)
            {
                var pdb_id_item = pdb_id_list[i];

                pdb_id_list[i] = default;

                var task = Task.Run(() =>
                {
                    List<protein_subsequence_info> psi_list = new List<protein_subsequence_info>();

                    if (cmd_params.class_name == "standard_coil")
                    {
                        psi_list = dataset_gen_coils.find_coils(pdb_id_item.dimer_type, pdb_id_item.pdb_id, cmd_params.class_id,
                            cmd_params.class_name, cmd_params.use_dssp3);
                    }
                    else if (cmd_params.class_name == "dimorphic_coil")
                    {
                        psi_list.Add(dataset_gen_dimorphic.get_dhc_item(cmd_params.class_id, cmd_params.class_name,
                            cmd_params.use_dssp3, true, false, pdb_id_item));
                    }
                    else
                    {
                        throw new Exception();
                    }

                    return psi_list;
                });

                tasks1.Add(task);

                //wait_tasks(tasks1.ToArray<Task>(), nameof(program), nameof(part1));
            }

            pdb_id_list.Clear();
            pdb_id_list = null;
            wait_tasks(tasks1.ToArray<Task>(), nameof(program), nameof(part1));
            var psi_list = tasks1.SelectMany(a => a.Result).ToList();

            psi_list = psi_list.Where(a => a.aa_subsequence.Length >= cmd_params.min_sequence_length).ToList();

            protein_subsequence_info.save(cache_filename, psi_list);

            return psi_list;
        }
    }
}

