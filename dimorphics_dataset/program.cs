using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Loader;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class program
    {
        // todo: tortuosity for NH, Protein. tort for CA, CB, NH, other types of atoms.
        // todo: average distance between AAs in interface, NH, Protein.
        // todo: distance between different types of atoms...
        // todo: protr

        public static string data_root_folder = $@"C:\betastrands_dataset\";

        //public class class_info
        //{
        //    public int class_id;
        //    public substructure_type substructure_type;
        //    public int min_subsequence_length;
        //    public string class_name;
        //}

        public static List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> get_dataset_pdb_id_list(string class_name_in_file = "Single")
        {
            //List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)>

            var dimorphics_data1 = io.ReadAllLines(Path.Combine(program.data_root_folder, $@"csv", $@"distinct dimorphics list.csv"))
                .Skip(1)
                .Where(a => !string.IsNullOrWhiteSpace(a.Replace(",", "")))
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

        public static subsequence_classification_data get_subsequence_classificiation_data(protein_subsequence_info psi, subsequence_classification_data.feature_types feature_types)
        {
            // todo: check if pdb_folder should be pdb or pdb_repair?

            var pdb_atoms = atom.load_atoms_pdb
                (
                    pdb_id: psi.pdb_id,
                    new atom.load_atoms_pdb_options()
                    {
                        find_intramolecular = true,
                        find_intermolecular = false,
                        load_rsa_data = true,
                        load_dssp_data = true,
                        load_stride_data = true,
                        load_ring_data = true,
                        load_mpsa_sec_struct_predictions = true,
                        load_blast_pssms = true,
                        load_iup_data = true,
                        load_sable = true,
                        load_ala_scan = true,
                        load_dna_binding_vars = true,
                    }
                    )
                .Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();

            var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == psi.chain_id).ToList();
            var pdb_chain_master_atoms = atom.select_amino_acid_master_atoms(psi.pdb_id, pdb_chain_atoms);

            var subsequence_atoms = psi.res_ids.SelectMany(a => pdb_chain_atoms.Where(b => a.res_id == b.residue_index && a.amino_acid == b.amino_acid && a.i_code == b.i_code).ToList()).ToList();
            var subsequence_master_atoms = psi.res_ids.Select(a => pdb_chain_master_atoms.First(b => a.res_id == b.residue_index && a.amino_acid == b.amino_acid && a.i_code == b.i_code)).ToList();

            if (subsequence_master_atoms.Count != psi.res_ids.Count) throw new Exception();


            var scd = new subsequence_classification_data
            {
                dimer_type = psi.dimer_type ?? "",
                symmetry_mode = psi.symmetry_mode ?? "",
                parallelism = psi.parallelism ?? "",

                class_id = psi.class_id,
                class_name = psi.class_name,

                pdb_id = psi.pdb_id,
                chain_id = psi.chain_id,
                res_ids = psi.res_ids,
                
                aa_subsequence = string.Join("", subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                
                dssp_multimer_subsequence = string.Join("", subsequence_master_atoms.Select(a => a.multimer_dssp).ToList()),
                dssp_monomer_subsequence = string.Join("", subsequence_master_atoms.Select(a => a.monomer_dssp).ToList()),
                
                stride_multimer_subsequence = string.Join("", subsequence_master_atoms.Select(a => a.multimer_stride).ToList()),
                stride_monomer_subsequence = string.Join("", subsequence_master_atoms.Select(a => a.monomer_stride).ToList()),
                
                pdb_chain_atoms = pdb_chain_atoms,
                pdb_chain_master_atoms = pdb_chain_master_atoms,
                subsequence_atoms = subsequence_atoms,
                subsequence_master_atoms = subsequence_master_atoms,
                
                foldx_energy_differences = null,
                neighbourhood_1d = null,
                neighbourhood_3d = null,
                
                parent = null,
                protein_1d = null,
                protein_3d = null
            };

            scd.foldx_energy_differences = info_foldx.load_calc_energy_differences(scd.pdb_id, scd.chain_id, scd.res_ids, false, protein_data_sources.subsequence_3d);

            if (feature_types.feature_types_neighbourhood_1d != null &&
                feature_types.feature_types_neighbourhood_1d.AsArray().Any(a => a.value))
            {
                scd.neighbourhood_1d = atom.get_intramolecular_neighbourhood_1d(scd);
            }

            if (feature_types.feature_types_protein_1d != null &&
                feature_types.feature_types_protein_1d.AsArray().Any(a => a.value))
            {
                scd.protein_1d = atom.get_intramolecular_protein_1d(scd);
            }

            if (feature_types.feature_types_neighbourhood_3d != null &&
                feature_types.feature_types_neighbourhood_3d.AsArray().Any(a => a.value))
            {
                scd.neighbourhood_3d = atom.get_intramolecular_neighbourhood_3d(scd);
            }

            if (feature_types.feature_types_protein_3d != null &&
                feature_types.feature_types_protein_3d.AsArray().Any(a => a.value))
            {
                scd.protein_3d = atom.get_intramolecular_protein_3d(scd);
            }

            return scd;
        }

        public static void wait_tasks(Task[] tasks)
        {
            do
            {
                var incomplete_tasks = tasks.Where(a => !a.IsCompleted).ToList();

                if (incomplete_tasks.Count > 0)
                {
                    Task.WaitAny(incomplete_tasks.ToArray<Task>());

                    var incomplete = tasks.Count(a => !a.IsCompleted);
                    var complete = tasks.Length - incomplete;
                    var pct = (double)complete / (double)tasks.Length;
                    io.WriteLine($@"{complete} / {tasks.Length} ( {pct:0.00} % )", nameof(program), nameof(Main));
                }
            } while (tasks.Any(a => !a.IsCompleted));

            Task.WaitAll(tasks.ToArray<Task>());
        }

        public static void Main(string[] args)
        {
            Console.CancelKeyPress += (sender, eventArgs) =>
            {
                io.WriteLine($@"Console.CancelKeyPress", nameof(program), nameof(Main));
            };
            AssemblyLoadContext.Default.Unloading += context =>
            {
                io.WriteLine($@"AssemblyLoadContext.Default.Unloading", nameof(program), nameof(Main));
            };
            AppDomain.CurrentDomain.ProcessExit += (sender, eventArgs) =>
            {
                io.WriteLine($@"AppDomain.CurrentDomain.ProcessExit", nameof(program), nameof(Main));
            };

            args = new string[] { "2i", "true", "+1", "dimorphic_coil", "3", "100", $@"e:\dataset\test\" };

            io.WriteLine(string.Join(" ", args), nameof(program), nameof(Main));
            // dimorphics_dataset.exe  [2i,2n,2p,3i,3n,3p] [use_dssp3=true|false] [class_id=-1|+1] [class_name] [min_subseq_len=3] [max_features=100] [output_folder]
            //
            // coils : dimorphics_dataset.exe  2i,2n,2p,3i,3n,3p True 3 0 ctl coils
            //
            // dimorphics : dimorphics_dataset.exe  2i,2n,2p,3i,3n,3p True 3 0 ctl dimorphics
            //


            if (args == null || args.Length == 0)
            {
                throw new Exception($@"No {nameof(args)} specified.");

                //args = new string[] { "2i,2n,2p,3i,3n,3p", "true", "3", "ctl", "func", "args", "0" };
            }

            var do_2d_interface = false;
            var do_2d_nh = false;
            var do_2d_protein = false;

            var do_3d_interface = false;
            var do_3d_nh = false;
            var do_3d_protein = false;


            var arg_index = 0;


            var feature_opts = args[arg_index++].Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries).Select(a => a.ToLowerInvariant()).ToList();
            do_2d_interface = feature_opts.Any(a => a.StartsWith("2i", StringComparison.InvariantCultureIgnoreCase));
            do_2d_nh = feature_opts.Any(a => a.StartsWith("2n", StringComparison.InvariantCultureIgnoreCase));
            do_2d_protein = feature_opts.Any(a => a.StartsWith("2p", StringComparison.InvariantCultureIgnoreCase));
            do_3d_interface = feature_opts.Any(a => a.StartsWith("3i", StringComparison.InvariantCultureIgnoreCase));
            do_3d_nh = feature_opts.Any(a => a.StartsWith("3n", StringComparison.InvariantCultureIgnoreCase));
            do_3d_protein = feature_opts.Any(a => a.StartsWith("3p", StringComparison.InvariantCultureIgnoreCase));
            feature_opts = feature_opts.Except(new string[] { "2i", "2n", "2p", "3i", "3n", "3p" }).ToList();
            if (feature_opts.Count > 0) throw new Exception("Unknown args: " + string.Join(", ", feature_opts));

            var dataset_name = string.Join("_", new string[] { (do_2d_interface ? "2i" : ""), (do_2d_nh ? "2n" : ""), (do_2d_protein ? "2p" : ""), (do_3d_interface ? "3i" : ""), (do_3d_nh ? "3n" : ""), (do_3d_protein ? "3p" : ""), }.Where(a => !string.IsNullOrWhiteSpace(a)).ToList());
            var use_dssp3 = bool.Parse(args[arg_index++]);
            var class_id = int.Parse(args[arg_index++], NumberStyles.Integer, CultureInfo.InvariantCulture);
            var class_name = args[arg_index++].ToLowerInvariant();
            var min_subsequence_length = int.Parse(args[arg_index++], NumberStyles.Integer, CultureInfo.InvariantCulture);
            var max_features = int.Parse(args[arg_index++], NumberStyles.Integer, CultureInfo.InvariantCulture);
            var output_folder = args[arg_index++];

            
            //var mode = args[arg_index++].ToLowerInvariant();
            //var func = args[arg_index++].ToLowerInvariant();


            io.WriteLine($@"{nameof(dataset_name)} = {dataset_name}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(use_dssp3)} = {use_dssp3}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(class_id)} = {class_id}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(class_name)} = {class_name}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(min_subsequence_length)} = {min_subsequence_length}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(max_features)} = {max_features}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(output_folder)} = {output_folder}", nameof(program), nameof(Main));


            var sw1 = new Stopwatch();
            sw1.Start();

            var feature_types = new subsequence_classification_data.feature_types()
            {
                feature_types_subsequence_1d = !do_2d_interface ? null : new subsequence_classification_data.feature_types_1d()
                {
                    pse_aac_sequence_classification_data = do_2d_interface,
                    sequence_geometry_classification_data = do_2d_interface,
                    mpsa_classification_data_subsequence = do_2d_interface,
                    intrinsically_unordered_data = do_2d_interface,
                    aa_index_classification_data = do_2d_interface,
                    sable_classification_data = do_2d_interface,
                    dna_binding_prediction_data = false, //must be false - protein level only
                    blast_pssm_subsequence_classification_data = do_2d_interface,
                    r_peptides = do_2d_interface
                },
                feature_types_neighbourhood_1d = !do_2d_nh ? null : new subsequence_classification_data.feature_types_1d()
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
                },
                feature_types_protein_1d = !do_2d_protein ? null : new subsequence_classification_data.feature_types_1d()
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
                },
                feature_types_subsequence_3d = !do_3d_interface ? null : new subsequence_classification_data.feature_types_3d()
                {
                    sasa_classification_data = do_3d_interface,
                    intramolecular_classification_data = do_3d_interface,
                    aa_aa_distances = do_3d_interface,
                    foldx_classification_data = do_3d_interface,
                    tortuosity_classification_data = do_3d_interface,
                    ring_classification_data = do_3d_interface,
                    pse_ssc_dssp_classification_data = do_3d_interface,
                },
                feature_types_neighbourhood_3d = !do_3d_nh ? null : new subsequence_classification_data.feature_types_3d()
                {
                    sasa_classification_data = do_3d_nh,
                    intramolecular_classification_data = do_3d_nh,
                    aa_aa_distances = do_3d_nh,
                    foldx_classification_data = do_3d_nh,
                    tortuosity_classification_data = do_3d_nh,
                    ring_classification_data = do_3d_nh,
                    pse_ssc_dssp_classification_data = do_3d_nh,
                },
                feature_types_protein_3d = !do_3d_protein ? null : new subsequence_classification_data.feature_types_3d()
                {
                    sasa_classification_data = do_3d_protein,
                    intramolecular_classification_data = do_3d_protein,
                    aa_aa_distances = do_3d_protein,
                    foldx_classification_data = do_3d_protein,
                    tortuosity_classification_data = do_3d_protein,
                    ring_classification_data = do_3d_protein,
                    pse_ssc_dssp_classification_data = do_3d_protein,
                }
            };

            io.WriteLine($@"{nameof(feature_types.feature_types_subsequence_1d)} = {feature_types.feature_types_subsequence_1d}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_1d)} = {feature_types.feature_types_neighbourhood_1d}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(feature_types.feature_types_protein_1d)} = {feature_types.feature_types_protein_1d}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(feature_types.feature_types_subsequence_3d)} = {feature_types.feature_types_subsequence_3d}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}", nameof(program), nameof(Main));
            io.WriteLine($@"{nameof(feature_types.feature_types_protein_3d)} = {feature_types.feature_types_protein_3d}", nameof(program), nameof(Main));





            var pdb_id_list = get_dataset_pdb_id_list();

            var atom_types = pdb_id_list.Select(a => a.pdb_id).Distinct().Select(a => atom.load_atoms_pdb(a,
                new atom.load_atoms_pdb_options()
                {
                    load_rsa_data = false,
                    load_dssp_data = false,
                    load_iup_data = false,
                    load_blast_pssms = false,
                    find_intramolecular = false,
                    load_ala_scan = false,
                    select_first_icode = false,
                    load_stride_data = false,
                    load_sable = false,
                    load_ring_data = false,
                    load_mpsa_sec_struct_predictions = false,
                    find_intermolecular = false,
                    load_dna_binding_vars = false,
                    first_model_only = false
                })).SelectMany(a => a.SelectMany(b => b.pdb_model_chain_atoms.Select(c => c.atom_type).ToList()).ToList()).ToList();
            var atom_types_count = atom_types.Distinct().Select(a => (type: a, count: atom_types.Count(b => b == a))).Select(a=> (type:a.type, count:a.count, pct: (double)a.count / (double)atom_types.Count))
                .OrderByDescending(a => a.Item2).ToList();
            
            atom_types_count.ForEach(a=> io.WriteLine($"{a.type} {a.count} {a.pct:0.00}"));

            pdb_id_list = pdb_id_list.Take(3).ToList();

            // 1. find subsequence details from analysis of pdb files
            io.WriteLine($@"{class_id} {class_name}: Loading pdb info...", nameof(program), nameof(Main));

            var tasks1 = new List<Task<List<protein_subsequence_info>>>();

            for (var i = 0; i < pdb_id_list.Count; i++)
            {
                var pdb_id_item = pdb_id_list[i];

                var task = Task.Run(() =>
                {
                    List<protein_subsequence_info> psi_list = new List<protein_subsequence_info>(); 
                    
                    if (class_name == "standard_coil") psi_list = dataset_gen_coils.find_coils(pdb_id_item.dimer_type, pdb_id_item.pdb_id, class_id, class_name, use_dssp3);
                    else if (class_name == "dimorphic_coil") psi_list.Add(dataset_gen_dimorphic.get_dhc_item(class_id, class_name, use_dssp3, true, false, pdb_id_item));
                    else throw new Exception();

                    return psi_list;
                });

                tasks1.Add(task);
            }

            wait_tasks(tasks1.ToArray<Task>());
            var psi_list = tasks1.SelectMany(a => a.Result).ToList();

            psi_list = psi_list.Where(a => a.aa_subsequence.Length >= min_subsequence_length).ToList();

            
            // 2. load available data sources
            io.WriteLine($@"{class_id} {class_name}: 2. Loading available data...", nameof(program), nameof(Main));
            
            var tasks2 = new List<Task<subsequence_classification_data>>();

            for (var i = 0; i < psi_list.Count; i++)
            {
                var psi = psi_list[i];

                var task = Task.Run(() =>
                {
                    var classificiation_data = get_subsequence_classificiation_data(psi, feature_types);

                    return classificiation_data;
                });

                tasks2.Add(task);
            }

            wait_tasks(tasks2.ToArray<Task>());
            var class_data_list = tasks2.Select(a => a.Result).ToList();


            // 3. encode the data as svm classification features
            io.WriteLine($@"{class_id} {class_name}: 3. Encoding data...", nameof(program), nameof(Main));

            var tasks3 =
                new List<Task<(instance_meta_data instance_meta_data,
                    List<subsequence_classification_data.feature_info> feature_info)>>();

            for (var i = 0; i < class_data_list.Count; i++)
            {
                var class_data_item = class_data_list[i];

                var task = Task.Run(() =>
                {
                    var encoded_class_data_item = subsequence_classification_data.encode_subsequence_classification_data_row(class_data_item, max_features, feature_types);

                    return encoded_class_data_item;
                });

                tasks3.Add(task);
            }

            wait_tasks(tasks3.ToArray<Task>());
            var data_encoded_list = tasks3.Select(a => a.Result).ToList();


            // 4. Save to file
            io.WriteLine($@"{class_id} {class_name}: 4. Saving encoded data to file...", nameof(program), nameof(Main));

            // get header row indexes in csv format
            var row_feature_header_csv = string.Join(",", Enumerable.Range(0, data_encoded_list.First().feature_info.Count));

            // get comments file header in csv format
            var row_comments_header_csv = string.Join(",", data_controller.get_row_comments_headers(data_encoded_list.First()));

            // get list of the feature headers in csv format
            var feature_headers = data_controller.get_feature_headers_lines_csv(data_encoded_list.First());
            var fn_headers = Path.Combine(output_folder, $"h__[{class_name}].csv");
            io.WriteAllLines(fn_headers, feature_headers, nameof(program), nameof(Main));

            var tasks4 = new List<Task<(string row_feature_values_csv, string row_comments_csv)>>();

            for (var i = 0; i < data_encoded_list.Count; i++)
            {
                var row_index = i;
                var data_encoded_item = data_encoded_list[i];
                
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
            }

            wait_tasks(tasks4.ToArray<Task>());
            

            var comments_lines = tasks4.Select(a => a.Result.row_comments_csv).ToList();
            comments_lines.Insert(0, row_comments_header_csv);
            var fn_comments = Path.Combine(output_folder, $"c__[{class_name}].csv");
            io.WriteAllLines(fn_comments, comments_lines, nameof(program), nameof(Main));

            var features_lines = tasks4.Select(a => a.Result.row_feature_values_csv).ToList();
            features_lines.Insert(0, row_feature_header_csv);
            var fn_features = Path.Combine(output_folder, $"f__[{class_name}].csv");
            io.WriteAllLines(fn_features, features_lines, nameof(program), nameof(Main));
            

            sw1.Stop();

            io.WriteLine($@"{class_id} {class_name}: Finished: ({sw1.Elapsed:dd\\:hh\\:mm\\:ss})", nameof(program), nameof(Main));
        }
    }
}

