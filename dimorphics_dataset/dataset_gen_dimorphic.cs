using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class dataset_gen_dimorphic
    {
        internal static List<protein_subsequence_info> run_dhc_dataset_maker(enum_substructure_type strand_type, int class_id, string class_name, bool use_dssp3 = true)
        {
            var class_name_in_file = "";
            var include_host_coil = false;
            var full_protein_seq = false;

            switch (strand_type)
            {
                case enum_substructure_type.dimorphic:
                    class_name_in_file = "Single";
                    break;

                case enum_substructure_type.dimorphic_coil:
                    class_name_in_file = "Single";
                    include_host_coil = true;
                    break;

                case enum_substructure_type.standard_strand:
                    class_name_in_file = "Multiple";
                    break;

                case enum_substructure_type.standard_strand_coil:
                    class_name_in_file = "Multiple";
                    include_host_coil = true;
                    break;

                case enum_substructure_type.dimorphic_full_protein_sequence:
                    class_name_in_file = "Single";
                    full_protein_seq = true;
                    break;

                case enum_substructure_type.standard_strand_full_protein_sequence:
                    class_name_in_file = "Multiple";
                    full_protein_seq = true;
                    break;

                case enum_substructure_type.standard_coil:
                    throw new Exception();

                default:
                    throw new Exception();
            }

            var dimorphics_data_all = program.get_dataset_pdb_id_list(null);

            var dimorphics_data = dimorphics_data_all.Where(a => a.class_name == class_name_in_file).ToList();

            var get_dhc_ret = get_dhc_list(class_id, class_name, use_dssp3, include_host_coil, full_protein_seq, dimorphics_data);

            return get_dhc_ret;
        }

        
        internal static List<protein_subsequence_info> get_dhc_list(int class_id, string class_name, bool use_dssp3, bool include_host_coil, bool full_protein_seq,
            List<(string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index)> dimorphics_data)
        {

            if (dimorphics_data == null)
            {
                throw new ArgumentNullException(nameof(dimorphics_data));
            }

            var tasks = new List<Task<protein_subsequence_info>>();

            for (var i = 0; i < dimorphics_data.Count; i++)
            {
                var dimorphics_interface = dimorphics_data[i];

                var task = Task.Run(() =>
                {
                    return get_dhc_item(class_id, class_name, use_dssp3, include_host_coil, full_protein_seq, dimorphics_interface);
                });

                tasks.Add(task);
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), nameof(dataset_gen_dimorphic), nameof(get_dhc_list));

            var psi_list = tasks.Select(a => a.Result).ToList();

            return psi_list;
        }

        internal static (int strand_array_index, List<int> strand_res_ids) find_subseq_index(string strand_seq, string strand_protein, List<int> strand_protein_res_ids, int optional_res_index = -1)
        {
            var strand_array_index = strand_protein?.IndexOf(strand_seq, StringComparison.InvariantCulture) ?? -1;

            if (optional_res_index > -1)
            {
                strand_array_index = strand_protein_res_ids.IndexOf(optional_res_index);

                if (strand_protein.Substring(strand_array_index, strand_seq.Length) != strand_seq)
                {
                    var indexes = strand_protein.Select((a, i) => i + strand_seq.Length < strand_protein.Length && strand_protein.Substring(i, strand_seq.Length) == strand_seq ? i : -1).Where(a => a > -1).ToList();
                    strand_array_index = indexes.Select(a => (index1: a, index2: Math.Abs(strand_array_index - a))).OrderBy(a => a.index2).First().index1;
                }
            }

            if (strand_protein?.Substring(strand_array_index, strand_seq.Length) != strand_seq)
            {
                throw new Exception("Not correct position");
            }

            var strand_res_ids = strand_protein_res_ids.Skip(strand_array_index).Take(strand_seq?.Length ?? 0).ToList();

            return (strand_array_index, strand_res_ids);
        }

        internal static char get_chain_id_from_chain_number(List<atom> pdb_atoms, int chain_number)
        {
            var chain_list_default_order = pdb_atoms.First(a => a.chain_order_in_pdb_file != null && a.chain_order_in_pdb_file.Count > 0).chain_order_in_pdb_file; //pdb_atoms.Select(a => a.chain_id).Distinct().ToList();
            var pdb_master_atoms = atom.select_amino_acid_master_atoms(null, pdb_atoms);
            var sequences = pdb_master_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, sequence: string.Join("", a.Select(b => b.amino_acid).ToList()))).ToList();
            var chain_id = chain_list_default_order[chain_number].chain_id;

            return chain_id;
        }

        internal static /*subsequence_classification_data*/protein_subsequence_info get_dhc_item(
            int class_id,
            string class_name,
            bool use_dssp3, 
            bool include_host_coil,
            bool full_protein_seq,
            (string pdb_id, string dimer_type, string class_name, string symmetry_mode, string parallelism, int chain_number, string strand_seq, string optional_res_index) dimorphics_interface
            )
        {
            var pdb_atoms = atom.load_atoms_pdb(dimorphics_interface.pdb_id, new atom.load_atoms_pdb_options()
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
                
            }).Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();


            var pdb_id = dimorphics_interface.pdb_id;
            var chain_number = dimorphics_interface.chain_number;

            var dimer_type = dimorphics_interface.dimer_type;
            var class_name2 = dimorphics_interface.class_name;
            var symmetry_mode = dimorphics_interface.symmetry_mode;
            var parallelism = dimorphics_interface.parallelism;
            var strand_seq = dimorphics_interface.strand_seq;
            var optional_res_index = dimorphics_interface.optional_res_index;

            //var chain_list_default_order = pdb_atoms.First(a => a.chain_order_in_pdb_file != null && a.chain_order_in_pdb_file.Count > 0).chain_order_in_pdb_file; //pdb_atoms.Select(a => a.chain_id).Distinct().ToList();
            var pdb_master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_atoms);
            var sequences = pdb_master_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, sequence: string.Join("", a.Select(b => b.amino_acid).ToList()))).ToList();
            //var chain_id = chain_list_default_order[chain_number].chain_id;

            var chain_id = get_chain_id_from_chain_number(pdb_atoms, chain_number);

            pdb_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
            var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
            var pdb_chain_master_atoms = pdb_master_atoms.Where(a => a.chain_id == chain_id).ToList();
            

            var strand_protein = sequences.First(a => a.chain_id == chain_id).sequence;
            //var strand_array_index = strand_protein.IndexOf(strand_seq, StringComparison.InvariantCulture);

            //var res_ids = pdb_chain_master_atoms.Select(a => a.residue_index).ToList();


            //if (!string.IsNullOrWhiteSpace(optional_res_index))
            //{
            //    strand_array_index = res_ids.IndexOf(int.Parse(optional_res_index, NumberStyles.Integer, CultureInfo.InvariantCulture));

            //    if (strand_protein.Substring(strand_array_index, strand_seq.Length) != strand_seq)
            //    {
            //        var indexes = strand_protein.Select((a, i) => i + strand_seq.Length < strand_protein.Length && strand_protein.Substring(i, strand_seq.Length) == strand_seq ? i : -1).Where(a => a > -1).ToList();
            //        strand_array_index = indexes.Select(a => (index1: a, index2: Math.Abs(strand_array_index - a))).OrderBy(a => a.index2).First().index1;
            //    }
            //}

            //if (strand_protein.Substring(strand_array_index, strand_seq.Length) != strand_seq)
            //{
            //    throw new Exception("Not correct position");
            //}

            var find_subseq_index_ret = find_subseq_index(strand_seq, strand_protein,
                pdb_chain_master_atoms.Select(a => a.residue_index).ToList(),
                int.TryParse(optional_res_index, NumberStyles.Integer, CultureInfo.InvariantCulture, out var opt_res_ix)
                    ? opt_res_ix
                    : -1);
            var strand_array_index = find_subseq_index_ret.strand_array_index;

            // get list of dimorphics with same pdbid and chain
            /*var other_interfaces_in_chain = dimorphics_data_all.Where(a => / * a != di && * / a.pdb_id == dimorphics_interface.pdb_id && a.chain_number == dimorphics_interface.chain_number).ToList();


            foreach (var other_interface in other_interfaces_in_chain)
            {
                // here, set interface marker properties for atoms

                var other_pdb_id = other_interface.pdb_id;
                var other_chain_number = other_interface.chain_number;
                var other_dimer_type = other_interface.dimer_type;
                var other_class_name = other_interface.class_name;
                var other_symmetry_mode = other_interface.symmetry_mode;
                var other_parallelism = other_interface.parallelism;
                var other_strand_seq = other_interface.strand_seq;
                var other_optional_res_index = other_interface.optional_res_index;

                var other_chain_list_default_order = pdb_atoms.First(a => a.chain_order_in_pdb_file != null && a.chain_order_in_pdb_file.Count > 0).chain_order_in_pdb_file; //pdb_atoms.Select(a => a.chain_id).Distinct().ToList();
                var other_pdb_master_atoms = atom.select_amino_acid_master_atoms(other_pdb_id, pdb_atoms);
                var other_sequences = other_pdb_master_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, sequence: string.Join("", a.Select(b => b.amino_acid).ToList()))).ToList();
                var other_chain_id = other_chain_list_default_order[other_chain_number].chain_id;

                var other_pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == other_chain_id).ToList();
                var other_pdb_chain_master_atoms = other_pdb_master_atoms.Where(a => a.chain_id == other_chain_id).ToList();

                var other_strand_protein = other_sequences.First(a => a.chain_id == other_chain_id).sequence;
                var other_strand_array_index = other_strand_protein.IndexOf(other_strand_seq, StringComparison.InvariantCulture);

                var other_res_ids = other_pdb_chain_master_atoms.Select(a => a.residue_index).ToList();

                if (!string.IsNullOrWhiteSpace(other_optional_res_index))
                {
                    other_strand_array_index = other_res_ids.IndexOf(int.Parse(other_optional_res_index));

                    if (other_strand_protein.Substring(other_strand_array_index, other_strand_seq.Length) != other_strand_seq)
                    {
                        var other_indexes = other_strand_protein.Select((a, i) => i + other_strand_seq.Length < other_strand_protein.Length && other_strand_protein.Substring(i, other_strand_seq.Length) == other_strand_seq ? i : -1).Where(a => a > -1).ToList();
                        other_strand_array_index = other_indexes.Select(a => (index1: a, index2: Math.Abs(other_strand_array_index - a))).OrderBy(a => a.index2).First().index1;
                    }
                }

                if (other_strand_protein.Substring(other_strand_array_index, other_strand_seq.Length) != other_strand_seq)
                {
                    throw new Exception("Not correct position");
                }

                var dimorphic_dssp3 = "";

                for (var i = other_strand_array_index; i < other_strand_array_index + other_strand_seq.Length; i++)
                {
                    other_pdb_chain_master_atoms[i].parallelism = other_parallelism.Trim().First();
                    other_pdb_chain_master_atoms[i].symmetry_mode = other_symmetry_mode.Trim().First();
                    other_pdb_chain_master_atoms[i].strand_interface_type = other_class_name.Trim().First();

                    other_pdb_chain_master_atoms[i].is_strand_interface_atom = 1;
                    other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_strand_interface_atom = 1);

                    if (string.Equals(other_class_name, "Single", StringComparison.InvariantCultureIgnoreCase) || string.Equals(other_class_name, "Dimorphic", StringComparison.InvariantCultureIgnoreCase))
                    {
                        other_pdb_chain_master_atoms[i].is_dimorphic_strand_interface_atom = 1;
                        other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_dimorphic_strand_interface_atom = 1);
                        dimorphic_dssp3 += other_pdb_chain_master_atoms[i].monomer_dssp3;
                    }
                    else if (string.Equals(other_class_name, "Multiple", StringComparison.InvariantCultureIgnoreCase) || string.Equals(other_class_name, "Standard", StringComparison.InvariantCultureIgnoreCase))
                    {
                        other_pdb_chain_master_atoms[i].is_standard_strand_interface_atom = 1;
                        other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_standard_strand_interface_atom = 1);
                    }
                    else if (string.Equals(other_class_name, "Hybrid", StringComparison.InvariantCultureIgnoreCase))
                    {
                        other_pdb_chain_master_atoms[i].is_hybrid_strand_interface_atom = 1;
                        other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_hybrid_strand_interface_atom = 1);
                    }
                }

                //if (!string.IsNullOrWhiteSpace(dimorphic_dssp3)) io_proxy.WriteLine(dimorphic_dssp3);
            }*/

            /*
             //use this code to generate the spreadsheets for dssp/stride to ss prediction comparisons in q3 format:

            foreach (var g in pdb_atoms.GroupBy(a => (a.pdb_id, a.chain_id)))
            {
                var p = g.Key.pdb_id;
                var c = g.Key.chain_id;
                var d = g.ToList();

                atom.compare_dssp_to_mpsa(p, c, d);
            }
            */

            //for (var i = strand_array_index; i < strand_array_index + strand_seq.Length; i++)
            //{
            //    pdb_chain_master_atoms[i].is_strand_interface_atom = 1;
            //    pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a=>a.is_strand_interface_atom=1);
            //}

            var strand_sequence_first_index = strand_array_index;
            var strand_sequence_last_index = (strand_array_index + strand_seq.Length) - 1;
            //var dimorphic_sequence_centre = (double)(dimorphic_sequence_first + dimorphic_sequence_last) / (double)2;
            //var dimorphic_sequence_length = dimorphic_seq.Length;

            var shc_sequence_before = "";
            var shc_start_array_index = strand_sequence_first_index;

            if (include_host_coil)
            {
                var continuous_strand1 = (use_dssp3 ? pdb_chain_master_atoms[strand_sequence_first_index].dssp3_multimer : pdb_chain_master_atoms[strand_sequence_first_index].dssp_multimer) == 'E';

                for (var i = strand_sequence_first_index - 1; i >= 0; i--)
                {
                    var index_dssp_code = use_dssp3 ? pdb_chain_master_atoms[i].dssp3_multimer : pdb_chain_master_atoms[i].dssp_multimer;

                    continuous_strand1 = continuous_strand1 && index_dssp_code == 'E';

                    if (!continuous_strand1 && index_dssp_code != 'C') break;

                    shc_start_array_index = i;
                    shc_sequence_before = pdb_chain_master_atoms[i].amino_acid + shc_sequence_before;
                }
            }
            //if (dhc_start_array_index == -1) dhc_start_array_index = 0;


            var shc_sequence_after = "";

            var shc_end_array_index = strand_sequence_last_index;

            if (include_host_coil)
            {
                var continuous_strand2 = (use_dssp3 ? pdb_chain_master_atoms[strand_sequence_last_index].dssp3_multimer : pdb_chain_master_atoms[strand_sequence_last_index].dssp_multimer) == 'E';

                for (var i = strand_sequence_last_index + 1; i < pdb_chain_master_atoms.Count; i++)
                {
                    var index_dssp_code = use_dssp3 ? pdb_chain_master_atoms[i].dssp3_multimer : pdb_chain_master_atoms[i].dssp_multimer;

                    continuous_strand2 = continuous_strand2 && index_dssp_code == 'E';

                    if (!continuous_strand2 && index_dssp_code != 'C') break;

                    shc_end_array_index = i;
                    shc_sequence_after += pdb_chain_master_atoms[i].amino_acid;
                }
            }


            //if (dhc_end_array_index == -1) dhc_end_array_index = pdb_chain_master_atoms.Count() - 1;

            var shc_subsequence_master_atoms = pdb_chain_master_atoms.Skip(shc_start_array_index).Take((shc_end_array_index - shc_start_array_index) + 1).ToList();

            if (full_protein_seq)
            {
                shc_subsequence_master_atoms = pdb_chain_master_atoms;
                shc_start_array_index = 0;
                shc_end_array_index = pdb_chain_master_atoms.Count - 1;
                strand_seq = string.Join("", shc_subsequence_master_atoms.Select(a => a.amino_acid).ToList());
                shc_sequence_before = "";
                shc_sequence_after = "";
            }

            //var shc_subsequence_atoms = pdb_chain_atoms.Where(a => shc_subsequence_master_atoms.Any(b => a.pdb_id == b.pdb_id && a.chain_id == b.chain_id && a.residue_index == b.residue_index && a.i_code == b.i_code)).ToList();

            var shc_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.amino_acid).ToList());

            //var shc_dssp_multimer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.dssp_multimer).ToList());
            //var shc_dssp_monomer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.dssp_monomer).ToList());

            //var shc_stride_multimer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.stride_multimer).ToList());
            //var shc_stride_monomer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.stride_monomer).ToList());

            var shc_length = (shc_end_array_index - shc_start_array_index) + 1;

            if (shc_sequence.Length != shc_length) throw new Exception("shc length error");
            if (shc_sequence_before + strand_seq + shc_sequence_after != shc_sequence) throw new Exception("shc seq error");

            //if (shc_sequence.Length < min_strand_host_coil_length) continue;

            var psi = new protein_subsequence_info()
            {
                class_id = class_id,
                class_name = class_name,
                pdb_id = pdb_id,
                chain_id = chain_id,
                dimer_type = dimer_type,
                parallelism = parallelism,
                symmetry_mode = symmetry_mode,
                res_ids = shc_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
                aa_subsequence = shc_sequence,
                //aa_before = "",
                //aa_after = "",
                //aa_protein = strand_protein
            };


            return psi;
            //var scd = new subsequence_classification_data
            //{
            //    dimer_type = dimer_type,
            //    parallelism = parallelism,
            //    symmetry_mode = symmetry_mode,
            //    class_id = class_id,
            //    class_name = class_name,

            //    pdb_id = pdb_id,
            //    chain_id = chain_id,
            //    res_ids = shc_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
            //    aa_subsequence = shc_sequence,
            //    dssp_multimer_subsequence = shc_dssp_multimer_sequence,
            //    dssp_monomer_subsequence = shc_dssp_monomer_sequence,
            //    stride_multimer_subsequence = shc_stride_multimer_sequence,
            //    stride_monomer_subsequence = shc_stride_monomer_sequence,
            //    pdb_chain_atoms = pdb_chain_atoms,
            //    pdb_chain_master_atoms = pdb_chain_master_atoms,
            //    subsequence_atoms = shc_subsequence_atoms,
            //    subsequence_master_atoms = shc_subsequence_master_atoms,
            //};

            //var fx = info_foldx.load_calc_energy_differences(scd.pdb_id, scd.chain_id, scd.res_ids, false, enum_protein_data_source.subsequence_3d);
            //scd.foldx_energy_differences = fx;

            //return scd;
        }
    }
}