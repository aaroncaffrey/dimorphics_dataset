﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class dataset_gen_coils
    {
        public const string module_name = nameof(dataset_gen_coils);

        internal static /*List<subsequence_classification_data>*/ List<protein_subsequence_info> find_coils(string dimer_type, string pdb_id, int chain_number, int class_id, string class_name, bool use_dssp3 = true, bool use_multimer_dssp = false, CancellationTokenSource cts = null)
        {
            const string method_name = nameof(find_coils);

            using var i_cts = new CancellationTokenSource();
            if (cts == null) cts = i_cts;

            // find the coils within a given pdb file
            //var use_multimer_dssp = true;

            var pdb_atoms = atom.load_atoms_pdb
                (
                pdb_id:pdb_id,
                new load_atoms_pdb_options()
                { 
                    find_3d_intramolecular=false,

                    load_1d_blast_pssms = false,
                    load_1d_iup_data = false,
                    load_1d_sable = false,
                    load_1d_dna_binding = false,

                    load_2d_mpsa_sec_struct_predictions = false,

                    load_3d_rsa_data =false,
                    load_3d_dssp_data=true,
                    load_3d_stride_data=false,
                    load_3d_ring_data=false,
                    load_3d_foldx_ala_scan = false,
                }
                )
                .Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();

            var chain_id_filter = dataset_gen_dimorphic.get_chain_id_from_chain_number(pdb_atoms, chain_number);

            pdb_atoms = pdb_atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.OrdinalIgnoreCase)).ToList();
            
            pdb_atoms = pdb_atoms.Where(a => a.chain_id == chain_id_filter).ToList();


            //pdb_atoms = pdb_atoms.Where(a => a.Multimer_DSSP3 == 'C').ToList(); // only need info about coils

            var chain_ids = pdb_atoms.Select(a => a.chain_id).Distinct().ToList();

            //var tasks = new List<Task<List<subsequence_classification_data>>>();
            var tasks = new List<Task<List<protein_subsequence_info>>>();

            var tasks_start_time = DateTime.Now;
            
            foreach (var loop_chain_id in chain_ids)
            {
                var chain_id = loop_chain_id;

                var task = Task.Run(() =>
                {
                    if (cts != null && cts.IsCancellationRequested) return default;

                    //var pdb_chain_coils = new List<subsequence_classification_data>();
                    var pdb_chain_coils = new List<protein_subsequence_info>();
                    
                    var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
                    var pdb_chain_master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_chain_atoms);
                    var pdb_chain_aa_sequence = string.Join(/*program.string_debug*/($@""), pdb_chain_master_atoms.Select(a => a.amino_acid).ToList());

                    var coil_subsequence_atoms = new List<atom>();
                    var coil_subsequence_master_atoms = new List<atom>();


                    for (var pdb_chain_master_atoms_index = 0; pdb_chain_master_atoms_index < pdb_chain_master_atoms.Count; pdb_chain_master_atoms_index++)
                    {
                        if (cts != null && cts.IsCancellationRequested) return default;

                        var master_atom = pdb_chain_master_atoms[pdb_chain_master_atoms_index];

                        var master_atom_ss = ' ';
                        if (use_dssp3 && use_multimer_dssp) {master_atom_ss = master_atom.dssp3_multimer;}
                        else if (use_dssp3 && !use_multimer_dssp) {master_atom_ss = master_atom.dssp3_monomer;}
                        else if (!use_dssp3 && use_multimer_dssp){ master_atom_ss = master_atom.dssp_multimer;}
                        else if (!use_dssp3 && !use_multimer_dssp) {master_atom_ss = master_atom.dssp_monomer;}
                        else {throw new Exception();}

                        var this_res_id = master_atom.residue_index;
                        var last_res_id = pdb_chain_master_atoms_index == 0 ? this_res_id - 1 : pdb_chain_master_atoms[pdb_chain_master_atoms_index - 1].residue_index;
                        var next_res_id = pdb_chain_master_atoms_index == pdb_chain_master_atoms.Count - 1 ? this_res_id + 1 : pdb_chain_master_atoms[pdb_chain_master_atoms_index + 1].residue_index;

                        var res_id_consecutive_to_last = this_res_id == last_res_id + 1;
                        var res_id_consecutive_to_next = this_res_id == next_res_id - 1;

                        //if ((use_dssp3 && master_atom.multimer_dssp3 == 'C') || (!use_dssp3 && master_atom.multimer_dssp == 'C'))
                        if (master_atom_ss == 'C')
                        {
                            coil_subsequence_master_atoms.Add(master_atom);
                            coil_subsequence_atoms.AddRange(master_atom.amino_acid_atoms);
                        }

                        if (coil_subsequence_master_atoms.Count > 0 && 
                            //( ((use_dssp3 && master_atom.multimer_dssp3 != 'C') || (!use_dssp3 && master_atom.multimer_dssp != 'C')) || !res_id_consecutive_to_next || pdb_chain_master_atoms_index == pdb_chain_master_atoms.Count - 1))
                            (master_atom_ss != 'C' || !res_id_consecutive_to_next || pdb_chain_master_atoms_index == pdb_chain_master_atoms.Count - 1))
                        {
                            if (coil_subsequence_master_atoms.Count == pdb_chain_master_atoms.Count)
                            {
                                break; // whole chain is a coil - not useful
                            }

                            //if (coil_subsequence_master_atoms.Count >= min_coil_length_required)
                            //{
                                var psi = new protein_subsequence_info()
                                {
                                    dimer_type = dimer_type,
                                    symmetry_mode = /*program.string_debug*/($@""),
                                    parallelism = /*program.string_debug*/($@""),
                                    class_id = class_id,
                                    class_name = class_name,
                                    pdb_id = pdb_id,
                                    chain_id = chain_id,
                                    res_ids = coil_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
                                    aa_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                                    //aa_before = /*program.string_debug*/($@""),
                                    //aa_after = /*program.string_debug*/($@""),
                                    //aa_protein = pdb_chain_aa_sequence,
                                };

                                pdb_chain_coils.Add(psi);
                                //var scd = new subsequence_classification_data
                                //{
                                //    dimer_type = dimer_type,
                                //    class_id = coils_class_id,
                                //    pdb_id = pdb_id,
                                //    chain_id = chain_id,
                                //    res_ids = coil_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
                                //    aa_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                                //    dssp_multimer_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.multimer_dssp).ToList()),
                                //    dssp_monomer_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.monomer_dssp).ToList()),
                                //    stride_multimer_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.stride_multimer).ToList()),
                                //    stride_monomer_subsequence = string.Join(/*program.string_debug*/($@""), coil_subsequence_master_atoms.Select(a => a.stride_monomer).ToList()),
                                //    pdb_chain_atoms = pdb_chain_atoms,
                                //    pdb_chain_master_atoms = pdb_chain_master_atoms,
                                //    subsequence_atoms = coil_subsequence_atoms,
                                //    subsequence_master_atoms = coil_subsequence_master_atoms,
                                //};
                                //var fx = info_foldx.load_calc_energy_differences(scd.pdb_id, scd.chain_id, scd.res_ids, false, enum_protein_data_source.subsequence_3d);
                                //scd.foldx_energy_differences = fx;
                                //
                                //pdb_chain_coils.Add(scd);
                            //}

                            coil_subsequence_master_atoms = new List<atom>();
                            coil_subsequence_atoms = new List<atom>();
                        }
                    }

                    return !cts.IsCancellationRequested ? pdb_chain_coils : default;
                }, cts.Token);

                tasks.Add(task);

                program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, -1, module_name, method_name, cts);
                if (cts != null && cts.IsCancellationRequested) return default;
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, 0, module_name, method_name, cts);
            if (cts != null && cts.IsCancellationRequested) return default;

            var pdb_coils = tasks.SelectMany(a => a.Result).ToList();

            if (chain_ids != null && chain_ids.Count > 1)
            {
                var chain_seqs = atom.amino_acid_sequence(null, pdb_atoms);

                if (chain_seqs.Any(a => chain_seqs.Any(b => a != b && (a.aa_sequence.Contains(b.aa_sequence, StringComparison.Ordinal) || b.aa_sequence.Contains(a.aa_sequence, StringComparison.Ordinal)))))
                {
                    pdb_coils = pdb_coils.GroupBy(scd => (scd.pdb_id, /*scd.chain_id,*/ scd.aa_subsequence))
                        .Select(group => group.First()).Distinct().ToList();
                }
            }

            return pdb_coils;
        }
    }
}
