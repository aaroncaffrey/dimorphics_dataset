using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class dataset_gen_coils
    {
        public static /*List<subsequence_classification_data>*/ List<protein_subsequence_info> find_coils(string dimer_type, string pdb_id, int class_id, string class_name, bool use_dssp3 = true, bool use_multimer_dssp = false)
        {
            // find the coils within a given pdb file
            //var use_multimer_dssp = true;

            var pdb_atoms = atom.load_atoms_pdb
                (
                pdb_id:pdb_id,
                new atom.load_atoms_pdb_options()
                { 
                    find_3d_intramolecular=false, 
                    
                    load_3d_rsa_data=false,
                    load_3d_dssp_data=true,
                    load_3d_stride_data=false,
                    load_3d_ring_data=false,
                    load_2d_mpsa_sec_struct_predictions=false,
                    load_2d_blast_pssms=false,
                    load_2d_iup_data=false,
                    load_2d_sable=false,
                    load_3d_foldx_ala_scan=false,
                    load_2d_dna_binding=false,
                }
                )
                .Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();

            pdb_atoms = pdb_atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.InvariantCultureIgnoreCase)).ToList();

            //pdb_atoms = pdb_atoms.Where(a => a.Multimer_DSSP3 == 'C').ToList(); // only need info about coils

            var chain_ids = pdb_atoms.Select(a => a.chain_id).Distinct().ToList();

            //var tasks = new List<Task<List<subsequence_classification_data>>>();
            var tasks = new List<Task<List<protein_subsequence_info>>>();

            foreach (var loop_chain_id in chain_ids)
            {
                var chain_id = loop_chain_id;

                var task = Task.Run(() =>
                {
                    //var pdb_chain_coils = new List<subsequence_classification_data>();
                    var pdb_chain_coils = new List<protein_subsequence_info>();
                    
                    var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
                    var pdb_chain_master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_chain_atoms);
                    var pdb_chain_aa_sequence = string.Join("", pdb_chain_master_atoms.Select(a => a.amino_acid).ToList());

                    var coil_subsequence_atoms = new List<atom>();
                    var coil_subsequence_master_atoms = new List<atom>();


                    for (var pdb_chain_master_atoms_index = 0; pdb_chain_master_atoms_index < pdb_chain_master_atoms.Count; pdb_chain_master_atoms_index++)
                    {
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
                                    symmetry_mode = "",
                                    parallelism = "",
                                    class_id = class_id,
                                    class_name = class_name,
                                    pdb_id = pdb_id,
                                    chain_id = chain_id,
                                    res_ids = coil_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
                                    aa_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                                    //aa_before = "",
                                    //aa_after = "",
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
                                //    aa_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                                //    dssp_multimer_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.multimer_dssp).ToList()),
                                //    dssp_monomer_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.monomer_dssp).ToList()),
                                //    stride_multimer_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.stride_multimer).ToList()),
                                //    stride_monomer_subsequence = string.Join("", coil_subsequence_master_atoms.Select(a => a.stride_monomer).ToList()),
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

                    return pdb_chain_coils;
                });

                tasks.Add(task);
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), nameof(dataset_gen_coils), nameof(find_coils));


            var pdb_coils = tasks.SelectMany(a => a.Result).ToList();

            pdb_coils = pdb_coils.GroupBy(scd => (scd.pdb_id, scd.aa_subsequence)).Select(group => group.First()).Distinct().ToList();

            return pdb_coils;
        }
    }
}
