using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class subsequence_classification_data_region
    {
        internal string pdb_id;
        internal char chain_id;
        internal List<(int residue_index, char i_code, char amino_acid)> res_ids;
        internal string aa_sequence;
        internal string dssp_monomer;
        internal string stride_monomer;
        internal string dssp_multimer;
        internal string stride_multimer;
        internal string dssp3_monomer;
        internal string stride3_monomer;
        internal string dssp3_multimer;
        internal string stride3_multimer;
        internal List<atom> atoms;
        internal List<atom> master_atoms;
        internal List<(string format, string prediction)> ss_predictions;
        internal foldx_energy_differences foldx_energy_differences;

        internal string unique_id()
        {
            return string.Join(/*program.string_debug*/($@"_"), new string[]
            {
                /*program.string_debug*/($@"{pdb_id}"),
                /*program.string_debug*/($@"{chain_id}"),
                string.Join(/*program.string_debug*/($@"_"), res_ids?
                    .Select(a => 
                        /*program.string_debug*/($@"{a.amino_acid}{a.residue_index}{(a.i_code != default && !char.IsWhiteSpace(a.i_code) ? /*program.string_debug*/($@"a.i_code") : /*program.string_debug*/($@""))}")
                    ).ToArray() ?? Array.Empty<string>())
            });
        }
        
        internal subsequence_classification_data_region(subsequence_classification_data scd, List<atom> region_atoms, bool load_ss_predictions = true, bool load_foldx_energy = true)
        {
            pdb_id = scd.pdb_id;
            chain_id = scd.chain_id;

            atoms = region_atoms;
            master_atoms = atom.select_amino_acid_master_atoms(null, atoms);
            res_ids = master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();
            aa_sequence = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.amino_acid).ToList());

            dssp_multimer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp_multimer).ToList());
            stride_multimer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride_multimer).ToList());
            dssp_monomer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp_monomer).ToList());
            stride_monomer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride_monomer).ToList());

            dssp3_multimer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp3_multimer).ToList());
            stride3_multimer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride3_multimer).ToList());
            dssp3_monomer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp3_monomer).ToList());
            stride3_monomer = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride3_monomer).ToList());

            if (load_ss_predictions)
            {
                if (master_atoms == null || master_atoms.Count == 0)
                {
                    var template_region = scd.get_regions().First(a => a.region != null && a.region.master_atoms != null && a.region.master_atoms.Count > 0).region.master_atoms;
                    ss_predictions = atom.get_dssp_and_mpsa_subsequences(template_region);
                    ss_predictions = ss_predictions.Select(a => (a.format, /*program.string_debug*/($@""))).ToList();
                }
                else
                {
                    ss_predictions = atom.get_dssp_and_mpsa_subsequences(master_atoms);
                }
            }

            if (load_foldx_energy)
            {
                foldx_energy_differences = info_foldx.load_calc_energy_differences(scd.pdb_id, scd.chain_id, res_ids, false);
            }
        }
    }
}
