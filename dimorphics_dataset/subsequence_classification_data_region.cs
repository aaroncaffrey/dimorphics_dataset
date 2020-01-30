using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Sockets;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using Newtonsoft.Json;

namespace dimorphics_dataset
{
    public class subsequence_classification_data_region
    {
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
        internal info_foldx.energy_differences foldx_energy_differences;

        public subsequence_classification_data_region(subsequence_classification_data scd, List<atom> region_atoms)
        {
            atoms = region_atoms;
            master_atoms = atom.select_amino_acid_master_atoms(null, atoms);
            res_ids = master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList();
            aa_sequence = string.Join("", master_atoms.Select(a => a.amino_acid).ToList());
            
            dssp_multimer = string.Join("", master_atoms.Select(a => a.dssp_multimer).ToList());
            stride_multimer = string.Join("", master_atoms.Select(a => a.stride_multimer).ToList());
            dssp_monomer = string.Join("", master_atoms.Select(a => a.dssp_monomer).ToList());
            stride_monomer = string.Join("", master_atoms.Select(a => a.stride_monomer).ToList());

            dssp3_multimer = string.Join("", master_atoms.Select(a => a.dssp3_multimer).ToList());
            stride3_multimer = string.Join("", master_atoms.Select(a => a.stride3_multimer).ToList());
            dssp3_monomer = string.Join("", master_atoms.Select(a => a.dssp3_monomer).ToList());
            stride3_monomer = string.Join("", master_atoms.Select(a => a.stride3_monomer).ToList());
        
            ss_predictions = null;
            foldx_energy_differences = info_foldx.load_calc_energy_differences(scd.pdb_id, scd.chain_id, res_ids, false);
        }
    }
}
