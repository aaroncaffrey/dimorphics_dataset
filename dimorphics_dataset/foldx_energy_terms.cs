using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal class foldx_energy_terms
    {
        internal string pdb_id;
        internal char chain_id;
        internal List<(int residue_index, char i_code, char amino_acid)> res_ids;

        internal int line_index;

        //internal string pdb_filename;
        //internal string wait_filename;
        //internal bool repaired;
        internal string Pdb;
        internal double SD;
        internal double total_energy;
        internal double Backbone_Hbond;
        internal double Sidechain_Hbond;
        internal double Van_der_Waals;
        internal double Electrostatics;
        internal double Solvation_Polar;
        internal double Solvation_Hydrophobic;
        internal double Van_der_Waals_clashes;
        internal double entropy_sidechain;
        internal double entropy_mainchain;
        internal double sloop_entropy;
        internal double mloop_entropy;
        internal double cis_bond;
        internal double torsional_clash;
        internal double backbone_clash;
        internal double helix_dipole;
        internal double water_bridge;
        internal double disulfide;
        internal double electrostatic_kon;
        internal double partial_covalent_bonds;
        internal double energy_Ionisation;
        internal double Entropy_Complex;

        internal (string name, double value)[] properties()
        {
            return new (string name, double value)[]
                {
                        (nameof(total_energy),              total_energy),
                        (nameof(Backbone_Hbond),            Backbone_Hbond),
                        (nameof(Sidechain_Hbond),           Sidechain_Hbond),
                        (nameof(Van_der_Waals),             Van_der_Waals),
                        (nameof(Electrostatics),            Electrostatics),
                        (nameof(Solvation_Polar),           Solvation_Polar),
                        (nameof(Solvation_Hydrophobic),     Solvation_Hydrophobic),
                        (nameof(Van_der_Waals_clashes),     Van_der_Waals_clashes),
                        (nameof(entropy_sidechain),         entropy_sidechain),
                        (nameof(entropy_mainchain),         entropy_mainchain),
                        (nameof(sloop_entropy),             sloop_entropy),
                        (nameof(mloop_entropy),             mloop_entropy),
                        (nameof(cis_bond),                  cis_bond),
                        (nameof(torsional_clash),           torsional_clash),
                        (nameof(backbone_clash),            backbone_clash),
                        (nameof(helix_dipole),              helix_dipole),
                        (nameof(water_bridge),              water_bridge),
                        (nameof(disulfide),                 disulfide),
                        (nameof(electrostatic_kon),         electrostatic_kon),
                        (nameof(partial_covalent_bonds),    partial_covalent_bonds),
                        (nameof(energy_Ionisation),         energy_Ionisation),
                        (nameof(Entropy_Complex),           Entropy_Complex),
                };
        }
    }
}
