namespace dimorphics_dataset
{
    internal class load_atoms_pdb_options
    {
        internal bool first_model_only = true;
        internal bool first_icode_only = true;

        internal bool load_uniprot;

        // 1d data
        internal bool load_1d_blast_pssms = true;
        internal bool load_1d_iup_data = true;
        internal bool load_1d_sable = true;
        internal bool load_1d_dna_binding = true;

        // 2d data
        internal bool load_2d_mpsa_sec_struct_predictions = true;

        // 3d data
        internal bool find_3d_intramolecular = true;
        //internal bool find_3d_intermolecular = true;
        internal bool load_3d_rsa_data = true;
        internal bool load_3d_dssp_data = true;
        internal bool load_3d_stride_data = true;
        internal bool load_3d_ring_data = true;
        internal bool load_3d_foldx_ala_scan = true;
    }
}
