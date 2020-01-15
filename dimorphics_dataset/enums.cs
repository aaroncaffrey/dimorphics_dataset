using System;

namespace dimorphics_dataset
{
    public enum enum_substructure_type
    {
        dimorphic,
        dimorphic_coil,
        dimorphic_full_protein_sequence,

        standard_strand,
        standard_strand_coil,
        standard_strand_full_protein_sequence,

        standard_coil,
    }

    public enum enum_protein_data_source
    {
        subsequence_1d,
        neighbourhood_1d,
        protein_1d,
        subsequence_3d,
        neighbourhood_3d,
        protein_3d
    }

    public enum enum_pssm_normalisation_method
    {

        norm_none, norm_whole_pssm, norm_subseq, norm_encoded_parts, norm_encoded_vector
    }


    public enum enum_pssm_value_type
    {
        standard,
        distances,
        intervals
    }

    public enum enum_ss_type
        { 
            DSSP3, 
            DSSP,
            STRIDE3,
            STRIDE
        }

        [Flags]
        public enum enum_get_dssp_and_mpsa_subsequences_params
        {
            none = 0b_0000_0000_0000,

            aa_seq = 0b_0000_0000_0001,

            monomer_dssp_seq = 0b_0000_0000_0010,
            monomer_stride_seq = 0b_0000_0000_0100,
            monomer_dssp3_seq = 0b_0000_0000_1000,
            monomer_stride3_seq = 0b_0000_0001_0000,

            multimer_dssp_seq = 0b_0000_0010_0000,
            multimer_stride_seq = 0b_0000_0100_0000,
            multimer_dssp3_seq = 0b_0000_1000_0000,
            multimer_stride3_seq = 0b_0001_0000_0000,
        }

        public enum enum_seq_type
        {
            amino_acid_sequence,
            secondary_structure_sequence
        }

        public enum enum_structure_oligomisation
        {
            monomer,
            multimer
        }


}
