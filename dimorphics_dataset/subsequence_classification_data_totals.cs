using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal static class subsequence_classification_data_totals
    {
        internal static int _peptides_call_count;
        internal static readonly object _peptides_call_count_lock = new object();
        internal static readonly object _peptides_cache_lock = new object();
        internal static readonly List<(string sequence, List<feature_info> features)> _peptides_cache = new List<(string sequence, List<feature_info> features)>();
        internal static int _protr_call_count;
        internal static readonly object _protr_call_count_lock = new object();
        internal static readonly object _protr_cache_lock = new object();
        internal static readonly List<(string sequence, List<feature_info> features)> _protr_cache = new List<(string sequence, List<feature_info> features)>();


        // 1d
        internal static int total_subsequence_1d_classification_data = -1;
        internal static int total_neighbourhood_1d_classification_data = -1;
        internal static int total_protein_1d_classification_data = -1;

        // 2d
        internal static int total_subsequence_2d_classification_data = -1;
        internal static int total_neighbourhood_2d_classification_data = -1;
        internal static int total_protein_2d_classification_data = -1;

        // 3d
        internal static int total_subsequence_3d_classification_data = -1;
        internal static int total_neighbourhood_3d_classification_data = -1;
        internal static int total_protein_3d_classification_data = -1;


        internal static int total_pse_aac_sequence_classification_data = -1;
        internal static int total_sable_sequence_classification_data = -1;
        internal static int total_mpsa_classification_data = -1;

        internal static int total_blast_pssm_subsequence_classification_data = -1;

        //internal static int total_blast_pssm_protein_classification_data = -1;
        internal static int total_aa_index_classification_data = -1;
        internal static int total_sequence_geometry_classification_data = -1;
        internal static int total_intrinsically_unordered_data = -1;
        internal static int total_dna_binding_prediction_data = -1;
        internal static int total_r_peptides_prediction_data = -1;
        internal static int total_r_protr_prediction_data = -1;

        //internal static int total_protein_pse_ssc_dssp_classification_data = -1;
        internal static int total_pse_ssc_dssp_classification_data = -1;
        internal static int total_protein_dssp_dist_classification_data = -1;

        internal static int total_foldx_classification_subsequence_3d_data = -1;
        internal static int total_foldx_classification_neighbourhood_3d_data = -1;
        internal static int total_foldx_classification_protein_3d_data = -1;

        internal static int total_ring_classification_data = -1;
        internal static int total_sasa_classification_data = -1;
        internal static int total_tortuosity_classification_data = -1;
        internal static int total_intramolecular_classification_data = -1;
        internal static int total_atom_distance_classification_data = -1;
        internal static int total_aa_aa_distances_classification_data = -1;

        
        internal static int total_atom_distances_classification_data = -1;
    }
}
