using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal static class subsequence_classification_data_templates
    {
        internal static readonly List<(string name, List<info_aaindex.info_aaindex_entry> list)> _aaindex_subset_templates = subsequence_classification_data_methods.aaindex_subset_templates_search();

        internal static subsequence_classification_data _template_scd;

        internal static List<feature_info> _calculate_sable_classification_data_template;
        internal static List<feature_info> _calculate_mpsa_classification_data_template;
        internal static List<feature_info> _calculate_ring_classification_data_template;
        internal static List<feature_info> _calculate_foldx_classification_data_subsequence_3d_template;
        internal static List<feature_info> _calculate_foldx_classification_data_neighbourhood_2d_template;
        internal static List<feature_info> _calculate_foldx_classification_data_neighbourhood_3d_template;
        internal static List<feature_info> _calculate_foldx_classification_data_protein_3d_template;
        internal static List<feature_info> _calculate_sequence_geometry_classification_data_template;
        internal static List<feature_info> _calculate_aa_index_classification_data_template;
        internal static List<feature_info> _calculate_dna_binding_prediction_data_template;
        internal static List<feature_info> _calculate_intrinsically_unordered_data_template;
        internal static List<feature_info> _calculate_blast_pssm_classification_data_template;
        internal static List<feature_info> _calculate_sasa_classification_data_template;
        internal static List<feature_info> _calculate_tortuosity_classification_data_template;
        internal static List<feature_info> _calculate_atom_distances_data_template;
        internal static List<feature_info> _calculate_aa_or_ss_sequence_classification_data_aa_template;
        internal static List<feature_info> _calculate_aa_or_ss_sequence_classification_data_ss_template;
        internal static List<feature_info> _peptides_data_template;
        internal static List<feature_info> _protr_data_template;
    }
}
