using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal class descriptive_stats_encoding_options
    {
        // defaults:
        // internal static readonly descriptive_stats_encoding_options options_all = new descriptive_stats_encoding_options("all", true);
        internal static readonly descriptive_stats_encoding_options options_default = new descriptive_stats_encoding_options("mean", false)
        {
            count = false,
            count_zero_values = false,
            count_non_zero_values = false,
            count_distinct_values = false,
            sum = false,
            mean_arithmetic = true, //
            mean_geometric = false,
            mean_harmonic = false,
            min = false,
            max = false,
            range = false,
            mid_range = false,
            variance = false,
            dev_standard = false,
            root_mean_square = false,
            skewness = false,
            kurtosis = false,
            interquartile_range = false,
            median_q1 = false,
            median_q2 = false,
            median_q3 = false,
            mad_mean_arithmetic = false,
            mad_mean_harmonic = false,
            mad_mean_geometric = false,
            mad_median_q1 = false,
            mad_median_q2 = false,
            mad_median_q3 = false,
            mad_mid_range = false,
        };

        internal const bool intervals_default = false;
        internal const bool distances_default = false;
        internal const bool interquartile_default = false;
        internal const bool abs_default = false;
        internal const bool rescale_default = false;

        // feature parameters:
        internal static descriptive_stats_encoding_options[] dse_options_aa_index = new[] { options_default };
        internal static bool dse_intervals_aa_index = intervals_default;
        internal static bool dse_distances_aa_index = distances_default;
        internal static bool dse_interquartile_aa_index = interquartile_default;
        internal static bool dse_abs_aa_index = abs_default;
        internal static bool dse_rescale_aa_index = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_entropy = new[] { options_default };
        internal static bool dse_intervals_sable_ds_entropy = intervals_default;
        internal static bool dse_distances_sable_ds_entropy = distances_default;
        internal static bool dse_interquartile_sable_ds_entropy = interquartile_default;
        internal static bool dse_abs_sable_ds_entropy = abs_default;
        internal static bool dse_rescale_sable_ds_entropy = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_burial_abs = new[] { options_default };
        internal static bool dse_intervals_sable_ds_burial_abs = intervals_default;
        internal static bool dse_distances_sable_ds_burial_abs = distances_default;
        internal static bool dse_interquartile_sable_ds_burial_abs = interquartile_default;
        internal static bool dse_abs_sable_ds_burial_abs = abs_default;
        internal static bool dse_rescale_sable_ds_burial_abs = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_burial_rel = new[] { options_default };
        internal static bool dse_intervals_sable_ds_burial_rel = intervals_default;
        internal static bool dse_distances_sable_ds_burial_rel = distances_default;
        internal static bool dse_interquartile_sable_ds_burial_rel = interquartile_default;
        internal static bool dse_abs_sable_ds_burial_rel = abs_default;
        internal static bool dse_rescale_sable_ds_burial_rel = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_distances = new[] { options_default };
        internal static bool dse_intervals_ring_distances = intervals_default;
        internal static bool dse_distances_ring_distances = distances_default;
        internal static bool dse_interquartile_ring_distances = interquartile_default;
        internal static bool dse_abs_ring_distances = abs_default;
        internal static bool dse_rescale_ring_distances = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_angles = new[] { options_default };
        internal static bool dse_intervals_ring_angles = intervals_default;
        internal static bool dse_distances_ring_angles = distances_default;
        internal static bool dse_interquartile_ring_angles = interquartile_default;
        internal static bool dse_abs_ring_angles = abs_default;
        internal static bool dse_rescale_ring_angles = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_energies = new[] { options_default };
        internal static bool dse_intervals_ring_energies = intervals_default;
        internal static bool dse_distances_ring_energies = distances_default;
        internal static bool dse_interquartile_ring_energies = interquartile_default;
        internal static bool dse_abs_ring_energies = abs_default;
        internal static bool dse_rescale_ring_energies = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_degrees = new[] { options_default };
        internal static bool dse_intervals_ring_degrees = intervals_default;
        internal static bool dse_distances_ring_degrees = distances_default;
        internal static bool dse_interquartile_ring_degrees = interquartile_default;
        internal static bool dse_abs_ring_degrees = abs_default;
        internal static bool dse_rescale_ring_degrees = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_rapdf = new[] { options_default };
        internal static bool dse_intervals_ring_rapdf = intervals_default;
        internal static bool dse_distances_ring_rapdf = distances_default;
        internal static bool dse_interquartile_ring_rapdf = interquartile_default;
        internal static bool dse_abs_ring_rapdf = abs_default;
        internal static bool dse_rescale_ring_rapdf = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_atom_distances = new[] { options_default };
        internal static bool dse_intervals_atom_distances = intervals_default;
        internal static bool dse_distances_atom_distances = distances_default;
        internal static bool dse_interquartile_atom_distances = interquartile_default;
        internal static bool dse_abs_atom_distances = abs_default;
        internal static bool dse_rescale_atom_distances = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_curves = new[] { options_default };
        internal static bool dse_intervals_tortuosity2_curves = intervals_default;
        internal static bool dse_distances_tortuosity2_curves = distances_default;
        internal static bool dse_interquartile_tortuosity2_curves = interquartile_default;
        internal static bool dse_abs_tortuosity2_curves = abs_default;
        internal static bool dse_rescale_tortuosity2_curves = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_displacement = new[] { options_default };
        internal static bool dse_intervals_tortuosity2_displacement = intervals_default;
        internal static bool dse_distances_tortuosity2_displacement = distances_default;
        internal static bool dse_interquartile_tortuosity2_displacement = interquartile_default;
        internal static bool dse_abs_tortuosity2_displacement = abs_default;
        internal static bool dse_rescale_tortuosity2_displacement = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_stat_values = new[] { options_default };
        internal static bool dse_intervals_tortuosity2_stat_values = intervals_default;
        internal static bool dse_distances_tortuosity2_stat_values = distances_default;
        internal static bool dse_interquartile_tortuosity2_stat_values = interquartile_default;
        internal static bool dse_abs_tortuosity2_stat_values = abs_default;
        internal static bool dse_rescale_tortuosity2_stat_values = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_sasa_all = new[] { options_default };
        internal static bool dse_intervals_sasa_all = intervals_default;
        internal static bool dse_distances_sasa_all = distances_default;
        internal static bool dse_interquartile_sasa_all = interquartile_default;
        internal static bool dse_abs_sasa_all = abs_default;
        internal static bool dse_rescale_sasa_all = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_ring_tap = new[] { options_default };
        internal static bool dse_intervals_ring_tap = intervals_default;
        internal static bool dse_distances_ring_tap = distances_default;
        internal static bool dse_interquartile_ring_tap = interquartile_default;
        internal static bool dse_abs_ring_tap = abs_default;
        internal static bool dse_rescale_ring_tap = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_iud_short = new[] { options_default };
        internal static bool dse_intervals_iud_short = intervals_default;
        internal static bool dse_distances_iud_short = distances_default;
        internal static bool dse_interquartile_iud_short = interquartile_default;
        internal static bool dse_abs_iud_short = abs_default;
        internal static bool dse_rescale_iud_short = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_iud_long = new[] { options_default };
        internal static bool dse_intervals_iud_long = intervals_default;
        internal static bool dse_distances_iud_long = distances_default;
        internal static bool dse_interquartile_iud_long = interquartile_default;
        internal static bool dse_abs_iud_long = abs_default;
        internal static bool dse_rescale_iud_long = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_iud_glob = new[] { options_default };
        internal static bool dse_intervals_iud_glob = intervals_default;
        internal static bool dse_distances_iud_glob = distances_default;
        internal static bool dse_interquartile_iud_glob = interquartile_default;
        internal static bool dse_abs_iud_glob = abs_default;
        internal static bool dse_rescale_iud_glob = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_anchor2 = new[] { options_default };
        internal static bool dse_intervals_anchor2 = intervals_default;
        internal static bool dse_distances_anchor2 = distances_default;
        internal static bool dse_interquartile_anchor2 = interquartile_default;
        internal static bool dse_abs_anchor2 = abs_default;
        internal static bool dse_rescale_anchor2 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_ala = new[] { options_default };
        internal static bool dse_intervals_foldx_ala = intervals_default;
        internal static bool dse_distances_foldx_ala = distances_default;
        internal static bool dse_interquartile_foldx_ala = interquartile_default;
        internal static bool dse_abs_foldx_ala = abs_default;
        internal static bool dse_rescale_foldx_ala = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan1 = new[] { options_default };
        internal static bool dse_intervals_foldx_pos_scan1 = intervals_default;
        internal static bool dse_distances_foldx_pos_scan1 = distances_default;
        internal static bool dse_interquartile_foldx_pos_scan1 = interquartile_default;
        internal static bool dse_abs_foldx_pos_scan1 = abs_default;
        internal static bool dse_rescale_foldx_pos_scan1 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan2 = new[] { options_default };
        internal static bool dse_intervals_foldx_pos_scan2 = intervals_default;
        internal static bool dse_distances_foldx_pos_scan2 = distances_default;
        internal static bool dse_interquartile_foldx_pos_scan2 = interquartile_default;
        internal static bool dse_abs_foldx_pos_scan2 = abs_default;
        internal static bool dse_rescale_foldx_pos_scan2 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan3 = new[] { options_default };
        internal static bool dse_intervals_foldx_pos_scan3 = intervals_default;
        internal static bool dse_distances_foldx_pos_scan3 = distances_default;
        internal static bool dse_interquartile_foldx_pos_scan3 = interquartile_default;
        internal static bool dse_abs_foldx_pos_scan3 = abs_default;
        internal static bool dse_rescale_foldx_pos_scan3 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan1 = new[] { options_default };
        internal static bool dse_intervals_foldx_bm_pos_scan1 = intervals_default;
        internal static bool dse_distances_foldx_bm_pos_scan1 = distances_default;
        internal static bool dse_interquartile_foldx_bm_pos_scan1 = interquartile_default;
        internal static bool dse_abs_foldx_bm_pos_scan1 = abs_default;
        internal static bool dse_rescale_foldx_bm_pos_scan1 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan2 = new[] { options_default };
        internal static bool dse_intervals_foldx_bm_pos_scan2 = intervals_default;
        internal static bool dse_distances_foldx_bm_pos_scan2 = distances_default;
        internal static bool dse_interquartile_foldx_bm_pos_scan2 = interquartile_default;
        internal static bool dse_abs_foldx_bm_pos_scan2 = abs_default;
        internal static bool dse_rescale_foldx_bm_pos_scan2 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan3 = new[] { options_default };
        internal static bool dse_intervals_foldx_bm_pos_scan3 = intervals_default;
        internal static bool dse_distances_foldx_bm_pos_scan3 = distances_default;
        internal static bool dse_interquartile_foldx_bm_pos_scan3 = interquartile_default;
        internal static bool dse_abs_foldx_bm_pos_scan3 = abs_default;
        internal static bool dse_rescale_foldx_bm_pos_scan3 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_sr1 = new[] { options_default };
        internal static bool dse_intervals_foldx_bm_sr1 = intervals_default;
        internal static bool dse_distances_foldx_bm_sr1 = distances_default;
        internal static bool dse_interquartile_foldx_bm_sr1 = interquartile_default;
        internal static bool dse_abs_foldx_bm_sr1 = abs_default;
        internal static bool dse_rescale_foldx_bm_sr1 = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm1_ds = new[] { options_default };
        internal static bool dse_intervals_pssm1_ds = intervals_default;
        internal static bool dse_distances_pssm1_ds = distances_default;
        internal static bool dse_interquartile_pssm1_ds = interquartile_default;
        internal static bool dse_abs_pssm1_ds = abs_default;
        internal static bool dse_rescale_pssm1_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm20col_ds = new[] { options_default };
        internal static bool dse_intervals_pssm20col_ds = intervals_default;
        internal static bool dse_distances_pssm20col_ds = distances_default;
        internal static bool dse_interquartile_pssm20col_ds = interquartile_default;
        internal static bool dse_abs_pssm20col_ds = abs_default;
        internal static bool dse_rescale_pssm20col_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm20row_ds = new[] { options_default };
        internal static bool dse_intervals_pssm20row_ds = intervals_default;
        internal static bool dse_distances_pssm20row_ds = distances_default;
        internal static bool dse_interquartile_pssm20row_ds = interquartile_default;
        internal static bool dse_abs_pssm20row_ds = abs_default;
        internal static bool dse_rescale_pssm20row_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm210_ds = new[] { options_default };
        internal static bool dse_intervals_pssm210_ds = intervals_default;
        internal static bool dse_distances_pssm210_ds = distances_default;
        internal static bool dse_interquartile_pssm210_ds = interquartile_default;
        internal static bool dse_abs_pssm210_ds = abs_default;
        internal static bool dse_rescale_pssm210_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm400_ds = new[] { options_default };
        internal static bool dse_intervals_pssm400_ds = intervals_default;
        internal static bool dse_distances_pssm400_ds = distances_default;
        internal static bool dse_interquartile_pssm400_ds = interquartile_default;
        internal static bool dse_abs_pssm400_ds = abs_default;
        internal static bool dse_rescale_pssm400_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm20colDT_ds = new[] { options_default };
        internal static bool dse_intervals_pssm20colDT_ds = intervals_default;
        internal static bool dse_distances_pssm20colDT_ds = distances_default;
        internal static bool dse_interquartile_pssm20colDT_ds = interquartile_default;
        internal static bool dse_abs_pssm20colDT_ds = abs_default;
        internal static bool dse_rescale_pssm20colDT_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm210DT_ds = new[] { options_default };
        internal static bool dse_intervals_pssm210DT_ds = intervals_default;
        internal static bool dse_distances_pssm210DT_ds = distances_default;
        internal static bool dse_interquartile_pssm210DT_ds = interquartile_default;
        internal static bool dse_abs_pssm210DT_ds = abs_default;
        internal static bool dse_rescale_pssm210DT_ds = rescale_default;

        internal static descriptive_stats_encoding_options[] dse_options_pssm400DT_ds = new[] { options_default };
        internal static bool dse_intervals_pssm400DT_ds = intervals_default;
        internal static bool dse_distances_pssm400DT_ds = distances_default;
        internal static bool dse_interquartile_pssm400DT_ds = interquartile_default;
        internal static bool dse_abs_pssm400DT_ds = abs_default;
        internal static bool dse_rescale_pssm400DT_ds = rescale_default;





        internal string options_name = $@"";

        internal bool count;
        internal bool count_zero_values;
        internal bool count_non_zero_values;
        internal bool count_distinct_values;
        internal bool sum;
        internal bool mean_arithmetic;
        internal bool mean_geometric;
        internal bool mean_harmonic;
        internal bool min;
        internal bool max;
        internal bool range;
        internal bool mid_range;
        internal bool variance;
        internal bool dev_standard;
        internal bool root_mean_square;
        internal bool skewness;
        internal bool kurtosis;
        internal bool interquartile_range;
        internal bool median_q1;
        internal bool median_q2;
        internal bool median_q3;
        internal bool mad_mean_arithmetic;
        internal bool mad_mean_harmonic;
        internal bool mad_mean_geometric;
        internal bool mad_median_q1;
        internal bool mad_median_q2;
        internal bool mad_median_q3;
        internal bool mad_mid_range;

        internal descriptive_stats_encoding_options(string name, bool enable = false)
        {
            options_name = name;

            count = enable;
            count_zero_values = enable;
            count_non_zero_values = enable;
            count_distinct_values = enable;
            sum = enable;
            mean_arithmetic = enable;
            mean_geometric = enable;
            mean_harmonic = enable;
            min = enable;
            max = enable;
            range = enable;
            mid_range = enable;
            variance = enable;
            dev_standard = enable;
            root_mean_square = enable;
            skewness = enable;
            kurtosis = enable;
            interquartile_range = enable;
            median_q1 = enable;
            median_q2 = enable;
            median_q3 = enable;
            mad_mean_arithmetic = enable;
            mad_mean_harmonic = enable;
            mad_mean_geometric = enable;
            mad_median_q1 = enable;
            mad_median_q2 = enable;
            mad_median_q3 = enable;
            mad_mid_range = enable;
        }

        //internal static descriptive_stats_encoding_options options_average()
        //{
        //    return new descriptive_stats_encoding_options("mean", false)
        //    {
        //        mean_arithmetic = true,
        //    };
        //}

        //internal static descriptive_stats_encoding_options options_average_sd()
        //{
        //    return new descriptive_stats_encoding_options("mean", false)
        //    {
        //        mean_arithmetic = true,
        //        dev_standard = true,
        //    };
        //}
    }
}
