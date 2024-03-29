﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class descriptive_stats_encoding_options
    {
        // defaults:
        // internal static readonly descriptive_stats_encoding_options options_all = new descriptive_stats_encoding_options(/*program.string_debug*/($@"all", true);
        internal static readonly descriptive_stats_encoding_options options_mean_arithmetic = new descriptive_stats_encoding_options(nameof(options_mean_arithmetic), false)
        {
            mean_arithmetic = true,
            dev_standard = true
        };

        //internal static readonly descriptive_stats_encoding_options options_mean_geometric = new descriptive_stats_encoding_options(nameof(options_mean_geometric), false)
        //{
            //mean_geometric_corrected = true,
            //mad_mean_geometric_corrected = true
            //mean_geometric_nonzero = true
        //};

        //internal static readonly descriptive_stats_encoding_options options_mean_harmonic = new descriptive_stats_encoding_options(nameof(options_mean_harmonic), false)
        //{
            //mean_harmonic_corrected = true,
            //mad_mean_harmonic_corrected = true
            //mean_harmonic_nonzero = true
        //};

        internal static readonly descriptive_stats_encoding_options options_median_q2 = new descriptive_stats_encoding_options(nameof(options_median_q2), false)
        {
            median_q2 = true,
            mad_median_q2 = true
            //mean_harmonic_nonzero = true
        };


        internal static readonly descriptive_stats_encoding_options options_all = new descriptive_stats_encoding_options(nameof(options_all), true)
        {
            count = true,
            count_zero_values = true,
            count_non_zero_values = true,
            count_distinct_values = true,
            sum = true,
            mean_arithmetic = true,
            //mean_geometric_corrected = true,
            //mean_geometric_nonzero = true,
            //mean_harmonic_corrected = true,
            //mean_harmonic_nonzero = true,
            min = true,
            max = true,
            range = true,
            mid_range = true,
            variance = true,
            dev_standard = true,
            root_mean_square = true,
            skewness = true,
            kurtosis = true,
            interquartile_range = true,
            median_q1 = true,
            median_q2 = true,
            median_q3 = true,
            mad_mean_arithmetic = true,
            //mad_mean_harmonic_corrected = true,
            //mad_mean_geometric_corrected = true,
            mad_median_q1 = true,
            mad_median_q2 = true,
            mad_median_q3 = true,
            mad_mid_range = true,

            intervals = null,
            distances = null,
            interquartile = null,
            abs = null,
            rescale = null
        };

        internal static readonly descriptive_stats_encoding_options options_all_plus = new descriptive_stats_encoding_options(options_all,nameof(options_all_plus), true)
        {
            intervals = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(intervals)}")),
            distances = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(distances)}")),
            interquartile = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(interquartile)}")),
            abs = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(abs)}")),
            rescale = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(rescale)}")),
        };

        //internal static readonly descriptive_stats_encoding_options options_test = new descriptive_stats_encoding_options(options_all, nameof(options_test), false)
        //{
        //    intervals = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(intervals)}")),
        //    //distances = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(distances)}")),
        //    //interquartile = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(interquartile)}")),
        //    //abs = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(abs)}")),
        //    //rescale = new descriptive_stats_encoding_options(options_all, /*program.string_debug*/($@"all_plus_{nameof(rescale)}")),
        //};



        // feature parameters:
        internal static descriptive_stats_encoding_options[] dse_options_aa_index = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_entropy = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_burial_abs = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_sable_ds_burial_rel = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_distances = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_angles = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_energies = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_degrees = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_rapdf = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_atom_distances = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_curves = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_displacement = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_tortuosity2_stat_values = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_sasa_all = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_ring_tap = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_iud_short = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_iud_long = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_iud_glob = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_anchor2 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_ala = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan1 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan2 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_pos_scan3 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan1 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan2 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_pos_scan3 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_foldx_bm_sr1 = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm1_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm20col_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm20row_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm210_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm400_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm20colDT_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm210DT_ds = new[] { new descriptive_stats_encoding_options(options_all) };
        internal static descriptive_stats_encoding_options[] dse_options_pssm400DT_ds = new[] { new descriptive_stats_encoding_options(options_all) };



        internal string options_name;

        internal descriptive_stats_encoding_options intervals;
        internal descriptive_stats_encoding_options distances;
        internal descriptive_stats_encoding_options interquartile;
        internal descriptive_stats_encoding_options abs;
        internal descriptive_stats_encoding_options rescale;

        internal bool count;
        internal bool count_zero_values;
        internal bool count_non_zero_values;
        internal bool count_distinct_values;
        internal bool sum;
        internal bool mean_arithmetic;
        //internal bool mean_geometric_corrected;
        //internal bool mean_geometric_nonzero;
        //internal bool mean_harmonic_corrected;
        //internal bool mean_harmonic_nonzero;
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
        //internal bool mad_mean_harmonic_corrected;
        //internal bool mad_mean_geometric_corrected;
        internal bool mad_median_q1;
        internal bool mad_median_q2;
        internal bool mad_median_q3;
        internal bool mad_mid_range;

        internal descriptive_stats_encoding_options(string name, bool? enable = null, string[] to_enable = null, string[] to_disable = null)
        {
            if (name != null)
            {
                options_name = /*program.string_debug*/(name);
            }

            if (enable != null)
            {
                set_enable(enable.Value);
            }

            if (to_enable != null && to_enable.Length > 0)
            {
                set_enable(to_enable, true);
            }

            if (to_disable != null && to_disable.Length > 0)
            {
                set_enable(to_disable, false);
            }
        }

        internal void set_enable(bool enable)
        {
            this.count = enable;
            this.count_zero_values = enable;
            this.count_non_zero_values = enable;
            this.count_distinct_values = enable;
            this.sum = enable;
            this.mean_arithmetic = enable;
            //this.mean_geometric_corrected = enable;
            //this.mean_geometric_nonzero = enable;
            //this.mean_harmonic_corrected = enable;
            //this.mean_harmonic_nonzero = enable;
            this.min = enable;
            this.max = enable;
            this.range = enable;
            this.mid_range = enable;
            this.variance = enable;
            this.dev_standard = enable;
            this.root_mean_square = enable;
            this.skewness = enable;
            this.kurtosis = enable;
            this.interquartile_range = enable;
            this.median_q1 = enable;
            this.median_q2 = enable;
            this.median_q3 = enable;
            this.mad_mean_arithmetic = enable;
            //this.mad_mean_harmonic_corrected = enable;
            //this.mad_mean_geometric_corrected = enable;
            this.mad_median_q1 = enable;
            this.mad_median_q2 = enable;
            this.mad_median_q3 = enable;
            this.mad_mid_range = enable;
        }


        internal void set_enable(string[] names, bool enable)
        {
            if (names == null || names.Length == 0)
            {
                return;
            }

            for (var i = 0; i < names.Length; i++)
            {
                if (string.Equals(names[i], nameof(count), StringComparison.OrdinalIgnoreCase)) count = enable;
                else if (string.Equals(names[i], nameof(count_zero_values), StringComparison.OrdinalIgnoreCase)) count_zero_values = enable;
                else if (string.Equals(names[i], nameof(count_non_zero_values), StringComparison.OrdinalIgnoreCase)) count_non_zero_values = enable;
                else if (string.Equals(names[i], nameof(count_distinct_values), StringComparison.OrdinalIgnoreCase)) count_distinct_values = enable;
                else if (string.Equals(names[i], nameof(sum), StringComparison.OrdinalIgnoreCase)) sum = enable;
                else if (string.Equals(names[i], nameof(mean_arithmetic), StringComparison.OrdinalIgnoreCase)) mean_arithmetic = enable;
                //else if (string.Equals(names[i], nameof(mean_geometric_corrected), StringComparison.OrdinalIgnoreCase)) mean_geometric_corrected = enable;
                //else if (string.Equals(names[i], nameof(mean_geometric_nonzero), StringComparison.OrdinalIgnoreCase)) mean_geometric_nonzero = enable;
                //else if (string.Equals(names[i], nameof(mean_harmonic_corrected), StringComparison.OrdinalIgnoreCase)) mean_harmonic_corrected = enable;
                //else if (string.Equals(names[i], nameof(mean_harmonic_nonzero), StringComparison.OrdinalIgnoreCase)) mean_harmonic_nonzero = enable;
                else if (string.Equals(names[i], nameof(min), StringComparison.OrdinalIgnoreCase)) min = enable;
                else if (string.Equals(names[i], nameof(max), StringComparison.OrdinalIgnoreCase)) max = enable;
                else if (string.Equals(names[i], nameof(range), StringComparison.OrdinalIgnoreCase)) range = enable;
                else if (string.Equals(names[i], nameof(mid_range), StringComparison.OrdinalIgnoreCase)) mid_range = enable;
                else if (string.Equals(names[i], nameof(variance), StringComparison.OrdinalIgnoreCase)) variance = enable;
                else if (string.Equals(names[i], nameof(dev_standard), StringComparison.OrdinalIgnoreCase)) dev_standard = enable;
                else if (string.Equals(names[i], nameof(root_mean_square), StringComparison.OrdinalIgnoreCase)) root_mean_square = enable;
                else if (string.Equals(names[i], nameof(skewness), StringComparison.OrdinalIgnoreCase)) skewness = enable;
                else if (string.Equals(names[i], nameof(kurtosis), StringComparison.OrdinalIgnoreCase)) kurtosis = enable;
                else if (string.Equals(names[i], nameof(interquartile_range), StringComparison.OrdinalIgnoreCase)) interquartile_range = enable;
                else if (string.Equals(names[i], nameof(median_q1), StringComparison.OrdinalIgnoreCase)) median_q1 = enable;
                else if (string.Equals(names[i], nameof(median_q2), StringComparison.OrdinalIgnoreCase)) median_q2 = enable;
                else if (string.Equals(names[i], nameof(median_q3), StringComparison.OrdinalIgnoreCase)) median_q3 = enable;
                else if (string.Equals(names[i], nameof(mad_mean_arithmetic), StringComparison.OrdinalIgnoreCase)) mad_mean_arithmetic = enable;
                //else if (string.Equals(names[i], nameof(mad_mean_harmonic_corrected), StringComparison.OrdinalIgnoreCase)) mad_mean_harmonic_corrected = enable;
                //else if (string.Equals(names[i], nameof(mad_mean_geometric_corrected), StringComparison.OrdinalIgnoreCase)) mad_mean_geometric_corrected = enable;
                else if (string.Equals(names[i], nameof(mad_median_q1), StringComparison.OrdinalIgnoreCase)) mad_median_q1 = enable;
                else if (string.Equals(names[i], nameof(mad_median_q2), StringComparison.OrdinalIgnoreCase)) mad_median_q2 = enable;
                else if (string.Equals(names[i], nameof(mad_median_q3), StringComparison.OrdinalIgnoreCase)) mad_median_q3 = enable;
                else if (string.Equals(names[i], nameof(mad_mid_range), StringComparison.OrdinalIgnoreCase)) mad_mid_range = enable;
            }
        }


        internal descriptive_stats_encoding_options(descriptive_stats_encoding_options dse_options, string name = null, bool? enable = null)
        {
            if (dse_options != null)
            {
                this.options_name = dse_options.options_name != null ? /*program.string_debug*/(dse_options.options_name) : null;

                this.intervals = dse_options.intervals != null ? new descriptive_stats_encoding_options(dse_options.intervals) : null;
                this.distances = dse_options.distances != null ? new descriptive_stats_encoding_options(dse_options.distances) : null;
                this.interquartile = dse_options.interquartile != null ? new descriptive_stats_encoding_options(dse_options.interquartile) : null;
                this.abs = dse_options.abs != null ? new descriptive_stats_encoding_options(dse_options.abs) : null;
                this.rescale = dse_options.rescale != null ? new descriptive_stats_encoding_options(dse_options.rescale) : null;

                this.count = dse_options.count;
                this.count_zero_values = dse_options.count_zero_values;
                this.count_non_zero_values = dse_options.count_non_zero_values;
                this.count_distinct_values = dse_options.count_distinct_values;
                this.sum = dse_options.sum;
                this.mean_arithmetic = dse_options.mean_arithmetic;
                //this.mean_geometric_corrected = dse_options.mean_geometric_corrected;
                //this.mean_geometric_nonzero = dse_options.mean_geometric_nonzero;
                //this.mean_harmonic_corrected = dse_options.mean_harmonic_corrected;
                //this.mean_harmonic_nonzero = dse_options.mean_harmonic_nonzero;
                this.min = dse_options.min;
                this.max = dse_options.max;
                this.range = dse_options.range;
                this.mid_range = dse_options.mid_range;
                this.variance = dse_options.variance;
                this.dev_standard = dse_options.dev_standard;
                this.root_mean_square = dse_options.root_mean_square;
                this.skewness = dse_options.skewness;
                this.kurtosis = dse_options.kurtosis;
                this.interquartile_range = dse_options.interquartile_range;
                this.median_q1 = dse_options.median_q1;
                this.median_q2 = dse_options.median_q2;
                this.median_q3 = dse_options.median_q3;
                this.mad_mean_arithmetic = dse_options.mad_mean_arithmetic;
                //this.mad_mean_harmonic_corrected = dse_options.mad_mean_harmonic_corrected;
                //this.mad_mean_geometric_corrected = dse_options.mad_mean_geometric_corrected;
                this.mad_median_q1 = dse_options.mad_median_q1;
                this.mad_median_q2 = dse_options.mad_median_q2;
                this.mad_median_q3 = dse_options.mad_median_q3;
                this.mad_mid_range = dse_options.mad_mid_range;
            }

            if (name != null)
            {
                this.options_name = /*program.string_debug*/(name);
            }

            if (enable != null)
            {
                set_enable(enable.Value);
            }
        }

        public static string[] keys = new string[]
        {
            nameof(count),
            nameof(count_zero_values),
            nameof(count_non_zero_values),
            nameof(count_distinct_values),
            nameof(sum),
            nameof(mean_arithmetic),
            //nameof(mean_geometric_corrected),
            //nameof(mean_geometric_nonzero),
            //nameof(mean_harmonic_corrected),
            //nameof(mean_harmonic_nonzero),
            nameof(min),
            nameof(max),
            nameof(range),
            nameof(mid_range),
            nameof(variance),
            nameof(dev_standard),
            nameof(root_mean_square),
            nameof(skewness),
            nameof(kurtosis),
            nameof(interquartile_range),
            nameof(median_q1),
            nameof(median_q2),
            nameof(median_q3),
            nameof(mad_mean_arithmetic),
            //nameof(mad_mean_harmonic_corrected),
            //nameof(mad_mean_geometric_corrected),
            nameof(mad_median_q1),
            nameof(mad_median_q2),
            nameof(mad_median_q3),
            nameof(mad_mid_range),
        };

        public List<(string key, bool value)> key_value_list()
        {
            var list = new List<(string key, bool value)> {
                (nameof(count),                         count),
                (nameof(count_zero_values),             count_zero_values),
                (nameof(count_non_zero_values),         count_non_zero_values),
                (nameof(count_distinct_values),         count_distinct_values),
                (nameof(sum),                           sum),
                (nameof(mean_arithmetic),               mean_arithmetic),
                //(nameof(mean_geometric_corrected),      mean_geometric_corrected),
                //(nameof(mean_geometric_nonzero),      mean_geometric_nonzero),
                //(nameof(mean_harmonic_corrected),       mean_harmonic_corrected),
                //(nameof(mean_harmonic_nonzero),       mean_harmonic_nonzero),
                (nameof(min),                           min),
                (nameof(max),                           max),
                (nameof(range),                         range),
                (nameof(mid_range),                     mid_range),
                (nameof(variance),                      variance),
                (nameof(dev_standard),                  dev_standard),
                (nameof(root_mean_square),              root_mean_square),
                (nameof(skewness),                      skewness),
                (nameof(kurtosis),                      kurtosis),
                (nameof(interquartile_range),           interquartile_range),
                (nameof(median_q1),                     median_q1),
                (nameof(median_q2),                     median_q2),
                (nameof(median_q3),                     median_q3),
                (nameof(mad_mean_arithmetic),           mad_mean_arithmetic),
                //(nameof(mad_mean_harmonic_corrected),   mad_mean_harmonic_corrected),
                //(nameof(mad_mean_geometric_corrected),  mad_mean_geometric_corrected),
                (nameof(mad_median_q1),                 mad_median_q1),
                (nameof(mad_median_q2),                 mad_median_q2),
                (nameof(mad_median_q3),                 mad_median_q3),
                (nameof(mad_mid_range),                 mad_mid_range),
            };

            if (intervals != null) list.AddRange(intervals.key_value_list().Select(a => (/*program.string_debug*/($"{nameof(intervals)}.{a.key}"), a.value)).ToArray());
            if (distances != null) list.AddRange(distances.key_value_list().Select(a => (/*program.string_debug*/($"{nameof(distances)}.{a.key}"), a.value)).ToArray());
            if (interquartile != null) list.AddRange(interquartile.key_value_list().Select(a => (/*program.string_debug*/($"{nameof(interquartile)}.{a.key}"), a.value)).ToArray());
            if (abs != null) list.AddRange(abs.key_value_list().Select(a => (/*program.string_debug*/($"{nameof(abs)}.{a.key}"), a.value)).ToArray());
            if (rescale != null) list.AddRange(rescale.key_value_list().Select(a => (/*program.string_debug*/($"{nameof(rescale)}.{a.key}"), a.value)).ToArray());

            return list;
        }
    }
}
