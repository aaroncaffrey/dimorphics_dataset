using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    public class descriptive_stats
    {
        internal string ds_group_name;
        internal string ds_member_name;

        internal uint count;
        internal uint count_zero_values;
        internal uint count_non_zero_values;
        internal uint count_distinct_values;
        internal double sum;
        internal double mean_arithmetic;
        //internal double mean_geometric_nonzeros;
        internal double mean_geometric_corrected;
        //internal double mean_harmonic_nonzeros;
        internal double mean_harmonic_corrected;
        internal double min;
        internal double max;
        internal double range;
        internal double mid_range;
        internal double variance;
        internal double dev_standard;
        internal double root_mean_square;
        internal double skewness;
        internal double kurtosis;
        internal double interquartile_range;
        internal double median_q1;
        internal double median_q2;
        internal double median_q3;
        internal double mad_mean_arithmetic;
        internal double mad_mean_harmonic_corrected;
        internal double mad_mean_geometric_corrected;
        internal double mad_median_q1;
        internal double mad_median_q2;
        internal double mad_median_q3;
        internal double mad_mid_range;

        //internal double mode;
        //internal double mad_mode;

        internal descriptive_stats intervals_descriptive_stats;
        internal descriptive_stats distances_descriptive_stats;
        internal descriptive_stats interquartile_range_descriptive_stats;
        internal descriptive_stats abs_descriptive_stats;
        internal descriptive_stats rescaled_descriptive_stats;

        internal static (double zeros, double nonzeros, double corrected) geometric_mean(double[] values)
        {
            if (values == null || values.Length == 0 || values.All(a => a == 0)) return (0.0, 0.0, 0.0);

            var use_log = false;
            var values_nonzero = values.Where(a => a != 0).ToArray();
            var num_zeros = values.Length - values_nonzero.Length;
            var sum = use_log ? 0.0 : 1.0;

            if (!use_log)
            {
                for (var i = 0; i < values_nonzero.Length; i++)
                {
                    sum *= values_nonzero[i];

                    if (double.IsInfinity(sum) || double.IsNaN(sum) || sum >= double.MaxValue)
                    {
                        use_log = true;
                        sum = 0.0;
                        break;
                    }
                }
            }

            if (use_log)
            {
                for (var i = 0; i < values_nonzero.Length; i++)
                {
                    sum += Math.Log(values_nonzero[i]);
                }
            }

            var result_nonzeros = use_log ? Math.Exp(sum / values_nonzero.Length) : Math.Pow(sum, 1.0 / values_nonzero.Length);
            var correction = ((double)values.Length - (double)num_zeros) / (double)values.Length;
            var result_corrected = result_nonzeros * correction;
            var result_zeros = num_zeros == 0 ? result_nonzeros : 0.0;

            return (result_zeros, result_nonzeros, result_corrected);
        }

        public static (double zeros, double nonzeros, double corrected) harmonic_mean(double[] values)
        {
            if (values == null || values.Length == 0 || values.All(a => a == 0)) return (0.0, 0.0, 0.0);

            //  hm = values.Length / values.Sum(i => 1.0 / i);
            var values_nonzero = values.Where(a => a != 0).ToArray();
            var num_zeros = values.Length - values_nonzero.Length;
            var result_nonzeros = values_nonzero.Length / values_nonzero.Sum(i => 1.0 / i); // same as (1 / values_nonzero.Select(i => 1.0 / i).Average());
            var correction = ((double)values.Length - (double)num_zeros) / (double)values.Length;
            var result_corrected = result_nonzeros * correction;
            var result_zeros = num_zeros == 0 ? result_nonzeros : 0.0;

            return (result_zeros, result_nonzeros, result_corrected);
        }

        public static double sample_variance(double[] samples)
        {
            if (samples == null || samples.Length <= 1)
            {
                return 0;
            }

            var variance = 0.0;
            var t = samples[0];
            for (int i = 1; i < samples.Length; i++)
            {
                t += samples[i];
                var diff = ((i + 1) * samples[i]) - t;
                variance += (diff * diff) / ((i + 1.0) * i);
            }

            return variance / (samples.Length - 1);
        }

        public static (double variance, double stdev) sample_standard_deviation(double[] samples)
        {
            if (samples == null || samples.Length <= 1)
            {
                return (0, 0);
            }

            var variance = sample_variance(samples);
            var stdev = Math.Sqrt(variance);

            return (variance, stdev);
        }

        public static double rms(double[] data)
        {
            if (data == null || data.Length == 0)
            {
                return double.NaN;
            }

            double mean = 0;
            ulong m = 0;
            for (int i = 0; i < data.Length; i++)
            {
                mean += (data[i] * data[i] - mean) / ++m;
            }

            return Math.Sqrt(mean);
        }

        public static (double skewness, double kurtosis, double mean, double variance, double stdev) shape(double[] data)
        {
            if (data == null || data.Length == 0 || data.All(a=> a == 0))
            {
                return (0.0, 0.0, 0.0, 0.0, 0.0);
            }

            var mean = 0.0;
            var error = 0.0;
            var skewness = 0.0;
            var kurtosis = 0.0;
            long n = 0;

            for (var index = 0; index < data.Length; index++)
            {
                var delta = data[index] - mean;
                var scaleDelta = delta / ++n;
                var scaleDeltaSqr = scaleDelta * scaleDelta;
                var tmpDelta = delta * (n - 1);

                mean += scaleDelta;

                kurtosis += tmpDelta * scaleDelta * scaleDeltaSqr * (n * n - 3 * n + 3) + 6 * scaleDeltaSqr * error - 4 * scaleDelta * skewness;

                skewness += tmpDelta * scaleDeltaSqr * (n - 2) - 3 * scaleDelta * error;
                error += tmpDelta * scaleDelta;
            }

            var variance = n > 1 ? error / (n - 1) : 0;
            var stdev = variance != 0 ? Math.Sqrt(variance) : 0;
            skewness = variance != 0 && n > 2 ? (double)n / ((n - 1) * (n - 2)) * (skewness / (variance * stdev)) : 0;
            kurtosis = variance != 0 && n > 3 ? ((double)n * n - 1) / ((n - 2) * (n - 3)) * (n * kurtosis / (error * error) - 3 + 6.0 / (n + 1)) : 0;

            return (skewness, kurtosis, mean, variance, stdev);
        }


        internal List<(string group_id, string member_id, string perspective_id, double perspective_value)> encode
        (
                descriptive_stats_encoding_options dse_options,
                /*bool same_group = true,*/
                /*bool individual = true,*/
                int level = 0
        )
        {
            return encode(
                this,
                dse_options,
                /*same_group,*/
                /*individual,*/
                level
            );
        }

        private static List<(string group_id, string member_id, string perspective_id, double perspective_value)> encode(
            descriptive_stats stats,
            descriptive_stats_encoding_options dse_options,
            /*bool same_group = true,*/
            /*bool individual = true,*/
            int level = 0)
        {
            const string module_name = nameof(descriptive_stats);
            const string method_name = nameof(encode);

            var result = new List<(string group_id, string member_id, string perspective_id, double perspective_value)>();

            if (stats == null)
            {
                stats = new descriptive_stats(null, dse_options, ds_group_name: /*program.string_debug*/($@""), ds_member_name: /*program.string_debug*/($@""), presorted: true);
            }

            /*if (!same_group && !individual)
            {
                throw new Exception(/*program.string_debug* /($@"{module_name}.{method_name}: {nameof(same_group)} and {nameof(individual)} both false."));
            }*/

            if (dse_options == null)
            {
                //descriptive_stats_encoding_options = descriptive_stats_encoding_options.options_mean_arithmetic;
                throw new ArgumentOutOfRangeException(nameof(dse_options), /*program.string_debug*/($@""));
            }

            if (dse_options?.intervals != null && stats.intervals_descriptive_stats != null && dse_options.intervals.key_value_list().Any(a => a.value))
            {
                var encoded_intervals_descriptive_stats = stats.intervals_descriptive_stats.encode(dse_options.intervals, /*same_group, individual,*/ level + 1);

                result.AddRange(encoded_intervals_descriptive_stats);
            }

            if (dse_options?.distances != null && stats.distances_descriptive_stats != null && dse_options.distances.key_value_list().Any(a => a.value))
            {
                var encoded_distances_descriptive_stats = stats.distances_descriptive_stats.encode(dse_options.distances, /*same_group, individual,*/ level + 1);

                result.AddRange(encoded_distances_descriptive_stats);
            }

            if (dse_options?.interquartile != null && stats.interquartile_range_descriptive_stats != null && dse_options.interquartile.key_value_list().Any(a => a.value))
            {
                var encoded_interquartile_range_descriptive_stats = stats.interquartile_range_descriptive_stats.encode(dse_options.interquartile, /*same_group, individual,*/ level + 1);

                result.AddRange(encoded_interquartile_range_descriptive_stats);
            }

            if (dse_options?.abs != null && stats.abs_descriptive_stats != null && dse_options.abs.key_value_list().Any(a => a.value))
            {
                var encoded_abs_descriptive_stats = stats.abs_descriptive_stats.encode(dse_options.abs, /*same_group, individual,*/ level + 1);

                result.AddRange(encoded_abs_descriptive_stats);
            }


            if (dse_options?.rescale != null && stats.rescaled_descriptive_stats != null && dse_options.rescale.key_value_list().Any(a => a.value))
            {
                var encoded_abs_descriptive_stats = stats.rescaled_descriptive_stats.encode(dse_options.rescale, /*same_group, individual,*/ level + 1);

                result.AddRange(encoded_abs_descriptive_stats);
            }


            /*if (same_group || individual)*/
            {
                var z = new List<(string group_id, string member_id, string perspective_id, double perspective_value)>();

                if (dse_options.count) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.count)}"), (double)stats.count)); }
                if (dse_options.count_zero_values) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.count_zero_values)}"), (double)stats.count_zero_values)); }
                if (dse_options.count_non_zero_values) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.count_non_zero_values)}"), (double)stats.count_non_zero_values)); }
                if (dse_options.count_distinct_values) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.count_distinct_values)}"), (double)stats.count_distinct_values)); }
                if (dse_options.min) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.min)}"), stats.min)); }
                if (dse_options.max) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.max)}"), stats.max)); }
                if (dse_options.range) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.range)}"), stats.range)); }
                if (dse_options.sum) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.sum)}"), stats.sum)); }
                if (dse_options.mid_range) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mid_range)}"), stats.mid_range)); }
                if (dse_options.median_q1) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.median_q1)}"), stats.median_q1)); }
                if (dse_options.median_q2) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.median_q2)}"), stats.median_q2)); }
                if (dse_options.median_q3) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.median_q3)}"), stats.median_q3)); }
                if (dse_options.root_mean_square) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.root_mean_square)}"), stats.root_mean_square)); }
                if (dse_options.mean_arithmetic) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mean_arithmetic)}"), stats.mean_arithmetic)); }

                if (dse_options.mean_geometric_corrected) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mean_geometric_corrected)}"), stats.mean_geometric_corrected)); }
                //if (dse_options.mean_geometric_nonzero) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mean_geometric_nonzeros)}"), stats.mean_geometric_nonzeros)); }

                if (dse_options.mean_harmonic_corrected) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mean_harmonic_corrected)}"), stats.mean_harmonic_corrected)); }
                //if (dse_options.mean_harmonic_nonzero) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mean_harmonic_nonzeros)}"), stats.mean_harmonic_nonzeros)); }

                if (dse_options.variance) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.variance)}"), stats.variance)); }
                if (dse_options.dev_standard) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.dev_standard)}"), stats.dev_standard)); }
                if (dse_options.mad_mean_arithmetic) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_mean_arithmetic)}"), stats.mad_mean_arithmetic)); }
                if (dse_options.mad_mean_harmonic_corrected) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_mean_harmonic_corrected)}"), stats.mad_mean_harmonic_corrected)); }
                if (dse_options.mad_mean_geometric_corrected) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_mean_geometric_corrected)}"), stats.mad_mean_geometric_corrected)); }
                if (dse_options.mad_median_q1) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_median_q1)}"), stats.mad_median_q1)); }
                if (dse_options.mad_median_q2) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_median_q2)}"), stats.mad_median_q2)); }
                if (dse_options.mad_median_q3) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_median_q3)}"), stats.mad_median_q3)); }
                if (dse_options.mad_mid_range) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_mid_range)}"), stats.mad_mid_range)); }
                if (dse_options.interquartile_range) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.interquartile_range)}"), stats.interquartile_range)); }
                if (dse_options.skewness) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.skewness)}"), stats.skewness)); }
                if (dse_options.kurtosis) { z.Add((group_id: /*program.string_debug*/($@"{stats.ds_group_name}"), member_id: /*program.string_debug*/($@"{stats.ds_member_name}"), perspective_id: /*program.string_debug*/($@"{nameof(stats.kurtosis)}"), stats.kurtosis)); }
                //if (descriptive_stats_encoding_options.mode) { z.Add((group_id: /*program.string_debug*/($@"{stats.group_id_name}", member_id: /*program.string_debug*/($@"{stats.member_id_name}", perspective_id: /*program.string_debug*/($@"{nameof(stats.mode)}", stats.mode)); }
                //if (descriptive_stats_encoding_options.mad_mode) { z.Add((group_id: /*program.string_debug*/($@"{stats.group_id_name}", member_id: /*program.string_debug*/($@"{stats.member_id_name}", perspective_id: /*program.string_debug*/($@"{nameof(stats.mad_mode)}", stats.mad_mode)); }

                if (z.Count > 0)
                {
                    /*if (same_group)*/
                    {
                        result.AddRange(z);
                    }

                    /*if (individual)
                    {
                        var z2 = z.Select(a => (string.Join(/*program.string_debug* /($@"_"), new[] {a.group_id, a.member_id, a.perspective_id}.Where(c => !string.IsNullOrWhiteSpace(c)).Distinct().ToArray()), a.member_id, a.perspective_id, a.perspective_value)).ToList();

                        result.AddRange(z2);
                    }*/
                }
            }


            if (level == 0 && result.Count == 0)
            {
                throw new ArgumentOutOfRangeException(nameof(dse_options), /*program.string_debug*/($@"{module_name}.{method_name}: no features are enabled in {nameof(dse_options)}."));
            }

            return result;
        }

        internal static descriptive_stats get_stat_values(
            double[] data,
            descriptive_stats_encoding_options dse_options,
            string group_id_name,
            string member_id_name,
            bool presorted
        )
        {
            return new descriptive_stats(data, dse_options, group_id_name, member_id_name, presorted);
        }

        internal descriptive_stats(
            double[] data,
            descriptive_stats_encoding_options dse_options,
            string ds_group_name,
            string ds_member_name,
            bool presorted
        )
        {

            this.ds_group_name = /*program.string_debug*/(ds_group_name);
            this.ds_member_name = /*program.string_debug*/(ds_member_name);


            if (data == null || data.Length == 0 || data.All(a => a == 0))
            {

                if (dse_options?.intervals != null)
                {
                    var interval_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.intervals), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    var interval_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.intervals), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    intervals_descriptive_stats = new descriptive_stats(Array.Empty<double>(), dse_options?.intervals, interval_group_name, interval_member_name, presorted: true);
                }

                if (dse_options?.distances != null)
                {
                    var distance_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.distances), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    var distance_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.distances), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    distances_descriptive_stats = new descriptive_stats(Array.Empty<double>(), dse_options?.distances, distance_group_name, distance_member_name, presorted: true);
                }

                if (dse_options?.interquartile != null)
                {
                    var interquartile_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.interquartile), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    var interquartile_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.interquartile), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    interquartile_range_descriptive_stats = new descriptive_stats(Array.Empty<double>(), dse_options?.interquartile, interquartile_group_name, interquartile_member_name, presorted: true);
                }

                if (dse_options?.abs != null)
                {
                    var abs_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.abs), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    var abs_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.abs), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    abs_descriptive_stats = new descriptive_stats(Array.Empty<double>(), dse_options?.abs, abs_group_name, abs_member_name, presorted: true);
                }

                if (dse_options?.rescale != null)
                {
                    var rescale_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.rescale), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    var rescale_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.rescale), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                    rescaled_descriptive_stats = new descriptive_stats(Array.Empty<double>(), dse_options?.rescale, rescale_group_name, rescale_member_name, presorted: true);
                }

                return;
            }

            if (data.Any(a => double.IsInfinity(a) || double.IsNaN(a)))
            {
                throw new ArgumentOutOfRangeException(nameof(data), /*program.string_debug*/($@""));
            }



            var sorted_data = presorted ? data : data.OrderBy(a => a).ToArray();

            if (dse_options.count || dse_options.sum || dse_options.mean_arithmetic || dse_options.dev_standard || dse_options.mad_mean_arithmetic)
            {
                count = (uint) sorted_data.Length;
            }

            if (dse_options.count_distinct_values)
            {
                count_distinct_values = (uint)sorted_data.Distinct().Count();
            }

            if (dse_options.count_non_zero_values || dse_options.count_zero_values)
            {
                count_non_zero_values = (uint)sorted_data.Count(a => a != 0);
            }

            if (dse_options.count_zero_values)
            {
                count_zero_values = (uint)sorted_data.Length - count_non_zero_values;
            }

            if (dse_options.sum || dse_options.mean_arithmetic || dse_options.dev_standard || dse_options.mad_mean_arithmetic)
            {
                sum = sorted_data.Sum();
            }

            if (dse_options.root_mean_square)
            {
                root_mean_square = rms(sorted_data);
            }

            if (dse_options.mean_arithmetic || dse_options.dev_standard || dse_options.mad_mean_arithmetic)
            {
                mean_arithmetic = count != 0 ? sum / count : 0d;
            }

            if (dse_options.mean_harmonic_corrected || /*dse_options.mean_harmonic_nonzero ||*/ dse_options.mad_mean_harmonic_corrected)
            {
                var hm = harmonic_mean(sorted_data);
                mean_harmonic_corrected = hm.corrected;
                //mean_harmonic_nonzeros = hm.nonzeros;
            }

            if (dse_options.mean_geometric_corrected || /*dse_options.mean_geometric_nonzero ||*/ dse_options.mad_mean_geometric_corrected)
            {
                var gm = geometric_mean(sorted_data);
                mean_geometric_corrected = gm.corrected;
                //mean_geometric_nonzeros = gm.nonzeros;
            }


            if (dse_options.variance || dse_options.dev_standard || dse_options.kurtosis || dse_options.skewness)
            {
                var stat = shape(sorted_data);
                variance = stat.variance;
                dev_standard = stat.stdev;
                kurtosis = stat.kurtosis;
                skewness = stat.skewness;
            }

            if (dse_options.min || dse_options.range || dse_options.mid_range || dse_options.mad_mid_range)
            {
                min = sorted_data[0];
            }

            if (dse_options.max || dse_options.range || dse_options.mid_range || dse_options.mad_mid_range)
            {
                max = sorted_data[^1];
            }

            if (dse_options.range)
            {
                range = max - min;
            }

            if (dse_options.mid_range || dse_options.mad_mid_range)
            {
                mid_range = (max + min) / 2.0;
            }

            if (dse_options.median_q1 || dse_options.mad_median_q1 || dse_options.interquartile_range)
            {
                median_q1 = percentile(sorted_data, 25);
            }

            if (dse_options.median_q2 || dse_options.mad_median_q2 || dse_options.interquartile_range)
            {
                median_q2 = percentile(sorted_data, 50);
            }

            if (dse_options.median_q3 || dse_options.mad_median_q3 || dse_options.interquartile_range)
            {
                median_q3 = percentile(sorted_data, 75);
            }

            if (dse_options.interquartile_range)
            {
                interquartile_range = Math.Abs(median_q3 - median_q1);
            }

            //var sorted_data_groups = sorted_data.GroupBy(x => x).ToArray();
            //var sorted_data_max_group_count = sorted_data_groups.Max(g => g.Count());
            //var modes = sorted_data_groups.Where(a => a.Count() == sorted_data_max_group_count).ToList();
            //mode = modes.Select(a => a.Key).DefaultIfEmpty(0).Average();

            if (dse_options.mad_mean_arithmetic) { mad_mean_arithmetic = mad(sorted_data, mean_arithmetic); }
            if (dse_options.mad_mean_geometric_corrected) { mad_mean_geometric_corrected = mad(sorted_data, mean_geometric_corrected); }
            if (dse_options.mad_mean_harmonic_corrected) { mad_mean_harmonic_corrected = mad(sorted_data, mean_harmonic_corrected); }
            if (dse_options.mad_median_q1) { mad_median_q1 = mad(sorted_data, median_q1); }
            if (dse_options.mad_median_q2) { mad_median_q2 = mad(sorted_data, median_q2); }
            if (dse_options.mad_median_q3) { mad_median_q3 = mad(sorted_data, median_q3); }
            if (dse_options.mad_mid_range) { mad_mid_range = mad(sorted_data, mid_range); }
            //mad_mode = mad(sorted_data, mode);


            fix_double();


            if (dse_options?.intervals != null && dse_options.intervals.key_value_list().Any(a => a.value))
            {
                var interval_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.intervals), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var interval_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.intervals), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var data_point_intervals = new double[sorted_data.Length - 1];

                for (var i = 1; i < sorted_data.Length; i++)
                {
                    data_point_intervals[i - 1] = Math.Abs(sorted_data[i] - sorted_data[i - 1]);
                }
                Array.Sort(data_point_intervals);

                intervals_descriptive_stats = new descriptive_stats(data_point_intervals, dse_options.intervals, interval_group_name, interval_member_name, presorted: true);
            }

            if (dse_options?.distances != null && dse_options.distances.key_value_list().Any(a => a.value))
            {
                var distance_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.distances), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var distance_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.distances), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var data_point_distances = new double[((sorted_data.Length - 1) * sorted_data.Length) / 2];
                var k = 0;

                for (var i = 0; i < sorted_data.Length; i++)
                {
                    for (var j = 0; j < sorted_data.Length; j++)
                    {
                        if (i <= j) continue;

                        data_point_distances[k] = Math.Abs(sorted_data[i] - sorted_data[j]);

                        k++;
                    }
                }

                Array.Sort(data_point_distances);
                distances_descriptive_stats = new descriptive_stats(data_point_distances, dse_options.distances, distance_group_name, distance_member_name, presorted: true);
            }

            if (dse_options?.interquartile != null && dse_options.interquartile.key_value_list().Any(a => a.value))
            {
                var interquartile_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.interquartile), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var interquartile_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.interquartile), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var interquartile_data = sorted_data.Where(a => a >= this.median_q1 && a <= this.median_q3).ToArray();
                interquartile_range_descriptive_stats = new descriptive_stats(interquartile_data, dse_options.interquartile, interquartile_group_name, interquartile_member_name, presorted: true);
            }

            if (dse_options?.abs != null && dse_options.abs.key_value_list().Any(a => a.value))
            {
                var abs_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.abs), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var abs_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.abs), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var has_neg = sorted_data.Any(a => a < 0);
                var sorted_data_abs = has_neg ? sorted_data.Select(Math.Abs).OrderBy(a => a).ToArray() : sorted_data;
                abs_descriptive_stats = new descriptive_stats(sorted_data_abs, dse_options.abs, abs_group_name, abs_member_name, presorted: true);
            }

            if (dse_options?.rescale != null && dse_options.rescale.key_value_list().Any(a => a.value))
            {
                var rescale_group_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.rescale), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));
                var rescale_member_name = /*program.string_debug*/(string.Join(/*program.string_debug*/($@"_"), new string[] { nameof(dse_options.rescale), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray()));

                var s = new scaling(sorted_data)
                {
                    rescale_scale_min = 1, // avoid zeros whilst keeping data on the same linear scale.. only applicable with single arrays, not for multi array comparison if the value size has meaning...
                    rescale_scale_max = 2
                };

                var rescaled_sorted_data = s.scale(sorted_data, scaling.scale_function.rescale);
                rescaled_descriptive_stats = new descriptive_stats(rescaled_sorted_data, dse_options.rescale, rescale_group_name, rescale_member_name, presorted: true);
            }
        }

        public void fix_double()
        {
            fix_double(ref sum);
            fix_double(ref mean_arithmetic);
            fix_double(ref mean_geometric_corrected);
            //fix_double(ref mean_geometric_nonzeros);
            fix_double(ref mean_harmonic_corrected);
            //fix_double(ref mean_harmonic_nonzeros);
            fix_double(ref min);
            fix_double(ref max);
            fix_double(ref range);
            fix_double(ref mid_range);
            fix_double(ref variance);
            fix_double(ref dev_standard);
            fix_double(ref root_mean_square);
            fix_double(ref skewness);
            fix_double(ref kurtosis);
            fix_double(ref interquartile_range);
            fix_double(ref median_q1);
            fix_double(ref median_q2);
            fix_double(ref median_q3);
            //fix_double(ref mode);
            fix_double(ref mad_mean_arithmetic);
            fix_double(ref mad_mean_harmonic_corrected);
            fix_double(ref mad_mean_geometric_corrected);
            fix_double(ref mad_median_q1);
            fix_double(ref mad_median_q2);
            fix_double(ref mad_median_q3);
            //fix_double(ref mad_mode);
            fix_double(ref mad_mid_range);
        }

        public static void fix_double(/*string name,*/ ref double value)
        {
            const double c_double_max = 1.79769e+308;
            const double c_double_min = -c_double_max;
            const double double_zero = 0.0;

            if (value == 0 || (value >= c_double_min && value <= c_double_max)) return;

            if (double.IsPositiveInfinity(value) || value >= c_double_max || value >= double.MaxValue)
            {
                //if (output) io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(fix_double)}: {name} = {value:G17} is positive infinity.");

                value = c_double_max;
            }
            else if (double.IsNegativeInfinity(value) || value <= c_double_min || value <= double.MinValue)
            {
                //if (output) io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(fix_double)}: {name} = {value:G17} is negative infinity.");

                value = c_double_min;
            }
            else if (double.IsNaN(value))
            {
                //if (output) io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(fix_double)}: {name} = {value:G17} is not a number.");

                value = double_zero;
            }
        }

        public static double fix_double(double value)
        {
            const double c_double_max = 1.79769e+308;
            const double c_double_min = -c_double_max;
            const double double_zero = 0.0;

            if (value == 0 || (value >= c_double_min && value <= c_double_max)) return value;

            if (double.IsPositiveInfinity(value) || value >= c_double_max || value >= double.MaxValue)
            {
                value = c_double_max;
            }
            else if (double.IsNegativeInfinity(value) || value <= c_double_min || value <= double.MinValue)
            {
                value = c_double_min;
            }
            else if (double.IsNaN(value))
            {
                value = double_zero;
            }

            return value;
        }

        internal static double mad(double[] values, double? centre)
        {
            if (values == null || values.Length == 0)
            {
                return 0;
            }

            if (centre == null) centre = values.Average();

            var mad = values.Sum(a => Math.Abs(centre.Value - a)) / (double)values.Length;

            return mad;
        }

        internal static double percentile(double[] sorted_data, double p)
        {
            if (sorted_data == null || sorted_data.Length == 0)
            {
                return 0;
            }

            if (sorted_data.Length == 1)
            {
                return sorted_data[0];
            }

            if (p >= 100.0)
            {
                return sorted_data[^1];
            }

            var position = (sorted_data.Length + 1) * p / (100.0);
            double left_number;
            double right_number;

            double n = p / (100.0) * (sorted_data.Length - 1) + (1.0);

            if (position >= 1)
            {
                left_number = sorted_data[(int)Math.Floor(n) - 1];
                right_number = sorted_data[(int)Math.Floor(n)];
            }
            else
            {
                left_number = sorted_data[0];
                right_number = sorted_data[1];
            }

            if (left_number == right_number)
            {
                return left_number;
            }
            else
            {
                double part = n - Math.Floor(n);
                return left_number + part * (right_number - left_number);
            }
        }
    }
}
