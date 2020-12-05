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
        internal double mean_geometric;
        internal double mean_harmonic;
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
        internal double mad_mean_harmonic;
        internal double mad_mean_geometric;
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

        internal static double geometric_mean(double[] values, bool ignore_zeros = true, bool add_one = false)
        {
            // note: it only makes sense to ignore zeros if they are non-responses/not-applicable.
            // other strategies exist, such as adding 1 to all values or changing 0 to 1

            if (values == null || values.Length == 0) return 0;

            if (add_one)
            {
                values = values.Select(a => a + 1).ToArray();
            }

            double sum = 1.0;
            int zeros = 0;
            var use_log = false;

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] == 0 && ignore_zeros)
                {
                    zeros++;
                    continue;
                }

                sum *= values[i];

                if (double.IsInfinity(sum) || double.IsNaN(sum))
                {
                    use_log = true;

                    break;
                }
            }

            if (!use_log)
            {
                if (zeros == values.Length) return 0;

                return Math.Pow(sum, 1.0 / (values.Length - zeros));
            }
            else
            {
                zeros = 0;
                sum = 0;

                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] == 0 && ignore_zeros)
                    {
                        zeros++;
                        continue;
                    }
                    sum += Math.Log(values[i]);
                }

                if (zeros == values.Length) return 0;

                return Math.Exp(sum / (values.Length - zeros));
            }
        }

        public static double harmonic_mean(double[] values, bool correct_zeros = true, bool ignore_zeros = true, bool add_one = false)
        {
            if (values == null || values.Length == 0) return 0;

            if (add_one)
            {
                values = values.Select(a => a + 1).ToArray();
            }

            var nonzero = values.Where(a => a != 0).ToArray();
            var zeros = values.Length - nonzero.Length;

            if (correct_zeros)
            {
                var nt = (double)values.Length;
                var n0 = (double)zeros;
                var correction = (nt - n0) / nt;
                var hm = (1 / nonzero.Select(i => 1.0 / i).Average()) * correction;
                return hm;
            }

            if (!ignore_zeros && zeros > 0) return 0;

            return values.Length / values.Sum(i => i != 0 ? 1.0 / i : 0);
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
            if (data == null || data.Length == 0)
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


        internal List<(string group_id, string member_id, string perspective_id, double perspective_value)> encode(
                descriptive_stats_encoding_options descriptive_stats_encoding_options)
        {
            return encode(
                this,
                descriptive_stats_encoding_options
            );
        }

        private static List<(string group_id, string member_id, string perspective_id, double perspective_value)> encode(
            descriptive_stats stats,
            descriptive_stats_encoding_options descriptive_stats_encoding_options)
        {
            var result = new List<(string group_id, string member_id, string perspective_id, double perspective_value)>();

            if (stats == null)
            {
                stats = new descriptive_stats(data: null, descriptive_stats_encoding_options, ds_group_name: $@"", ds_member_name: $@"", presorted: true);
            }

            if (descriptive_stats_encoding_options == null)
            {
                descriptive_stats_encoding_options = descriptive_stats_encoding_options.options_default;
            }

            if (descriptive_stats_encoding_options?.intervals != null && stats.intervals_descriptive_stats != null && descriptive_stats_encoding_options.intervals.key_value_list().Any(a => a.value))
            {
                var encoded_intervals_descriptive_stats = encode(stats.intervals_descriptive_stats, descriptive_stats_encoding_options.intervals);

                result.AddRange(encoded_intervals_descriptive_stats);
            }

            if (descriptive_stats_encoding_options?.distances != null && stats.distances_descriptive_stats != null && descriptive_stats_encoding_options.distances.key_value_list().Any(a => a.value))
            {
                var encoded_distances_descriptive_stats = encode(stats.distances_descriptive_stats, descriptive_stats_encoding_options.distances);

                result.AddRange(encoded_distances_descriptive_stats);
            }

            if (descriptive_stats_encoding_options?.interquartile != null && stats.interquartile_range_descriptive_stats != null && descriptive_stats_encoding_options.interquartile.key_value_list().Any(a => a.value))
            {
                var encoded_interquartile_range_descriptive_stats = encode(stats.interquartile_range_descriptive_stats, descriptive_stats_encoding_options.interquartile);

                result.AddRange(encoded_interquartile_range_descriptive_stats);
            }

            if (descriptive_stats_encoding_options?.abs != null && stats.abs_descriptive_stats != null && descriptive_stats_encoding_options.abs.key_value_list().Any(a => a.value))
            {
                var encoded_abs_descriptive_stats = encode(stats.abs_descriptive_stats, descriptive_stats_encoding_options.abs);

                result.AddRange(encoded_abs_descriptive_stats);
            }


            if (descriptive_stats_encoding_options?.rescale != null && stats.rescaled_descriptive_stats != null && descriptive_stats_encoding_options.rescale.key_value_list().Any(a => a.value))
            {
                var encoded_abs_descriptive_stats = encode(stats.rescaled_descriptive_stats, descriptive_stats_encoding_options.rescale);

                result.AddRange(encoded_abs_descriptive_stats);
            }

            var z = new List<(string group_id, string member_id, string perspective_id, double perspective_value)>();
            if (descriptive_stats_encoding_options.count) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.count)}", (double)stats.count)); }
            if (descriptive_stats_encoding_options.count_zero_values) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.count_zero_values)}", (double)stats.count_zero_values)); }
            if (descriptive_stats_encoding_options.count_non_zero_values) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.count_non_zero_values)}", (double)stats.count_non_zero_values)); }
            if (descriptive_stats_encoding_options.count_distinct_values) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.count_distinct_values)}", (double)stats.count_distinct_values)); }
            if (descriptive_stats_encoding_options.min) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.min)}", stats.min)); }
            if (descriptive_stats_encoding_options.max) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.max)}", stats.max)); }
            if (descriptive_stats_encoding_options.range) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.range)}", stats.range)); }
            if (descriptive_stats_encoding_options.sum) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.sum)}", stats.sum)); }
            if (descriptive_stats_encoding_options.mid_range) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mid_range)}", stats.mid_range)); }
            //if (descriptive_stats_encoding_options.mode) { z.Add((group_id: $@"{stats.group_id_name}", member_id: $@"{stats.member_id_name}", perspective_id: $@"{nameof(stats.mode)}", stats.mode)); }
            if (descriptive_stats_encoding_options.median_q1) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.median_q1)}", stats.median_q1)); }
            if (descriptive_stats_encoding_options.median_q2) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.median_q2)}", stats.median_q2)); }
            if (descriptive_stats_encoding_options.median_q3) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.median_q3)}", stats.median_q3)); }

            if (descriptive_stats_encoding_options.root_mean_square) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.root_mean_square)}", stats.root_mean_square)); }
            if (descriptive_stats_encoding_options.mean_arithmetic) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mean_arithmetic)}", stats.mean_arithmetic)); }
            if (descriptive_stats_encoding_options.mean_harmonic) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mean_harmonic)}", stats.mean_harmonic)); }
            if (descriptive_stats_encoding_options.mean_geometric) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mean_geometric)}", stats.mean_geometric)); }

            if (descriptive_stats_encoding_options.variance) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.variance)}", stats.variance)); }
            if (descriptive_stats_encoding_options.dev_standard) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.dev_standard)}", stats.dev_standard)); }

            if (descriptive_stats_encoding_options.mad_mean_arithmetic) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_mean_arithmetic)}", stats.mad_mean_arithmetic)); }
            if (descriptive_stats_encoding_options.mad_mean_harmonic) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_mean_harmonic)}", stats.mad_mean_harmonic)); }
            if (descriptive_stats_encoding_options.mad_mean_geometric) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_mean_geometric)}", stats.mad_mean_geometric)); }

            if (descriptive_stats_encoding_options.mad_median_q1) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_median_q1)}", stats.mad_median_q1)); }
            if (descriptive_stats_encoding_options.mad_median_q2) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_median_q2)}", stats.mad_median_q2)); }
            if (descriptive_stats_encoding_options.mad_median_q3) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_median_q3)}", stats.mad_median_q3)); }
            //if (descriptive_stats_encoding_options.mad_mode) { z.Add((group_id: $@"{stats.group_id_name}", member_id: $@"{stats.member_id_name}", perspective_id: $@"{nameof(stats.mad_mode)}", stats.mad_mode)); }
            if (descriptive_stats_encoding_options.mad_mid_range) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.mad_mid_range)}", stats.mad_mid_range)); }

            if (descriptive_stats_encoding_options.interquartile_range) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.interquartile_range)}", stats.interquartile_range)); }
            if (descriptive_stats_encoding_options.skewness) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.skewness)}", stats.skewness)); }
            if (descriptive_stats_encoding_options.kurtosis) { z.Add((group_id: $@"{stats.ds_group_name}", member_id: $@"{stats.ds_member_name}", perspective_id: $@"{nameof(stats.kurtosis)}", stats.kurtosis)); }
            result.AddRange(z);


            return result;
        }

        internal static descriptive_stats get_stat_values(
            double[] data,
            descriptive_stats_encoding_options descriptive_stats_encoding_options,
            string group_id_name,
            string member_id_name,
            bool presorted
        )
        {
            return new descriptive_stats(data, descriptive_stats_encoding_options, group_id_name, member_id_name, presorted);
        }

        internal descriptive_stats(
            double[] data,
            descriptive_stats_encoding_options descriptive_stats_encoding_options,
            string ds_group_name,
            string ds_member_name,
            bool presorted
        )
        {

            this.ds_group_name = ds_group_name;
            this.ds_member_name = ds_member_name;


            if (data == null || data.Length == 0 || data.All(a => a == 0))
            {

                if (descriptive_stats_encoding_options?.intervals != null)
                {
                    var interval_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.intervals), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    var interval_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.intervals), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    intervals_descriptive_stats = new descriptive_stats(Array.Empty<double>(), descriptive_stats_encoding_options?.intervals, interval_group_name, interval_member_name, presorted: true);
                }

                if (descriptive_stats_encoding_options?.distances != null)
                {
                    var distance_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.distances), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    var distance_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.distances), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    distances_descriptive_stats = new descriptive_stats(Array.Empty<double>(), descriptive_stats_encoding_options?.distances, distance_group_name, distance_member_name, presorted: true);
                }

                if (descriptive_stats_encoding_options?.interquartile != null)
                {
                    var interquartile_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.interquartile), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    var interquartile_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.interquartile), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    interquartile_range_descriptive_stats = new descriptive_stats(Array.Empty<double>(), descriptive_stats_encoding_options?.interquartile, interquartile_group_name, interquartile_member_name, presorted: true);
                }

                if (descriptive_stats_encoding_options?.abs != null)
                {
                    var abs_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.abs), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    var abs_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.abs), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    abs_descriptive_stats = new descriptive_stats(Array.Empty<double>(), descriptive_stats_encoding_options?.abs, abs_group_name, abs_member_name, presorted: true);
                }

                if (descriptive_stats_encoding_options?.rescale != null)
                {
                    var rescale_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.rescale), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    var rescale_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.rescale), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                    rescaled_descriptive_stats = new descriptive_stats(Array.Empty<double>(), descriptive_stats_encoding_options?.rescale, rescale_group_name, rescale_member_name, presorted: true);
                }

                return;
            }

            if (data.Any(a => double.IsInfinity(a) || double.IsNaN(a)))
            {
                throw new ArgumentOutOfRangeException(nameof(data), $@"");
            }

         

            var sorted_data = presorted ? data : data.OrderBy(a => a).ToArray();

            count = (uint)sorted_data.Length;
            count_distinct_values = (uint)sorted_data.Distinct().Count();
            count_non_zero_values = (uint)sorted_data.Count(a => a != 0);
            count_zero_values = (uint)sorted_data.Length - count_non_zero_values;

            sum = sorted_data.Sum();
            root_mean_square = rms(sorted_data);
            mean_arithmetic = count != 0 ? sum / count : 0d;
            mean_harmonic = harmonic_mean(sorted_data);
            mean_geometric = geometric_mean(sorted_data);



            var stat = shape(sorted_data);


            variance = stat.variance;
            dev_standard = stat.stdev;
            kurtosis = stat.kurtosis;
            skewness = stat.skewness;


            min = sorted_data[0];
            max = sorted_data[^1];
            range = max - min;
            mid_range = (max + min) / 2.0;

            median_q1 = percentile(sorted_data, 25);
            median_q3 = percentile(sorted_data, 75);
            median_q2 = percentile(sorted_data, 50);
            interquartile_range = Math.Abs(median_q3 - median_q1);

            //var sorted_data_groups = sorted_data.GroupBy(x => x).ToArray();
            //var sorted_data_max_group_count = sorted_data_groups.Max(g => g.Count());
            //var modes = sorted_data_groups.Where(a => a.Count() == sorted_data_max_group_count).ToList();
            //mode = modes.Select(a => a.Key).DefaultIfEmpty(0).Average();

            mad_mean_arithmetic = mad(sorted_data, mean_arithmetic);
            mad_mean_geometric = mad(sorted_data, mean_geometric);
            mad_mean_harmonic = mad(sorted_data, mean_harmonic);
            mad_median_q1 = mad(sorted_data, median_q1);
            mad_median_q2 = mad(sorted_data, median_q2);
            mad_median_q3 = mad(sorted_data, median_q3);
            mad_mid_range = mad(sorted_data, mid_range);
            //mad_mode = mad(sorted_data, mode);


            fix_double();


            if (descriptive_stats_encoding_options?.intervals != null && descriptive_stats_encoding_options.intervals.key_value_list().Any(a => a.value))
            {
                var interval_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.intervals), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var interval_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.intervals), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var data_point_intervals = new double[sorted_data.Length - 1];

                for (var i = 1; i < sorted_data.Length; i++)
                {
                    data_point_intervals[i - 1] = Math.Abs(sorted_data[i] - sorted_data[i - 1]);
                }
                Array.Sort(data_point_intervals);

                intervals_descriptive_stats = new descriptive_stats(data_point_intervals, descriptive_stats_encoding_options.intervals, interval_group_name, interval_member_name, presorted: true);
            }

            if (descriptive_stats_encoding_options?.distances != null && descriptive_stats_encoding_options.distances.key_value_list().Any(a => a.value))
            {
                var distance_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.distances), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var distance_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.distances), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
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
                distances_descriptive_stats = new descriptive_stats(data_point_distances, descriptive_stats_encoding_options.distances, distance_group_name, distance_member_name, presorted: true);
            }

            if (descriptive_stats_encoding_options?.interquartile != null && descriptive_stats_encoding_options.interquartile.key_value_list().Any(a => a.value))
            {
                var interquartile_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.interquartile), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var interquartile_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.interquartile), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var interquartile_data = sorted_data.Where(a => a >= this.median_q1 && a <= this.median_q3).ToArray();
                interquartile_range_descriptive_stats = new descriptive_stats(interquartile_data, descriptive_stats_encoding_options.interquartile, interquartile_group_name, interquartile_member_name, presorted: true);
            }

            if (descriptive_stats_encoding_options?.abs != null && descriptive_stats_encoding_options.abs.key_value_list().Any(a => a.value))
            {
                var abs_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.abs), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var abs_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.abs), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var has_neg = sorted_data.Any(a => a < 0);
                var sorted_data_abs = has_neg ? sorted_data.Select(Math.Abs).OrderBy(a => a).ToArray() : sorted_data;
                abs_descriptive_stats = new descriptive_stats(sorted_data_abs, descriptive_stats_encoding_options.abs, abs_group_name, abs_member_name, presorted: true);
            }

            if (descriptive_stats_encoding_options?.rescale != null && descriptive_stats_encoding_options.rescale.key_value_list().Any(a => a.value))
            {
                var rescale_group_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.rescale), this.ds_group_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());
                var rescale_member_name = string.Join($@"_", new string[] { nameof(descriptive_stats_encoding_options.rescale), this.ds_member_name }.Where(a => !string.IsNullOrWhiteSpace(a)).ToArray());

                var s = new scaling(sorted_data)
                {
                    rescale_scale_min = 1, // avoid zeros whilst keeping data on the same scale.. only applicable with single arrays, not for multi array comparison if the value size has meaning...
                    rescale_scale_max = 2
                };
                
                var rescaled_sorted_data = s.scale(sorted_data, scaling.scale_function.rescale);
                rescaled_descriptive_stats = new descriptive_stats(rescaled_sorted_data, descriptive_stats_encoding_options.rescale, rescale_group_name, rescale_member_name, presorted: true);
            }
        }

        public void fix_double()
        {
            fix_double(ref sum);
            fix_double(ref mean_arithmetic);
            fix_double(ref mean_geometric);
            fix_double(ref mean_harmonic);
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
            fix_double(ref mad_mean_harmonic);
            fix_double(ref mad_mean_geometric);
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

            //var output = false;

            if (double.IsPositiveInfinity(value) || value >= c_double_max || value >= double.MaxValue)
            {
                //if (output) io_proxy.WriteLine($@"{nameof(fix_double)}: {name} = {value:G17} is positive infinity.");

                value = c_double_max;
            }
            else if (double.IsNegativeInfinity(value) || value <= c_double_min || value <= double.MinValue)
            {
                //if (output) io_proxy.WriteLine($@"{nameof(fix_double)}: {name} = {value:G17} is negative infinity.");

                value = c_double_min;
            }
            else if (double.IsNaN(value))
            {
                //if (output) io_proxy.WriteLine($@"{nameof(fix_double)}: {name} = {value:G17} is not a number.");

                value = double_zero;
            }
        }

        public static double fix_double(double value)
        {
            const double c_double_max = (double)1.79769e+308;
            const double c_double_min = (double)-c_double_max;
            const double double_zero = (double)0;

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
