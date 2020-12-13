using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace dimorphics_dataset
{
    internal static class info_sable
    {
        //internal static List<sable_item> sable_data = sable.load(Path.Combine(program.data_root_folder,"sable\OUT_SABLE_RES.txt");

        internal static List<info_sable_item> load(string filename)// = Path.Combine(program.data_root_folder,"sable\OUT_SABLE_RES.txt")
        {
            var lines = io_proxy.ReadAllLines(filename, nameof(info_sable), nameof(load)).ToList();

            var queries = new List<List<string>>();

            var temp_query = new List<string>();
            for (var i = 0; i <= lines.Count; i++)
            {

                if (i == lines.Count || lines[i].StartsWith(/*program.string_debug*/($@"Query:"), StringComparison.Ordinal))
                {
                    if (temp_query != null && temp_query.Count > 0)
                    {
                        queries.Add(temp_query);
                    }

                    temp_query = new List<string>();

                    if (i == lines.Count) break;

                }

                if (!string.IsNullOrWhiteSpace(lines[i]))
                {
                    temp_query.Add(lines[i]);
                }
            }

            var result = new List<info_sable_item>();

            foreach (var query in queries)
            {
                var query_name = query.First(a => a.StartsWith(/*program.string_debug*/($@"Query:"), StringComparison.Ordinal)).Split().Last();

                // 3 sectionss
                var SECTION_SS = query.Skip(query.FindIndex(a => string.Equals(a, /*program.string_debug*/($@"SECTION_SS"), StringComparison.Ordinal))+1).TakeWhile(a => !string.Equals(a, /*program.string_debug*/($@"END_SECTION"), StringComparison.Ordinal)).SkipWhile(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Where(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Select(a => a.Trim()).ToList();

                var SECTION_SS_joined = (
                    query_seq: string.Join(/*program.string_debug*/($@""), SECTION_SS.Where((a, i) =>          i % 3 == 0).ToList()).ToCharArray(),
                    predicted_ss_seq: string.Join(/*program.string_debug*/($@""), SECTION_SS.Where((a, i) =>   i % 3 == 1).ToList()).ToCharArray(),
                    confidence_level: string.Join(/*program.string_debug*/($@""), SECTION_SS.Where((a, i) =>   i % 3 == 2).ToList()).Select(a=>double.Parse(a.ToString(CultureInfo.InvariantCulture), NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 10).ToList()
                    );

                // 4 sections H-> E-> C->
                var SECTION_SS_PROBABILITIES = query.Skip(query.FindIndex(a => string.Equals(a, /*program.string_debug*/($@"SECTION_SS_PROBABILITIES"), StringComparison.Ordinal))+1).TakeWhile(a => !string.Equals(a, /*program.string_debug*/($@"END_SECTION"), StringComparison.Ordinal)).SkipWhile(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Where(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Select(a => a.Trim()).ToList();

                var SECTION_SS_PROBABILITIES_joined = (
                    query_seq: string.Join(/*program.string_debug*/($@""), SECTION_SS_PROBABILITIES.Where((a, i) =>   i % 4 == 0).ToList()).Replace(/*program.string_debug*/($@" "),/*program.string_debug*/($@""), StringComparison.Ordinal).ToCharArray(),
                    prob_h: string.Join(/*program.string_debug*/($@""), SECTION_SS_PROBABILITIES.Where((a, i) =>   i % 4 == 1).ToList()).Replace(/*program.string_debug*/($@"H->"), /*program.string_debug*/($@""), StringComparison.Ordinal).Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 100).ToList(),
                    prob_e: string.Join(/*program.string_debug*/($@""), SECTION_SS_PROBABILITIES.Where((a, i) => i % 4 == 2).ToList()).Replace(/*program.string_debug*/($@"E->"), /*program.string_debug*/($@""), StringComparison.Ordinal).Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 100).ToList(),
                    prob_c: string.Join(/*program.string_debug*/($@""), SECTION_SS_PROBABILITIES.Where((a, i) =>  i % 4 == 3).ToList()).Replace(/*program.string_debug*/($@"C->"), /*program.string_debug*/($@""), StringComparison.Ordinal).Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 100).ToList()
                    );

                // 3 sections
                var SECTION_SA = query.Skip(query.FindIndex(a => string.Equals(a, /*program.string_debug*/($@"SECTION_SA"), StringComparison.Ordinal))+1).TakeWhile(a => !string.Equals(a, /*program.string_debug*/($@"END_SECTION"), StringComparison.Ordinal)).SkipWhile(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Where(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Select(a => a.Trim()).ToList();

                var SECTION_SA_joined = (
                    query_seq: string.Join(/*program.string_debug*/($@""), SECTION_SA.Where((a, i) =>   i % 3 == 0).ToList()).ToCharArray(),
                    burial_level: string.Join(/*program.string_debug*/($@""), SECTION_SA.Where((a, i) =>   i % 3 == 1).ToList()).Select(a => double.Parse(a.ToString(CultureInfo.InvariantCulture), NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 10).ToList(), //0=buried, 9=exposed
                    confidence_level: string.Join(/*program.string_debug*/($@""), SECTION_SA.Where((a, i) => i % 3 == 2).ToList()).Select(a => double.Parse(a.ToString(CultureInfo.InvariantCulture), NumberStyles.Float, NumberFormatInfo.InvariantInfo) / 10).ToList()
                );

                //every 2 lines
                var SECTION_SA_ABSOLUTE = query.Skip(query.FindIndex(a => string.Equals(a, /*program.string_debug*/($@"SECTION_SA_ABSOLUTE"), StringComparison.Ordinal))+1).TakeWhile(a => !string.Equals(a, /*program.string_debug*/($@"END_SECTION"), StringComparison.Ordinal)).SkipWhile(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Where(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Select(a => a.Trim()).ToList();

                var SECTION_SA_ABSOLUTE_joined = (
                    seq: string.Join(/*program.string_debug*/($@""), SECTION_SA_ABSOLUTE.Where((a, i) => i % 2 == 0).ToList()).Replace(/*program.string_debug*/($@" "),/*program.string_debug*/($@""), StringComparison.Ordinal).ToCharArray(),
                    absolute_burial_value: string.Join(/*program.string_debug*/($@" "), SECTION_SA_ABSOLUTE.Where((a, i) => i % 2 == 1).ToList()).Split(new []{' '},StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)/100).ToList()
                );

                // every 2 lines /*program.string_debug*/($@"ENTROPY->"
                var SECTION_ENTROPY = query.Skip(query.FindIndex(a => string.Equals(a, /*program.string_debug*/($@"SECTION_ENTROPY"), StringComparison.Ordinal))+1).TakeWhile(a => !string.Equals(a, /*program.string_debug*/($@"END_SECTION"), StringComparison.Ordinal)).SkipWhile(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Where(a => !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal)).Select(a => a.Trim()).ToList();

                var SECTION_ENTROPY_joined = (
                    seq: string.Join(/*program.string_debug*/($@""), SECTION_ENTROPY.Where((a, i) => i % 2 == 0).ToList()).Replace(/*program.string_debug*/($@" "), /*program.string_debug*/($@""), StringComparison.Ordinal).ToCharArray(),
                    entropy: string.Join(/*program.string_debug*/($@""), SECTION_ENTROPY.Where((a, i) => i % 2 == 1).ToList()).Replace(/*program.string_debug*/($@"ENTROPY->"), /*program.string_debug*/($@""), StringComparison.Ordinal).Split(new char[]{' '},StringSplitOptions.RemoveEmptyEntries).Select(a=>double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList()
                );

                var joined = SECTION_ENTROPY_joined.seq.Select((a, i) => new info_sable_item()
                {
                    seq_name = query_name,
                    seq_index = i,
                    amino_acid = a,
                    absolute_burial_value = SECTION_SA_ABSOLUTE_joined.absolute_burial_value[i],
                    entropy_value = SECTION_ENTROPY_joined.entropy[i],
                    predicted_ss = SECTION_SS_joined.predicted_ss_seq[i],
                    predicted_ss_confidence = SECTION_SS_joined.confidence_level[i],
                    prob_h = SECTION_SS_PROBABILITIES_joined.prob_h[i],
                    prob_e = SECTION_SS_PROBABILITIES_joined.prob_e[i],
                    prob_c = SECTION_SS_PROBABILITIES_joined.prob_c[i],
                    relative_burial_confidence = SECTION_SA_joined.confidence_level[i],
                    relative_burial_value = SECTION_SA_joined.burial_level[i],

                }).ToList();

                var normalse_ss_prob = true;
                var normalise_entropy = true;
                var normalise_burial_rel = true;
                var normalise_burial_abs = true;

                if (normalse_ss_prob)
                {
                    //var min = joined.Min(a => new double[] { a.prob_h, a.prob_e, a.prob_c }.Min());
                    //var max = joined.Min(a => new double[] { a.prob_h, a.prob_e, a.prob_c }.Max());

                    

                    foreach (var item in joined)
                    {
                        var prob_total = item.prob_h + item.prob_e + item.prob_c;
                        if (prob_total != 0)
                        {
                            item.prob_h /= prob_total;
                            item.prob_e /= prob_total;
                            item.prob_c /= prob_total;
                        }

                       
                    }
                }

                if (normalise_entropy)
                {
                    var min = joined.Min(a => a.entropy_value);
                    var max = joined.Max(a => a.entropy_value);

                }

                if (normalise_burial_rel)
                {
                    var min = joined.Min(a => a.relative_burial_value);
                    var max = joined.Max(a => a.relative_burial_value);

                }

                if (normalise_burial_abs)
                {
                    var min = joined.Min(a => a.absolute_burial_value);
                    var max = joined.Max(a => a.absolute_burial_value);

                }

                result.AddRange(joined);
            }

            return result;

        }
    }
}
