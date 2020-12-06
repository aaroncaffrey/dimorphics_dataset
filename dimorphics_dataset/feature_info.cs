using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace dimorphics_dataset
{
    public class feature_info
    {
        public string alphabet;
        public string stats;
        public int dimension;
        public string category;
        public string source;
        public string group;
        public string member;
        public string perspective;
        public double feature_value;

        public feature_info()
        {

        }

        public feature_info(string csv_line, int csv_offset = 0)
        {
            if (csv_line == null) throw new ArgumentNullException(nameof(csv_line));

            var x = csv_line.Split(',');

            var ix = csv_offset;
            alphabet = x[ix++];
            stats = x[ix++];
            dimension = int.Parse(x[ix++], NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            category = x[ix++];
            source = x[ix++];
            @group = x[ix++];
            member = x[ix++];
            perspective = x[ix++];
            feature_value = double.Parse(x[ix++], NumberStyles.Float, NumberFormatInfo.InvariantInfo);

            //"alphabet,dim?,1,r_peptides,,r_peptides_aIndex,r_peptides_aIndex,default,0"
            //var j = 0;

            //var start = 0;
            //var len = 0;
            //for (var i = 0; i <= csv_line.Length; i++)
            //{
            //    if (i == csv_line.Length || csv_line[i] == ',')
            //    {
            //        if (len > 0)
            //        {
            //            var txt = csv_line.Substring(start, len);
            //
            //            if (j==0 + csv_offset) alphabet = txt;
            //            else if (j==1+csv_offset) dimension = int.Parse(txt, NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            //            else if (j==2+csv_offset) category = txt;
            //            else if (j==3+csv_offset) source = txt;
            //            else if (j==4+csv_offset) @group = txt;
            //            else if (j==5+csv_offset) member = txt;
            //            else if (j==6+csv_offset) perspective = txt;
            //            else if (j==7+csv_offset) feature_value = double.Parse(txt, NumberStyles.Float, NumberFormatInfo.InvariantInfo);
            //        }
            //
            //        j++;
            //        start = i + 1;
            //        len = 0;
            //    }
            //    else
            //    {
            //        len++;
            //    }
            //}


        }

        public feature_info(feature_info feature_info)
        {
            if (feature_info == null)
            {
                throw new ArgumentNullException(nameof(feature_info));
            }

            this.alphabet = feature_info.alphabet;
            this.stats = feature_info.stats;
            this.dimension = feature_info.dimension;
            this.category = feature_info.category;
            this.source = feature_info.source;
            this.@group = feature_info.@group;
            this.member = feature_info.member;
            this.perspective = feature_info.perspective;
            this.feature_value = feature_info.feature_value;
        }

        public bool validate()
        {
            var ok = (
                !string.IsNullOrWhiteSpace(alphabet) &&
                /*!string.IsNullOrWhiteSpace(stats) &&*/
                (dimension >= 0 && dimension <= 3) &&
                !string.IsNullOrWhiteSpace(source) &&
                !string.IsNullOrWhiteSpace(group) &&
                !string.IsNullOrWhiteSpace(member) &&
                !string.IsNullOrWhiteSpace(perspective) &&
                (!double.IsNaN(feature_value) && !double.IsInfinity(feature_value)));

            return ok;
        }

        public (string key, string value)[] key_value_list(bool include_value = false)
        {
            if (!include_value)
            {
                return new (string key, string value)[]
                {
                    (nameof(alphabet), alphabet),
                    (nameof(stats), stats),
                    (nameof(dimension), dimension.ToString(CultureInfo.InvariantCulture)),
                    (nameof(category), category),
                    (nameof(source), source),
                    (nameof(@group), @group),
                    (nameof(member), member),
                    (nameof(perspective), perspective)
                };
            }
            else
            {
                return new (string key, string value)[]
                {
                    (nameof(alphabet), alphabet),
                    (nameof(stats), stats),
                    (nameof(dimension), dimension.ToString(CultureInfo.InvariantCulture)),
                    (nameof(category), category),
                    (nameof(source), source),
                    (nameof(@group), @group),
                    (nameof(member), member),
                    (nameof(perspective), perspective),
                    (nameof(feature_value), $@"{feature_value:G17}")
                };
            }
        }


        public string key_value_list_hr()
        {
            var data = key_value_list();

            return string.Join($@", ", data);
        }

        public override string ToString()
        {
            return key_value_list_hr();
        }

        public static List<string> get_feature_headers_lines_csv(List<feature_info> feature_list)
        {
            //var header_list = feature_list.Select(a => new feature_info(a) { feature_value = 0 }).ToList();

            List<string> headers = new List<string> { $@"fid,{nameof(feature_info.alphabet)},{nameof(feature_info.stats)},{nameof(feature_info.dimension)},{nameof(feature_info.category)},{nameof(feature_info.source)},{nameof(feature_info.@group)},{nameof(feature_info.member)},{nameof(feature_info.perspective)}" };
            headers.AddRange(feature_list.Select((a, fid) => $@"{fid},{a.alphabet},{a.stats},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList());

            return headers;

            // check for duplicately named features
            //var header_list_str = header_list.Select((a, fid) => $@"{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();
            //var dupes = find_duplicate_strings(header_list_str);
            //if (dupes != null && dupes.Count > 0)
            //{
            //    throw new Exception($@"{module_name}.{method_name}: Duplicate headers found: " + string.Join($@", ", dupes));
            //}

            //var header_list_str_c = new List<string>();
            //header_list_str_c.Add($@"fid,{nameof(subsequence_classification_data.feature_info.alphabet)},{nameof(subsequence_classification_data.feature_info.dimension)},{nameof(subsequence_classification_data.feature_info.category)},{nameof(subsequence_classification_data.feature_info.source)},{nameof(subsequence_classification_data.feature_info.@group)},{nameof(subsequence_classification_data.feature_info.member)},{nameof(subsequence_classification_data.feature_info.perspective)}");
            //header_list_str_c.AddRange(header_list.Select((a, fid) => $@"{fid.ToString(CultureInfo.InvariantCulture)},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList());

            //return header_list_str_c;
        }
    }
}
