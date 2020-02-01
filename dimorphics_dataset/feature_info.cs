using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace dimorphics_dataset
{
    public class feature_info
    {
        public string alphabet;
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

        public feature_info(feature_info feature_info)
        {
            if (feature_info == null)
            {
                throw new ArgumentNullException(nameof(feature_info));
            }

            this.alphabet = feature_info.alphabet;
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
            var ok = (!string.IsNullOrWhiteSpace(alphabet) &&
                (dimension >= 0 && dimension <= 3) &&
                !string.IsNullOrWhiteSpace(source) &&
                !string.IsNullOrWhiteSpace(group) &&
                !string.IsNullOrWhiteSpace(member) &&
                !string.IsNullOrWhiteSpace(perspective) &&
                (!double.IsNaN(feature_value) && !double.IsInfinity(feature_value)));

            return ok;
        }



        public override string ToString()
        {
            //var header_list_str_c = header_list.Select((a, fid) => $"{fid.ToString(CultureInfo.InvariantCulture)},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();


            var data = new List<(string key, string value)>
                {
                    (nameof(alphabet), alphabet),
                    (nameof(dimension), dimension.ToString(CultureInfo.InvariantCulture)),
                    (nameof(category), category),
                    (nameof(source), source),
                    (nameof(@group), @group),
                    (nameof(member), member),
                    (nameof(perspective), perspective)
                };

            //data.Add((nameof(feature_value), feature_value.ToString("G17", CultureInfo.InvariantCulture)));

            return string.Join(", ", data);
        }

        public static List<string> get_feature_headers_lines_csv(List<feature_info> feature_list)
        {
            //var header_list = feature_list.Select(a => new feature_info(a) { feature_value = 0 }).ToList();

            List<string> headers = new List<string>();
            headers.Add($@"fid,{nameof(feature_info.alphabet)},{nameof(feature_info.dimension)},{nameof(feature_info.category)},{nameof(feature_info.source)},{nameof(feature_info.@group)},{nameof(feature_info.member)},{nameof(feature_info.perspective)}");
            headers.AddRange(feature_list.Select((a, fid) => $"{fid},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList());

            return headers;

            // check for duplicately named features
            //var header_list_str = header_list.Select((a, fid) => $"{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();
            //var dupes = find_duplicate_strings(header_list_str);
            //if (dupes != null && dupes.Count > 0)
            //{
            //    throw new Exception("Duplicate headers found: " + string.Join(", ", dupes));
            //}

            //var header_list_str_c = new List<string>();
            //header_list_str_c.Add($@"fid,{nameof(subsequence_classification_data.feature_info.alphabet)},{nameof(subsequence_classification_data.feature_info.dimension)},{nameof(subsequence_classification_data.feature_info.category)},{nameof(subsequence_classification_data.feature_info.source)},{nameof(subsequence_classification_data.feature_info.@group)},{nameof(subsequence_classification_data.feature_info.member)},{nameof(subsequence_classification_data.feature_info.perspective)}");
            //header_list_str_c.AddRange(header_list.Select((a, fid) => $"{fid.ToString(CultureInfo.InvariantCulture)},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList());

            //return header_list_str_c;
        }
    }
}
