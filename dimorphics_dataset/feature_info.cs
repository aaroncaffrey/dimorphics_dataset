using System;
using System.Collections.Generic;
using System.Globalization;
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
    }
}
