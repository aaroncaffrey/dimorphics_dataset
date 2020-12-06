using System;
using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal static class info_dssp
    {
        internal static List<info_dssp_item> Load(string file)
        {
            const string SsDataMarker1 = @"  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA ";
            const string SsDataMarker2 = @"  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN";

            var lines = io_proxy.ReadAllLines(file, nameof(info_dssp), nameof(Load));
            
            var dssp_record_list = new List<info_dssp_item>();

            var marker_found = false;

            foreach (var line in lines)
            {
                if (!marker_found && (string.Equals(line, SsDataMarker1, StringComparison.Ordinal) || string.Equals(line, SsDataMarker2, StringComparison.Ordinal)))
                {
                    marker_found = true;
                    continue;
                }

                if (marker_found)
                {
                    var dr = new info_dssp_item(line);
                    dssp_record_list.Add(dr);
                }
            }

            return dssp_record_list;
        }
    }
}