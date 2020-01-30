using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class data_controller
    {
        //public static List<subsequence_classification_data> extract_features_list
        //      (
        //      string dataset_name,
        //      program.class_info class_info,
        //      List<subsequence_classification_data> data_list,
        //      subsequence_classification_data.feature_types feature_types,
        //      int max_features
        //  )
        //{
        //    if (string.IsNullOrWhiteSpace(dataset_name))
        //    {
        //        throw new Exception();
        //    }

        //    var output_folder = $@"F:\dataset\{dataset_name}\";


        //    io_proxy.WriteLine($"{class_info.class_name}: Started: encoding and saving class #{class_info.class_id} {class_info.class_name}");
        //    io_proxy.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");

        //    var all_output_data = new List<string>();
        //    var comment_rows = new List<string>();
        //    var sw2 = new Stopwatch();
        //    sw2.Start();



        //    var results = new List<(int row_index, program.class_info class_info, string features_filename, string headers_list_filename, string comments_filename)>();
        //    var wait_files = new List<string>();
        //    for (var index = 0; index < data_list.Count; index++)
        //    {
        //        var data_list_row = data_list[index];

        //        var files = extract_features_item(output_folder, index, class_info, feature_types, max_features, data_list_row);

        //        results.Add((index, class_info, files.features_filename, files.headers_list_filename, files.comments_filename));

        //        wait_files.Add(files.features_filename);
        //        wait_files.Add(files.comments_filename);
        //        wait_files.Add(files.headers_list_filename);
        //    }

        //    sw2.Stop();
        //    var ts2 = sw2.Elapsed;

        //    wait_for_files(wait_files);

        //    // combine features files
        //    var feature_files = results.OrderBy(a => a.row_index).Select(a => a.features_filename).ToList();
        //    var feature_data = feature_files.Select((a, i) => File.ReadAllLines(a).Skip(i == 0 ? 0 : 1).ToList()).ToList();
        //    if (feature_data.Any(a => a.Count < 1)) throw new Exception();
        //    var all_features_file = Path.Combine(output_folder, $"f__[{class_info.class_name}].csv");
        //    io.WriteAllLines(all_features_file, all_output_data, nameof(program), nameof(extract_features_list));
        //    io.WriteLine("Saved: " + all_features_file);


        //    // combine comments files
        //    var comments_files = results.OrderBy(a => a.row_index).Select(a => a.comments_filename).ToList();
        //    var comments_data = comments_files.Select((a, i) => File.ReadAllLines(a).Skip(i == 0 ? 0 : 1).ToList()).ToList();
        //    if (comments_data.Any(a => a.Count < 1)) throw new Exception();
        //    var comment_file = Path.Combine(output_folder, $"c__[{class_info.class_name}].csv");
        //    io.WriteAllLines(comment_file, comment_rows, nameof(program), nameof(extract_features_list));
        //    io.WriteLine("Saved: " + comment_file);


        //    // combine headers files
        //    var headers_files = results.OrderBy(a => a.row_index).Select(a => a.headers_list_filename).ToList();
        //    var headers = headers_files.Select(a => File.ReadAllLines(a)).ToList();
        //    for (var i = 0; i < headers.Count; i++) for (var j = 0; j < headers.Count; j++) if (i < j) if (!headers[i].SequenceEqual(headers[j])) throw new Exception();
        //    var header_data = headers.First();


        //    io_proxy.WriteLine($"{class_info.class_name}: ({ts2.Days}d {ts2.Hours}h {ts2.Minutes}m {ts2.Seconds}s)");
        //    io_proxy.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");
        //    io_proxy.WriteLine($"{class_info.class_name}: Finished: encoding and saving class #{class_info.class_id} {class_info.class_name}");
        //    io_proxy.WriteLine();


        //    return data_list;
        //}

        //public static void wait_for_files(List<string> filenames)
        //{

        //}

        //public static (string features_filename, string headers_list_filename, string comments_filename) extract_features_item
        //(
        //    string output_folder,
        //    int row_index,
        //    program.class_info class_info,
        //    subsequence_classification_data.feature_types feature_types,
        //    int max_features,
        //    subsequence_classification_data row
        //)
        //{
        //    var features_filename = Path.Combine(output_folder, $@"features_header_{class_info.class_name}_{row_index}.csv");
        //    var headers_list_filename = Path.Combine(output_folder, $@"features_header_{class_info.class_name}_{row_index}.csv");
        //    var comments_filename = Path.Combine(output_folder, $@"features_header_{class_info.class_name}_{row_index}.csv");

        //    if (!File.Exists(features_filename) || !File.Exists(headers_list_filename) || !File.Exists(comments_filename))
        //    {
        //        // save features (incl. header row of fids)
        //        var row_encoded = subsequence_classification_data.encode_subsequence_classification_data_row(row, max_features, feature_types);


        //        var row_feature_values = row_encoded.feature_info.Select((a, fid) => a.feature_value.ToString("G17", CultureInfo.InvariantCulture)).ToList();
        //        var row_feature_header_str = string.Join(",", Enumerable.Range(0, row_feature_values.Count));
        //        var row_feature_values_str = string.Join(",", row_feature_values);
        //        var row_features_file_data = new List<string>();
        //        row_features_file_data.Add(row_feature_header_str);
        //        row_features_file_data.Add(row_feature_values_str);

        //        io.WriteAllLines(features_filename, row_features_file_data, nameof(program), nameof(extract_features_item));



        //        // save feature list (incl. header row)
        //        var feature_headers_file_data = get_feature_list_csv(row_encoded);

        //        io.WriteAllLines(headers_list_filename, feature_headers_file_data, nameof(program), nameof(extract_features_item));



        //        // save comments (incl. header row)
        //        var row_comments_file_data = get_row_comments(row_index, row_encoded);

        //        io.WriteAllLines(comments_filename, row_comments_file_data, nameof(program), nameof(extract_features_item));

        //    }

        //    return (features_filename, headers_list_filename, comments_filename);
        //    //return row_encoded;
        //}

        

        //public static List<string> find_duplicate_strings(List<string> query)
        //{
        //    if (query == null || query.Count == 0) return null;

        //    var query_distinct = query.Distinct().ToList();

        //    if (query.Count == query_distinct.Count) return null;

        //    var ret = query_distinct.Where(a => query.Count(b => a == b) > 1).ToList();

        //    return ret;
        //}

        public static List<string> get_feature_headers_lines_csv((subsequence_classification_data instance_meta_data, List<feature_info> feature_info) row_encoded)
        {
            var header_list = row_encoded.feature_info.Select(a => new feature_info(a) { feature_value = 0 }).ToList();

            List<string> headers = new List<string>();
            headers.Add($@"fid,{nameof(feature_info.alphabet)},{nameof(feature_info.dimension)},{nameof(feature_info.category)},{nameof(feature_info.source)},{nameof(feature_info.@group)},{nameof(feature_info.member)},{nameof(feature_info.perspective)}");
            headers.AddRange(header_list.Select((a, fid) => $"{fid},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList());
            
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



        public static List<string> get_row_comments_headers((subsequence_classification_data instance_meta_data, List<feature_info> feature_info) row_encoded)
        {
            var comment_headers = new List<string>();

            comment_headers.Add("row_index");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.pdb_id)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.chain_id)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.dimer_type)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.class_id)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.class_name)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.parallelism)}");
            comment_headers.Add($"{nameof(row_encoded.instance_meta_data.symmetry_mode)}");

            foreach (var region in row_encoded.instance_meta_data.get_regions())
            {
                comment_headers.Add($"{region.region_name}_{nameof(region.region.aa_sequence.Length)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.aa_sequence)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.res_ids)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.dssp3_monomer)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.dssp3_multimer)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.dssp_monomer)}");
                comment_headers.Add($"{region.region_name}_{nameof(region.region.dssp_multimer)}");
                comment_headers.AddRange(region.region.ss_predictions.Select(a => $@"{region.region_name}_{a.format}").ToList());
            }

            return comment_headers;
        }

        public static List<string> get_row_comments(int row_index, (subsequence_classification_data instance_meta_data, List<feature_info> feature_info) row_encoded)
        {
            var row_comments = new List<string>();

            row_comments.Add(row_index.ToString(CultureInfo.InvariantCulture));
            row_comments.Add($"{(row_encoded.instance_meta_data.pdb_id)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.chain_id)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.dimer_type)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.class_id)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.class_name)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.parallelism)}");
            row_comments.Add($"{(row_encoded.instance_meta_data.symmetry_mode)}");

            foreach (var region in row_encoded.instance_meta_data.get_regions())
            {
                row_comments.Add($"{(region.region.aa_sequence.Length)}");
                row_comments.Add($"{(region.region.aa_sequence)}");
                row_comments.Add($"{(string.Join(" ", region.region.res_ids.Select(a => $@"{a.amino_acid}{a.residue_index}{(a.i_code != default && a.i_code != ' ' ? a.i_code.ToString(CultureInfo.InvariantCulture) : "")}").ToList()))}");
                row_comments.Add($"{(region.region.dssp3_monomer)}");
                row_comments.Add($"{(region.region.dssp3_multimer)}");
                row_comments.Add($"{(region.region.dssp_monomer)}");
                row_comments.Add($"{(region.region.dssp_multimer)}");
                row_comments.AddRange(region.region.ss_predictions.Select(a => a.prediction).ToList());
            }
            
            return row_comments;
        }

    }
}
