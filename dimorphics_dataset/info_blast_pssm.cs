using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class info_blast_pssm
    {
        internal static string run_psi_blast_get_pssm(string seq_filename, string pssm_filename, string blast_database, string database_folder, int num_iterations = 1, string evalue = @"0.1", string inclusion_ethresh = @"0.1", bool remote = true, bool run = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(run_psi_blast_get_pssm);

            //run_psi_blast_get_pssm: stdout: Database: Non - redundant UniProtKB / SwissProt sequences
            //run_psi_blast_get_pssm: stdout: 468,734 sequences; 176,603,381 total letters
            //run_psi_blast_get_pssm: stdout:
            //run_psi_blast_get_pssm: stdout:

            //run_psi_blast_get_pssm: stderr: Warning: [psiblast] Error initializing remote BLAST database data loader: Protein BLAST database 'db/swissprot nr' does not exist in the NCBI servers
            //run_psi_blast_get_pssm: stderr: BLAST query/ options error: GI / accession / sequence mismatch: protein input required but nucleotide provided
            //run_psi_blast_get_pssm: stderr: Please refer to the BLAST + user manual.


            //psiblast_args.Add($@"-query query/{Path.GetFileName(seq_filename)}");
            //psiblast_args.Add($@"-db $@""c:\blast\bin\db\{blast_database}""");


            // get pssm for each globin
            var psiblast_exe_filename = @"c:\blast\psiblast.exe";

            var num_threads = Environment.ProcessorCount;

            var psiblast_args = new List<string>
            {
                $@"-query $@""{seq_filename}""",
                $@"-db {(remote ? blast_database : Path.Combine(database_folder, blast_database))}",
                $@"-evalue={evalue}",
                $@"-inclusion_ethresh={inclusion_ethresh}",
                $@"-out_ascii_pssm=""{pssm_filename}"""
            };
            

            if (remote)
            {
                psiblast_args.Add($@"-remote");
            }
            else
            {
                psiblast_args.Add($@"-num_iterations={num_iterations}");
                psiblast_args.Add($@"-num_threads={num_threads}");
            }




            var start = new ProcessStartInfo
            {
                FileName = psiblast_exe_filename,
                WorkingDirectory = Path.GetDirectoryName(psiblast_exe_filename) ?? $@"",
                Arguments = string.Join($@" ", psiblast_args),
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                WindowStyle = ProcessWindowStyle.Hidden
            };

            if (run)
            {

                io_proxy.WriteLine($@"{module_name}.{method_name}: run: ""{start.FileName}"" {start.Arguments}", module_name, method_name);

                using var process = Process.Start(start);

                if (process == null)
                {
                    throw new Exception($@"{module_name}.{method_name}: {nameof(process)} is null");
                }


                using (var reader = process.StandardOutput)
                {
                    var stdout = reader.ReadToEnd();
                    if (!string.IsNullOrWhiteSpace(stdout))
                    {
                        stdout = $"\r\n{stdout}".Replace($"\r\n", $"\r\n{module_name}.{method_name}: {nameof(stdout)}: ", StringComparison.Ordinal);
                        io_proxy.WriteLine(stdout, module_name, method_name);
                    }
                }

                using (var reader = process.StandardError)
                {
                    var stderr = reader.ReadToEnd();

                    if (!string.IsNullOrWhiteSpace(stderr))
                    {
                        stderr = $"\r\n{stderr}".Replace($"\r\n", $"\r\n{module_name}.{method_name}: {nameof(stderr)}: ", StringComparison.Ordinal);
                        io_proxy.WriteLine(stderr, module_name, method_name);
                    }
                }

                process.WaitForExit();
            }

            return $@"{start.FileName} {start.Arguments}";
        }

        internal static List<info_blast_pssm_entry> normalise_pssm(List<info_blast_pssm_entry> pssm)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(normalise_pssm);

            if (pssm == null || pssm.Count == 0)
            {
                return new List<info_blast_pssm_entry>();
            }

            pssm = pssm.Select(a => new info_blast_pssm_entry(a)).ToList();

            var scores = pssm.Select(a => a.score).Distinct().ToList();
            var min_score = scores.Min();
            var max_score = scores.Max();

            for (var i = 0; i < pssm.Count; i++)
            {
                if (max_score - min_score == 0) pssm[i].score = 0;
                else pssm[i].score = (pssm[i].score - min_score) / (max_score - min_score);
            }

            return pssm;
        }
        internal static double[] normalise_array(double[] pssm)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(normalise_array);

            if (pssm == null || pssm.Length == 0) return Array.Empty<double>();

            pssm = pssm.ToArray();

            var scores = pssm.Select(a => a).Distinct().ToList();
            var min_score = scores.Min();
            var max_score = scores.Max();

            for (var i = 0; i < pssm.Length; i++)
            {
                if (max_score - min_score == 0) pssm[i] = 0;
                else pssm[i] = (pssm[i] - min_score) / (max_score - min_score);
            }

            return pssm;
        }

        internal static double[] normalise_array(double[] pssm, double min_score, double max_score)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(normalise_array);

            if (pssm == null || pssm.Length == 0) return Array.Empty<double>();

            pssm = pssm.ToArray();

            for (var i = 0; i < pssm.Length; i++)
            {
                if (max_score - min_score == 0) pssm[i] = 0;
                else pssm[i] = (pssm[i] - min_score) / (max_score - min_score);
            }

            return pssm;
        }

        internal static double[] intervals(double[] pssm, bool normalise = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(intervals);
            // note: data always input in the same meaningful order - so considered pre-ordered

            if (pssm == null || pssm.Length == 0)
            {
                throw new ArgumentNullException(nameof(pssm));
            }

            var x = new double[pssm.Length - 1];

            for (var i = 1; i < pssm.Length; i++)
            {
                x[i - 1] = Math.Abs(pssm[i] - pssm[i - 1]);
            }

            if (normalise) x = normalise_array(x);

            return x;
        }

        internal static double[] distances(double[] pssm, bool normalise = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(distances);

            //var x = pssm.SelectMany((a, i) => pssm.Where((b, j) => i < j).Select((b, j) => Math.Abs(a - b)).ToArray()).ToArray();

            // return the distances between the values

            var x = pssm.SelectMany((a, i) => pssm.Where((b, j) => i < j).Select((b, j) => Math.Abs(a - b)).ToArray()).ToArray();

            if (normalise) x = normalise_array(x);

            return x;
        }

        internal static List<info_blast_pssm_entry> load_psi_blast_pssm(string pssm_filename, bool normalise_pssm = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(load_psi_blast_pssm);


            var pssm = new List<info_blast_pssm_entry>();

            if (!File.Exists(pssm_filename) || new FileInfo(pssm_filename).Length == 0)
            {
                return new List<info_blast_pssm_entry>();
            }

            var line_list = io_proxy.ReadAllLines(pssm_filename, module_name, method_name).Skip(2).ToList();
            var pssm_end_index = line_list.IndexOf($@"");
            if (pssm_end_index > -1) line_list = line_list.GetRange(0, pssm_end_index);

            line_list[0] = $@"ID AA {line_list[0]} W1 W2";



            //lineList = lineList.Skip(1).ToList();

            var line_list_split = line_list.Select(a => a.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).ToList()).ToList();

            line_list_split = line_list_split.Select(a => a.Take(22).ToList()).ToList();

            for (var row = 1; row < line_list_split.Count; row++)
            {
                for (var col = 2; col < line_list_split[row].Count; col++)
                {
                    if (!double.TryParse(line_list_split[row][col], out var score))
                    {
                        score = 0;
                    }

                    var col_aa_provided = line_list_split[0][col][0];
                    var row_aa_provided = line_list_split[row][1][0];
                    var row_index_provided = int.Parse(line_list_split[row][0], NumberStyles.Integer, NumberFormatInfo.InvariantInfo);

                    var e = new info_blast_pssm_entry()
                    {
                        query_sequence_aa_pos = row_index_provided,
                        query_sequence_aa = row_aa_provided,
                        matrix_row_index = row - 1,
                        matrix_column_index = col - 2,
                        position_aa = col_aa_provided,
                        score = score
                    };

                    pssm.Add(e);
                }
            }

            if (normalise_pssm)
            {
                pssm = dimorphics_dataset.info_blast_pssm.normalise_pssm(pssm);
            }

            return pssm;
        }

        //internal static List<(int id, string name, List<string> groups)> alphabets = new List<(int id, string name, List<string> groups)>()
        //{
        //    (0, $@"Normal",new List<string>(){
        //        $@"A", $@"R", $@"N", $@"D", $@"C", $@"Q", $@"E", $@"G", $@"H", $@"I", $@"L", $@"K", $@"M", $@"F", $@"P", $@"S", $@"T", $@"W", $@"Y", $@"V"
        //    }),
        //    (1, $@"Physicochemical",new List<string>(){
        //        $@"AVFPMILW", $@"DE", $@"RK", $@"STYHCNGQ"
        //    }),
        //    (2, $@"Hydrophobicity",new List<string>(){
        //        $@"AGILPV", $@"FYW", $@"DENQRHSTK", $@"CM"
        //    }),
        //    (3, $@"UniProtKb",new List<string>(){
        //        $@"LAGVIP", $@"DE", $@"ST", $@"RKH", $@"FYW", $@"NQ", $@"CM"
        //    }),
        //    (4, $@"PdbSum",new List<string>(){
        //        $@"HKR", $@"DE", $@"STNQ", $@"AVLIM", $@"FYW", $@"PG", $@"C"
        //    }),
        //    (5, $@"Venn",new List<string>(){
        //        $@"ILV", $@"TS", $@"AGCS", $@"VPAGCSTDN", $@"NQ", $@"HKR", $@"DE", $@"CSTNQDEHKRYW", $@"ILVAMFGP", $@"DEHKR", $@"ILVACTMFYWHK", $@"FYWH", $@"MC"
        //    })
        //};

     

        internal static double[] pssm_to_vector1(List<info_blast_pssm_entry> pssm, enum_pssm_value_type pssm_value_type, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector1);

            var x = (pssm == null || pssm.Count == 0) ? new[] { 0d } : pssm.Select(a => a.score).ToArray();

            if (x == null || x.Length == 0) x = new[] { 0d };

            if (normalise_all) x = normalise_array(x);
            //var sum = scores.Sum();
            //var average = scores.Average();

            if (pssm_value_type == enum_pssm_value_type.distances)
            {
                x = distances(x, normalise_all);
            }
            else if (pssm_value_type == enum_pssm_value_type.intervals)
            {
                x = intervals(x, normalise_all);
            }

            return x; //(scores, sum, average);
        }

        internal static List<(string alphabet, List<(string col_aa, double[] values)> x)> pssm_to_vector20col(List<info_blast_pssm_entry> pssm, enum_pssm_value_type pssm_value_type, bool normalise_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector20col);

            // join all cols with same amino acid, calc average
            //var cols = pssm.Select(a => a.position_aa).Distinct().ToList();
            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string col_aa, double[] values)> x)>();

            foreach (var alphabet in feature_calcs.aa_alphabets)
            {
                var alphabet_groups = alphabet.groups;

                var vector = new List<(string col_aa, double[] values)>();

                foreach (var c1 in alphabet_groups)
                {
                    double[] x = (pssm == null || pssm.Count == 0) ? new[] { 0d } : pssm.Where(a => c1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal)).Select(a => a.score).ToArray();

                    if (x == null || x.Length == 0) x = new[] { 0d };

                    if (normalise_col) x = normalise_array(x);
                    if (x == null || x.Length == 0) x = new[] { 0d };


                    if (pssm_value_type == enum_pssm_value_type.distances)
                    {
                        x = distances(x, normalise_col);
                    }
                    else if (pssm_value_type == enum_pssm_value_type.intervals)
                    {
                        x = intervals(x, normalise_col);
                    }
                    if (x == null || x.Length == 0) x = new[] { 0d };

                    //var sum = x.Sum();
                    //var average = x.Average();
                    vector.Add((c1.group_amino_acids, x));
                }

                if (normalise_all)
                {
                    var values = vector.SelectMany(a => a.values).Distinct().ToList();
                    var min = values.Min();
                    var max = values.Max();
                    vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                }


                result.Add((alphabet.name, vector));
            }

            return result;
        }

        internal static List<(string alphabet, List<(string col_aa, int lag, double[] values)> x)> pssm_to_vector20col_DT(List<info_blast_pssm_entry> pssm, int max_lag, enum_pssm_value_type pssm_value_type, bool normalise_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector20col_DT);

            // join all cols with same amino acid, calc average
            //var cols = pssm.Select(a => a.position_aa).Distinct().ToList();
            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string col_aa, int lag, double[] values)> x)>();

            var tasks = new List<Task<(string name, List<(string col_aa, int lag, double[] values)> vector)>>();
            var tasks_start_time = DateTime.Now;
            
            foreach (var l_alphabet in feature_calcs.aa_alphabets)
            {
                var alphabet = l_alphabet;

                for (var l_lag = 1; l_lag <= max_lag; l_lag++)
                {
                    var lag = l_lag;

                    Task<(string name, List<(string col_aa, int lag, double[] values)> vector)> task = Task.Run(() =>
                    {
                        var alphabet_groups = alphabet.groups;

                        var vector = new List<(string col_aa, int lag, double[] values)>();

                        for (var r1_index = 0; r1_index < alphabet_groups.Count; r1_index++)
                        {
                            var r1 = alphabet_groups[r1_index];

                            double[] x = new[] {0d};

                            if (pssm != null && pssm.Count > 0)
                            {
                                var col_pssm = pssm.Where(a => r1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal)).ToList();
                                // z = list of column 

                                //for (var r1_index = 0; r1_index < alphabet_groups.Count; r1_index++)
                                //{
                                //    var r1 = alphabet_groups[r1_index];

                                for (var r2_index = 0; r2_index < alphabet_groups.Count; r2_index++)
                                {
                                    var r2 = alphabet_groups[r2_index];

                                    if (r2_index < r1_index) continue;

                                    var r1_match = col_pssm.Where(a => r1.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)).ToList();
                                    var r2_match = col_pssm.Where(a => r2.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)).ToList();
                                    var join = r1_match.Union(r2_match).ToList();
                                    var z = @join.Where(a => @join.Any(b => Math.Abs(a.matrix_row_index - b.matrix_row_index) == lag)).ToList();
                                    if (z.Count > 0)
                                    {
                                        x = z.Select(a => a.score).ToArray();
                                    }
                                }
                            }
                            //}

                            if (x == null || x.Length == 0) x = new[] {0d};

                            if (normalise_col) x = normalise_array(x);

                            if (x == null || x.Length == 0) x = new[] {0d};

                            if (pssm_value_type == enum_pssm_value_type.distances)
                            {
                                x = distances(x, normalise_col);
                            }
                            else if (pssm_value_type == enum_pssm_value_type.intervals)
                            {
                                x = intervals(x, normalise_col);
                            }

                            if (x == null || x.Length == 0) x = new[] {0d};

                            //var sum = x.Sum();
                            //var average = x.Average();
                            vector.Add((r1.group_amino_acids, lag, x));
                        }

                        if (normalise_all)
                        {
                            var values = vector.SelectMany(a => a.values).Distinct().ToList();
                            var min = values.Min();
                            var max = values.Max();
                            vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                        }


                        //result.Add((alphabet.name, vector));
                        return (alphabet.name, vector);

                    });
                    tasks.Add(task);
                }
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, module_name, method_name);


            tasks.ForEach(a=>result.Add(a.Result));

            return result;
        }

        internal static List<(string alphabet, List<(string row_aa, double[] values)> x)> pssm_to_vector20row(List<info_blast_pssm_entry> pssm, enum_pssm_value_type pssm_value_type, bool normalise_row = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector20row);

            // join all rows with the same amino acid, calculate average
            // note: it doesn't make sense to have DT for 20row

            //var cols = pssm.Select(a => a.position_aa).Distinct().ToList();
            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string row_aa, double[] values)> x)>();

            foreach (var alphabet in feature_calcs.aa_alphabets)
            {
                var alphabet_groups = alphabet.groups;

                var vector = new List<(string row_aa, double[] values)>();

                foreach (var c1 in alphabet_groups)
                {
                    double[] x = (pssm == null || pssm.Count == 0) ? new[] { 0d } : pssm.Where(a => c1.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)).Select(a => a.score).ToArray();
                    if (x == null || x.Length == 0) x = new[] { 0d };
                    if (normalise_row) x = normalise_array(x);
                    if (x == null || x.Length == 0) x = new[] { 0d };

                    if (pssm_value_type == enum_pssm_value_type.distances)
                    {
                        x = distances(x, normalise_row);
                    }
                    else if (pssm_value_type == enum_pssm_value_type.intervals)
                    {
                        x = intervals(x, normalise_row);
                    }

                    if (x == null || x.Length == 0) x = new[] { 0d };


                    //var sum = x.Sum();
                    //var average = x.Average();
                    vector.Add((c1.group_amino_acids, x));
                }

                if (normalise_all)
                {
                    var values = vector.SelectMany(a => a.values).Distinct().ToList();
                    var min = values.Min();
                    var max = values.Max();
                    vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                }

                result.Add((alphabet.name, vector));
            }

            return result;
        }

        internal static List<(string alphabet, List<(string row_aa, string col_aa, double[] values)> x)> pssm_to_vector210(List<info_blast_pssm_entry> pssm, enum_pssm_value_type pssm_value_type, bool normalise_row_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector210);
            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string row_aa, string col_aa, double[] values)> x)>();

            foreach (var alphabet in feature_calcs.aa_alphabets)
            {
                var vector = new List<(string row_aa, string col_aa, double[] values)>();


                for (var c1_index = 0; c1_index < alphabet.groups.Count; c1_index++)
                {
                    var c1 = alphabet.groups[c1_index];
                    for (var c2_index = 0; c2_index < alphabet.groups.Count; c2_index++)
                    {
                        var c2 = alphabet.groups[c2_index];
                        if (c2_index < c1_index) continue;

                        double[] x = (pssm == null || pssm.Count == 0) ? new[] { 0d } : pssm.Where(a => (c1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c2.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)) || (c2.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c1.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal))).Select(a => a.score).ToArray();
                        if (x == null || x.Length == 0) x = new[] { 0d };
                        if (normalise_row_col) x = normalise_array(x);

                        if (x == null || x.Length == 0) x = new[] { 0d };


                        if (pssm_value_type == enum_pssm_value_type.distances)
                        {
                            x = distances(x, normalise_row_col);
                        }
                        else if (pssm_value_type == enum_pssm_value_type.intervals)
                        {
                            x = intervals(x, normalise_row_col);
                        }

                        if (x == null || x.Length == 0) x = new[] { 0d };

                        vector.Add((c2.group_amino_acids, c1.group_amino_acids, x));
                    }
                }


                if (normalise_all)
                {
                    var values = vector.SelectMany(a => a.values).Distinct().ToList();
                    var min = values.Min();
                    var max = values.Max();
                    vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                }

                result.Add((alphabet.name, vector));
            }

            return result;
        }

        internal static List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)> x)> pssm_to_vector210_DT(List<info_blast_pssm_entry> pssm, int max_lag, enum_pssm_value_type pssm_value_type, bool normalise_row_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector210_DT);

            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)> x)>();

            var tasks = new List<Task<(string name, List<(string row_aa, string col_aa, int lag, double[] values)> vector)>>();
            var tasks_start_time = DateTime.Now;
            
            foreach (var l_alphabet in feature_calcs.aa_alphabets)
            {
                var alphabet = l_alphabet;

                for (var l_lag = 1; l_lag <= max_lag; l_lag++)
                {
                    var lag = l_lag;

                    Task<(string name, List<(string row_aa, string col_aa, int lag, double[] values)> vector)> task = Task.Run(() =>
                    {
                        var vector = new List<(string row_aa, string col_aa, int lag, double[] values)>();


                        for (var c1_index = 0; c1_index < alphabet.groups.Count; c1_index++)
                        {
                            var c1 = alphabet.groups[c1_index];
                            for (var c2_index = 0; c2_index < alphabet.groups.Count; c2_index++)
                            {
                                var c2 = alphabet.groups[c2_index];
                                if (c2_index < c1_index) continue;




                                double[] x = new[] {0d};

                                if (pssm != null && pssm.Count > 0)
                                {
                                    var z = pssm.Where(a => (c1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c2.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)) || (c2.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c1.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal))).ToList();
                                    z = z.OrderBy(a => a.matrix_row_index).ToList();
                                    z = z.Where(a => z.Any(b => Math.Abs(a.matrix_row_index - b.matrix_row_index) == lag)).ToList();
                                    if (z.Count > 0)
                                    {
                                        x = z.Select(a => a.score).ToArray();
                                    }
                                }

                                if (x == null || x.Length == 0) x = new[] {0d};

                                if (normalise_row_col) x = normalise_array(x);

                                if (x == null || x.Length == 0) x = new[] {0d};

                                if (pssm_value_type == enum_pssm_value_type.distances)
                                {
                                    x = distances(x, normalise_row_col);
                                }
                                else if (pssm_value_type == enum_pssm_value_type.intervals)
                                {
                                    x = intervals(x, normalise_row_col);
                                }

                                if (x == null || x.Length == 0) x = new[] {0d};


                                vector.Add((c2.group_amino_acids, c1.group_amino_acids, lag, x));
                            }
                        }


                        if (normalise_all)
                        {
                            var values = vector.SelectMany(a => a.values).Distinct().ToList();
                            var min = values.Min();
                            var max = values.Max();
                            vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                        }

                       //result.Add((alphabet.name, vector));

                        return (alphabet.name, vector);
                    });
                    tasks.Add(task);
                }
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, module_name, method_name);


            tasks.ForEach(a=>result.Add(a.Result));

            return result;
        }

        internal static List<(string alphabet, List<(string row_aa, string col_aa, double[] values)> x)> pssm_to_vector400(List<info_blast_pssm_entry> pssm, enum_pssm_value_type pssm_value_type, bool normalise_row_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector400);

            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string row_aa, string col_aa, double[] values)>)>();

            foreach (var alphabet in feature_calcs.aa_alphabets)
            {
                
                var vector = new List<(string row_aa, string col_aa, double[] values)>();


                for (var c1_index = 0; c1_index < alphabet.groups.Count; c1_index++)
                {
                    var c1 = alphabet.groups[c1_index];

                    for (var c2_index = 0; c2_index < alphabet.groups.Count; c2_index++)
                    {
                        //if (c2_index > c1_index) continue;
                        
                        var c2 = alphabet.groups[c2_index];
                        double[] x = (pssm == null || pssm.Count == 0) ? new[] {0d} : pssm.Where(a => c1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c2.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)).Select(a => a.score).ToArray();
                        if (x == null || x.Length == 0) x = new[] {0d};
                        if (normalise_row_col) x = normalise_array(x);

                        if (x == null || x.Length == 0) x = new[] {0d};


                        if (pssm_value_type == enum_pssm_value_type.distances)
                        {
                            x = distances(x, normalise_row_col);
                        }
                        else if (pssm_value_type == enum_pssm_value_type.intervals)
                        {
                            x = intervals(x, normalise_row_col);
                        }

                        if (x == null || x.Length == 0) x = new[] {0d};


                        vector.Add((c2.group_amino_acids, c1.group_amino_acids, x));
                    }
                }

                if (normalise_all)
                {
                    var values = vector.SelectMany(a => a.values).Distinct().ToList();
                    var min = values.Min();
                    var max = values.Max();
                    vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                }

                result.Add((alphabet.name, vector));
            }

            return result;
        }


        internal static List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)> x)> pssm_to_vector400_DT(List<info_blast_pssm_entry> pssm, int max_lag, enum_pssm_value_type pssm_value_type, bool normalise_row_col = false, bool normalise_all = false)
        {
            const string module_name = nameof(info_blast_pssm);
            const string method_name = nameof(pssm_to_vector400_DT);

            //var cols = $@"ARNDCQEGHILKMFPSTWYV".ToCharArray(); //pssm.Select(a => a.position_aa).Distinct().ToList();

            var result = new List<(string alphabet, List<(string row_aa, string col_aa, int lag, double[] values)>)>();

            var tasks = new List<Task<(string name, List<(string row_aa, string col_aa, int lag, double[] values)> vector)>>();
            var tasks_start_time = DateTime.Now;
            
            foreach (var l_alphabet in feature_calcs.aa_alphabets)
            {
                var alphabet = l_alphabet;

                for (var l_lag = 1; l_lag <= max_lag; l_lag++)
                {
                    var lag = l_lag;

                    Task<(string name, List<(string row_aa, string col_aa, int lag, double[] values)> vector)> task = Task.Run(() =>
                    {
                        var vector = new List<(string row_aa, string col_aa, int lag, double[] values)>();


                        foreach (var c1 in alphabet.groups)
                        {
                            foreach (var c2 in alphabet.groups)
                            {

                                double[] x = new[] {0d};

                                if (pssm != null && pssm.Count > 0)
                                {
                                    var z = pssm.Where(a => c1.group_amino_acids.Contains(a.position_aa, StringComparison.Ordinal) && c2.group_amino_acids.Contains(a.query_sequence_aa, StringComparison.Ordinal)).ToList();
                                    z = z.OrderBy(a => a.matrix_row_index).ToList();
                                    z = z.Where(a => z.Any(b => Math.Abs(a.matrix_row_index - b.matrix_row_index) == lag)).ToList();
                                    if (z.Count > 0)
                                    {
                                        x = z.Select(a => a.score).ToArray();
                                    }
                                }


                                if (x == null || x.Length == 0) x = new[] {0d};

                                if (normalise_row_col) x = normalise_array(x);

                                if (x == null || x.Length == 0) x = new[] {0d};


                                if (pssm_value_type == enum_pssm_value_type.distances)
                                {
                                    x = distances(x, normalise_row_col);
                                }
                                else if (pssm_value_type == enum_pssm_value_type.intervals)
                                {
                                    x = intervals(x, normalise_row_col);
                                }

                                if (x == null || x.Length == 0) x = new[] {0d};


                                vector.Add((c2.group_amino_acids, c1.group_amino_acids, lag, x));

                            }
                        }

                        if (normalise_all)
                        {
                            var values = vector.SelectMany(a => a.values).Distinct().ToList();
                            var min = values.Min();
                            var max = values.Max();
                            vector.ForEach(a => a.values = normalise_array(a.values, min, max));
                        }

                        //result.Add((alphabet.name, vector));
                        return (alphabet.name, vector);
                    });
                    tasks.Add(task);
                }
            }

            //Task.WaitAll(tasks.ToArray<Task>());
            program.wait_tasks(tasks.ToArray<Task>(), tasks_start_time, module_name, method_name);

            tasks.ForEach(a=>result.Add(a.Result));

            return result;
        }



    }
}
