using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class info_foldx
    {
        internal static string foldx_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"));
        internal static string pdb_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"pdb"));






        private static readonly object file_write_lock = new object();

        internal static foldx_energy_differences load_calc_energy_differences(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run, string source = @"all", bool write_bat = false)//int nh_first_res_id, int nh_last_res_id)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(load_calc_energy_differences);

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"{nameof(calc_energy_differences)}(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run);");
#endif

            //var nh_first_res_id = res_ids.Min(a => a.residue_index);
            //var nh_last_res_id = res_ids.Max(a => a.residue_index);

            var result = new foldx_energy_differences();
            // fragments should be first residue to last residue (even non-interacting residues), since we want to mutate only some of them

            // returns energy differences upon mutation of a neighbourhood
            // can be viewed as: difference over whole monomer, difference over fragment only, and also, repaired/unrepaired, but we only do repaired.

            //var pdb_folder = Path.Combine(program.data_root_folder,"pdb\";
            //var foldx_folder = Path.Combine(program.data_root_folder,"foldx\";
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            // make monomer from dimer
            var monomer_file = atom.extract_split_pdb_chains(pdb_id, chain_id).First();

            // repair
            var monomer_file_repair = foldx_repair_pdb(Path.GetFileNameWithoutExtension(monomer_file), run);

            var repair_res_ids = io_proxy.ReadAllLines(monomer_file_repair, module_name, method_name).Where(a => a.StartsWith(/*program.string_debug*/($@"ATOM"), StringComparison.Ordinal)).Select(a => int.Parse(a.Substring(22, 4), NumberStyles.Integer, NumberFormatInfo.InvariantInfo)).Distinct().OrderBy(a => a).ToList();

            // filter res ids to remove res ides which were in the original pdb structure file but removed by the foldx repair
            res_ids = res_ids.Where(a => repair_res_ids.Contains(a.residue_index)).ToList();

            // calculate energy differences before/after mutation

            //result.foldx_ala_scanning_result_protein = foldx_caller.load_foldx_ala_scanning(monomer_file_repair, chain_id, null, run);
            result.foldx_ala_scanning_result_subsequence = info_foldx.load_foldx_ala_scanning(monomer_file_repair, chain_id, res_ids, run);
            result.foldx_position_scanning_result_subsequence = info_foldx.load_foldx_position_scanning((monomer_file_repair, chain_id, res_ids), run);
            result.foldx_buildmodel_position_scan_result_subsequence = info_foldx.load_foldx_buildmodel_position_scan((monomer_file_repair, chain_id, res_ids), run);
            result.foldx_buildmodel_subsequence_mutant_result_subsequence = info_foldx.load_foldx_buildmodel_subsequence_mutant((monomer_file_repair, chain_id, res_ids), run);

            if (write_bat)
            {
                lock (file_write_lock)
                {
                    var fn1 = Path.Combine(/*program.string_debug*/($@"{foldx_folder}"), /*program.string_debug*/($@"foldx_calc_ala_scanning_{source}.bat.skip"));
                    io_proxy.AppendAllLines(fn1, new[] { /*program.string_debug*/($@"if not exist ""{result.foldx_ala_scanning_result_subsequence.wait_filename}"" {result.foldx_ala_scanning_result_subsequence.cmd_line}") }, module_name, method_name);

                    var fn2 = Path.Combine(/*program.string_debug*/($@"{foldx_folder}"), /*program.string_debug*/($@"foldx_calc_position_scanning_{source}.bat"));
                    io_proxy.AppendAllLines(fn2, new[] { /*program.string_debug*/($@"if not exist ""{result.foldx_position_scanning_result_subsequence.wait_filename}"" {result.foldx_position_scanning_result_subsequence.cmd_line}") }, module_name, method_name);

                    var fn3 = Path.Combine(/*program.string_debug*/($@"{foldx_folder}"), /*program.string_debug*/($@"foldx_calc_buildmodel_position_scan_{source}.bat"));
                    io_proxy.AppendAllLines(fn3, new[] { /*program.string_debug*/($@"if not exist ""{result.foldx_buildmodel_position_scan_result_subsequence.wait_filename}"" {result.foldx_buildmodel_position_scan_result_subsequence.cmd_line}") }, module_name, method_name);

                    var fn4 = Path.Combine(/*program.string_debug*/($@"{foldx_folder}"), /*program.string_debug*/($@"foldx_calc_buildmodel_subsequence_mutant_{source}.bat"));
                    io_proxy.AppendAllLines(fn4, new[] { /*program.string_debug*/($@"if not exist ""{result.foldx_buildmodel_subsequence_mutant_result_subsequence.wait_filename}"" {result.foldx_buildmodel_subsequence_mutant_result_subsequence.cmd_line}") }, module_name, method_name);
                }
            }

            return result;
        }






        internal static readonly
            List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>
            foldx_residues =
                new
                    List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>()
                    { // add MSE/MSA ?? the special MET substitute

                        (00, /*program.string_debug*/($@"glycine"), /*program.string_debug*/($@"GLY"), 'G', /*program.string_debug*/($@"GLY"), 'G', true, /*program.string_debug*/($@"standard aa")),
                        (01, /*program.string_debug*/($@"alanine"), /*program.string_debug*/($@"ALA"), 'A', /*program.string_debug*/($@"ALA"), 'A', true, /*program.string_debug*/($@"standard aa")),
                        (02, /*program.string_debug*/($@"leucine"), /*program.string_debug*/($@"LEU"), 'L', /*program.string_debug*/($@"LEU"), 'L', true, /*program.string_debug*/($@"standard aa")),
                        (03, /*program.string_debug*/($@"valine"), /*program.string_debug*/($@"VAL"), 'V', /*program.string_debug*/($@"VAL"), 'V', true, /*program.string_debug*/($@"standard aa")),
                        (04, /*program.string_debug*/($@"isoleucine"), /*program.string_debug*/($@"ILE"), 'I', /*program.string_debug*/($@"ILE"), 'I', true, /*program.string_debug*/($@"standard aa")),
                        (05, /*program.string_debug*/($@"proline"), /*program.string_debug*/($@"PRO"), 'P', /*program.string_debug*/($@"PRO"), 'P', true, /*program.string_debug*/($@"standard aa")),
                        (06, /*program.string_debug*/($@"arginine"), /*program.string_debug*/($@"ARG"), 'R', /*program.string_debug*/($@"ARG"), 'R', true, /*program.string_debug*/($@"standard aa")),
                        (07, /*program.string_debug*/($@"threonine"), /*program.string_debug*/($@"THR"), 'T', /*program.string_debug*/($@"THR"), 'T', true, /*program.string_debug*/($@"standard aa")),
                        (08, /*program.string_debug*/($@"serine"), /*program.string_debug*/($@"SER"), 'S', /*program.string_debug*/($@"SER"), 'S', true, /*program.string_debug*/($@"standard aa")),
                        (09, /*program.string_debug*/($@"cysteine"), /*program.string_debug*/($@"CYS"), 'C', /*program.string_debug*/($@"CYS"), 'C', true, /*program.string_debug*/($@"standard aa")),
                        (10, /*program.string_debug*/($@"methionine"), /*program.string_debug*/($@"MET"), 'M', /*program.string_debug*/($@"MET"), 'M', true, /*program.string_debug*/($@"methionine")),
                        (11, /*program.string_debug*/($@"lysine"), /*program.string_debug*/($@"LYS"), 'K', /*program.string_debug*/($@"LYS"), 'K', true, /*program.string_debug*/($@"standard aa")),
                        (12, /*program.string_debug*/($@"glutamic"), /*program.string_debug*/($@"GLU"), 'E', /*program.string_debug*/($@"GLU"), 'E', true, /*program.string_debug*/($@"standard aa")),
                        (13, /*program.string_debug*/($@"glutamine"), /*program.string_debug*/($@"GLN"), 'Q', /*program.string_debug*/($@"GLN"), 'Q', true, /*program.string_debug*/($@"standard aa")),
                        (14, /*program.string_debug*/($@"aspartic"), /*program.string_debug*/($@"ASP"), 'D', /*program.string_debug*/($@"ASP"), 'D', true, /*program.string_debug*/($@"standard aa")),
                        (15, /*program.string_debug*/($@"asparagine"), /*program.string_debug*/($@"ASN"), 'N', /*program.string_debug*/($@"ASN"), 'N', true, /*program.string_debug*/($@"standard aa")),
                        (16, /*program.string_debug*/($@"tryptophane"), /*program.string_debug*/($@"TRP"), 'W', /*program.string_debug*/($@"TRP"), 'W', true, /*program.string_debug*/($@"standard aa")),
                        (17, /*program.string_debug*/($@"tyrosine"), /*program.string_debug*/($@"TYR"), 'Y', /*program.string_debug*/($@"TYR"), 'Y', true, /*program.string_debug*/($@"standard aa")),
                        (18, /*program.string_debug*/($@"phenylalanine"), /*program.string_debug*/($@"PHE"), 'F', /*program.string_debug*/($@"PHE"), 'F', true, /*program.string_debug*/($@"standard aa")),
                        (19, /*program.string_debug*/($@"histidine"), /*program.string_debug*/($@"HIS"), 'H', /*program.string_debug*/($@"HIS"), 'H', true, /*program.string_debug*/($@"standard aa")),
                        (20, /*program.string_debug*/($@"phoshoporylated threonine"), /*program.string_debug*/($@"PTR"), 'y', /*program.string_debug*/($@"THR"), 'T', true, /*program.string_debug*/($@"phoshoporylated threonine")),
                        (21, /*program.string_debug*/($@"phosphorylated tyrosine"), /*program.string_debug*/($@"TPO"), 'p',  /*program.string_debug*/($@"TYR"), 'Y', true, /*program.string_debug*/($@"phosphorylated tyrosine")),
                        (22, /*program.string_debug*/($@"phosphorylated serine"), /*program.string_debug*/($@"SEP"), 's', /*program.string_debug*/($@"SER"), 'S', true, /*program.string_debug*/($@"phosphorylated serine")),
                        (23, /*program.string_debug*/($@"hydroxiproline"), /*program.string_debug*/($@"HYP"), 'h', /*program.string_debug*/($@"PRO"), 'P', true, /*program.string_debug*/($@"hydroxiproline")),
                        (24, /*program.string_debug*/($@"sulfotyrosine"), /*program.string_debug*/($@"TYS"), 'z', /*program.string_debug*/($@"TYR"), 'Y', true, /*program.string_debug*/($@"sulfotyrosine")),
                        (25, /*program.string_debug*/($@"monomethylated lysine"), /*program.string_debug*/($@"MLZ"), 'k', /*program.string_debug*/($@"LYS"), 'K', true, /*program.string_debug*/($@"monomethylated lysine")),
                        (26, /*program.string_debug*/($@"dimethylated lysine"), /*program.string_debug*/($@"MLY"), 'm', /*program.string_debug*/($@"LYS"), 'K', true, /*program.string_debug*/($@"dimethylated lysine")),
                        (27, /*program.string_debug*/($@"trimethylated lysine"), /*program.string_debug*/($@"M3L"), 'l', /*program.string_debug*/($@"LYS"), 'K', true, /*program.string_debug*/($@"trimethylated lysine")),
                        (28, /*program.string_debug*/($@"charged ND1 histidine"), /*program.string_debug*/($@"H1S"), 'o', /*program.string_debug*/($@"HIS"), 'H', true, /*program.string_debug*/($@"charged ND1 histidine")),
                        (29, /*program.string_debug*/($@"charged NE2 histidine"), /*program.string_debug*/($@"H2S"), 'e', /*program.string_debug*/($@"HIS"), 'H', true, /*program.string_debug*/($@"charged NE2 histidine")),
                        (30, /*program.string_debug*/($@"neutral histidine"), /*program.string_debug*/($@"H3S"), 'f', /*program.string_debug*/($@"HIS"), 'H', true, /*program.string_debug*/($@"neutral histidine")),
                        (31, /*program.string_debug*/($@"adenosine"), /*program.string_debug*/($@"A"), 'a', /*program.string_debug*/($@"A"), 'a', true, /*program.string_debug*/($@"adenosine")),
                        (32, /*program.string_debug*/($@"guanosine"), /*program.string_debug*/($@"G"), 'g', /*program.string_debug*/($@"G"), 'g', true, /*program.string_debug*/($@"guanosine")),
                        (33, /*program.string_debug*/($@"cytosine"), /*program.string_debug*/($@"C"), 'c', /*program.string_debug*/($@"C"), 'c', true, /*program.string_debug*/($@"cytosine")),
                        (34, /*program.string_debug*/($@"thymidine"), /*program.string_debug*/($@"T"), 't',/*program.string_debug*/($@"T"), 't', true, /*program.string_debug*/($@"thymidine")),
                        (35, /*program.string_debug*/($@"6-methylated adenosine"), /*program.string_debug*/($@"6MA"), 'b', /*program.string_debug*/($@"6MA"), 'b', true, /*program.string_debug*/($@"6-methylated adenosine")),
                        (36, /*program.string_debug*/($@"5-methylated cytosine"), /*program.string_debug*/($@"5CM"), 'd', /*program.string_debug*/($@"5CM"), 'd', true, /*program.string_debug*/($@"5-methylated cytosine")),
                        (37, /*program.string_debug*/($@"adenosine triphosphate"), /*program.string_debug*/($@"ATP"), ' ', /*program.string_debug*/($@"ATP"), ' ', false, /*program.string_debug*/($@"adenosine triphosphate")),
                        (38, /*program.string_debug*/($@"adenosine diphosphate"), /*program.string_debug*/($@"ADP"), ' ', /*program.string_debug*/($@"ADP"), ' ', false, /*program.string_debug*/($@"adenosine diphosphate")),
                        (39, /*program.string_debug*/($@"guanosine triphosphate"), /*program.string_debug*/($@"GTP"), ' ', /*program.string_debug*/($@"GTP"), ' ', false, /*program.string_debug*/($@"guanosine triphosphate")),
                        (40, /*program.string_debug*/($@"guanosine diphosphate"), /*program.string_debug*/($@"GDP"), ' ',/*program.string_debug*/($@"GDP"), ' ', false, /*program.string_debug*/($@"guanosine diphosphate")),
                        (41, /*program.string_debug*/($@"dioctadecylglycerol-3-phosphatidyl-choline"), /*program.string_debug*/($@"LIP"), ' ',/*program.string_debug*/("LIP"), ' ', false, /*program.string_debug*/($@"dioctadecylglycerol-3-phosphatidyl-choline")),
                        (42, /*program.string_debug*/($@"calcium"), /*program.string_debug*/($@"CA"), ' ',/*program.string_debug*/($@"CA"), ' ', false, /*program.string_debug*/($@"calcium")),
                        (43, /*program.string_debug*/($@"magnesium"), /*program.string_debug*/($@"MG"), ' ',/*program.string_debug*/($@"MG"), ' ', false, /*program.string_debug*/($@"magnesium")),
                        (44, /*program.string_debug*/($@"manganese"), /*program.string_debug*/($@"MN"), ' ',/*program.string_debug*/($@"MN"), ' ', false, /*program.string_debug*/($@"manganese")),
                        (45, /*program.string_debug*/($@"sodium"), /*program.string_debug*/($@"NA"), ' ',/*program.string_debug*/($@"NA"), ' ', false, /*program.string_debug*/($@"sodium")),
                        (46, /*program.string_debug*/($@"zinc"), /*program.string_debug*/($@"ZN"), ' ',/*program.string_debug*/($@"ZN"), ' ', false, /*program.string_debug*/($@"zinc")),
                        (47, /*program.string_debug*/($@"iron"), /*program.string_debug*/($@"FE"), ' ',/*program.string_debug*/($@"FE"), ' ', false, /*program.string_debug*/($@"iron")),
                        (48, /*program.string_debug*/($@"copper"), /*program.string_debug*/($@"CU"), ' ',/*program.string_debug*/($@"CU"), ' ', false, /*program.string_debug*/($@"copper")),
                        (49, /*program.string_debug*/($@"cobalt"), /*program.string_debug*/($@"CO"), ' ',/*program.string_debug*/($@"CO"), ' ', false, /*program.string_debug*/($@"cobalt")),
                        (50, /*program.string_debug*/($@"potasium"), /*program.string_debug*/($@"K"), ' ',/*program.string_debug*/($@"K"), ' ', false, /*program.string_debug*/($@"potasium")),
                        (51, /*program.string_debug*/($@"water"), /*program.string_debug*/($@"HOH"), ' ',/*program.string_debug*/($@"HOH"), ' ', false, /*program.string_debug*/($@"water")),
                    };

        internal static readonly List<(int index, string full_name, string foldx_aa_code3, char foldx_aa_code1, string standard_aa_code3, char standard_aa_code1, bool is_mutable, string residue_type)>
            foldx_residues_aa_mutable = foldx_residues.Where(a => (a.index >= 0) && (a.index <= 30)).ToList();

        private static readonly List<string> call_foldx_lock = new List<string>();

        internal static string[] read_all_lines_until_success(string wait_file, int max_tries = int.MaxValue, int delay_ms = 10)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(read_all_lines_until_success);

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"{nameof(read_all_lines_until_success)}(string wait_file, int max_tries = int.MaxValue, int delay_ms = 10);");
#endif

            var file_accessed = false;
            var access_attempt_count = 0;
            string[] data = null;

            while (!file_accessed)
            {
                if (access_attempt_count >= max_tries) throw new IOException(/*program.string_debug*/($@"{module_name}.{method_name}: File could not be accessed ""{wait_file}"""));

                access_attempt_count++;

                try
                {

                    if (File.Exists(wait_file))// && new FileInfo(wait_file).Length > 0)
                    {
                        var f_len1 = new FileInfo(wait_file).Length;

                        data = io_proxy.ReadAllLines(wait_file, module_name, method_name);

                        var f_len2 = new FileInfo(wait_file).Length;

                        if (f_len1 != f_len2)
                        {
                            data = null;

                            Task.Delay(delay_ms).Wait();

                        }
                        else
                        {
                            file_accessed = true;

                        }
                    }
                    else
                    {
                        Task.Delay(delay_ms).Wait();
                    }

                }
                catch (IOException e)
                {
                    io_proxy.log_exception(e, $@"", module_name, method_name);

                    Task.Delay(delay_ms).Wait();
                }
            }

            return data;
        }

        internal static (string cmd_line, string[] data) call_foldx(string pdb_file_id, string foldx_command, string foldx_args, string wait_filename, string lock_code_local, bool run, string pdb_folder = null)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(call_foldx);

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"call_foldx(pdb_file_id={pdb_file_id}, foldx_command={foldx_command}, foldx_args={foldx_args}, wait_filename={wait_filename}, lock_code_local={lock_code_local}");
#endif
            if (string.IsNullOrWhiteSpace(pdb_folder))
            {
                pdb_folder = info_foldx.pdb_folder;
            }

            var foldx_exe = /*program.string_debug*/($@"{foldx_folder}foldx.exe");
            var pdb_file = /*program.string_debug*/($@"{pdb_folder}{Path.GetFileNameWithoutExtension(pdb_file_id)}.pdb");


            var lock_code = string.Join(/*program.string_debug*/($@"_"), new string[] { pdb_file_id, foldx_command, foldx_args, wait_filename, lock_code_local });


            var wait_for_complete = false;

            lock (call_foldx_lock)
            {
                if (call_foldx_lock.Contains(lock_code))
                {
                    wait_for_complete = true;
                }
                else
                {
                    call_foldx_lock.Add(lock_code);
                }
            }

            while (wait_for_complete)
            {
                Task.Delay(10).Wait();

                lock (call_foldx_lock)
                {
                    wait_for_complete = call_foldx_lock.Contains(lock_code);

                    if (!wait_for_complete)
                    {
                        call_foldx_lock.Add(lock_code);
                    }
                }
            }

            string[] data = null;


            if (string.IsNullOrWhiteSpace(foldx_args))
            {
                //throw new ArgumentNullException(nameof(foldx_args));
            }

            var foldx_args2 = foldx_args;
            if (!foldx_args2.Contains(/*program.string_debug*/($@"--command="), StringComparison.Ordinal)) { foldx_args2 += /*program.string_debug*/($@" --command={foldx_command}"); }
            if (!foldx_args2.Contains(/*program.string_debug*/($@"--pdb="), StringComparison.Ordinal)) { foldx_args2 += /*program.string_debug*/($@" --pdb={Path.GetFileName(pdb_file)}"); }
            if (!foldx_args2.Contains(/*program.string_debug*/($@"--pdb-dir="), StringComparison.Ordinal)) { foldx_args2 += /*program.string_debug*/($@" --pdb-dir=""{Path.GetDirectoryName(pdb_file) ?? /*program.string_debug*/($@".")}"""); }
            if (!foldx_args2.Contains(/*program.string_debug*/($@"--out-pdb="), StringComparison.Ordinal) && !string.Equals(foldx_command, /*program.string_debug*/($@"RepairPDB"), StringComparison.OrdinalIgnoreCase)) { foldx_args2 += /*program.string_debug*/($@" --out-pdb=false"); }

            var start = new ProcessStartInfo { FileName = foldx_exe, WorkingDirectory = Path.GetDirectoryName(foldx_exe) ?? /*program.string_debug*/($@""), Arguments = foldx_args2, UseShellExecute = false, CreateNoWindow = false, RedirectStandardOutput = true, RedirectStandardError = true };
            var cmd_line = /*program.string_debug*/($@"""{start.FileName}"" {start.Arguments}");

            if (string.IsNullOrWhiteSpace(wait_filename) || !File.Exists(wait_filename) || new FileInfo(wait_filename).Length <= 0)
            {

                if (run)
                {
                    io_proxy.WriteLine(/*program.string_debug*/($@"{module_name}.{method_name}: run: ""{start.FileName}"" {start.Arguments}"), module_name, method_name);

                    using (var process = Process.Start(start))
                    {
                        if (process == null) throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: {nameof(process)} is null"));


                        using (var reader = process.StandardOutput)
                        {
                            var stdout = reader.ReadToEnd();
                            stdout = stdout.Replace(/*program.string_debug*/($"\r\n"), /*program.string_debug*/($"\r\n{module_name}.{method_name}: {nameof(stdout)}: "), StringComparison.Ordinal);
                            if (!string.IsNullOrWhiteSpace(stdout))
                            {
                                io_proxy.WriteLine(/*program.string_debug*/($@"{module_name}.{method_name}: {nameof(stdout)}: {stdout}"), module_name, method_name);
                            }

                            var stderr = process.StandardError.ReadToEnd();
                            stderr = stderr.Replace(/*program.string_debug*/($"\r\n"), /*program.string_debug*/($"\r\n{module_name}.{method_name}: {nameof(stderr)}: "), StringComparison.Ordinal);
                            if (!string.IsNullOrWhiteSpace(stderr))
                            {
                                io_proxy.WriteLine(/*program.string_debug*/($@"{module_name}.{method_name}: {nameof(stderr)}: {stderr}"), module_name, method_name);
                            }
                        }

                        process.WaitForExit();
                    }

                    if (run)
                    {
                        data = read_all_lines_until_success(wait_filename);
                    }
                    else
                    {
                        data = File.Exists(wait_filename) && new FileInfo(wait_filename).Length > 0 ? io_proxy.ReadAllLines(wait_filename, module_name, method_name) : Array.Empty<string>();
                    }
                }
            }
            else
            {
                if (!string.IsNullOrWhiteSpace(wait_filename))
                {
                    if (run)
                    {
                        data = read_all_lines_until_success(wait_filename);
                    }
                    else
                    {
                        data = File.Exists(wait_filename) && new FileInfo(wait_filename).Length > 0 ? io_proxy.ReadAllLines(wait_filename, module_name, method_name) : Array.Empty<string>();
                    }
                }
            }


            lock (call_foldx_lock)
            {
                call_foldx_lock.Remove(lock_code);
            }

            return (cmd_line, data);
        }

        internal static string foldx_repair_pdb(string pdb_id, bool run, string pdb_folder = null, string repair_folder = null)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(foldx_repair_pdb);

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"foldx_repair_pdb(pdb_id = ""{pdb_id}"")");
#endif
            if (string.IsNullOrWhiteSpace(pdb_folder))
            {
                pdb_folder = info_foldx.pdb_folder;
            }

            if (string.IsNullOrWhiteSpace(repair_folder))
            {
                repair_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"foldx"), /*program.string_debug*/($@"pdb"));
            }

            var pdb_file = Path.Combine(/*program.string_debug*/($@"{pdb_folder}"), /*program.string_debug*/($@"{Path.GetFileNameWithoutExtension(pdb_id)}.pdb"));
            var pdb_file_repair = Path.Combine(/*program.string_debug*/($@"{repair_folder}"), /*program.string_debug*/($@"{Path.GetFileNameWithoutExtension(pdb_id)}_Repair.pdb"));

            if (File.Exists(pdb_file_repair) && new FileInfo(pdb_file_repair).Length > 0) return pdb_file_repair;

            var foldx_cmd = /*program.string_debug*/($@"RepairPDB");
            var foldx_args = /*program.string_debug*/($@"--output-dir=""{Path.GetDirectoryName(pdb_file) ?? /*program.string_debug*/($@".")}""");
            var call_foldx_result = call_foldx(pdb_id, foldx_cmd, foldx_args, pdb_file_repair, pdb_id, run, pdb_folder);

            return pdb_file_repair;
        }




        internal static (string cmd_line, string wait_filename, List<foldx_energy_terms_ps> data) load_foldx_buildmodel_position_scan((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run, bool write_list = false)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(load_foldx_buildmodel_position_scan);

            var (pdb_id, chain_id, res_ids) = interface_residues;

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"load_foldx_buildmodel_position_scan(pdb_id = ""{pdb_id}"", chain_id = ""{chain_id}"", res_ids = ""{res_ids}"")");
#endif

            if (res_ids == null || res_ids.Count == 0) return (/*program.string_debug*/($@""), /*program.string_debug*/($@""), null);

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            var lock_code = /*program.string_debug*/($@"bm_ps_{pdb_id}{chain_id}{string.Join(/*program.string_debug*/($@"_"), res_ids)}");

            var foldx_mutation_positions_data = new List<string>();

            var mutation_positions_data = new List<(
                char original_amino_acid,
                char chain_id,
                int residue_index,
                char mutant_foldx_amino_acid1,
                string mutant_foldx_amino_acid3,
                char mutant_standard_amino_acid1,
                string mutant_standard_amino_acid3)>();

            foreach (var master_atom in res_ids)
            {
                foreach (var mutant_amino_acid in foldx_residues_aa_mutable)//.Select(b => b.foldx_aa_code1).ToList())
                {

                    var mp_data = (
                        original_amino_acid: master_atom.amino_acid,
                        chain_id,
                        master_atom.residue_index,
                        mutant_foldx_amino_acid1: mutant_amino_acid.foldx_aa_code1,
                        mutant_foldx_amino_acid3: mutant_amino_acid.foldx_aa_code3,
                        mutant_standard_amino_acid1: mutant_amino_acid.standard_aa_code1,
                        mutant_standard_amino_acid3: mutant_amino_acid.standard_aa_code3);

                    mutation_positions_data.Add(mp_data);

                    var mp = /*program.string_debug*/($@"{mp_data.original_amino_acid}{mp_data.chain_id}{mp_data.residue_index}{mp_data.mutant_foldx_amino_acid1};");
                    foldx_mutation_positions_data.Add(mp);
                }
            }

            var first_amino_acid = /*program.string_debug*/($@"{res_ids.First().amino_acid}{res_ids.First().residue_index}");
            var last_amino_acid = /*program.string_debug*/($@"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}");
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();



            var mutant_list_file = Path.Combine(foldx_folder, /*program.string_debug*/($@"bm_ps"), /*program.string_debug*/($@"individual_list_bm_ps_{pdb_id}{chain_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}.txt"));

            if (write_list)
            {
                lock (file_write_lock)
                {
                    //Directory.CreateDirectory(Path.GetDirectoryName(mutant_list_file));
                    io_proxy.WriteAllLines(mutant_list_file, foldx_mutation_positions_data, module_name, method_name);
                }
            }

            var foldx_cmd = /*program.string_debug*/($@"BuildModel");
            var foldx_output_pdb = false;
            var foldx_number_of_runs = 1;
            var output_file_tag = /*program.string_debug*/($@"bm_ps_{pdb_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}");
            var wait_filename = Path.Combine(foldx_folder, /*program.string_debug*/($@"bm_ps"), /*program.string_debug*/($@"Dif_{output_file_tag}_{pdb_id}.fxout"));
            var foldx_args = /*program.string_debug*/($@"--mutant-file={mutant_list_file} --numberOfRuns={foldx_number_of_runs} --out-pdb={foldx_output_pdb} --output-file={output_file_tag}");

            var foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = foldx_result_1.cmd_line;
            var foldx_result = foldx_result_1.data;

            if (foldx_result == null || foldx_result.Length == 0)
            {
                return (cmd_line, wait_filename, new List<foldx_energy_terms_ps>());
            }

            var sd = false;
            var marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith(/*program.string_debug*/($"Pdb\ttotal energy"), StringComparison.Ordinal));

            if (marker_index < 0)
            {
                marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith(/*program.string_debug*/($"Pdb\tSD\ttotal energy"), StringComparison.Ordinal));

                if (marker_index > -1)
                {
                    sd = true;
                }
            }

            var foldx_result2 = foldx_result.Skip(marker_index + 1).ToList();

            var results = foldx_result2.Select((a, i) =>
            {
                var b = a.Split();

                var j = 0;

                var c = new foldx_energy_terms_ps()
                {
                    line_index = i,
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    res_ids = res_ids,

                    mutation_positions_data = mutation_positions_data[i],

                    Pdb = b[j++],
                    SD = sd ? double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    total_energy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Backbone_Hbond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Sidechain_Hbond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Van_der_Waals = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Electrostatics = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Solvation_Polar = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Solvation_Hydrophobic = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Van_der_Waals_clashes = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    entropy_sidechain = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    entropy_mainchain = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    sloop_entropy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    mloop_entropy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    cis_bond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    torsional_clash = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    backbone_clash = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    helix_dipole = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    water_bridge = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    disulfide = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    electrostatic_kon = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    partial_covalent_bonds = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    energy_Ionisation = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Entropy_Complex = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo)

                };

                return c;
            }).ToList();

            if (results.Count != foldx_mutation_positions_data.Count)
            {
                throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: Missing foldx model results"));
            }

            return (cmd_line, wait_filename, results);
        }

        internal static (string cmd_line, string wait_filename, List<foldx_energy_terms_sm> data) load_foldx_buildmodel_subsequence_mutant((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run, bool save_list = false)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(load_foldx_buildmodel_subsequence_mutant);

            var (pdb_id, chain_id, res_ids) = interface_residues;

            if (res_ids == null || res_ids.Count == 0) return (/*program.string_debug*/($@""), /*program.string_debug*/($@""), null);

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"load_foldx_buildmodel_subsequence_mutant(pdb_id = ""{pdb_id}"", chain_id = ""{chain_id}"", res_ids = ""{res_ids}"")");
#endif
            // FoldX --command=BuildModel --pdb=BM.pdb --mutant-file=individual_list.txt
            // individual_list.txt = e.g. FA39L,FB39L; (res_aa, chain, res_num, mutant_aa)


            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            var lock_code = /*program.string_debug*/($@"bm_if_subs_{pdb_id}{chain_id}{string.Join(/*program.string_debug*/($@"_"), res_ids)}");

            var mutation_positions_data = new List<List<(char original_amino_acid1, char chain_id, int residue_index, char mutant_foldx_amino_acid1, string mutant_foldx_amino_acid3, char mutant_standard_amino_acid1, string mutant_standard_amino_acid3)>>();
            var foldx_mutation_positions_data = new List<string>();

            foreach (var mutant_amino_acid in foldx_residues_aa_mutable)
            {
                var mp_data = res_ids.Select(a => (original_amino_acid1: a.amino_acid, chain_id: chain_id, residue_index: a.residue_index, mutant_foldx_amino_acid1: mutant_amino_acid.foldx_aa_code1, mutant_foldx_amino_acid3: mutant_amino_acid.foldx_aa_code3, mutant_standard_amino_acid1: mutant_amino_acid.standard_aa_code1, mutant_standard_amino_acid3: mutant_amino_acid.standard_aa_code3)).ToList();
                mutation_positions_data.Add(mp_data);

                var mp = /*program.string_debug*/($@"{string.Join(/*program.string_debug*/($@","), mp_data.Select(a => /*program.string_debug*/($@"{a.original_amino_acid1}{a.chain_id}{a.residue_index}{a.mutant_foldx_amino_acid1}")).ToList())};");

                foldx_mutation_positions_data.Add(mp);
            }

            var first_amino_acid = /*program.string_debug*/($@"{res_ids.First().amino_acid}{res_ids.First().residue_index}");
            var last_amino_acid = /*program.string_debug*/($@"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}");
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var mutant_list_file = Path.Combine(foldx_folder, /*program.string_debug*/($@"bm_if_subs"), /*program.string_debug*/($@"individual_list_bm_if_subs_{pdb_id}{chain_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}.txt"));


            if (save_list)
            {
                lock (file_write_lock)
                {
                    //Directory.CreateDirectory(Path.GetDirectoryName(mutant_list_file));
                    io_proxy.WriteAllLines(mutant_list_file, foldx_mutation_positions_data, module_name, method_name);
                }
            }

            var foldx_cmd = /*program.string_debug*/($@"BuildModel");
            var foldx_output_pdb = false;
            var foldx_number_of_runs = 1;
            var output_file_tag = /*program.string_debug*/($@"bm_if_subs_{pdb_id}_{first_amino_acid}_{last_amino_acid}_{reside_index_sum}");
            var wait_filename = Path.Combine(foldx_folder, /*program.string_debug*/($@"bm_if_subs"), /*program.string_debug*/($@"Dif_{output_file_tag}_{pdb_id}.fxout"));

            var foldx_args = /*program.string_debug*/($@"--mutant-file={mutant_list_file} --numberOfRuns={foldx_number_of_runs} --out-pdb={foldx_output_pdb} --output-file={output_file_tag}");

            var foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = foldx_result_1.cmd_line;

            var foldx_result = foldx_result_1.data;

            if (foldx_result == null || foldx_result.Length <= 0)
            {
                return (cmd_line, wait_filename, new List<foldx_energy_terms_sm>());
            }

            var sd = false;
            var marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith(/*program.string_debug*/($"Pdb\ttotal energy"), StringComparison.Ordinal));

            if (marker_index < 0)
            {
                marker_index = foldx_result.ToList().FindIndex(a => a.StartsWith(/*program.string_debug*/($"Pdb\tSD\ttotal energy"), StringComparison.Ordinal));

                if (marker_index > -1)
                {
                    sd = true;
                }
            }

            var foldx_result2 = foldx_result.Skip(marker_index + 1).ToList();

            var results = foldx_result2.Select((a, i) =>
            {
                var b = a.Split();

                var j = 0;

                var foldx_energy_terms = new foldx_energy_terms_sm
                {
                    mutation_positions_data = mutation_positions_data[i],
                    line_index = i,
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    res_ids = res_ids,

                    Pdb = b[j++],
                    SD = sd ? double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    total_energy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Backbone_Hbond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Sidechain_Hbond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Van_der_Waals = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Electrostatics = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Solvation_Polar = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Solvation_Hydrophobic = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Van_der_Waals_clashes = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    entropy_sidechain = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    entropy_mainchain = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    sloop_entropy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    mloop_entropy = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    cis_bond = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    torsional_clash = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    backbone_clash = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    helix_dipole = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    water_bridge = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    disulfide = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    electrostatic_kon = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    partial_covalent_bonds = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    energy_Ionisation = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo),
                    Entropy_Complex = double.Parse(b[j++], NumberStyles.Float, NumberFormatInfo.InvariantInfo)
                };

                return foldx_energy_terms;

            }).ToList();

            if (results.Count != foldx_mutation_positions_data.Count)
            {
                throw new Exception(/*program.string_debug*/($@"{module_name}.{method_name}: Missing foldx model results"));
            }


            return (cmd_line, wait_filename, results);
        }

        internal static string fix_non_standard_naming(string res_name, bool fix)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(fix_non_standard_naming);

            if (string.IsNullOrWhiteSpace(res_name))
            {
                throw new ArgumentNullException(nameof(res_name));
            }

            if (fix)
            {
                if (res_name.Length == 3)
                {
                    var ri = foldx_residues.FindIndex(a => string.Equals(a.foldx_aa_code3, res_name, StringComparison.Ordinal));
                    if (ri > -1)
                    {
                        var r = foldx_residues[ri];

                        if (!string.Equals(r.foldx_aa_code3, r.standard_aa_code3, StringComparison.Ordinal))
                        {
                            return r.standard_aa_code3;
                        }

                        return res_name;
                    }
                }

                if (res_name.Length == 1)
                {
                    var ri = foldx_residues.FindIndex(a => a.foldx_aa_code1 == res_name[0]);
                    if (ri > -1)
                    {
                        var r = foldx_residues[ri];

                        if (r.foldx_aa_code1 != r.standard_aa_code1)
                        {
                            return r.standard_aa_code1.ToString(CultureInfo.InvariantCulture);
                        }

                        return res_name;
                    }
                }
            }

            return res_name;
        }

        internal static (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) load_foldx_ala_scanning(string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids, bool run)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(load_foldx_ala_scanning);

            //    var (pdb_id, chain_id, res_ids) = interface_residues;
#if DEBUG
            ////if (program.verbose_debug) Program.WriteLine(/*program.string_debug*/($@"load_foldx_ala_scanning(pdb_id = ""{pdb_id}"", chain_id = ""{chain_id}"", res_ids = ""{res_ids}"")");
#endif
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            //if (res_ids == null || res_ids.Count == 0) return (/*program.string_debug*/($@""), /*program.string_debug*/($@""), null);

            var wait_filename = Path.Combine(foldx_folder, "ala_scan", /*program.string_debug*/($@"{pdb_id}_AS.fxout"));

            string[] file_data = null;



            //var cmd_line = /*program.string_debug*/($@"");


            //if (!File.Exists(wait_filename))
            //{

            //var first_amino_acid = /*program.string_debug*/($@"{res_ids.First().amino_acid}{res_ids.First().residue_index}";
            //var last_amino_acid = /*program.string_debug*/($@"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}";
            //var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var lock_code = /*program.string_debug*/($@"{pdb_id}{chain_id}{(res_ids != null ? string.Join(/*program.string_debug*/($@"_"), res_ids) : /*program.string_debug*/($@""))}");

            var foldx_cmd = /*program.string_debug*/($@"AlaScan");
            var foldx_args = /*program.string_debug*/($@""); ///*program.string_debug*/($@"--output-file={file_tag}";

            var call_foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var cmd_line = call_foldx_result_1.cmd_line;

            var call_foldx_result = call_foldx_result_1.data;

            file_data = call_foldx_result;

            //}
            //else
            //{
            //    file_data = program.ReadAllLines(wait_filename);
            //}

            if (file_data == null || file_data.Length == 0)
            {
                return (cmd_line, wait_filename, new List<foldx_ala_scanning_result>());
            }

            //if (file_data != null && file_data.Length > 0)
            // {
            var split_lines = file_data.Select(b => b.Split()).ToList();

            //var fix_nsn = false;

            var results = split_lines.Select(c =>
                new foldx_ala_scanning_result()
                {
                    pdb_id = pdb_id,
                    chain_id = chain_id,
                    residue_index = int.Parse(c[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo),

                    //res_name = fix_non_standard_naming(c[0], fix_nsn),

                    original_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[0], StringComparison.OrdinalIgnoreCase)).foldx_aa_code1,
                    original_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[0], StringComparison.OrdinalIgnoreCase)).standard_aa_code1,
                    original_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[0], StringComparison.OrdinalIgnoreCase)).foldx_aa_code3,
                    original_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[0], StringComparison.OrdinalIgnoreCase)).standard_aa_code3,

                    mutant_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[3], StringComparison.OrdinalIgnoreCase)).foldx_aa_code1,
                    mutant_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[3], StringComparison.OrdinalIgnoreCase)).standard_aa_code1,
                    mutant_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[3], StringComparison.OrdinalIgnoreCase)).foldx_aa_code3,
                    mutant_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, c[3], StringComparison.OrdinalIgnoreCase)).standard_aa_code3,


                    ddg = double.Parse(c[7], NumberStyles.Float, NumberFormatInfo.InvariantInfo)
                }).ToList();



            results = results.Where(a => string.Equals(pdb_id, a.pdb_id, StringComparison.OrdinalIgnoreCase) && chain_id == a.chain_id).ToList();

            if (res_ids != null && res_ids.Count > 0)
            {
                results = results.Where(a => res_ids.Any(b => b.residue_index == a.residue_index)).ToList();
            }


            return (cmd_line, wait_filename, results);
        }



        internal static (string cmd_line, string wait_filename, List<foldx_position_scanning_result> data) load_foldx_position_scanning((string pdb_id, char chain_id, List<(int residue_index, char i_code, char amino_acid)> res_ids) interface_residues, bool run)
        {
            const string module_name = nameof(info_foldx);
            const string method_name = nameof(load_foldx_position_scanning);

            var (pdb_id, chain_id, res_ids) = interface_residues;

#if DEBUG
            //if (program.verbose_debug) io.WriteLine(/*program.string_debug*/($@"load_foldx_position_scanning(pdb_id = ""{pdb_id}"", chain_id = ""{chain_id}"", res_ids = ""{res_ids}"")");
#endif

            if (res_ids == null || res_ids.Count == 0) return (/*program.string_debug*/($@""), /*program.string_debug*/($@""), null);

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var lock_code = /*program.string_debug*/($@"{pdb_id}{chain_id}{string.Join(/*program.string_debug*/($@"_"), res_ids)}");

            var first_amino_acid = /*program.string_debug*/($@"{res_ids.First().amino_acid}{res_ids.First().residue_index}");
            var last_amino_acid = /*program.string_debug*/($@"{res_ids.Last().amino_acid}{res_ids.Last().residue_index}");
            var reside_index_sum = res_ids.Select(a => a.residue_index).Sum();

            var mutation_code = 'd';
            var mutation_position_ids = res_ids.Select(a => /*program.string_debug*/($@"{a.amino_acid}{chain_id}{a.residue_index}{mutation_code}")).ToList();
            var mutation_positions = /*program.string_debug*/($@"{string.Join(/*program.string_debug*/($@","), mutation_position_ids)}");
            var file_tag = /*program.string_debug*/($@"{pdb_id}_{mutation_position_ids.First()}_{mutation_position_ids.Last()}_{reside_index_sum}");
            var wait_filename = Path.Combine(foldx_folder, /*program.string_debug*/($@"ps"), /*program.string_debug*/($@"PS_{file_tag}_scanning_output.txt"));
            var foldx_cmd = /*program.string_debug*/($@"PositionScan");
            var foldx_args = /*program.string_debug*/($@"--positions={mutation_positions} --output-file={file_tag}");

            var call_foldx_result_1 = call_foldx(pdb_id, foldx_cmd, foldx_args, wait_filename, lock_code, run);
            var call_foldx_result = call_foldx_result_1.data;

            List<foldx_position_scanning_result> results = null;

            //var fix_nsn = false;

            if (call_foldx_result != null && call_foldx_result.Length > 0)
            {
                var split_lines = call_foldx_result.Select(b => b.Split()).ToList();

                results = split_lines.Select(line =>
                {

                    return new foldx_position_scanning_result()
                    {
                        pdb_id = pdb_id,
                        chain_id = chain_id,


                        original_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.OrdinalIgnoreCase)).foldx_aa_code1,
                        original_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.OrdinalIgnoreCase)).standard_aa_code1,

                        original_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.OrdinalIgnoreCase)).foldx_aa_code3,
                        original_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => string.Equals(a.foldx_aa_code3, line[0].Substring(0, 3), StringComparison.OrdinalIgnoreCase)).standard_aa_code3,

                        mutant_foldx_amino_acid_1 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][^1])).foldx_aa_code1,
                        mutant_standard_amino_acid_1 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][^1])).standard_aa_code1,
                        mutant_foldx_amino_acid_3 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][^1])).foldx_aa_code3,
                        mutant_standard_amino_acid_3 = foldx_residues_aa_mutable.First(a => a.foldx_aa_code1.Equals(line[0][^1])).standard_aa_code3,

                        //residue_index = int.Parse(line[0].Substring(4, line[0].Length - 5), NumberStyles.Integer, NumberFormatInfo.InvariantInfo),
                        residue_index = int.Parse(line[0][4..^1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo),
                        //res_mutant_aa = line[0][line[0].Length - 1],
                        ddg = double.Parse(line[1], NumberStyles.Float, NumberFormatInfo.InvariantInfo)
                    };
                }).ToList();
            }

            return (call_foldx_result_1.cmd_line, wait_filename, results);
        }
    }
}
