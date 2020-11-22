using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace dimorphics_dataset
{

    internal class cmd_params
    {
        internal string[] area;
        internal bool use_dssp3;
        internal int class_id;
        internal string class_name;
        internal int min_sequence_length;
        internal int max_features;
        internal string output_folder;
        internal int? first_index;
        internal int? last_index;
        internal bool? verbose;
        //internal string tag;
        internal bool? use_children;

        internal cmd_params()
        {

        }

        internal cmd_params(cmd_params p)
        {
            area = p.area;
            use_dssp3 = p.use_dssp3;
            class_id = p.class_id;
            class_name = p.class_name;
            min_sequence_length = p.min_sequence_length;
            max_features = p.max_features;
            output_folder = p.output_folder;
            first_index = p.first_index;
            last_index = p.last_index;
            verbose = p.verbose;
            use_children = p.use_children;
            //tag = p.tag;
        }

        internal static string get_param(string name, List<string> args)
        {
            if (args == null)
            {
                return null;
            }

            var ix = args.FindIndex(a => string.Equals(a, $"-{name}", StringComparison.InvariantCultureIgnoreCase));

            if (ix == -1) return null;

            var value = (ix == args.Count - 1) ? "" : args[ix + 1];

            if (ix < args.Count - 1) args.RemoveAt(ix + 1);

            args.RemoveAt(ix);

            return value;
        }

        internal const string area1i = "1i";
        internal const string area1n = "1n";
        internal const string area1p = "1p";
        internal const string area2i = "2i";
        internal const string area2n = "2n";
        internal const string area2p = "2p";
        internal const string area3i = "3i";
        internal const string area3n = "3n";
        internal const string area3p = "3p";

        internal static readonly string[] defined_areas = {area1i, area1n, area1p, area2i, area2n, area2p, area3i, area3n, area3p};

        internal static cmd_params get_params(string[] args)
        {
            const string standard_coil = "standard_coil";
            const string dimorphic_coil = "dimorphic_coil";

            

            if (args == null || args.Length == 0) // || (args.Length == 1 && args[0] == "-test"))
            {
                var test = args?.Any(a => a == "-test") ?? false;

                
                var classes = new[] {(class_id: +1, class_name: dimorphic_coil), (class_id: -1, class_name: standard_coil)};
                var output_folder = @"e:\dataset\";
                program.verbose = true;
                var exe = Path.GetFileName(Process.GetCurrentProcess().MainModule.FileName);
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} usage:");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{exe} -area=[{string.Join(",", defined_areas)}] -use_dssp3=[true|false] -class_id=[-1|+1] -class_name=[class_name] -min_sequence_length=[3] -max_features=[100] -output_folder=[path] [-verbose=true|false] [-tag=tag] [-first_index=[i] -last_index=[j]]");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} examples:");
                io_proxy.WriteLine($@"");

                foreach (var defined_area in defined_areas)
                {
                    io_proxy.WriteLine($@"{defined_area[0]}d {(defined_area[1] == 'i' || defined_area[1] == 'I' ? "interface subsequence" : "")}{(defined_area[1] == 'n' || defined_area[1] == 'N' ? "neighbourhood" : "")}{(defined_area[1] == 'p' || defined_area[1] == 'P' ? "protein chain" : "")} area:");
                    //io_proxy.WriteLine($@"{exe} -area={a} -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{(test ? $@"test\" : "")} -verbose=true 1> e:\dataset\{(test ? $@"test\" : "")}stdout_{a}_dimorphic_coil.txt 2> e:\dataset\{(test ? $@"test\" : "")}stderr_{a}_dimorphic_coil.txt");
                    //io_proxy.WriteLine($@"{exe} -area={a} -use_dssp3=true -class_id=-1 -class_name=standard_coil  -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{(test ? $@"test\" : "")} -verbose=true 1> e:\dataset\{(test ? $@"test\" : "")}stdout_{a}_standard_coil.txt  2> e:\dataset\{(test ? $@"test\" : "")}stderr_{a}_standard_coil.txt");

                    foreach (var b in classes)
                    {
                        var stdout_file = Path.Combine(output_folder, $@"stdout_{defined_area}_({(b.class_id > 0 ? "+" : "")}{b.class_id})_({b.class_name}).txt");
                        var stderr_file = Path.Combine(output_folder, $@"stderr_{defined_area}_({(b.class_id > 0 ? "+" : "")}{b.class_id})_({b.class_name}).txt");

                        io_proxy.WriteLine($@"{exe} -area={defined_area} -use_dssp3=true -class_id={(b.class_id > 0 ? "+" : "")}{b.class_id} -class_name={b.class_name} -min_sequence_length=3 -max_features=100 -output_folder={output_folder} -use_children=true -verbose=true 1> {stdout_file} 2> {stderr_file}");
                    }

                    io_proxy.WriteLine($@"");
                }

                //Environment.Exit(0);
                //return default;
                //throw new Exception($@"No {nameof(args)} specified.");
                return null;
            }

            io_proxy.WriteLine(string.Join(" ", args), nameof(program), nameof(get_params));

            var args_list = args.SelectMany(a => a.Split(new char[] {'='}, StringSplitOptions.RemoveEmptyEntries)).ToList();

            var param_names = new[]
            {
                nameof(area),
                nameof(use_dssp3),
                nameof(class_id),
                nameof(class_name),
                nameof(min_sequence_length),
                nameof(max_features),
                nameof(output_folder),
                nameof(first_index),
                nameof(last_index),
                nameof(verbose),
                nameof(use_children),
            };

            var param_list = param_names.Select(a => (name: a, value: get_param(a, args_list))).ToList();

            var ret = new cmd_params()
            {
                area = param_list.First(z=>z.name==nameof(area)).value?.ToLowerInvariant().Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries).SelectMany(a => a.Select((b, i) => i % 2 == 0 && i < a.Length - 1 ? a.Substring(i, 2) : "").Where(b => !string.IsNullOrWhiteSpace(b)).ToList()).OrderBy(a => a).ToArray() ?? null,
                use_dssp3 = bool.Parse(param_list.First(z=>z.name==nameof(use_dssp3)).value),
                class_id = int.Parse(param_list.First(z=>z.name==nameof(class_id)).value, NumberStyles.Integer, CultureInfo.InvariantCulture),
                class_name = param_list.First(z=>z.name==nameof(class_name)).value?.ToLowerInvariant(),
                min_sequence_length = int.Parse(param_list.First(z=>z.name==nameof(min_sequence_length)).value, NumberStyles.Integer, CultureInfo.InvariantCulture),
                max_features = int.Parse(param_list.First(z=>z.name==nameof(max_features)).value, NumberStyles.Integer, CultureInfo.InvariantCulture),
                output_folder = param_list.First(z=>z.name==nameof(output_folder)).value,
                first_index = int.TryParse(param_list.First(z=>z.name==nameof(first_index)).value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_first_index) ? tp_first_index : (int?)null,
                last_index = int.TryParse(param_list.First(z=>z.name==nameof(last_index)).value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_last_index) ? tp_last_index : (int?)null,
                verbose = bool.TryParse(param_list.First(z=>z.name==nameof(verbose)).value, out var tp_verbose) ? tp_verbose : (bool?)null,
                use_children = bool.TryParse(param_list.First(z=>z.name==nameof(use_children)).value, out var tp_use_children) ? tp_use_children : (bool?)null,
            };

            ret.verbose ??= (ret.first_index == null && ret.last_index == null);
            ret.first_index ??= ret.last_index;
            ret.last_index ??= ret.first_index;
            ret.use_children ??= true;

            if (args_list.Any())
            {
                io_proxy.WriteLine($@"Args unknown: {string.Join(" ", args_list)}", nameof(program),
                    nameof(get_params));

                Environment.Exit(0);
            }

            if (ret.area == null || ret.area.Any(string.IsNullOrWhiteSpace)) throw new ArgumentNullException(nameof(args), nameof(ret.area));
            if (ret.area.Except(defined_areas).Any()) throw new ArgumentOutOfRangeException(nameof(args), nameof(ret.area));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.name == nameof(use_dssp3)).value)) throw new ArgumentNullException(nameof(args), nameof(ret.use_dssp3));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.name == nameof(class_id)).value)) throw new ArgumentNullException(nameof(args), nameof(ret.class_id));
            if (string.IsNullOrWhiteSpace(ret.class_name)) throw new ArgumentNullException(nameof(args), nameof(ret.class_name));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.name == nameof(min_sequence_length)).value)) throw new ArgumentNullException(nameof(args), nameof(ret.min_sequence_length));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.name == nameof(max_features)).value)) throw new ArgumentNullException(nameof(args), nameof(ret.max_features));
            if (string.IsNullOrWhiteSpace(ret.output_folder)) throw new ArgumentNullException(nameof(args), nameof(ret.output_folder));


            if (ret.first_index == null && ret.last_index == null)
            {
                io_proxy.WriteLine($@"{nameof(ret.area)} = ""{string.Join(", ", ret.area)}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.use_dssp3)} = ""{ret.use_dssp3}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.class_id)} = ""{ret.class_id}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.class_name)} = ""{ret.class_name}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.min_sequence_length)} = ""{ret.min_sequence_length}""",
                    nameof(program), nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.max_features)} = ""{ret.max_features}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.output_folder)} = ""{ret.output_folder}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.first_index)} = ""{ret.first_index}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.last_index)} = ""{ret.last_index}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.verbose)} = ""{ret.verbose}""", nameof(program),
                    nameof(get_params));
                io_proxy.WriteLine($@"{nameof(ret.use_children)} = ""{ret.use_children}""", nameof(program),
                    nameof(get_params));
                //io_proxy.WriteLine($@"{nameof(ret.tag)} = ""{ret.tag}""", nameof(program), nameof(get_params));
            }

            return ret;
        }
    }

}
