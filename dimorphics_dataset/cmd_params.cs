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
        internal string tag;

        public cmd_params()
        {
            
        }

        public cmd_params(cmd_params p)
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
            tag = p.tag;
        }

        public static string get_param(string name, List<string> args)
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

        public static cmd_params get_params(string[] args)
        {
            if (args == null || args.Length == 0 || (args.Length == 1 && args[0] == "-test"))
            {
                var test = args?.Any(a => a == "-test") ?? false;

                program.verbose = true;
                var exe = Path.GetFileName(Process.GetCurrentProcess().MainModule.FileName);
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} usage:");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine(
                    $@"{exe} -area=[2i,2n,2p,3i,3n,3p] -use_dssp3=[true|false] -class_id=[-1|+1] -class_name=[class_name] -min_sequence_length=[3] -max_features=[100] -output_folder=[path] [-verbose=true|false] [-tag=tag] [-first_index=[i] -last_index=[j]]");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} examples:");
                io_proxy.WriteLine($@"");

                foreach (var a in new[] { "2i", "2n", "2p", "3i", "3n", "3p" })
                {
                    io_proxy.WriteLine(
                        $@"{a[0]}d {(a[1] == 'i' ? "interface subsequence" : "")}{(a[1] == 'n' ? "neighbourhood" : "")}{(a[1] == 'p' ? "protein" : "")} area:");
                    io_proxy.WriteLine(
                        $@"{exe} -area={a} -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{(test ? $@"test\" : "")}{a}\dimorphic_coil\ -verbose=true 1> e:\dataset\{(test ? $@"test\" : "")}{a}\dimorphic_coil\stdout.txt 2> e:\dataset\{(test ? $@"test\" : "")}{a}\dimorphic_coil\stderr.txt");
                    io_proxy.WriteLine(
                        $@"{exe} -area={a} -use_dssp3=true -class_id=-1 -class_name=standard_coil  -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\{(test ? $@"test\" : "")}{a}\standard_coil\ -verbose=true 1> e:\dataset\{(test ? $@"test\" : "")}{a}\standard_coil\stdout.txt 2> e:\dataset\{(test ? $@"test\" : "")}{a}\standard_coil\stderr.txt");
                    io_proxy.WriteLine($@"");
                }

                Environment.Exit(0);
                //return default;
                //throw new Exception($@"No {nameof(args)} specified.");
            }

            io_proxy.WriteLine(string.Join(" ", args), nameof(program), nameof(get_params));

            var args_list = args.SelectMany(a => a.Split(new char[] { '=' }, StringSplitOptions.RemoveEmptyEntries))
                .ToList();

            var ret = new cmd_params()
            {
                area = get_param(nameof(cmd_params.area), args_list)?.ToLowerInvariant()
                    .Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries).OrderBy(a => a)
                    .ToArray() ?? null,
                use_dssp3 = bool.Parse(get_param(nameof(cmd_params.use_dssp3), args_list)),
                class_id = int.Parse(get_param(nameof(cmd_params.class_id), args_list), NumberStyles.Integer,
                    CultureInfo.InvariantCulture),
                class_name = get_param(nameof(cmd_params.class_name), args_list)?.ToLowerInvariant(),
                min_sequence_length = int.Parse(get_param(nameof(cmd_params.min_sequence_length), args_list), NumberStyles.Integer,
                    CultureInfo.InvariantCulture),
                max_features = int.Parse(get_param(nameof(cmd_params.max_features), args_list), NumberStyles.Integer,
                    CultureInfo.InvariantCulture),
                output_folder = get_param(nameof(cmd_params.output_folder), args_list),
                first_index = int.TryParse(get_param(nameof(cmd_params.first_index), args_list), NumberStyles.Integer,
                    CultureInfo.InvariantCulture, out var tp_first_index)
                    ? tp_first_index
                    : (int?)null,
                last_index = int.TryParse(get_param(nameof(cmd_params.last_index), args_list), NumberStyles.Integer,
                    CultureInfo.InvariantCulture, out var tp_last_index)
                    ? tp_last_index
                    : (int?)null,
                verbose = (bool.TryParse(get_param(nameof(cmd_params.verbose), args_list), out var tp_verbose)
                    ? tp_verbose
                    : (bool?)null),
                tag = get_param(nameof(cmd_params.tag), args_list)
            };

            ret.verbose ??= (ret.first_index == null && ret.last_index == null);
            ret.first_index ??= ret.last_index;
            ret.last_index ??= ret.first_index;

            if (args_list.Any())
            {
                io_proxy.WriteLine($@"Args unknown: {string.Join(" ", args_list)}", nameof(program),
                    nameof(get_params));

                Environment.Exit(0);
            }

            if (ret.area == null || ret.area.Any(a => string.IsNullOrWhiteSpace(a)))
                throw new ArgumentNullException(nameof(ret.area));
            if (ret.area.Except(new string[] { "2i", "2n", "2p", "3i", "3n", "3p" }).Any())
                throw new ArgumentOutOfRangeException(nameof(ret.area));
            if (ret.use_dssp3 == null) throw new ArgumentNullException(nameof(ret.use_dssp3));
            if (ret.class_id == null) throw new ArgumentNullException(nameof(ret.class_id));
            if (string.IsNullOrWhiteSpace(ret.class_name)) throw new ArgumentNullException(nameof(ret.class_name));
            if (ret.min_sequence_length == null) throw new ArgumentNullException(nameof(ret.min_sequence_length));
            if (ret.max_features == null) throw new ArgumentNullException(nameof(ret.max_features));
            if (string.IsNullOrWhiteSpace(ret.output_folder))
                throw new ArgumentNullException(nameof(ret.output_folder));


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
                io_proxy.WriteLine($@"{nameof(ret.tag)} = ""{ret.tag}""", nameof(program), nameof(get_params));
            }

            return ret;
        }
    }

}
