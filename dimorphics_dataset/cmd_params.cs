using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;

namespace dimorphics_dataset
{

    internal class cmd_params
    {
        internal bool parse_ok;

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

        internal string get_output_file_tag()
        {
            var x = new List<string>();

            // 1d tags
            if (feature_types_1d_interface?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1i);
            }

            if (feature_types_1d_neighbourhood?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1n);
            }

            if (feature_types_1d_chain?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1p);
            }

            // 2d tags
            if (feature_types_2d_interface?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2i);
            }

            if (feature_types_2d_neighbourhood?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2n);
            }

            if (feature_types_2d_chain?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2p);
            }

            // 3d tags
            if (feature_types_3d_interface?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3i);
            }

            if (feature_types_3d_neighbourhood?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3n);
            }

            if (feature_types_3d_chain?.key_value_list()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3p);
            }

            return String.Join("", x);
        }


        internal static List<(string key, string value)> get_params(string[] args)
        {
            var x = new List<(string key, string value)>();

            if (args == null || args.Length == 0)
            {
                return x;
            }

            args = args.SelectMany(a => a.Split(new char[] { '=' }, StringSplitOptions.RemoveEmptyEntries)).Where(a => !string.IsNullOrWhiteSpace(a)).ToArray();

            var name = "";
            var value = "";

            for (var i = 0; i < args.Length; i++)
            {
                var arg = args[i];

                var starts_dash = arg[0] == '-';
                var is_num = double.TryParse(arg, NumberStyles.Float, NumberFormatInfo.InvariantInfo, out var val_out);

                var is_name = starts_dash && (!is_num || string.IsNullOrEmpty(name)) && arg.Length > 1;
                var is_value = !is_name;
                var is_final_index = i == args.Length - 1;

                if (is_name)
                {
                    if (!string.IsNullOrWhiteSpace(name) || !string.IsNullOrWhiteSpace(value))
                    {
                        x.Add((name, value));
                    };

                    name = arg.Substring(1);
                    value = "";
                }

                if (is_value) value += (value.Length > 0 ? " " : "") + arg;

                if (is_final_index)
                {
                    x.Add((name, value));
                }
            }

            return x;
        }

        internal const string area1i = @"1i";
        internal const string area1n = @"1n";
        internal const string area1p = @"1p";

        internal const string area2i = @"2i";
        internal const string area2n = @"2n";
        internal const string area2p = @"2p";

        internal const string area3i = @"3i";
        internal const string area3n = @"3n";
        internal const string area3p = @"3p";

        internal static readonly string[] defined_areas = { area1i, area1n, area1p, area2i, area2n, area2p, area3i, area3n, area3p };

        internal feature_types_1d feature_types_1d_interface = new feature_types_1d(enum_protein_data_source.interface_1d, false);
        internal feature_types_1d feature_types_1d_neighbourhood = new feature_types_1d(enum_protein_data_source.neighbourhood_1d, false);
        internal feature_types_1d feature_types_1d_chain = new feature_types_1d(enum_protein_data_source.chain_1d, false);

        internal feature_types_2d feature_types_2d_interface = new feature_types_2d(enum_protein_data_source.interface_2d, false);
        internal feature_types_2d feature_types_2d_neighbourhood = new feature_types_2d(enum_protein_data_source.neighbourhood_2d, false);
        internal feature_types_2d feature_types_2d_chain = new feature_types_2d(enum_protein_data_source.chain_2d, false);
        
        internal feature_types_3d feature_types_3d_interface = new feature_types_3d(enum_protein_data_source.interface_3d, false);
        internal feature_types_3d feature_types_3d_neighbourhood = new feature_types_3d(enum_protein_data_source.neighbourhood_3d, false);
        internal feature_types_3d feature_types_3d_chain = new feature_types_3d(enum_protein_data_source.chain_3d, false);
        


        internal cmd_params(string[] args)
        {
            const string module_name = nameof(cmd_params);
            const string method_name = nameof(cmd_params);

            const string standard_coil = @"standard_coil";
            const string dimorphic_coil = @"dimorphic_coil";

            var param_list = get_params(args);
            io_proxy.WriteLine(string.Join($@", ", param_list), nameof(program), nameof(cmd_params));

          
            if (param_list.Count == 0) 
            {
                var classes = new[] { (class_id: +1, class_name: dimorphic_coil), (class_id: -1, class_name: standard_coil) };
                var output_folder = @"e:\dataset\";
                program.verbose = true;
                var exe = Path.GetFileName(Process.GetCurrentProcess().MainModule.FileName);
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} usage:");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{exe} -area=[{string.Join(",", defined_areas)}] -{nameof(use_dssp3)}=[true|false] -{nameof(class_id)}=[-1|+1] -{nameof(class_name)}=[class_name] -{nameof(min_sequence_length)}=[3] -{nameof(max_features)}=[100] -{nameof(output_folder)}=[path] [-{nameof(verbose)}=true|false] [-{nameof(first_index)}=[i] -{nameof(last_index)}=[j]]");
                io_proxy.WriteLine($@"");
                io_proxy.WriteLine($@"{nameof(dimorphics_dataset)} examples:");
                io_proxy.WriteLine($@"");

                foreach (var defined_area in defined_areas)
                {
                    io_proxy.WriteLine($@"{defined_area[0]}d {(defined_area[1] == 'i' || defined_area[1] == 'I' ? $@"interface subsequence" : $@"")}{(defined_area[1] == 'n' || defined_area[1] == 'N' ? $@"neighbourhood" : $@"")}{(defined_area[1] == 'p' || defined_area[1] == 'P' ? $@"protein chain" : $@"")} area:");

                    foreach (var b in classes)
                    {
                        var stdout_file = Path.Combine(output_folder, $@"stdout_{defined_area}_({(b.class_id > 0 ? $@"+" : $@"")}{b.class_id})_({b.class_name}).txt");
                        var stderr_file = Path.Combine(output_folder, $@"stderr_{defined_area}_({(b.class_id > 0 ? $@"+" : $@"")}{b.class_id})_({b.class_name}).txt");

                        io_proxy.WriteLine($@"{exe} -area={defined_area} -{nameof(use_dssp3)}=true -{nameof(class_id)}={(b.class_id > 0 ? $@"+" : $@"")}{b.class_id} -{nameof(class_name)}={b.class_name} -{nameof(min_sequence_length)}=3 -{nameof(max_features)}=100 -{nameof(output_folder)}={output_folder} -{nameof(use_children)}=true -{nameof(verbose)}=true 1> {stdout_file} 2> {stderr_file}");
                    }

                    io_proxy.WriteLine($@"");
                }

                Environment.Exit(0);
                //return default;
                //throw new Exception($@"No {nameof(args)} specified.");
            }

            if (param_list.Any(a => string.Equals(a.key, area1i, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1i, StringComparison.InvariantCultureIgnoreCase));
                feature_types_1d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area1n, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1n, StringComparison.InvariantCultureIgnoreCase));
                feature_types_1d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area1p, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1p, StringComparison.InvariantCultureIgnoreCase));
                feature_types_1d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2i, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2i, StringComparison.InvariantCultureIgnoreCase));
                feature_types_2d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2n, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2n, StringComparison.InvariantCultureIgnoreCase));
                feature_types_2d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2p, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2p, StringComparison.InvariantCultureIgnoreCase));
                feature_types_2d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3i, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3i, StringComparison.InvariantCultureIgnoreCase));
                feature_types_3d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3n, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3n, StringComparison.InvariantCultureIgnoreCase));
                feature_types_3d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3p, StringComparison.InvariantCultureIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3p, StringComparison.InvariantCultureIgnoreCase));
                feature_types_3d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(v.value, "true", StringComparison.InvariantCultureIgnoreCase));
            }

            feature_types_1d_interface
              .set_enable(param_list
                  .Where(a => a.key.StartsWith($"{area1i}.", StringComparison.InvariantCultureIgnoreCase) &&
                              feature_types_1d.keys.Contains(a.key.Substring($"{area1i}.".Length))
                  )
                  .Select(a =>
                      (
                          enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                          name: a.key.Substring($"{area1i}.".Length)
                      )
                  )
                  .ToArray()
              );

            feature_types_1d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area1n}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_1d.keys.Contains(a.key.Substring($"{area1n}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area1n}.".Length)
                        )
                    )
                    .ToArray()
                );


            feature_types_1d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area1p}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_1d.keys.Contains(a.key.Substring($"{area1p}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area1p}.".Length)
                        )
                    )
                    .ToArray()
                );

            feature_types_2d_interface
                     .set_enable(param_list
                         .Where(a => a.key.StartsWith($"{area2i}.", StringComparison.InvariantCultureIgnoreCase) &&
                                     feature_types_2d.keys.Contains(a.key.Substring($"{area2i}.".Length))
                         )
                         .Select(a =>
                             (
                                 enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                                 name: a.key.Substring($"{area2i}.".Length)
                             )
                         )
                         .ToArray()
                     );

            feature_types_2d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area2n}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_2d.keys.Contains(a.key.Substring($"{area2n}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area2n}.".Length)
                        )
                    )
                    .ToArray()
                );


            feature_types_2d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area2p}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_2d.keys.Contains(a.key.Substring($"{area2p}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area2p}.".Length)
                        )
                    )
                    .ToArray()
                );


            feature_types_3d_interface
                     .set_enable(param_list
                         .Where(a => a.key.StartsWith($"{area3i}.", StringComparison.InvariantCultureIgnoreCase) &&
                                     feature_types_3d.keys.Contains(a.key.Substring($"{area3i}.".Length))
                         )
                         .Select(a =>
                             (
                                 enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                                 name: a.key.Substring($"{area3i}.".Length)
                             )
                         )
                         .ToArray()
                     );

            feature_types_3d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area3n}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_3d.keys.Contains(a.key.Substring($"{area3n}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area3n}.".Length)
                        )
                    )
                    .ToArray()
                );


            feature_types_3d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith($"{area3p}.", StringComparison.InvariantCultureIgnoreCase) &&
                                feature_types_3d.keys.Contains(a.key.Substring($"{area3p}.".Length))
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.InvariantCultureIgnoreCase) || string.Equals(a.value, "true", StringComparison.InvariantCultureIgnoreCase)),
                            name: a.key.Substring($"{area3p}.".Length)
                        )
                    )
                    .ToArray()
                );

            


            if (param_list.Any(a => string.Equals(a.key, nameof(area), StringComparison.InvariantCultureIgnoreCase))) area = param_list.FirstOrDefault(z => z.key == nameof(area)).value?.ToLowerInvariant().Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries).SelectMany(a => a.Select((b, i) => i % 2 == 0 && i < a.Length - 1 ? a.Substring(i, 2) : $@"").Where(b => !string.IsNullOrWhiteSpace(b)).ToList()).OrderBy(a => a).ToArray() ?? null;
            if (param_list.Any(a => string.Equals(a.key, nameof(use_dssp3), StringComparison.InvariantCultureIgnoreCase))) use_dssp3 = bool.Parse(param_list.FirstOrDefault(z => z.key == nameof(use_dssp3)).value);
            if (param_list.Any(a => string.Equals(a.key, nameof(class_id), StringComparison.InvariantCultureIgnoreCase))) class_id = int.Parse(param_list.FirstOrDefault(z => z.key == nameof(class_id)).value, NumberStyles.Integer, CultureInfo.InvariantCulture);
            if (param_list.Any(a => string.Equals(a.key, nameof(class_name), StringComparison.InvariantCultureIgnoreCase))) class_name = param_list.FirstOrDefault(z => z.key == nameof(class_name)).value?.ToLowerInvariant();
            if (param_list.Any(a => string.Equals(a.key, nameof(min_sequence_length), StringComparison.InvariantCultureIgnoreCase))) min_sequence_length = int.Parse(param_list.FirstOrDefault(z => z.key == nameof(min_sequence_length)).value, NumberStyles.Integer, CultureInfo.InvariantCulture);
            if (param_list.Any(a => string.Equals(a.key, nameof(max_features), StringComparison.InvariantCultureIgnoreCase))) max_features = int.Parse(param_list.FirstOrDefault(z => z.key == nameof(max_features)).value, NumberStyles.Integer, CultureInfo.InvariantCulture);
            if (param_list.Any(a => string.Equals(a.key, nameof(output_folder), StringComparison.InvariantCultureIgnoreCase))) output_folder = param_list.FirstOrDefault(z => z.key == nameof(output_folder)).value;
            if (param_list.Any(a => string.Equals(a.key, nameof(first_index), StringComparison.InvariantCultureIgnoreCase))) first_index = int.TryParse(param_list.FirstOrDefault(z => z.key == nameof(first_index)).value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_first_index) ? tp_first_index : (int?)null;
            if (param_list.Any(a => string.Equals(a.key, nameof(last_index), StringComparison.InvariantCultureIgnoreCase))) last_index = int.TryParse(param_list.FirstOrDefault(z => z.key == nameof(last_index)).value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var tp_last_index) ? tp_last_index : (int?)null;
            if (param_list.Any(a => string.Equals(a.key, nameof(verbose), StringComparison.InvariantCultureIgnoreCase))) verbose = bool.TryParse(param_list.FirstOrDefault(z => z.key == nameof(verbose)).value, out var tp_verbose) ? tp_verbose : (bool?)null;
            if (param_list.Any(a => string.Equals(a.key, nameof(use_children), StringComparison.InvariantCultureIgnoreCase))) use_children = bool.TryParse(param_list.FirstOrDefault(z => z.key == nameof(use_children)).value, out var tp_use_children) ? tp_use_children : (bool?)null;


            verbose ??= (first_index == null && last_index == null);
            first_index ??= last_index;
            last_index ??= first_index;
            use_children ??= true;

            //if (args_list.Any())
            //{
            //    io_proxy.WriteLine($@"Args unknown: {string.Join($@" ", args_list)}", nameof(program), nameof(cmd_params));
            //    Environment.Exit(0);
            //}


            if (first_index == null && last_index == null)
            {
                io_proxy.WriteLine();

                io_proxy.WriteLine($@"{nameof(area)} = ""{string.Join($@", ", area ?? Array.Empty<string>())}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(use_dssp3)} = ""{use_dssp3}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(class_id)} = ""{class_id}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(class_name)} = ""{class_name}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(min_sequence_length)} = ""{min_sequence_length}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(max_features)} = ""{max_features}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(output_folder)} = ""{output_folder}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(first_index)} = ""{first_index}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(last_index)} = ""{last_index}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(verbose)} = ""{verbose}""", module_name, method_name);
                io_proxy.WriteLine($@"{nameof(use_children)} = ""{use_children}""", module_name, method_name);

                io_proxy.WriteLine();

                feature_types_1d_interface.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_1d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_1d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_1d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));
                
                io_proxy.WriteLine();

                feature_types_1d_chain.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_1d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_2d_interface.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_2d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_2d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_2d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_2d_chain.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_2d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_3d_interface.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_3d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_3d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_3d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));

                io_proxy.WriteLine();

                feature_types_3d_chain.key_value_list().ForEach(a => io_proxy.WriteLine($@"{nameof(feature_types_3d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));
            }

            if (area == null || area.Any(string.IsNullOrWhiteSpace)) throw new ArgumentNullException(nameof(args), nameof(area));
            if (area.Except(defined_areas).Any()) throw new ArgumentOutOfRangeException(nameof(args), nameof(area));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.key == nameof(use_dssp3)).value)) throw new ArgumentNullException(nameof(args), nameof(use_dssp3));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.key == nameof(class_id)).value)) throw new ArgumentNullException(nameof(args), nameof(class_id));
            if (string.IsNullOrWhiteSpace(class_name)) throw new ArgumentNullException(nameof(args), nameof(class_name));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.key == nameof(min_sequence_length)).value)) throw new ArgumentNullException(nameof(args), nameof(min_sequence_length));
            if (string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => a.key == nameof(max_features)).value)) throw new ArgumentNullException(nameof(args), nameof(max_features));
            if (string.IsNullOrWhiteSpace(output_folder)) throw new ArgumentNullException(nameof(args), nameof(output_folder));


            parse_ok = true;
        }
    }

}
/*internal static feature_types feature_types_params(string[] area)
{
    var module_name = nameof(dimorphics_dataset.feature_types);
    var method_name = nameof(feature_types_params);

    if (area == null || area.Length == 0)
    {
        return null;
    }

    area = area.SelectMany(a => a.Split(new char[] { ' ', ';', ',', '-', '|' }, StringSplitOptions.RemoveEmptyEntries)).Where(a => !string.IsNullOrWhiteSpace(a)).ToArray();


    var do_1d_interface = area?.Any(a => string.Equals(a, cmd_params.area1i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_1d_nh = area?.Any(a => string.Equals(a, cmd_params.area1n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_1d_protein = area?.Any(a => string.Equals(a, cmd_params.area1p, StringComparison.InvariantCultureIgnoreCase)) ?? false;

    var do_2d_interface = area?.Any(a => string.Equals(a, cmd_params.area2i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_2d_nh = area?.Any(a => string.Equals(a, cmd_params.area2n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_2d_protein = area?.Any(a => string.Equals(a, cmd_params.area2p, StringComparison.InvariantCultureIgnoreCase)) ?? false;

    var do_3d_interface = area?.Any(a => string.Equals(a, cmd_params.area3i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_3d_nh = area?.Any(a => string.Equals(a, cmd_params.area3n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
    var do_3d_protein = area?.Any(a => string.Equals(a, cmd_params.area3p, StringComparison.InvariantCultureIgnoreCase)) ?? false;


    var feature_types = new feature_types()
    {
        // 1d
        feature_types_1d_interface = !do_1d_interface
            ? null
            : new feature_types_1d(enum_protein_data_source.interface_1d)
            {
                pse_aac = do_1d_interface,
                geometry = do_1d_interface,
                        //mpsa_classification_data_subsequence = do_1d_interface,
                        iupred2a = do_1d_interface,
                aaindex = do_1d_interface,
                sable = do_1d_interface,
                stackdppred = false, //must be false - 1d protein level only
                        blast_pssm = do_1d_interface,
                r_peptides = do_1d_interface,
                r_protr = do_1d_interface,
            },
        feature_types_neighbourhood_1d = !do_1d_nh
            ? null
            : new feature_types_1d(enum_protein_data_source.neighbourhood_1d)
            {
                pse_aac = do_1d_nh,
                geometry = do_1d_nh,
                        //mpsa_classification_data_subsequence = do_1d_nh,
                        iupred2a = do_1d_nh,
                aaindex = do_1d_nh,
                sable = do_1d_nh,
                stackdppred = false, //must be false - 1d protein level only
                        blast_pssm = do_1d_nh,
                r_peptides = do_1d_nh,
                r_protr = do_1d_nh,
            },
        feature_types_chain_1d = !do_1d_protein
            ? null
            : new feature_types_1d(enum_protein_data_source.chain_1d)
            {
                pse_aac = do_1d_protein,
                geometry = do_1d_protein,
                        //mpsa_classification_data_subsequence = do_1d_protein,
                        iupred2a = do_1d_protein,
                aaindex = do_1d_protein,
                sable = do_1d_protein,
                stackdppred = do_1d_protein,
                blast_pssm = do_1d_protein,
                r_peptides = do_1d_protein,
                r_protr = do_1d_protein,
            },

        // 2d
        feature_types_2d_interface = !do_2d_interface
            ? null
            : new feature_types_2d(enum_protein_data_source.interface_2d)
            {
                mpsa = do_2d_interface,
            },
        feature_types_neighbourhood_2d = !do_2d_nh
            ? null
            : new feature_types_2d(enum_protein_data_source.neighbourhood_2d)
            {

                mpsa = do_2d_nh,

            },
        feature_types_chain_2d = !do_2d_protein
            ? null
            : new feature_types_2d(enum_protein_data_source.chain_2d)
            {
                mpsa = do_2d_protein,

            },


        // 3d

        feature_types_3d_interface = !do_3d_interface
            ? null
            : new feature_types_3d(enum_protein_data_source.interface_3d)
            {
                sasa = do_3d_interface,
                intramolecular = do_3d_interface,

                foldx = do_3d_interface,
                tortuosity = do_3d_interface,
                ring = do_3d_interface,
                pse_ssc_dssp = do_3d_interface,
            },
        feature_types_neighbourhood_3d = !do_3d_nh
            ? null
            : new feature_types_3d(enum_protein_data_source.neighbourhood_3d)
            {
                sasa = do_3d_nh,
                intramolecular = do_3d_nh,

                foldx = do_3d_nh,
                tortuosity = do_3d_nh,
                ring = do_3d_nh,
                pse_ssc_dssp = do_3d_nh,
            },
        feature_types_chain_3d = !do_3d_protein
            ? null
            : new feature_types_3d(enum_protein_data_source.chain_3d)
            {
                sasa = do_3d_protein,
                intramolecular = do_3d_protein,

                foldx = do_3d_protein,
                tortuosity = do_3d_protein,
                ring = do_3d_protein,
                pse_ssc_dssp = do_3d_protein,
            }
    };

    // 1d
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_1d_interface)} = {feature_types.feature_types_1d_interface}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_1d)} = {feature_types.feature_types_neighbourhood_1d}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_1d)} = {feature_types.feature_types_chain_1d}", module_name, method_name);

    //2d
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_2d_interface)} = {feature_types.feature_types_2d_interface}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_2d)} = {feature_types.feature_types_neighbourhood_2d}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_2d)} = {feature_types.feature_types_chain_2d}", module_name, method_name);

    //3d
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_3d_interface)} = {feature_types.feature_types_3d_interface}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}", module_name, method_name);
    io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_3d)} = {feature_types.feature_types_chain_3d}", module_name, method_name);
    return feature_types;
}*/